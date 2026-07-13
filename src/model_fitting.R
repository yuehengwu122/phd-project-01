# ============================================================
# Stan Model Compilation and Fitting
# ============================================================

# --- Model cache (compile once, reuse) ---
if (!exists(".compiled_models", envir = globalenv(), inherits = FALSE)) {
  .compiled_models <- new.env(parent = emptyenv())
}


#' Get a compiled Stan model (with caching)
#'
#' @param type One of "full", "linear", "simulation"
#' @param model_path Directory containing Stan files
#' @return Compiled stanmodel object
get_model <- function(
  type = c("full", "linear", "simulation"),
  model_path = "models/"
) {
  type <- match.arg(type)

  stan_file <- switch(type,
    full       = "fit_full_model.stan",
    linear     = "fit_linear_model.stan",
    simulation = "simulation.stan"
  )

  path <- file.path(model_path, stan_file)
  if (!file.exists(path)) stop("Stan file not found: ", path)
  path_norm <- normalizePath(path, mustWork = TRUE)
  path_mtime <- file.info(path_norm)$mtime

  entry <- .compiled_models[[type]]
  if (is.list(entry) &&
      !is.null(entry$model) &&
      identical(entry$path, path_norm) &&
      identical(entry$mtime, path_mtime)) {
    return(entry$model)
  }

  message("Compiling Stan model: ", type, " ...")
  model_obj <- rstan::stan_model(file = path_norm)
  .compiled_models[[type]] <- list(
    model = model_obj,
    path = path_norm,
    mtime = path_mtime
  )
  model_obj
}


#' Fit GP Regression Model
#'
#' Main fitting function. Accepts raw (unscaled) X and y,
#' rescales X to [-0.5, 0.5] internally, fits the Stan model,
#' and returns a structured result object.
#'
#' @param X Predictor matrix (N x P), original scale
#' @param y Response vector (length N)
#' @param model Model type: "full" (GP + linear) or "linear" (linear only)
#' @param d_order Moment prior order for delta (0 = local, 1-3 = non-local)
#' @param chains Number of MCMC chains (default: 2)
#' @param iter Total iterations per chain (default: 2000)
#' @param warmup Warmup iterations per chain (default: 500)
#' @param adapt_delta Target acceptance rate (default: 0.95)
#' @param max_treedepth Maximum tree depth (default: 12)
#' @param seed Random seed (default: 42)
#' @param refresh Print frequency (default: 500 = print every 500 iterations)
#' @param model_path Directory containing Stan files
#' @param ... Additional arguments passed to rstan::sampling()
#' @return A gp_fit object (list) with components:
#'   fit, model_type, d_order, stan_data, orig, coef_orig, predictor_names
fit_model <- function(
  X,
  y,
  model = c("full", "linear"),
  d_order = 0L,
  chains = 2L,
  iter = 2000L,
  warmup = 500L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  seed = 42L,
  refresh = 500L,
  model_path = "models/",
  ...
) {
  model <- match.arg(model)
  d_order <- as.integer(d_order)
  stopifnot(d_order %in% 0:3)

  # --- Prepare X ---
  X_orig <- as.matrix(X)
  if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)
  y <- as.numeric(y)
  stopifnot(nrow(X_orig) == length(y))

  N <- nrow(X_orig)
  P <- ncol(X_orig)

  # Preserve predictor names
  pred_names <- colnames(X_orig)
  if (is.null(pred_names)) pred_names <- paste0("X", seq_len(P))

  # Rescale X to [-0.5, 0.5]
  X_scaled <- apply(X_orig, 2, rescale_to_unit)
  if (is.vector(X_scaled)) X_scaled <- matrix(X_scaled, ncol = 1)
  colnames(X_scaled) <- pred_names

  # --- Build Stan data ---
  stan_data <- list(N = N, P = P, X = X_scaled, y = y)

  if (model == "full") {
    stan_data <- c(stan_data, moment_prior_stan_data(d_order))
  }

  # --- Compile & fit ---
  stan_mod <- get_model(model, model_path = model_path)

  fit <- rstan::sampling(
    object = stan_mod,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    seed = seed,
    refresh = refresh,
    verbose = FALSE,
    open_progress = FALSE,
    ...
  )

  # --- Back-transform linear coefficients to original scale ---
  post <- rstan::extract(fit)
  beta1_draws <- post$beta1

  y_bar <- mean(y)
  x_bar <- colMeans(X_orig)
  x_min <- apply(X_orig, 2, min)
  x_max <- apply(X_orig, 2, max)
  rng <- x_max - x_min

  # beta1_orig[j] = beta1_scaled[j] / (x_max[j] - x_min[j])
  if (is.matrix(beta1_draws)) {
    beta1_draws_orig <- sweep(beta1_draws, 2, rng, "/")
  } else {
    beta1_draws_orig <- matrix(beta1_draws / rng[1], ncol = 1)
  }

  # beta0 = y_bar - x_bar' * beta1_orig
  if (is.matrix(beta1_draws_orig)) {
    beta0_draws_orig <- as.vector(y_bar - beta1_draws_orig %*% x_bar)
  } else {
    beta0_draws_orig <- y_bar - beta1_draws_orig * x_bar
  }

  # --- Assemble result ---
  result <- list(
    fit = fit,
    model_type = model,
    d_order = d_order,
    stan_data = stan_data,
    orig = list(X = X_orig, y = y),
    coef_orig = list(
      beta0_draws = beta0_draws_orig,
      beta1_draws = beta1_draws_orig
    ),
    predictor_names = pred_names
  )
  class(result) <- "gp_fit"
  result
}


#' Generate GP Data from Simulation Model
#'
#' @param delta True effect size parameter
#' @param lambda GP lengthscale parameter
#' @param sigma Noise standard deviation
#' @param N Sample size
#' @param x Covariate values (vector of length N)
#' @param n_samples Number of datasets to generate
#' @param seed Random seed
#' @param model_path Directory containing Stan files
#' @return Matrix of simulated y values (n_samples x N)
generate_gp_data <- function(
  delta, lambda, sigma, N, x,
  n_samples = 1L,
  seed = 30L,
  model_path = "models/"
) {
  simu_mod <- get_model("simulation", model_path = model_path)

  simu_data <- list(
    delta = delta,
    lambda = lambda,
    sigma = sigma,
    N = N,
    x = as.numeric(x)
  )

  simu_fit <- rstan::sampling(
    simu_mod,
    data = simu_data,
    warmup = 0,
    iter = n_samples,
    chains = 1,
    algorithm = "Fixed_param",
    seed = seed,
    refresh = 0
  )

  rstan::extract(simu_fit)$y
}


#' Simulate GP draws in R with optional orthogonalization
#'
#' This mirrors the model definition in fit_full_model.stan when
#' orthogonalize = TRUE: the latent GP draw is projected onto the
#' orthogonal complement of span{1, x_c}, so the resulting nonlinear
#' component contains no intercept or linear trend.
#'
#' @param x Covariate values
#' @param delta Nonlinear magnitude parameter
#' @param lambda GP lengthscale parameter
#' @param sigma Residual scale parameter used in alpha = delta * sigma
#' @param n_draws Number of latent function draws
#' @param orthogonalize If TRUE, project draws to the nonlinear subspace
#' @param beta0 Linear intercept added to the returned mean function
#' @param beta1 Linear slope(s); scalar for univariate x, vector otherwise
#' @param add_noise If TRUE, add Gaussian noise with sd = sigma to y
#' @param seed Optional random seed
#' @param nugget Small diagonal jitter for numerical stability
#' @return List with x, mu, f_raw, f, and y
simulate_gp_draws <- function(
  x,
  delta,
  lambda,
  sigma,
  n_draws = 1L,
  orthogonalize = FALSE,
  beta0 = 0,
  beta1 = 0,
  add_noise = FALSE,
  seed = NULL,
  nugget = 1e-10
) {
  x_mat <- as.matrix(x)
  if (is.vector(x_mat)) x_mat <- matrix(x_mat, ncol = 1)

  N <- nrow(x_mat)
  P <- ncol(x_mat)

  if (!is.null(seed)) set.seed(seed)

  x_c <- scale(x_mat, center = TRUE, scale = FALSE)
  x_vec <- as.numeric(x_mat[, 1])
  alpha <- delta * sigma

  D2 <- outer(x_vec, x_vec, FUN = function(a, b) (a - b)^2)
  K <- alpha^2 * exp(-0.5 * D2 / lambda^2)
  diag(K) <- diag(K) + nugget

  R <- chol(K)
  z <- matrix(stats::rnorm(n_draws * N), nrow = n_draws, ncol = N)
  f_raw <- z %*% R

  f_nl <- f_raw
  if (orthogonalize) {
    Z_aug <- cbind(1, x_c)
    Q_Z <- qr.Q(qr(Z_aug))
    f_nl <- f_raw - (f_raw %*% Q_Z) %*% t(Q_Z)
  }

  beta1 <- rep(beta1, length.out = P)
  mu_vec <- as.numeric(beta0 + x_c %*% beta1)
  mu <- matrix(mu_vec, nrow = n_draws, ncol = N, byrow = TRUE)

  y <- mu + f_nl
  if (add_noise) {
    y <- y + matrix(stats::rnorm(n_draws * N, sd = sigma), nrow = n_draws, ncol = N)
  }

  list(
    x = x_mat,
    mu = mu,
    f_raw = f_raw,
    f = f_nl,
    y = y
  )
}


#' Generate orthogonalized GP data
#'
#' Convenience wrapper for the univariate case used in prior-elicitation
#' plots: beta0 = 0 and beta1 = 0, with the nonlinear component projected
#' onto the orthogonal complement of span{1, x_c}.
#'
#' @param delta Nonlinear magnitude parameter
#' @param lambda GP lengthscale parameter
#' @param sigma Noise standard deviation
#' @param N Sample size
#' @param x Covariate values (vector of length N)
#' @param n_samples Number of datasets to generate
#' @param seed Random seed
#' @param add_noise If TRUE, add Gaussian noise with sd = sigma
#' @param return_latent If TRUE, return latent components as well
#' @return Matrix of y values, or a full list when return_latent = TRUE
generate_orthogonalized_gp_data <- function(
  delta, lambda, sigma, N, x,
  n_samples = 1L,
  seed = 30L,
  add_noise = TRUE,
  return_latent = FALSE
) {
  sim <- simulate_gp_draws(
    x = x,
    delta = delta,
    lambda = lambda,
    sigma = sigma,
    n_draws = n_samples,
    orthogonalize = TRUE,
    beta0 = 0,
    beta1 = 0,
    add_noise = add_noise,
    seed = seed
  )

  if (return_latent) return(sim)
  sim$y
}


#' Simulate from the projected-and-standardized GP parameterization
#'
#' This matches fit_full_model_projected_scaled.stan. The unit-amplitude
#' kernel is first projected to the orthogonal complement of span{1, x_c},
#' then rescaled so that trace(K_tilde_std) / N = 1. Under this
#' parameterization, delta controls the projected RMS on the latent scale.
#'
#' @param x Covariate values
#' @param delta Projected nonlinear magnitude parameter
#' @param lambda GP lengthscale parameter
#' @param sigma Residual scale parameter
#' @param n_draws Number of latent function draws
#' @param beta0 Linear intercept
#' @param beta1 Linear slope(s)
#' @param add_noise If TRUE, add Gaussian noise with sd = sigma to y
#' @param seed Optional random seed
#' @param nugget Small diagonal jitter
#' @return List with x, mu, f, y, and scale_factor
simulate_projected_scaled_gp_draws <- function(
  x,
  delta,
  lambda,
  sigma,
  n_draws = 1L,
  beta0 = 0,
  beta1 = 0,
  add_noise = FALSE,
  seed = NULL,
  nugget = 1e-10
) {
  x_mat <- as.matrix(x)
  if (is.vector(x_mat)) x_mat <- matrix(x_mat, ncol = 1)

  N <- nrow(x_mat)
  P <- ncol(x_mat)

  if (!is.null(seed)) set.seed(seed)

  x_c <- scale(x_mat, center = TRUE, scale = FALSE)
  x_vec <- as.numeric(x_mat[, 1])
  Z_aug <- cbind(1, x_c)
  Q_Z <- qr.Q(qr(Z_aug))
  P_perp <- diag(N) - Q_Z %*% t(Q_Z)

  D2 <- outer(x_vec, x_vec, FUN = function(a, b) (a - b)^2)
  R_lambda <- exp(-0.5 * D2 / lambda^2)
  diag(R_lambda) <- diag(R_lambda) + nugget

  K_proj <- P_perp %*% R_lambda %*% P_perp
  scale_factor <- sum(diag(K_proj)) / N
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    stop("Projected kernel scale factor is non-positive.")
  }

  K_std <- K_proj / scale_factor
  K_latent <- (sigma * delta)^2 * K_std
  K_latent <- 0.5 * (K_latent + t(K_latent))

  R <- chol(K_latent + diag(1e-12, N))
  z <- matrix(stats::rnorm(n_draws * N), nrow = n_draws, ncol = N)
  f <- z %*% R

  beta1 <- rep(beta1, length.out = P)
  mu_vec <- as.numeric(beta0 + x_c %*% beta1)
  mu <- matrix(mu_vec, nrow = n_draws, ncol = N, byrow = TRUE)

  y <- mu + f
  if (add_noise) {
    y <- y + matrix(stats::rnorm(n_draws * N, sd = sigma), nrow = n_draws, ncol = N)
  }

  list(
    x = x_mat,
    mu = mu,
    f = f,
    y = y,
    scale_factor = scale_factor
  )
}

# ============================================================
# fit_model_select: GP regression with per-predictor selection
#
# Allows excluding predictors from the linear part, the nonlinear
# (GP) part, or both. Returns a gp_fit-compatible object so that
# downstream utilities (bridge sampling, summaries, etc.) work
# without modification.
#
# Source this file AFTER src/_setup.R.
# ============================================================

#' Fit GP Regression Model with Predictor Selection
#'
#' @param X                 Predictor matrix (N x P), original scale
#' @param y                 Response vector (length N)
#' @param include_linear    Integer/logical vector of length P (1 = include)
#' @param include_nonlinear Integer/logical vector of length P (1 = include)
#' @param d_order           Moment-prior order on delta (0-3)
#' @param chains,iter,warmup,adapt_delta,max_treedepth,seed,refresh  MCMC controls
#' @param stan_file         Path to the selection Stan file
#' @param ...               Passed to rstan::sampling()
#' @return gp_fit object (same structure as fit_model())
fit_model_select <- function(
  X,
  y,
  include_linear,
  include_nonlinear,
  d_order = 0L,
  chains = 2L,
  iter = 2000L,
  warmup = 500L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  seed = 42L,
  refresh = 500L,
  stan_file = "models/fit_full_model_select.stan",
  ...
) {
  d_order <- as.integer(d_order)
  stopifnot(d_order %in% 0:3)

  # --- Prepare X ---
  X_orig <- as.matrix(X)
  if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)
  y <- as.numeric(y)
  stopifnot(nrow(X_orig) == length(y))

  N <- nrow(X_orig)
  P <- ncol(X_orig)

  pred_names <- colnames(X_orig)
  if (is.null(pred_names)) pred_names <- paste0("X", seq_len(P))

  # --- Validate masks ---
  include_linear    <- as.integer(as.logical(include_linear))
  include_nonlinear <- as.integer(as.logical(include_nonlinear))
  stopifnot(length(include_linear)    == P)
  stopifnot(length(include_nonlinear) == P)

  # Rescale X to [-0.5, 0.5]
  X_scaled <- apply(X_orig, 2, rescale_to_unit)
  if (is.vector(X_scaled)) X_scaled <- matrix(X_scaled, ncol = 1)
  colnames(X_scaled) <- pred_names

  # --- Build Stan data ---
  stan_data <- list(
    N = N, P = P, X = X_scaled, y = y,
    include_linear    = include_linear,
    include_nonlinear = include_nonlinear
  )
  stan_data <- c(stan_data, moment_prior_stan_data(d_order))

  # --- Compile (cached per stan_file path) ---
  if (!exists(".compiled_models_select", envir = globalenv(), inherits = FALSE)) {
    .compiled_models_select <<- new.env(parent = emptyenv())
  }
  key <- normalizePath(stan_file, mustWork = TRUE)
  key_mtime <- file.info(key)$mtime
  entry <- .compiled_models_select[[key]]
  if (!is.list(entry) ||
      is.null(entry$model) ||
      !identical(entry$mtime, key_mtime)) {
    message("Compiling Stan model (select): ", stan_file)
    .compiled_models_select[[key]] <- list(
      model = rstan::stan_model(file = key),
      mtime = key_mtime
    )
  }
  stan_mod <- .compiled_models_select[[key]]$model

  fit <- rstan::sampling(
    object = stan_mod,
    data = stan_data,
    iter = iter,
    warmup = warmup,
    chains = chains,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    seed = seed,
    refresh = refresh,
    verbose = FALSE,
    open_progress = FALSE,
    ...
  )

  # --- Back-transform linear coefficients (full-length beta1 from generated quantities) ---
  post <- rstan::extract(fit)
  beta1_draws <- post$beta1   # draws x P, zeros for excluded-linear predictors

  y_bar <- mean(y)
  x_bar <- colMeans(X_orig)
  x_min <- apply(X_orig, 2, min)
  x_max <- apply(X_orig, 2, max)
  rng <- x_max - x_min

  beta1_draws_orig <- sweep(beta1_draws, 2, rng, "/")
  # Ensure zeros stay zeros (division by rng for excluded cols is harmless here
  # since those draws are identically 0).
  beta0_draws_orig <- as.vector(y_bar - beta1_draws_orig %*% x_bar)

  result <- list(
    fit = fit,
    model_type = "full_select",
    d_order = d_order,
    stan_data = stan_data,
    orig = list(X = X_orig, y = y),
    coef_orig = list(
      beta0_draws = beta0_draws_orig,
      beta1_draws = beta1_draws_orig
    ),
    predictor_names = pred_names,
    include_linear    = include_linear,
    include_nonlinear = include_nonlinear
  )
  class(result) <- c("gp_fit_select", "gp_fit")
  result
}


# ============================================================
# Helpers for downstream utilities
# ============================================================

#' Check whether a gp_fit object has GP terms available
#'
#' @param result A gp_fit-like object
#' @return TRUE for full GP fits, including selection fits
is_gp_full_fit <- function(result) {
  inherits(result, "gp_fit") &&
    isTRUE(result$model_type %in% c("full", "full_select"))
}


#' Get the active linear/nonlinear inclusion masks for a fit
#'
#' @param result A gp_fit object
#' @return Named list with integer vectors `linear` and `nonlinear`
get_gp_fit_masks <- function(result) {
  stopifnot(inherits(result, "gp_fit"))

  P <- result$stan_data$P

  linear_mask <- result$include_linear
  nonlinear_mask <- result$include_nonlinear

  if (is.null(linear_mask)) linear_mask <- rep(1L, P)
  if (is.null(nonlinear_mask)) nonlinear_mask <- rep(1L, P)

  list(
    linear = as.integer(linear_mask),
    nonlinear = as.integer(nonlinear_mask)
  )
}


#' Coerce posterior draws to a draws x P matrix
#'
#' @param x Posterior draws, matrix or vector
#' @param P Number of predictors
#' @return Numeric matrix with `P` columns
as_gp_draw_matrix <- function(x, P) {
  if (is.null(dim(x))) {
    matrix(as.numeric(x), ncol = P)
  } else {
    as.matrix(x)
  }
}


#' Extract full-length posterior draws for a model parameter
#'
#' For `full_select` fits, reduced parameters are expanded to length `P`
#' using generated quantities so downstream summaries can index by the
#' original predictor positions.
#'
#' @param result A gp_fit object
#' @param post Optional posterior list from `rstan::extract(result$fit)`
#' @param param One of "theta", "beta1", "delta", or "lambda"
#' @return Numeric matrix with one column per original predictor
extract_gp_draws_full <- function(result, post = NULL, param = c("theta", "beta1", "delta", "lambda")) {
  stopifnot(inherits(result, "gp_fit"))
  param <- match.arg(param)

  if (is.null(post)) post <- rstan::extract(result$fit)
  P <- result$stan_data$P

  if (identical(result$model_type, "full_select")) {
    if (param == "theta") {
      beta1_full <- as_gp_draw_matrix(post$beta1, P)
      return(sweep(beta1_full, 1, post$sigma, "/"))
    }
    if (param == "beta1") {
      return(as_gp_draw_matrix(post$beta1, P))
    }
    if (param == "delta") {
      return(as_gp_draw_matrix(post$delta_full, P))
    }
    if (param == "lambda") {
      return(as_gp_draw_matrix(post$lambda_full, P))
    }
  }

  if (param == "beta1" && !is.null(post$beta1)) {
    return(as_gp_draw_matrix(post$beta1, P))
  }

  as_gp_draw_matrix(post[[param]], P)
}
