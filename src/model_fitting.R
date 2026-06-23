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
