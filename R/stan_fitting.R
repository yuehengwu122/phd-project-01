#' Compile Stan Model
#' @param model_name Model name (e.g., "single_gp", "null")
#' @param model_path Path to Stan files directory (default: "stan/")
#' @return Compiled Stan model object
compile_model <- function(model_name, model_path = "stan/") {
  # Model file mapping table
  model_files <- list(
    null = "single_gp/null_4.stan",
    single_gp_p1 = "single_gp/fit_single_07.stan",
    single_gp_p2 = "single_gp/fit_single_07_p2.stan",
    simu = "simu.stan"
  )
  
  if (!model_name %in% names(model_files)) {
    stop(paste("Unknown model:", model_name))
  }
  
  file_path <- file.path(model_path, model_files[[model_name]])
  rstan::stan_model(file = file_path)
}

#' Fit Stan Model with Standard Settings (rstan)
#' @param model Compiled rstan model object
#' @param X Predictor matrix (N x P)
#' @param y Response vector (length N)
#' @param delta_prior Prior specification for delta ("p1" or "p2"). Default is "p1".
#' @param M Number of basis functions for HSGP approximation (default: 40)
#' @param chains Number of chains (default: 2)
#' @param iter Total iterations per chain (default: 1000)
#' @param warmup Warmup/burn-in iterations (default: iter/2)
#' @param seed Random seed for reproducibility
#' @return List containing model_name, data, orig, coef_orig, delta_prior, and fit object
fit_model <- function(model, X, y, delta_prior = "p1", M = 40,
                      chains = 2, iter = 1000, warmup = iter/2, seed = 42) {
  
  # Capture parameter name as model_name
  model_name <- deparse(substitute(model))
  
  # Store original X before rescaling
  X_orig <- as.matrix(X)
  if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)
  
  # Prepare data for Stan (rescale X to range [0,1])
  fit_data <- list(
    N = nrow(X_orig),
    y = as.numeric(y),
    P = ncol(X_orig),
    X = apply(X_orig, 2, rescale_to_range1),
    M = M
  )
  
  fit <- rstan::sampling(
    model, 
    data = fit_data, 
    iter = iter, 
    warmup = warmup,
    chains = chains, 
    seed = seed,
    refresh = 1000,
    verbose = FALSE, 
    open_progress = FALSE
  )
  
  # Compute coefficients on original scale
  post <- rstan::extract(fit)
  beta1_draws <- post$beta1
  
  # Data summaries on original scale
  y_bar <- mean(y)
  x_bar <- colMeans(X_orig)
  x_min <- apply(X_orig, 2, min)
  x_max <- apply(X_orig, 2, max)
  rng <- x_max - x_min
  
  # Back-transform slopes: beta1_orig[j] = beta1_c[j] / (x_max[j] - x_min[j])
  if (is.matrix(beta1_draws)) {
    beta1_draws_orig <- sweep(beta1_draws, 2, rng, "/")  # S x P
  } else {
    beta1_draws_orig <- matrix(beta1_draws / rng[1], ncol = 1)
  }
  
  # Back-transform intercept per draw: beta0 = y_bar - x_bar' * beta1_orig
  if (is.matrix(beta1_draws_orig)) {
    beta0_draws_orig <- as.vector(y_bar - beta1_draws_orig %*% x_bar)  # length S
  } else {
    beta0_draws_orig <- y_bar - beta1_draws_orig * x_bar
  }
  
  list(
    model_name = model_name,
    data = fit_data,
    orig = list(X = X_orig, y = y),
    coef_orig = list(beta0_draws = beta0_draws_orig, beta1_draws = beta1_draws_orig),
    delta_prior = delta_prior,
    fit = fit
  )
}

#' Generate GP Data from Simulation Model
#' @param delta True effect size parameter
#' @param lambda GP lengthscale parameter
#' @param sigma Noise standard deviation
#' @param N Sample size
#' @param x Covariate values
#' @param n_samples Number of datasets to generate
#' @return Matrix of simulated y values (n_samples × N)
generate_gp_data <- function(delta, lambda, sigma, N, x, n_samples = 1) {  
  simu_data <- list(
    delta = delta, 
    lambda = lambda, 
    sigma = sigma,
    N = N, 
    x = as.numeric(x)
  )
  
  simu_fit <- rstan::sampling(
    gp_simu, 
    data = simu_data,
    warmup = 0, 
    iter = n_samples,
    chains = 1, 
    algorithm = "Fixed_param", 
    seed = 30,
    refresh = 0
  )
  
  rstan::extract(simu_fit)$y
}


#' Reconstruct Result Object from Existing Stan Fit
#' 
#' Creates a result object compatible with summary_gp_fit() and plot_gp_trends()
#' from an existing Stan fit object without needing to refit.
#'
#' @param fit Stan fit object (rstan stanfit)
#' @param X_orig Original (unscaled) predictor matrix
#' @param y Response vector
#' @param model_name Model name string (optional)
#' @param delta_prior Prior specification ("p1" or "p2"). Default is "p1".
#' @param M Number of basis functions (default: 40)
#' @return List with same structure as fit_model() output
reconstruct_result <- function(fit, X_orig, y, model_name = NULL, 
                                delta_prior = "p1", M = 40) {
  
  X_orig <- as.matrix(X_orig)
  if (is.vector(X_orig)) X_orig <- matrix(X_orig, ncol = 1)
  y <- as.numeric(y)
  
  # Rescale X to range [0,1] as done during fitting
  X_scaled <- apply(X_orig, 2, rescale_to_range1)
  
  # Build data list matching fit_model structure
  fit_data <- list(
    N = nrow(X_orig),
    y = y,
    P = ncol(X_orig),
    X = X_scaled,
    M = M
  )
  
  # Compute coefficients on original scale
  post <- rstan::extract(fit)
  beta1_draws <- post$beta1
  
  y_bar <- mean(y)
  x_bar <- colMeans(X_orig)
  x_min <- apply(X_orig, 2, min)
  x_max <- apply(X_orig, 2, max)
  rng <- x_max - x_min
  
  if (is.matrix(beta1_draws)) {
    beta1_draws_orig <- sweep(beta1_draws, 2, rng, "/")
  } else {
    beta1_draws_orig <- matrix(beta1_draws / rng[1], ncol = 1)
  }
  
  if (is.matrix(beta1_draws_orig)) {
    beta0_draws_orig <- as.vector(y_bar - beta1_draws_orig %*% x_bar)
  } else {
    beta0_draws_orig <- y_bar - beta1_draws_orig * x_bar
  }
  
  list(
    model_name = model_name,
    data = fit_data,
    orig = list(X = X_orig, y = y),
    coef_orig = list(beta0_draws = beta0_draws_orig, beta1_draws = beta1_draws_orig),
    delta_prior = delta_prior,
    fit = fit
  )
}