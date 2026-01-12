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
#' @param data Data list
#' @param chains Number of chains (default: 2)
#' @param iter Total iterations per chain (default: 1000)
#' @param warmup Warmup/burn-in iterations (default: iter/2)
#' @param seed Random seed for reproducibility
#' @return List containing model_name, data, and fit object
fit_model <- function(model, data, chains = 2, iter = 1000, 
                      warmup = iter/2, seed = 42) {
  
  # Capture parameter name as model_name
  model_name <- deparse(substitute(model))
  
  fit <- rstan::sampling(
    model, 
    data = data, 
    iter = iter, 
    warmup = warmup,
    chains = chains, 
    seed = seed,
    refresh = 1000,
    verbose = FALSE, 
    open_progress = FALSE
  )

  list(model_name = model_name, data = data, fit = fit)
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
  stan_simu <- compile_model("simu")
  
  simu_data <- list(
    delta = delta, 
    lambda = lambda, 
    sigma = sigma,
    N = N, 
    x = as.numeric(x)
  )
  
  simu_fit <- rstan::sampling(
    stan_simu, 
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