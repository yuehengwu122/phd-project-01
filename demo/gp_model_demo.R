# =============================================================================
# GP Model Demo Script
# =============================================================================
# This script demonstrates Gaussian Process models with two dataset options
# and four different model variants for comparison.
#
# Models:
# - Exact GP vs HSGP approximation
# - Independent vs multivariate priors on coefficients
#
# Usage: Simply run this script section by section or all at once
# =============================================================================

# Load required functions and libraries
source("_common.R")

# =============================================================================
# DATA SELECTION
# =============================================================================
# Choose which dataset to analyze by setting this variable:
USE_SIMULATED_DATA <- TRUE  # Set to FALSE for diabetes data

if (USE_SIMULATED_DATA) {
  # Option 1: 2D Simulated Data (fast, good for quick testing)
  cat("Using 2D simulated data...\n")
  
  set.seed(2026)
  N_sim <- 100
  P_sim <- 2
  
  # Generate predictors
  X_sim <- matrix(runif(N_sim * P_sim, -1, 1), N_sim, P_sim)
  X_sim_scaled <- apply(X_sim, 2, rescale_to_range1)
  
  # Generate response with nonlinear functions
  f1 <- function(x) 0.5 * sin(3 * pi * x)
  f2 <- function(x) 0.3 * (2 * x^2 - 1)
  
  y_sim <- 0.2 + 0.1 * X_sim_scaled[,1] + 0.05 * X_sim_scaled[,2] + 
           f1(X_sim_scaled[,1]) + f2(X_sim_scaled[,2]) + 
           rnorm(N_sim, 0, 0.15)
  
  # Prepare data for Stan
  fit_data <- list(
    N = N_sim,
    y = y_sim,
    P = P_sim, 
    X = X_sim_scaled,
    M = 25  # Number of basis functions for HSGP
  )
  
} else {
  # Option 2: Diabetes Data (more complex, real application)
  cat("Using diabetes data...\n")
  
  library("care")
  data(efron2004)
  
  x_mat <- as.matrix(efron2004$x)   
  X <- x_mat[,-2]  # Remove "sex" (binary variable)
  y <- efron2004$y[,1]
  
  # Prepare data for Stan
  fit_data <- list(
    N = nrow(X),
    y = y,
    P = ncol(X), 
    X = apply(X, 2, rescale_to_range1),
    M = 40  # More basis functions for higher dimensional data
  )
}

cat("Data prepared: N =", fit_data$N, ", P =", fit_data$P, "\n\n")

# =============================================================================
# MODEL COMPILATION
# =============================================================================
cat("Compiling Stan models...\n")

# Four model variants to compare
exactgp_indep <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_13.stan")
exactgp_multi <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_15.stan")
hsgp_indep <- rstan::stan_model(file = "stan/multiple_gp/fit_add_hsgp_02.stan")
hsgp_multi <- rstan::stan_model(file = "stan/multiple_gp/fit_add_hsgp_01.stan")

cat("Models compiled successfully!\n\n")

# =============================================================================
# MODEL FITTING
# =============================================================================
# Fit all four models
# Note: Using fewer iterations for demo purposes (faster runtime)

cat("Fitting models...\n")

# Caution: If using diabetes dataset (USE_SIMULATED_DATA = FALSE), 
# fitting may take over 10 hours.

# cat("1. Exact GP with independent priors...\n")
# fit_exactgp_indep <- fit_model(
#   exactgp_indep,
#   data = fit_data,
#   chains = 2,
#   iter = 2000,
#   warmup = 500
# )

# cat("2. Exact GP with multivariate priors...\n")
# fit_exactgp_multi <- fit_model(
#   exactgp_multi,
#   data = fit_data,
#   chains = 2,
#   iter = 2000,
#   warmup = 500
# )

cat("3. HSGP with independent priors...\n")
fit_hsgp_indep <- fit_model(
  hsgp_indep,
  data = fit_data,
  chains = 2,
  iter = 2000,
  warmup = 500
)

cat("4. HSGP with multivariate priors...\n")
fit_hsgp_multi <- fit_model(
  hsgp_multi,
  data = fit_data,
  chains = 2,
  iter = 2000,
  warmup = 500
)

cat("All models fitted!\n\n")

# =============================================================================
# RESULTS & DIAGNOSTICS
# =============================================================================

cat("=== MODEL SUMMARIES ===\n")
print_quick_summary(fit_exactgp_indep)
print_quick_summary(fit_exactgp_multi)
print_quick_summary(fit_hsgp_indep)
print_quick_summary(fit_hsgp_multi)

cat("\n=== BAYES FACTORS (Savage-Dickey ratios) ===\n")
cat("Exact GP - Independent priors:\n")
compute_savage_dickey_ratios(fit_exactgp_indep)

cat("\nExact GP - Multivariate priors:\n")
compute_savage_dickey_ratios(fit_exactgp_multi)

cat("\nHSGP - Independent priors:\n")
compute_savage_dickey_ratios(fit_hsgp_indep)

cat("\nHSGP - Multivariate priors:\n")
compute_savage_dickey_ratios(fit_hsgp_multi)

# =============================================================================
# VISUALIZATION
# =============================================================================

cat("\n=== CREATING PLOTS ===\n")

cat("Plotting posterior distributions...\n")
plot_posterior_delta(fit_exactgp_indep, model_name = "Exact GP / Independent")
plot_posterior_delta(fit_exactgp_multi, model_name = "Exact GP / Multivariate") 
plot_posterior_delta(fit_hsgp_indep, model_name = "HSGP / Independent")
plot_posterior_delta(fit_hsgp_multi, model_name = "HSGP / Multivariate")

cat("Plotting GP trends...\n")
# GP trends (exact GP models)
plot_gp_trends(fit_exactgp_indep, model_name = "Exact GP / Independent")
plot_gp_trends(fit_exactgp_multi, model_name = "Exact GP / Multivariate")

# HSGP trends  
plot_hsgp_trends(fit_hsgp_indep, model_name = "HSGP / Independent")
plot_hsgp_trends(fit_hsgp_multi, model_name = "HSGP / Multivariate")

cat("Demo completed! Check the plots to compare model performance.\n")
