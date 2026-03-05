# HSGP Experiment: Testing Hilbert Space GP Approximation
# Goal: Compare HSGP vs exact GP on speed and accuracy
# Date: 2026-01-04

# Setup -------------------------------------------------------------------
source("_common.R")

set.seed(2026)


# Part 1: Simulated 2D Dataset --------------------------------------------
cat("=== Part 1: Testing HSGP on Simulated 2D Data ===\n\n")

# Generate simulated data with nonlinear effects
N_sim <- 100
P_sim <- 2

# Generate predictors
X_sim <- matrix(runif(N_sim * P_sim, -2, 2), N_sim, P_sim)

# True functions: one linear, one nonlinear
f1 <- function(x) 0.5 * x # linear
f2 <- function(x) sin(2 * x) # nonlinear

# Generate response
y_sim <- f1(X_sim[, 1]) + f2(X_sim[, 2]) + rnorm(N_sim, 0, 0.3)

# Plot true functions
par(mfrow = c(1, 2))
plot(X_sim[, 1], y_sim, pch = 16, xlab = "X1", ylab = "y")
lines(X_sim[, 1], f1(X_sim[, 1]), col = "blue", lwd = 2)
ord_X2 <- order(X_sim[, 2])
plot(X_sim[ord_X2, 2], y_sim[ord_X2], pch = 16, xlab = "X2", ylab = "y")
lines(X_sim[ord_X2, 2], f2(X_sim[ord_X2, 2]), col = "blue", lwd = 2)


# Rescale to [-0.5, 0.5] as expected by models
X_sim_scaled <- apply(X_sim, 2, rescale_to_range1)

# Compile models
cat("Compiling Stan models...\n")
gp_exact_01 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_13.stan")
gp_exact_02 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_p1_15.stan")
gp_hsgp_01 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_hsgp_01.stan")
gp_hsgp_02 <- rstan::stan_model(file = "stan/multiple_gp/fit_add_hsgp_02.stan")

# Prepare data

data_sim <- list(
  N = N_sim,
  P = P_sim,
  X = X_sim_scaled,
  y = y_sim,
  M = 25, # number of basis functions
  N_grid = 50 # grid size for predictions
)

# Fit exact GP
cat("\n--- Fitting Exact GP (Simulated Data) ---\n")
time_exact_sim <- system.time({
  fit_exact_01 <- fit_model(
    gp_exact_01,
    data = data_sim,
    chains = 2,
    iter = 4500,
    warmup = 500
  )
})

time_exact_sim_02 <- system.time({
  fit_exact_02 <- fit_model(
    gp_exact_02,
    data = data_sim,
    chains = 2,
    iter = 4500,
    warmup = 500
  )
})

cat("Exact GP time (simulated):", round(time_exact_sim[3], 2), "seconds\n")

# Fit HSGP
cat("\n--- Fitting HSGP (Simulated Data) ---\n")
time_hsgp_sim_01 <- system.time({
  fit_hsgp_01 <- fit_model(
    gp_hsgp_01,
    data = data_sim,
    chains = 2,
    iter = 4500,
    warmup = 500
  )
})

time_hsgp_sim_02 <- system.time({
  fit_hsgp_02 <- fit_model(
    gp_hsgp_02,
    data = data_sim,
    chains = 2,
    iter = 4500,
    warmup = 500
  )
})

print_quick_summary(fit_exact_01)
print_quick_summary(fit_exact_02)
print_quick_summary(fit_hsgp_01)
print_quick_summary(fit_hsgp_02)

# plot posterior of delta
plot_posterior(fit_exact_01, model_name = "exact-GP / marginal-θ")
plot_posterior(fit_exact_02, model_name = "exact-GP / multivariate-θ")
plot_posterior(fit_hsgp_01, model_name = "HSGP / multivariate-θ")
plot_posterior(fit_hsgp_02, model_name = "HSGP / marginal-θ")

# Plot fitted trends
plot_gp_trends(fit_exact_01, model_name = "exact-GP / marginal-θ")
plot_gp_trends(fit_exact_02, model_name = "exact-GP / multivariate-θ")
plot_hsgp_trends(fit_hsgp_01, model_name = "HSGP / multivariate-θ")
plot_hsgp_trends(fit_hsgp_02, model_name = "HSGP / marginal-θ")

# Compute Savage-Dickey ratios
compute_savage_dickey_ratios(fit_exact_01)
compute_savage_dickey_ratios(fit_exact_02)
compute_savage_dickey_ratios(fit_hsgp_01)
compute_savage_dickey_ratios(fit_hsgp_02)

bridgesampling::bridge_sampler(fit_exact_01$fit)
bridgesampling::bridge_sampler(fit_exact_02$fit)
bridgesampling::bridge_sampler(fit_hsgp_01$fit)
bridgesampling::bridge_sampler(fit_hsgp_02$fit)
