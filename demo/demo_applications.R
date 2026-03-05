# =============================================================================
# Demo: Fitted Model Objects for Multiple Application Cases
# =============================================================================
# This script demonstrates how to load and inspect fitted model objects
# from two application cases: A3 (Diabetes) and A4 (Boston Housing).
# =============================================================================

source("_common.R")

# -----------------------------------------------------------------------------
# Application A3: Diabetes Dataset
# -----------------------------------------------------------------------------

fit_a3_dep <- readRDS("analysis/application/fits/fit_a3_dep.rds")
fit_a3_dep$fit # View Stan fit object

# Convergence diagnostics and posterior summary
print_quick_summary(
  fit_a3_dep,
  key_params = c("delta", "lambda", "sigma", "theta")
)

# Posterior distributions with prior overlay
plot_posterior_delta(fit_a3_dep, bins = 50, plot_prior = TRUE, ylim = c(0, 0.8))
plot_posterior_theta(fit_a3_dep, plot_prior = TRUE, bins = 50, ylim = c(0, 1.6))

# Savage-Dickey ratios for variable selection
compute_savage_dickey_ratios(fit_a3_dep)


# -----------------------------------------------------------------------------
# Application A4: Boston Housing Dataset
# -----------------------------------------------------------------------------

fit_a4_dep <- readRDS("analysis/application/fits/fit_a4_dep.rds")
fit_a4_dep$fit # View Stan fit object

# Convergence diagnostics and posterior summary
print_quick_summary(
  fit_a4_dep,
  key_params = c("delta", "lambda", "sigma", "theta")
)

# Posterior distributions with prior overlay
plot_posterior_delta(fit_a4_dep, plot_prior = TRUE, bins = 50, ylim = c(0, 0.8))
plot_posterior_theta(fit_a4_dep, plot_prior = TRUE, bins = 50, ylim = c(0, 1.5))

# Savage-Dickey ratios for variable selection
compute_savage_dickey_ratios(fit_a4_dep)
