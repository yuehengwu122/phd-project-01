# ============================================================
# Application 1: Neuroscience
# Facebook Friends vs Gray Matter Density
#
# Dataset: Lock5Data::FacebookFriends
# - Predictor: sqrt(FBfriends), 1-dimensional
# - Response:  GMdensity (gray matter density)
# ============================================================

source("src/_setup.R")

# --- Data ---
library(Lock5Data)
data(FacebookFriends)

X <- matrix(sqrt(FacebookFriends$FBfriends), ncol = 1)
colnames(X) <- "sqrt_FBfriends"
y <- FacebookFriends$GMdensity

cat("N =", nrow(X), ", P =", ncol(X), "\n")

# --- Fit: full GP model (local prior, d_order = 0) ---
fit_d0 <- fit_model(X, y, model = "full", d_order = 0,
                     chains = 2, iter = 5500, warmup = 500)

# --- Fit: linear-only model (for bridge sampling BF) ---
fit_lin <- fit_model(X, y, model = "linear",
                      chains = 2, iter = 5500, warmup = 500)

# --- Diagnostics ---
print(fit_d0)
summary_gp_fit(fit_d0)

# --- Bayes factors ---
# Interval-null BF table (point null at delta = 0)
compute_bf_table(fit_d0)

# Bridge sampling BF (full vs linear)
cat("Bridge sampling log BF10:", compute_bf_bridge(fit_d0, fit_lin, log = TRUE), "\n")

# --- Plots ---
plot_posterior_delta(fit_d0, bins = 50, plot_prior = TRUE)
plot_posterior_theta(fit_d0, bins = 50, plot_prior = TRUE)
plot_gp_trends(fit_d0)

# --- Save ---
saveRDS(fit_d0,  "results/fit_01_neuro_d0.rds")
saveRDS(fit_lin, "results/fit_01_neuro_linear.rds")
