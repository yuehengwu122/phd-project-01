# ============================================================
# Application 4: Boston Housing
#
# Dataset: MASS::Boston
# - Predictors: 11 continuous features (excluding chas, rad)
# - Response:   medv (median home value)
# ============================================================

source("src/_setup.R")

# --- Data ---
library(MASS)
data(Boston)

X <- as.matrix(Boston[, -c(4, 9, 14)])   # drop chas (binary), rad (ordinal), medv (response)
y <- Boston$medv

cat("N =", nrow(X), ", P =", ncol(X), "\n")
cat("Predictors:", paste(colnames(X), collapse = ", "), "\n")

# --- Exploratory: marginal scatter + GAM smooth ---
df <- as.data.frame(cbind(y = y, X))
df |>
  tidyr::pivot_longer(-y, names_to = "predictor", values_to = "value") |>
  ggplot2::ggplot(ggplot2::aes(x = value, y = y)) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::stat_smooth(method = "gam", se = TRUE, col = "blue") +
  ggplot2::facet_wrap(~predictor, scales = "free") +
  ggplot2::ggtitle("Marginal scatter plots (with GAM smooth)")

# --- Fit: local prior (d_order = 0) ---
fit_d0 <- fit_model(X, y, model = "full", d_order = 0,
                     chains = 2, iter = 5500, warmup = 500)

# --- Fit: non-local prior order 2 ---
fit_d2 <- fit_model(X, y, model = "full", d_order = 2,
                     chains = 2, iter = 5500, warmup = 500)

# --- Fit: linear-only model ---
fit_lin <- fit_model(X, y, model = "linear",
                      chains = 2, iter = 5500, warmup = 500)

# --- Diagnostics & Summary ---
summary_gp_fit(fit_d0)
summary_gp_fit(fit_d2)

# --- Bayes factors ---
compute_bf_table(fit_d0)
compute_bf_table(fit_d2)

cat("Bridge sampling log BF10 (d0):", compute_bf_bridge(fit_d0, fit_lin, log = TRUE), "\n")
cat("Bridge sampling log BF10 (d2):", compute_bf_bridge(fit_d2, fit_lin, log = TRUE), "\n")

# --- Plots ---
plot_posterior_delta(fit_d0, bins = 50, plot_prior = TRUE, ylim = c(0, 0.8))
plot_posterior_theta(fit_d0, bins = 50, plot_prior = TRUE, ylim = c(0, 1.5))

plot_gp_trends(fit_d0, title = "GP Trends (local prior, d=0)", ylim = c(0, 50))
plot_gp_trends(fit_d2, title = "GP Trends (moment prior, d=2)", ylim = c(0, 50))

# --- Save ---
saveRDS(fit_d0,  "results/fit_04_boston_d0.rds")
saveRDS(fit_d2,  "results/fit_04_boston_d2.rds")
saveRDS(fit_lin, "results/fit_04_boston_linear.rds")
