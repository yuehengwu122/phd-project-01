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

fit_ig <- fit_model(X, y, model = "full_ig", d_order = 0,
                    ig_a = 2.0,     # ig_b defaults to 0.5 * N
                    chains = 2, iter = 5500, warmup = 500)
saveRDS(fit_ig, "results/fit_04_boston_ig.rds")


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







# --- Fit: GAM (mgcv) ---
predictors <- colnames(X)
data <- df

formula_mgcv <- as.formula(paste0(
  "y ~ ",
  paste(predictors, collapse = " + "),
  " + ",
  paste0("s(", predictors, ", bs='tp', m=1, k=10)", collapse = " + ")
))
gam_mgcv <- mgcv::gam(formula_mgcv, data = data)
summary(gam_mgcv)

# --- Extract and save results for the selected fit ---
selected_fit <- gam_mgcv
summary_gam <- summary(selected_fit)

# Predicted values (with standard errors) on the original data
pred_gam <- mgcv::predict.gam(selected_fit, newdata = data, se.fit = TRUE)
gam_results <- data.frame(
  y = data$y,
  fitted = pred_gam$fit,
  se = pred_gam$se.fit,
  residual = data$y - pred_gam$fit
)

# --- Plot estimated trend for each predictor using GAM predictions ---
trend_list <- lapply(predictors, function(pred) {
  grid <- seq(min(data[[pred]], na.rm = TRUE), max(data[[pred]], na.rm = TRUE), length.out = 120)
  medians <- apply(X, 2, median)
  newdata <- as.data.frame(matrix(nrow = length(grid), ncol = length(predictors)))
  colnames(newdata) <- predictors
  for (j in seq_along(predictors)) {
    newdata[[predictors[j]]] <- medians[j]
  }
  newdata[[pred]] <- grid
  p <- mgcv::predict.gam(selected_fit, newdata = newdata, se.fit = TRUE)
  data.frame(
    predictor = pred,
    value = grid,
    fitted = p$fit,
    se = p$se.fit,
    lower = p$fit - 2 * p$se.fit,
    upper = p$fit + 2 * p$se.fit
  )
})

trend_df <- do.call(rbind, trend_list)

# prepare long-form original data for overlaying scatter points
df_long <- data |>
  tidyr::pivot_longer(cols = tidyselect::all_of(predictors), names_to = "predictor", values_to = "value")

plot_trends <- ggplot2::ggplot(trend_df, ggplot2::aes(x = value, y = fitted)) +
  ggplot2::geom_point(data = df_long, ggplot2::aes(x = value, y = y), size = 0.7, alpha = 0.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.6) +
  ggplot2::geom_line(color = "blue") +
  ggplot2::facet_wrap(~predictor, scales = "free_x") +
  ggplot2::labs(title = "Estimated GAM trends (others at median)", x = "Predictor value", y = "Fitted y")

plot_trends

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
ggplot2::ggsave("results/gam_trends_04_boston.png", plot_trends, width = 11, height = 6)
