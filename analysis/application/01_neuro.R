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
fit_d0 <- fit_model(
  X,
  y,
  model = "full",
  d_order = 0,
  chains = 2,
  iter = 5500,
  warmup = 500
)

# --- Fit: linear-only model (for bridge sampling BF) ---
fit_lin <- fit_model(
  X,
  y,
  model = "linear",
  chains = 2,
  iter = 5500,
  warmup = 500
)

# --- Diagnostics ---
print(fit_d0)
summary_gp_fit(fit_d0)

# --- Bayes factors ---
# Interval-null BF table (point null at delta = 0)
compute_bf_table(fit_d0)

# Bridge sampling BF (full vs linear)
cat(
  "Bridge sampling log BF10:",
  compute_bf_bridge(fit_d0, fit_lin, log = TRUE),
  "\n"
)

# --- Plots ---
plot_posterior_delta(fit_d0, bins = 50, plot_prior = TRUE)
plot_posterior_theta(fit_d0, bins = 50, plot_prior = TRUE)

plot_neuro_trend <- plot_gp_trends(
  fit_d0,
  x_label = "Square root of number of Facebook friends",
  y_label = "Grey matter density",
  predictor_labels = c(sqrt_FBfriends = "Facebook friends"),
  show_legend = TRUE
)
plot_neuro_trend <- plot_neuro_trend +
  ggplot2::theme(
    panel.border = ggplot2::element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = ggplot2::element_blank()
  )
plot_neuro_trend

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}
ggplot2::ggsave(
  "results/fig_01_neuro_trend.png",
  plot_neuro_trend,
  width = 6.2,
  height = 4.6,
  dpi = 320
)
ggplot2::ggsave(
  "results/fig_01_neuro_trend.pdf",
  plot_neuro_trend,
  width = 6.2,
  height = 4.6
)

# --- Save ---
saveRDS(fit_d0, "results/fit_01_neuro_d0.rds")
saveRDS(fit_lin, "results/fit_01_neuro_linear.rds")
