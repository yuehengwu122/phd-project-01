source("src/_setup.R")

library(dplyr)
library(ggplot2)
library(readr)

# ============================================================
# Epsilon calibration figure for interval Bayes factors
#
# Purpose:
# Show what very small projected nonlinear effects look like under the
# projected-scaled delta definition, and record the corresponding
# expected projected RMS values.
# ============================================================

report_dir <- "experiment/output"
results_dir <- "experiment/output"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

x_grid <- seq(-0.5, 0.5, length.out = 200)
sigma_fixed <- 1
delta_vals <- c(0, 0.1, 0.2, 0.3, 0.5)
lambda_vals <- c(0.15, 0.25, 0.4)
n_draws <- 4L
seed_base <- 260616L

panel_seed <- function(seed_base, delta, lambda) {
  seed_base + as.integer(round(1000 * delta + 10000 * lambda))
}

expected_rms_tbl <- expand.grid(
  delta = delta_vals,
  lambda = lambda_vals,
  stringsAsFactors = FALSE
) |>
  mutate(
    sigma = sigma_fixed,
    expected_projected_rms = sigma * delta
  )

readr::write_csv(
  expected_rms_tbl,
  file.path(results_dir, "260616_epsilon_calibration_expected_rms.csv")
)

plot_df <- list()
row_id <- 1L

for (delta in delta_vals) {
  for (lambda in lambda_vals) {
    sim <- simulate_projected_scaled_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = n_draws,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base, delta, lambda)
    )

    for (draw_id in seq_len(n_draws)) {
      plot_df[[row_id]] <- data.frame(
        x = x_grid,
        value = sim$f[draw_id, ],
        draw = factor(draw_id),
        delta = sprintf("delta = %.1f", delta),
        lambda = paste0("lambda = ", lambda),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }
}

plot_df <- bind_rows(plot_df) |>
  mutate(
    delta = factor(delta, levels = sprintf("delta = %.1f", delta_vals)),
    lambda = factor(lambda, levels = paste0("lambda = ", lambda_vals))
  )

p <- ggplot(plot_df, aes(x = x, y = value, color = draw, group = draw)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_line(linewidth = 0.85, alpha = 0.92, show.legend = FALSE) +
  facet_grid(delta ~ lambda, scales = "fixed") +
  scale_color_manual(values = c("#0B6E4F", "#C75146", "#2E5AAC", "#8E6C88")) +
  labs(
    title = "Small projected nonlinear effects under the redefined delta",
    subtitle = "Under the projected-scaled parameterization, expected projected RMS is approximately sigma * delta",
    x = "x",
    y = "latent nonlinear function"
  ) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey94", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

out_path <- file.path(report_dir, "260616_epsilon_calibration_figure.png")
ggsave(out_path, p, width = 10.5, height = 10.5, dpi = 220, bg = "white")

message("Saved: ", out_path)
message("Saved: ", file.path(results_dir, "260616_epsilon_calibration_expected_rms.csv"))
