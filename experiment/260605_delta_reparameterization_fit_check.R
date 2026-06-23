source("src/_setup.R")

library(dplyr)
library(ggplot2)
library(readr)

# ============================================================
# Calibration fit check for delta reparameterization
#
# Produces:
# - results/simulated/260605_delta_calibration_fit.csv
# - results/simulated/260605_delta_calibration_fit_summary.csv
# - analysis/simulation/report/260605_projected_rms_vs_lambda.png
# - analysis/simulation/report/260605_delta_posterior_vs_projected_rms.png
# ============================================================

report_dir <- "experiment/output"
results_dir <- "experiment/output"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

N <- 80L
x_grid <- seq(-0.5, 0.5, length.out = N)
X <- matrix(x_grid, ncol = 1)

sigma_fixed <- 1
delta_vals <- c(0.5, 1, 2)
lambda_vals <- c(0.15, 0.3, 0.6)
n_reps <- 3L

fit_chains <- 1L
fit_iter <- 800L
fit_warmup <- 400L
fit_seed_base <- 260605L

projected_rms_theoretical <- function(x, delta, lambda, sigma) {
  x_mat <- matrix(x, ncol = 1)
  x_c <- scale(x_mat, center = TRUE, scale = FALSE)
  Z_aug <- cbind(1, x_c)
  Q_Z <- qr.Q(qr(Z_aug))
  P_perp <- diag(nrow(x_mat)) - Q_Z %*% t(Q_Z)

  D2 <- outer(x, x, FUN = function(a, b) (a - b)^2)
  R_lambda <- exp(-0.5 * D2 / lambda^2)
  K_tilde <- P_perp %*% R_lambda %*% P_perp

  sigma * delta * sqrt(sum(diag(K_tilde)) / nrow(x_mat))
}

rows <- list()
row_id <- 1L

for (delta in delta_vals) {
  for (lambda in lambda_vals) {
    expected_rms <- projected_rms_theoretical(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed
    )

    for (rep_id in seq_len(n_reps)) {
      sim_seed <- fit_seed_base + rep_id + as.integer(1000 * delta + 10000 * lambda)

      sim <- simulate_gp_draws(
        x = x_grid,
        delta = delta,
        lambda = lambda,
        sigma = sigma_fixed,
        n_draws = 1L,
        orthogonalize = TRUE,
        beta0 = 0,
        beta1 = 0,
        add_noise = TRUE,
        seed = sim_seed
      )

      y <- as.numeric(sim$y[1, ])
      f <- as.numeric(sim$f[1, ])

      fit <- fit_model(
        X = X,
        y = y,
        model = "full",
        d_order = 0L,
        chains = fit_chains,
        iter = fit_iter,
        warmup = fit_warmup,
        adapt_delta = 0.95,
        max_treedepth = 12L,
        seed = sim_seed,
        refresh = 0L
      )

      post <- rstan::extract(fit$fit, pars = c("delta", "lambda"))
      delta_draws <- as.numeric(post$delta)
      lambda_draws <- as.numeric(post$lambda)

      rows[[row_id]] <- data.frame(
        delta_input = delta,
        lambda_input = lambda,
        replicate = rep_id,
        projected_rms_expected = expected_rms,
        projected_rms_draw = sqrt(mean(f^2)),
        y_rms = sqrt(mean(y^2)),
        delta_post_mean = mean(delta_draws),
        delta_post_median = median(delta_draws),
        delta_post_q10 = as.numeric(quantile(delta_draws, 0.10)),
        delta_post_q90 = as.numeric(quantile(delta_draws, 0.90)),
        lambda_post_mean = mean(lambda_draws),
        lambda_post_median = median(lambda_draws),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L

      rm(fit, post, delta_draws, lambda_draws)
      gc()
    }
  }
}

res <- bind_rows(rows) |>
  mutate(
    delta_label = factor(paste0("delta input = ", delta_input),
                         levels = paste0("delta input = ", delta_vals)),
    lambda_label = factor(paste0("lambda input = ", lambda_input),
                          levels = paste0("lambda input = ", lambda_vals))
  )

summary_tbl <- res |>
  group_by(delta_input, lambda_input) |>
  summarise(
    projected_rms_expected = mean(projected_rms_expected),
    projected_rms_draw_mean = mean(projected_rms_draw),
    delta_post_median_mean = mean(delta_post_median),
    delta_post_median_sd = sd(delta_post_median),
    lambda_post_median_mean = mean(lambda_post_median),
    .groups = "drop"
  )

readr::write_csv(res, file.path(results_dir, "260605_delta_calibration_fit.csv"))
readr::write_csv(summary_tbl, file.path(results_dir, "260605_delta_calibration_fit_summary.csv"))

p_rms <- ggplot(summary_tbl, aes(x = lambda_input, y = projected_rms_expected,
                                 color = factor(delta_input), group = factor(delta_input))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_manual(values = c("#2E5AAC", "#0B6E4F", "#C75146"), name = "delta input") +
  labs(
    title = "Projected nonlinear RMS implied by the current parameterization",
    x = "lambda input",
    y = "Expected projected RMS"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_delta <- ggplot(summary_tbl, aes(x = projected_rms_draw_mean, y = delta_post_median_mean,
                                   color = factor(lambda_input))) +
  geom_point(size = 2.5) +
  geom_text(aes(label = paste0("delta=", delta_input)), nudge_y = 0.05, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("#8E6C88", "#2E5AAC", "#C75146"), name = "lambda input") +
  labs(
    title = "Posterior delta versus observed projected nonlinear RMS",
    x = "Mean draw-level projected RMS",
    y = "Mean posterior median delta"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(report_dir, "260605_projected_rms_vs_lambda.png"),
  p_rms, width = 7, height = 4.5, dpi = 220, bg = "white"
)
ggsave(
  file.path(report_dir, "260605_delta_posterior_vs_projected_rms.png"),
  p_delta, width = 7, height = 4.5, dpi = 220, bg = "white"
)

print(summary_tbl)
