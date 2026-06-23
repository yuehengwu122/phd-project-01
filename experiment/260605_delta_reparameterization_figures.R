source("src/_setup.R")

library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# Delta reparameterization figures
#
# Produces:
# - raw vs orthogonalized prior draws
# - raw-draw decomposition
# - RMS comparison under old vs projected-scaled delta
# - draw comparison under old vs projected-scaled delta
# - Figure-1-like projected-scaled prior draws
# ============================================================

report_dir <- "experiment/output"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

x_grid <- seq(-0.5, 0.5, length.out = 200)
sigma_fixed <- 1
n_draws <- 4L

delta_vals_old <- c(1, 2)
lambda_vals_old <- c(0.15, 0.3, 0.6)
delta_vals_rms <- c(0.5, 1, 2)
lambda_vals_rms <- c(0.15, 0.3, 0.6)
delta_vals_fig1 <- c(0, 1, 2, 3)
lambda_vals_fig1 <- c(0.15, 0.25, 0.4)

seed_base_old <- 260605L
seed_base_fig1 <- 260611L

theme_gp <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

panel_seed <- function(seed_base, delta, lambda) {
  seed_base + as.integer(100 * delta + 1000 * lambda)
}

projector_terms <- function(x) {
  x_mat <- matrix(x, ncol = 1)
  x_c <- scale(x_mat, center = TRUE, scale = FALSE)
  Z_aug <- cbind(1, x_c)
  Q_Z <- qr.Q(qr(Z_aug))
  list(x_c = x_c, Q_Z = Q_Z, P_perp = diag(nrow(x_mat)) - Q_Z %*% t(Q_Z))
}

expected_projected_rms_old <- function(x, delta, lambda, sigma) {
  proj <- projector_terms(x)
  D2 <- outer(x, x, FUN = function(a, b) (a - b)^2)
  R_lambda <- exp(-0.5 * D2 / lambda^2)
  K_tilde <- proj$P_perp %*% R_lambda %*% proj$P_perp
  sigma * delta * sqrt(sum(diag(K_tilde)) / length(x))
}

expected_projected_rms_new <- function(delta, sigma) {
  sigma * delta
}

# ============================================================
# 1) Raw vs orthogonalized prior draws
# ============================================================

plot_df <- list()
row_id <- 1L

for (delta in delta_vals_old) {
  for (lambda in lambda_vals_old) {
    sim_raw <- simulate_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = n_draws,
      orthogonalize = FALSE,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base_old, delta, lambda)
    )

    sim_orth <- simulate_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = n_draws,
      orthogonalize = TRUE,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base_old, delta, lambda)
    )

    for (draw_id in seq_len(n_draws)) {
      plot_df[[row_id]] <- data.frame(
        x = x_grid,
        value = sim_raw$f[draw_id, ],
        draw = factor(draw_id),
        type = "Raw GP",
        delta = paste0("delta = ", delta),
        lambda = paste0("lambda = ", lambda)
      )
      row_id <- row_id + 1L

      plot_df[[row_id]] <- data.frame(
        x = x_grid,
        value = sim_orth$f[draw_id, ],
        draw = factor(draw_id),
        type = "Orthogonalized GP",
        delta = paste0("delta = ", delta),
        lambda = paste0("lambda = ", lambda)
      )
      row_id <- row_id + 1L
    }
  }
}

plot_df <- bind_rows(plot_df) |>
  mutate(
    type = factor(type, levels = c("Raw GP", "Orthogonalized GP")),
    delta = factor(delta, levels = paste0("delta = ", delta_vals_old)),
    lambda = factor(lambda, levels = paste0("lambda = ", lambda_vals_old))
  )

p_raw_vs_orth <- ggplot(plot_df, aes(x = x, y = value, group = interaction(draw, type), color = draw)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_line(linewidth = 0.8, alpha = 0.9, show.legend = FALSE) +
  facet_grid(delta ~ lambda + type, scales = "fixed") +
  scale_color_manual(values = c("#0B6E4F", "#C75146", "#2E5AAC", "#8E6C88")) +
  labs(
    title = "Raw GP vs orthogonalized nonlinear GP draws",
    subtitle = "Intercept and slope fixed at zero; orthogonalized draws are projected onto span{1, x_c}^perp",
    x = "x",
    y = "latent function"
  ) +
  theme_gp()

ggsave(
  file.path(report_dir, "260605_raw_vs_orthogonalized_gp.png"),
  p_raw_vs_orth, width = 14, height = 6, dpi = 220, bg = "white"
)

# ============================================================
# 2) Raw-draw decomposition
# ============================================================

decomp_df <- list()
row_id <- 1L
proj_old <- projector_terms(x_grid)

for (delta in delta_vals_old) {
  for (lambda in lambda_vals_old) {
    sim_raw <- simulate_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = 1L,
      orthogonalize = FALSE,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base_old, delta, lambda)
    )

    f_raw <- as.numeric(sim_raw$f[1, ])
    f_linear <- as.numeric(proj_old$Q_Z %*% (crossprod(proj_old$Q_Z, f_raw)))
    f_orth <- f_raw - f_linear

    decomp_df[[row_id]] <- data.frame(
      x = x_grid,
      value = f_raw,
      component = "Raw draw",
      delta = paste0("delta = ", delta),
      lambda = paste0("lambda = ", lambda)
    )
    row_id <- row_id + 1L

    decomp_df[[row_id]] <- data.frame(
      x = x_grid,
      value = f_linear,
      component = "Projected linear part",
      delta = paste0("delta = ", delta),
      lambda = paste0("lambda = ", lambda)
    )
    row_id <- row_id + 1L

    decomp_df[[row_id]] <- data.frame(
      x = x_grid,
      value = f_orth,
      component = "Orthogonal residual",
      delta = paste0("delta = ", delta),
      lambda = paste0("lambda = ", lambda)
    )
    row_id <- row_id + 1L
  }
}

decomp_df <- bind_rows(decomp_df) |>
  mutate(
    component = factor(component, levels = c("Raw draw", "Projected linear part", "Orthogonal residual")),
    delta = factor(delta, levels = paste0("delta = ", delta_vals_old)),
    lambda = factor(lambda, levels = paste0("lambda = ", lambda_vals_old))
  )

p_decomp <- ggplot(decomp_df, aes(x = x, y = value, color = component)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_line(linewidth = 0.9) +
  facet_grid(delta ~ lambda) +
  scale_color_manual(values = c(
    "Raw draw" = "#404040",
    "Projected linear part" = "#C75146",
    "Orthogonal residual" = "#2E5AAC"
  )) +
  labs(
    title = "Decomposition of a raw GP draw",
    subtitle = "Each raw draw is split into its projection onto span{1, x_c} and its orthogonal residual",
    x = "x",
    y = "latent function",
    color = NULL
  ) +
  theme_gp() +
  theme(legend.position = "bottom")

ggsave(
  file.path(report_dir, "260605_raw_gp_decomposition.png"),
  p_decomp, width = 10, height = 6, dpi = 220, bg = "white"
)

# ============================================================
# 3) RMS comparison: old vs projected-scaled delta
# ============================================================

summary_df <- expand.grid(
  delta = delta_vals_rms,
  lambda = lambda_vals_rms,
  stringsAsFactors = FALSE
) |>
  mutate(
    old_rms = mapply(expected_projected_rms_old, MoreArgs = list(x = x_grid, sigma = sigma_fixed),
                     delta = delta, lambda = lambda),
    new_rms = mapply(expected_projected_rms_new, delta = delta, sigma = sigma_fixed)
  ) |>
  pivot_longer(
    cols = c(old_rms, new_rms),
    names_to = "parameterization",
    values_to = "expected_rms"
  ) |>
  mutate(
    parameterization = recode(
      parameterization,
      old_rms = "Current delta definition",
      new_rms = "Projected-scaled delta definition"
    ),
    delta_label = factor(paste0("delta = ", delta), levels = paste0("delta = ", delta_vals_rms))
  )

p_rms <- ggplot(summary_df, aes(x = lambda, y = expected_rms, color = parameterization, group = parameterization)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~delta_label, nrow = 1) +
  scale_color_manual(values = c(
    "Current delta definition" = "#C75146",
    "Projected-scaled delta definition" = "#2E5AAC"
  )) +
  labs(
    title = "Projected nonlinear RMS under old and new delta definitions",
    subtitle = "Under the projected-scaled parameterization, delta targets the final projected nonlinear magnitude",
    x = "lambda",
    y = "Expected projected RMS",
    color = NULL
  ) +
  theme_gp() +
  theme(legend.position = "bottom")

ggsave(
  file.path(report_dir, "260605_projected_scaled_rms_comparison.png"),
  p_rms, width = 10, height = 4.5, dpi = 220, bg = "white"
)

# ============================================================
# 4) Draw comparison: old vs projected-scaled delta
# ============================================================

draw_df <- list()
row_id <- 1L

for (delta in c(1, 2)) {
  for (lambda in lambda_vals_rms) {
    sim_old <- simulate_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = n_draws,
      orthogonalize = TRUE,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base_old, delta, lambda)
    )

    sim_new <- simulate_projected_scaled_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = n_draws,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base_old, delta, lambda)
    )

    for (draw_id in seq_len(n_draws)) {
      draw_df[[row_id]] <- data.frame(
        x = x_grid,
        value = sim_old$f[draw_id, ],
        draw = factor(draw_id),
        parameterization = "Current delta definition",
        delta = paste0("delta = ", delta),
        lambda = paste0("lambda = ", lambda)
      )
      row_id <- row_id + 1L

      draw_df[[row_id]] <- data.frame(
        x = x_grid,
        value = sim_new$f[draw_id, ],
        draw = factor(draw_id),
        parameterization = "Projected-scaled delta definition",
        delta = paste0("delta = ", delta),
        lambda = paste0("lambda = ", lambda)
      )
      row_id <- row_id + 1L
    }
  }
}

draw_df <- bind_rows(draw_df) |>
  mutate(
    parameterization = factor(parameterization,
                              levels = c("Current delta definition", "Projected-scaled delta definition")),
    delta = factor(delta, levels = paste0("delta = ", c(1, 2))),
    lambda = factor(lambda, levels = paste0("lambda = ", lambda_vals_rms))
  )

p_draws <- ggplot(draw_df, aes(x = x, y = value, color = draw, group = interaction(draw, parameterization))) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_line(linewidth = 0.8, alpha = 0.9, show.legend = FALSE) +
  facet_grid(delta ~ lambda + parameterization, scales = "fixed") +
  scale_color_manual(values = c("#0B6E4F", "#C75146", "#2E5AAC", "#8E6C88")) +
  labs(
    title = "Projected nonlinear draws under old and new delta definitions",
    subtitle = "Same delta and lambda inputs; the new definition stabilizes projected nonlinear magnitude across lambda",
    x = "x",
    y = "latent function"
  ) +
  theme_gp()

ggsave(
  file.path(report_dir, "260605_projected_scaled_draw_comparison.png"),
  p_draws, width = 14, height = 6, dpi = 220, bg = "white"
)

# ============================================================
# 5) Figure-1-like projected-scaled prior draws
# ============================================================

fig1_df <- list()
row_id <- 1L

for (delta in delta_vals_fig1) {
  for (lambda in lambda_vals_fig1) {
    sim <- simulate_projected_scaled_gp_draws(
      x = x_grid,
      delta = delta,
      lambda = lambda,
      sigma = sigma_fixed,
      n_draws = n_draws,
      beta0 = 0,
      beta1 = 0,
      add_noise = FALSE,
      seed = panel_seed(seed_base_fig1, delta, lambda)
    )

    for (draw_id in seq_len(n_draws)) {
      fig1_df[[row_id]] <- data.frame(
        x = x_grid,
        value = sim$f[draw_id, ],
        draw = factor(draw_id),
        delta = paste0("delta = ", delta),
        lambda = paste0("lambda = ", lambda)
      )
      row_id <- row_id + 1L
    }
  }
}

fig1_df <- bind_rows(fig1_df) |>
  mutate(
    delta = factor(delta, levels = paste0("delta = ", delta_vals_fig1)),
    lambda = factor(lambda, levels = paste0("lambda = ", lambda_vals_fig1))
  )

p_fig1 <- ggplot(fig1_df, aes(x = x, y = value, color = draw, group = draw)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
  geom_line(linewidth = 0.85, alpha = 0.92, show.legend = FALSE) +
  facet_grid(delta ~ lambda, scales = "fixed") +
  scale_color_manual(values = c("#0B6E4F", "#C75146", "#2E5AAC", "#8E6C88")) +
  labs(
    title = "Projected nonlinear prior draws under the redefined delta",
    subtitle = "Each panel shows latent draws after projection and standardization; delta now targets the final nonlinear magnitude",
    x = "x",
    y = "latent nonlinear function"
  ) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme_gp()

ggsave(
  file.path(report_dir, "260605_projected_scaled_figure1_like.png"),
  p_fig1, width = 10.5, height = 9.2, dpi = 220, bg = "white"
)
