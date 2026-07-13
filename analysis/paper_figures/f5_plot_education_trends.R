source("src/_setup.R")

output_dir <- "analysis/paper_figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

rescale_to_unit_like <- function(v, ref) {
  ref_min <- min(ref, na.rm = TRUE)
  ref_max <- max(ref, na.rm = TRUE)
  (v - ref_min) / (ref_max - ref_min) - 0.5
}

se_kern <- function(D2, alpha, l) {
  (alpha^2) * exp(-0.5 * D2 / (l^2))
}

ktilde <- function(K, Q) {
  U <- K %*% Q
  M <- crossprod(Q, U)
  Kt <- K - U %*% t(Q) - Q %*% t(U) + Q %*% M %*% t(Q)
  0.5 * (Kt + t(Kt))
}

ktilde_cross <- function(Kc, Q_pred, Q_train) {
  A <- Kc %*% Q_train
  B_Kc <- Q_pred %*% crossprod(Q_pred, Kc)
  Kc - A %*% t(Q_train) - B_Kc + Q_pred %*% crossprod(Q_pred, A) %*% t(Q_train)
}

chol_with_jitter <- function(Sigma, jitter_seq = c(0, 1e-10, 1e-8, 1e-6, 1e-4)) {
  Sigma <- 0.5 * (Sigma + t(Sigma))
  for (jit in jitter_seq) {
    Sigma_try <- Sigma
    if (jit > 0) diag(Sigma_try) <- diag(Sigma_try) + jit
    U <- try(chol(Sigma_try), silent = TRUE)
    if (!inherits(U, "try-error")) return(U)
  }
  stop("Cholesky decomposition failed for all jitter values.")
}

kidiq <- read.csv("data/kidiq_data.csv")
fit_int <- readRDS("results/fit_02_education_interaction.rds")
post <- rstan::extract(fit_int)

x_raw <- kidiq$mom_iq
y <- kidiq$kid_score
z_bin <- kidiq$mom_hs

x <- rescale_to_unit(x_raw)
z <- z_bin - 0.5

N_grid <- 120
S_all <- length(post$sigma)
nd <- min(400, S_all)
idx <- sample.int(S_all, nd)

y_bar <- mean(y)
x_bar <- mean(x)
y_c <- y - y_bar
x_c <- x - x_bar

Q_Z <- qr.Q(qr(cbind(1, x_c)))

x_grid_raw <- seq(min(x_raw), max(x_raw), length.out = N_grid)
x_grid_scaled <- rescale_to_unit_like(x_grid_raw, x_raw)
x_grid_c <- x_grid_scaled - x_bar
Q_star <- qr.Q(qr(cbind(1, x_grid_c)))

D2_train <- outer(x_c, x_c, "-")^2
D2_cross <- outer(x_grid_c, x_c, "-")^2

sigma_s <- post$sigma[idx]
beta1_s <- post$beta1[idx]
dm_s <- post$delta_main[idx]
lam_m_s <- post$lambda_main[idx]
di_s <- post$delta_int[idx]
lam_i_s <- post$lambda_int[idx]

z_star_values <- c(0, -0.5, 0.5)
curve_levels <- c("Common", "HS = 0", "HS = 1")
curve_labels <- setNames(curve_levels, z_star_values)
pred_mat <- lapply(z_star_values, function(zg) matrix(NA_real_, nd, N_grid))

for (s in seq_len(nd)) {
  sig <- sigma_s[s]
  b1 <- beta1_s[s]
  am <- sig * dm_s[s]
  ai <- sig * di_s[s]
  lam_m <- lam_m_s[s]
  lam_i <- lam_i_s[s]

  Km <- se_kern(D2_train, am, lam_m)
  Ki <- se_kern(D2_train, ai, lam_i)
  diag(Km) <- diag(Km) + 1e-10
  diag(Ki) <- diag(Ki) + 1e-10

  Ktm <- ktilde(Km, Q_Z)
  Kti <- ktilde(Ki, Q_Z)
  DKD <- z * t(z * t(Kti))

  Sigma <- Ktm + DKD
  diag(Sigma) <- diag(Sigma) + sig^2
  U <- chol_with_jitter(Sigma)

  resid <- y_c - x_c * b1
  a_vec <- backsolve(U, forwardsolve(t(U), resid))

  Kc_m_tilde <- ktilde_cross(se_kern(D2_cross, am, lam_m), Q_star, Q_Z)
  Kc_i_tilde <- ktilde_cross(se_kern(D2_cross, ai, lam_i), Q_star, Q_Z)

  f_main_star <- as.numeric(Kc_m_tilde %*% a_vec)
  f_int_star <- as.numeric(Kc_i_tilde %*% (z * a_vec))
  lin_star <- y_bar + b1 * x_grid_c

  for (gi in seq_along(z_star_values)) {
    pred_mat[[gi]][s, ] <- lin_star + f_main_star + z_star_values[gi] * f_int_star
  }
}

plot_df <- do.call(
  rbind,
  lapply(seq_along(z_star_values), function(gi) {
    mat <- pred_mat[[gi]]
    data.frame(
      x = x_grid_raw,
      mean = colMeans(mat),
      lo = apply(mat, 2, quantile, 0.025),
      hi = apply(mat, 2, quantile, 0.975),
      curve = curve_levels[gi]
    )
  })
)
plot_df$curve <- factor(plot_df$curve, levels = curve_levels)

plot_df_ribbon <- subset(plot_df, curve != "Common")

scatter_df <- data.frame(
  mom_iq = x_raw,
  kid_score = y,
  curve = factor(ifelse(z_bin == 0, "HS = 0", "HS = 1"), levels = curve_levels)
)

curve_colors <- c(
  "Common" = "black",
  "HS = 0" = "#2C7FB8",
  "HS = 1" = "#E69F00"
)
curve_fills <- c(
  "Common" = grDevices::adjustcolor("black", alpha.f = 0),
  "HS = 0" = "#2C7FB8",
  "HS = 1" = "#E69F00"
)
curve_linetypes <- c(
  "Common" = "22",
  "HS = 0" = "solid",
  "HS = 1" = "solid"
)

p <- ggplot() +
  geom_point(
    data = scatter_df,
    aes(x = mom_iq, y = kid_score, color = curve),
    size = 1.4,
    alpha = 0.34,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = plot_df_ribbon,
    aes(x = x, ymin = lo, ymax = hi, fill = curve),
    alpha = 0.18,
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(
    data = plot_df,
    aes(x = x, y = mean, color = curve, linetype = curve),
    linewidth = 1.05,
    show.legend = TRUE
  ) +
  scale_color_manual(values = curve_colors, breaks = curve_levels, name = NULL) +
  scale_fill_manual(values = curve_fills, breaks = curve_levels, name = NULL) +
  scale_linetype_manual(values = curve_linetypes, breaks = curve_levels, name = NULL) +
  labs(
    x = "Mother's IQ",
    y = "Child test score"
  ) +
  guides(
    fill = "none",
    color = guide_legend(order = 1, nrow = 1, byrow = TRUE, override.aes = list(alpha = 1)),
    linetype = guide_legend(order = 1, nrow = 1, byrow = TRUE)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = grid::unit(0.8, "lines"),
    plot.margin = margin(6, 8, 6, 8)
  )

ggsave(
  filename = file.path(output_dir, "fig_education_trends.png"),
  plot = p,
  width = 7.2,
  height = 4.8,
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "plot-app-2.png"),
  plot = p,
  width = 7.2,
  height = 4.8,
  dpi = 320,
  bg = "white"
)
