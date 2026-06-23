# ============================================================
# Application 2: Education
#
# Dataset: data/kidiq_data.csv
# - Predictor: mom_iq (mother's IQ), rescaled to [-0.5, 0.5]
# - Response:  kid_score (child test score, raw)
# - Group:     mom_hs (HS completion), deviation-coded {-0.5, +0.5}
#
# Model: fit_interaction_model.stan
#   y = beta_0 + beta_1 x + f_main(x) + z f_int(x) + eps
#   f_main ~ GP(0, sigma^2 delta_main^2 K~_{lambda_main})
#   f_int  ~ GP(0, sigma^2 delta_int^2  K~_{lambda_int})
# ============================================================

source("src/_setup.R")

# --- Data ---
kidiq <- read.csv("data/kidiq_data.csv")

x_raw <- kidiq$mom_iq      # mother's IQ (raw)
y     <- kidiq$kid_score   # child test score
z_bin <- kidiq$mom_hs      # 0/1: mother completed HS

# Rescale x to [-0.5, 0.5]; deviation-code z to {-0.5, +0.5}
x <- rescale_to_unit(x_raw)
z <- z_bin - 0.5

N <- length(y)
cat("N =", N, "\n")
cat("Group sizes: HS=0:", sum(z_bin == 0), "| HS=1:", sum(z_bin == 1), "\n")

# --- Stan data (d_order = 0, local prior) ---
stan_data <- c(
  list(N = N, x = x, y = y, z = z),
  moment_prior_stan_data(0L)
)

# --- Compile interaction model (cached) ---
if (!exists(".int_model", envir = globalenv(), inherits = FALSE)) {
  message("Compiling fit_interaction_model.stan ...")
  .int_model <<- rstan::stan_model(file = "models/fit_interaction_model.stan")
}

# --- Fit ---
fit_int <- rstan::sampling(
  .int_model,
  data    = stan_data,
  iter    = 5500,
  warmup  = 500,
  chains  = 2,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed    = 42,
  refresh = 500,
  verbose = FALSE,
  open_progress = FALSE
)

# --- Convergence ---
diag <- extract_diagnostics(fit_int)
max_rhat <- max(diag$rhat, na.rm = TRUE)
cat(
  "Max Rhat:",
  round(max_rhat, 3),
  "| Divergences:",
  diag$n_divergent,
  "| Treedepth hits:",
  diag$n_max_treedepth,
  "| Mean accept:",
  round(diag$mean_accept_stat, 3),
  "\n"
)
if (diag$n_divergent > 0) {
  cat("  WARNING:", diag$n_divergent, "divergent transitions\n")
}
if (max_rhat > 1.01) {
  cat("  WARNING: Max Rhat > 1.01\n")
}

# --- Posterior summary ---
pars_key <- c("beta1", "delta_main", "delta_int", "lambda_main", "lambda_int", "sigma")
print(
  rstan::summary(fit_int, pars = pars_key)$summary[,
    c("mean", "sd", "2.5%", "97.5%", "Rhat", "n_eff")
  ],
  digits = 3
)

# --- Bayes factors ---
post <- rstan::extract(fit_int)

bf_main <- bf10_interval_null_delta(post$delta_main, d_order = d_order, eps = 0)
cat("log BF10 (delta_main):", round(bf_main$log_bf10, 3), "\n")

bf_int <- bf10_interval_null_delta(post$delta_int, d_order = d_order, eps = 0)
cat("log BF10 (delta_int):", round(bf_int$log_bf10, 3), "\n")

X_mat <- matrix(x, ncol = 1)
prior0_th <- prior0_theta_analytical(X_mat)
post0_th <- post0_theta_logspline(post$theta)
log_bf10_linear <- log(prior0_th / post0_th)
cat("log BF10 (linear slope):", round(log_bf10_linear, 3), "\n")

# ============================================================
# Plot 3: Common trend + group-specific trends (Figure 7 parity)
#
# Posterior predictive mean of the latent signal at (x*, z*):
#   mu*(x*, z*) = y_bar + beta_1 x*_c
#                 + K~_main_cross * a
#                 + z* * (K~_int_cross * (diag(z) * a))
#
# with a = Sigma^{-1} (y_c - x_c beta_1),
# Sigma = K~_main + diag(z) K~_int diag(z) + sigma^2 I.
#
# Three curves:
#   z* = 0     common trend (dashed)
#   z* = -0.5  HS = 0
#   z* = +0.5  HS = 1
# ============================================================

N_grid <- 100

# --- Centering (mirrors Stan transformed data) ---
y_bar_s <- mean(y)
x_bar_s <- mean(x)
y_c_s <- y - y_bar_s
x_c_s <- x - x_bar_s

# --- Orthogonal projectors ---
Q_Z <- qr.Q(qr(cbind(1, x_c_s))) # training: N x 2

# Grid on the original IQ scale
x_grid_raw <- seq(min(x_raw), max(x_raw), length.out = N_grid)
x_grid_s <- rescale_to_unit_like(x_grid_raw, x_raw) # see helper below
x_grid_c <- x_grid_s - x_bar_s
Q_star <- qr.Q(qr(cbind(1, x_grid_c))) # grid:    N_grid x 2

# --- Squared-distance matrices ---
D2_train <- outer(x_c_s, x_c_s, "-")^2
D2_cross <- outer(x_grid_c, x_c_s, "-")^2

# --- Kernel helpers ---
se_kern <- function(D2, alpha, l) alpha^2 * exp(-0.5 * D2 / l^2)

ktilde <- function(K, Q) {
  U <- K %*% Q
  M <- crossprod(Q, U)
  Kt <- K - U %*% t(Q) - Q %*% t(U) + Q %*% M %*% t(Q)
  0.5 * (Kt + t(Kt))
}

# Cross-kernel tilde: P_perp_pred K P_perp_train
ktilde_cross <- function(Kc, Q_pred, Q_train) {
  A <- Kc %*% Q_train # N_grid x 2
  B_Kc <- Q_pred %*% crossprod(Q_pred, Kc) # N_grid x N
  Kc - A %*% t(Q_train) - B_Kc + Q_pred %*% crossprod(Q_pred, A) %*% t(Q_train)
}

# --- Posterior subsample ---
set.seed(1)
S_all <- length(post$sigma)
nd <- min(400, S_all)
idx <- sample.int(S_all, nd)

sigma_s <- post$sigma[idx]
beta1_s <- post$beta1[idx]
dm_s <- post$delta_main[idx]
lam_m_s <- post$lambda_main[idx]
di_s <- post$delta_int[idx]
lam_i_s <- post$lambda_int[idx]

# Storage: three groups -- common (z*=0), HS=0 (z*=-0.5), HS=1 (z*=+0.5)
z_stars <- c(0, -0.5, 0.5)
grp_labels <- c("Common", "No HS", "HS")
pred_mat <- lapply(z_stars, function(zg) matrix(NA_real_, nd, N_grid))

for (s in seq_len(nd)) {
  sig <- sigma_s[s]
  b1 <- beta1_s[s]
  am <- sig * dm_s[s]
  ai <- sig * di_s[s]
  lam_m <- lam_m_s[s]
  lam_i <- lam_i_s[s]

  # Training kernels and P_perp-projected versions
  Km <- se_kern(D2_train, am, lam_m)
  Ki <- se_kern(D2_train, ai, lam_i)
  diag(Km) <- diag(Km) + 1e-10 # jitter, matches Stan nugget
  diag(Ki) <- diag(Ki) + 1e-10

  Ktm <- ktilde(Km, Q_Z)
  Kti <- ktilde(Ki, Q_Z)

  # D_z K~_int D_z: scale rows and columns by z
  DKD <- z * t(z * t(Kti))

  # Sigma and Cholesky
  Sigma <- Ktm + DKD
  diag(Sigma) <- diag(Sigma) + sig^2
  Sigma <- 0.5 * (Sigma + t(Sigma))
  L <- chol(Sigma) # upper-triangular R s.t. R' R = Sigma

  # a = Sigma^{-1} (y_c - x_c beta_1)
  resid <- y_c_s - x_c_s * b1
  a_vec <- backsolve(L, forwardsolve(t(L), resid))

  # Tilde-projected cross-covariances
  Kc_m_tilde <- ktilde_cross(se_kern(D2_cross, am, lam_m), Q_star, Q_Z)
  Kc_i_tilde <- ktilde_cross(se_kern(D2_cross, ai, lam_i), Q_star, Q_Z)

  f_main_star <- as.numeric(Kc_m_tilde %*% a_vec)
  f_int_star <- as.numeric(Kc_i_tilde %*% (z * a_vec))
  lin_star <- y_bar_s + b1 * x_grid_c

  for (gi in seq_along(z_stars)) {
    pred_mat[[gi]][s, ] <- lin_star + f_main_star + z_stars[gi] * f_int_star
  }
}

# --- Assemble for plotting ---
plot_df <- do.call(
  rbind,
  lapply(seq_along(z_stars), function(gi) {
    mat <- pred_mat[[gi]]
    data.frame(
      x = x_grid_raw,
      mean = colMeans(mat),
      lo = apply(mat, 2, quantile, 0.025),
      hi = apply(mat, 2, quantile, 0.975),
      group = grp_labels[gi]
    )
  })
)
plot_df$group <- factor(plot_df$group, levels = grp_labels)

df_raw <- data.frame(
  mom_iq = x_raw,
  kid_score = y,
  group = factor(z_bin, levels = c(0, 1), labels = c("No HS", "HS"))
)

# Colors and linetypes to mirror draft Figure 7
grp_colors <- c("Common" = "black", "No HS" = "steelblue", "HS" = "darkorange")
grp_lty <- c("Common" = "dashed", "No HS" = "solid", "HS" = "solid")

# Show ribbons only for the two group-specific curves, not the common trend
plot_df_ribbon <- plot_df[plot_df$group != "Common", ]

ggplot2::ggplot() +
  ggplot2::geom_point(
    data = df_raw,
    ggplot2::aes(x = mom_iq, y = kid_score, color = group),
    size = 0.9,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  ggplot2::geom_ribbon(
    data = plot_df_ribbon,
    ggplot2::aes(x = x, ymin = lo, ymax = hi, fill = group),
    alpha = 0.20
  ) +
  ggplot2::geom_line(
    data = plot_df,
    ggplot2::aes(x = x, y = mean, color = group, linetype = group),
    linewidth = 1.0
  ) +
  ggplot2::scale_color_manual(values = grp_colors) +
  ggplot2::scale_fill_manual(values = grp_colors[c("No HS", "HS")]) +
  ggplot2::scale_linetype_manual(values = grp_lty) +
  ggplot2::labs(
    title = "Mother's IQ vs. Child Test Scores",
    subtitle = "Common nonlinear trend (dashed) and group-specific trends (95% CI)",
    x = "Mother's IQ",
    y = "Child test score",
    color = NULL,
    linetype = NULL,
    fill = NULL
  ) +
  .posterior_theme()

# --- Save ---
saveRDS(fit_int, "results/fit_02_education_interaction.rds")

# ============================================================
# BF10 method comparison for delta_main and delta_int
# Methods: logspline, raw (ecdf), laplace, jeffreys
# ============================================================

EPS     <- 0.2
d_order <- 0L
prior0  <- prior_prob_delta(d_order, EPS)

compute_group_bf <- function(x, eps, prior0) {
  x      <- x[is.finite(x) & x >= 0]
  n      <- length(x)
  x_sort <- sort(x)
  k      <- findInterval(eps, x_sort)
  data.frame(
    k              = k,
    n              = n,
    prior0         = prior0,
    post0_laplace  = (k + 1)   / (n + 2),
    post0_jeffreys = (k + 0.5) / (n + 1)
  )
}

delta_list  <- list(delta_main = post$delta_main, delta_int = post$delta_int)
bf_rows <- lapply(names(delta_list), function(nm) {
  x         <- delta_list[[nm]]
  post0_ls  <- post_prob_delta_logspline(x, EPS)
  post0_raw <- post_prob_delta_ecdf(x, EPS)$post0
  gb        <- compute_group_bf(x, EPS, prior0)
  data.frame(
    predictor          = nm,
    log_bf10_logspline = round(log(prior0 / post0_ls),             4),
    log_bf10_raw       = round(log(prior0 / post0_raw),            4),
    log_bf10_laplace   = round(log(prior0 / gb$post0_laplace),     4),
    log_bf10_jeffreys  = round(log(prior0 / gb$post0_jeffreys),    4)
  )
})

bf_tbl <- do.call(rbind, bf_rows)
cat(sprintf("\neps = %.1f,  prior0 = %.6f,  d_order = %d\n\n", EPS, prior0, d_order))
print(bf_tbl, row.names = FALSE)


# ============================================================
# Bridge-sampling validation of SDR Bayes factors
#
# Two comparisons:
#   1. Main nonlinear BF:
#      Numerator:   full interaction model (fit_int, already fitted above)
#      Denominator: interaction-only model (no f_main)
#                   -> fit_interaction_only_model.stan
#
#   2. Interaction nonlinear BF:
#      Numerator:   full interaction model (fit_int)
#      Denominator: main-only model (no f_int)
#                   -> fit_full_model.stan with P = 1
# ============================================================

cat("\n========== Bridge-Sampling Validation ==========\n")

# ----------------------------------------------------------
# Reduced model 1: interaction-only (no f_main)
# ----------------------------------------------------------
if (!exists(".int_only_model", envir = globalenv(), inherits = FALSE)) {
  message("Compiling fit_interaction_only_model.stan ...")
  .int_only_model <<- rstan::stan_model(
    file = "models/fit_interaction_only_model.stan"
  )
}

fit_int_only <- rstan::sampling(
  .int_only_model,
  data = stan_data, # same data as the full interaction model
  iter = 5500,
  warmup = 500,
  chains = 2,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 42,
  refresh = 500,
  verbose = FALSE,
  open_progress = FALSE
)

diag_io <- extract_diagnostics(fit_int_only)
cat(
  "Interaction-only model -- Max Rhat:",
  round(max(diag_io$rhat, na.rm = TRUE), 3),
  "| Divergences:",
  diag_io$n_divergent,
  "\n"
)

# ----------------------------------------------------------
# Reduced model 2: main-only (no f_int)
#   = fit_full_model.stan with P = 1, using x_raw as the single predictor.
#   fit_model() handles rescaling, moment-prior data, etc.
# ----------------------------------------------------------
fit_main_only <- fit_model(
  X = matrix(x_raw, ncol = 1),
  y = y,
  model = "full",
  d_order = d_order,
  iter = 5500,
  warmup = 500,
  chains = 2,
  adapt_delta = 0.95,
  max_treedepth = 12,
  seed = 42,
  refresh = 500
)

diag_mo <- extract_diagnostics(fit_main_only$fit)
cat(
  "Main-only model -- Max Rhat:",
  round(max(diag_mo$rhat, na.rm = TRUE), 3),
  "| Divergences:",
  diag_mo$n_divergent,
  "\n"
)

# ----------------------------------------------------------
# Bridge sampling
# ----------------------------------------------------------
bridge_full <- bridgesampling::bridge_sampler(fit_int, silent = TRUE)
bridge_int_only <- bridgesampling::bridge_sampler(fit_int_only, silent = TRUE)
bridge_main_only <- bridgesampling::bridge_sampler(fit_main_only$fit, silent = TRUE)

# BF10 for main nonlinear effect: full / interaction-only
log_bf10_main_bs <- bridgesampling::bf(bridge_full, bridge_int_only, log = TRUE)$bf

# BF10 for interaction nonlinear effect: full / main-only
log_bf10_int_bs <- bridgesampling::bf(bridge_full, bridge_main_only, log = TRUE)$bf

cat("\n--- Bridge-sampling BF10 ---\n")
cat("log BF10 (main nonlinear, BS):", round(log_bf10_main_bs, 3), "\n")
cat("log BF10 (interaction nonlinear, BS):", round(log_bf10_int_bs, 3), "\n")

cat("\n--- SDR BF10 (from above) ---\n")
cat("log BF10 (main nonlinear, SDR):", round(bf_main$log_bf10, 3), "\n")
cat("log BF10 (interaction nonlinear, SDR):", round(bf_int$log_bf10, 3), "\n")

cat("\n--- Comparison (BS - SDR) ---\n")
cat("Main:        ", round(log_bf10_main_bs - bf_main$log_bf10, 3), "\n")
cat("Interaction: ", round(log_bf10_int_bs - bf_int$log_bf10, 3), "\n")

# --- Save reduced fits ---
saveRDS(fit_int_only, "results/fit_02_education_int_only.rds")
saveRDS(fit_main_only, "results/fit_02_education_main_only.rds")



#### ===========================================================================
# Need stan_data to be loaded (run data prep section of the script first)
.int_model <- rstan::stan_model(file = "models/fit_interaction_model.stan")
fit_dummy <- rstan::sampling(
  .int_model, data = stan_data,
  iter = 2, chains = 1, warmup = 0,
  refresh = 0, open_progress = FALSE
)
# Transplant the C++ backend into the loaded fit
fit_int@.MISC <- fit_dummy@.MISC

# Restore fit_int_only
.int_only_model <- rstan::stan_model(file = "models/fit_interaction_only_model.stan")
fit_dummy_io <- rstan::sampling(
  .int_only_model, data = stan_data,
  iter = 2, chains = 1, warmup = 0,
  refresh = 0, open_progress = FALSE
)
fit_int_only@.MISC <- fit_dummy_io@.MISC

# Restore fit_main_only$fit
.full_model <- rstan::stan_model(file = "models/fit_full_model.stan")
fit_dummy_mo <- rstan::sampling(
  .full_model, data = fit_main_only$stan_data,
  iter = 2, chains = 1, warmup = 0,
  refresh = 0, open_progress = FALSE
)
fit_main_only$fit@.MISC <- fit_dummy_mo@.MISC

