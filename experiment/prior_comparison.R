# ============================================================
# Prior Comparison: Moment vs Inverse-Moment
#
# Group A: Moment priors (your current), d = 0, 1, 2, 3
#   π_M(δ) ∝ δ^d × g_κ(δ),  κ calibrated so median = 2
#
# Group B: Inverse-moment priors, d = 1, 2, 3, 4
#   π_I(δ) ∝ δ^{-(ν+1)} exp(-(τ/δ²)^d),  τ calibrated so median = 2
#
# Usage:
#   source("src/moment_prior.R")  # adjust path
#   source("prior_comparison.R")
#
# ============================================================
# USER SETTINGS
# ============================================================

NU_IMOMENT <- 2   # tail parameter for IM priors (try 1, 3, 5, 7, ...)
                   # ν=1 gives Cauchy-like tails (very heavy)
                   # ν=5 gives tails similar to half-t(6)
                   # larger ν = lighter tails

# ============================================================

library(ggplot2)

# ============================================================
# 1. Inverse-moment prior on δ > 0
# ============================================================

# J&R eq (10) adapted to the half-line δ > 0:
#
#   π_I(δ) = C × δ^{-(ν+1)} × exp(-(τ/δ²)^d)
#
# where C is the normalising constant (computed numerically).

#' Unnormalised inverse-moment density on δ > 0
#' @param delta  Positive values
#' @param d      Order parameter (1, 2, 3, ...)
#' @param tau    Scale parameter (> 0)
#' @param nu     Tail parameter
imoment_unnorm <- function(delta, d, tau, nu = NU_IMOMENT) {
  delta <- pmax(delta, 1e-300)
  log_dens <- -(nu + 1) * log(delta) - (tau / delta^2)^d
  exp(log_dens)
}

#' Normalising constant for inverse-moment prior on (0, ∞)
imoment_norm_const <- function(d, tau, nu = NU_IMOMENT) {
  integrate(function(x) imoment_unnorm(x, d, tau, nu),
            lower = 1e-10, upper = Inf,
            rel.tol = 1e-10, subdivisions = 1000)$value
}

#' Normalised inverse-moment density
imoment_density <- function(delta, d, tau, nu = NU_IMOMENT, C_norm = NULL) {
  if (is.null(C_norm)) C_norm <- imoment_norm_const(d, tau, nu)
  imoment_unnorm(delta, d, tau, nu) / C_norm
}

#' CDF of inverse-moment prior: Pr(δ < eps)
imoment_cdf <- function(eps, d, tau, nu = NU_IMOMENT, C_norm = NULL) {
  if (eps <= 0) return(0)
  if (is.null(C_norm)) C_norm <- imoment_norm_const(d, tau, nu)
  integrate(function(x) imoment_unnorm(x, d, tau, nu) / C_norm,
            lower = 1e-10, upper = eps,
            rel.tol = 1e-10, subdivisions = 1000)$value
}

#' Quantile of inverse-moment prior by numerical inversion
imoment_quantile <- function(p, d, tau, nu = NU_IMOMENT) {
  C_norm <- imoment_norm_const(d, tau, nu)
  upper <- 1
  while (imoment_cdf(upper, d, tau, nu, C_norm) < p) upper <- upper * 2
  uniroot(function(q) imoment_cdf(q, d, tau, nu, C_norm) - p,
          lower = 1e-10, upper = upper, tol = 1e-10)$root
}

#' Find tau so that the inverse-moment prior has median = target
find_tau_for_median <- function(d, target_median = 2, nu = NU_IMOMENT,
                                tau_range = c(1e-6, 1e4)) {
  objective <- function(log_tau) {
    tau <- exp(log_tau)
    med <- imoment_quantile(0.5, d, tau, nu)
    med - target_median
  }
  log_tau <- uniroot(objective,
                     interval = log(tau_range),
                     tol = 1e-8)$root
  exp(log_tau)
}


# ============================================================
# 2. Calibrate all priors
# ============================================================

cat("Calibrating moment priors (d = 0, 1, 2, 3)...\n")
moment_d <- 0:3
moment_kappa <- sapply(moment_d, get_kappa)
for (i in seq_along(moment_d)) {
  cat(sprintf("  Moment d=%d: kappa = %.4f\n", moment_d[i], moment_kappa[i]))
}

cat(sprintf("\nCalibrating inverse-moment priors (d = 1, 2, 3, 4, nu = %d, median = 2)...\n",
            NU_IMOMENT))
imoment_d <- 1:4
imoment_tau <- numeric(4)
for (i in seq_along(imoment_d)) {
  cat(sprintf("  Finding tau for d=%d ... ", imoment_d[i]))
  imoment_tau[i] <- find_tau_for_median(imoment_d[i], target_median = 2)
  med_check <- imoment_quantile(0.5, imoment_d[i], imoment_tau[i])
  cat(sprintf("tau = %.4f, median = %.4f\n", imoment_tau[i], med_check))
}


# ============================================================
# 3. Diagnostics table
# ============================================================

eps_vals <- c(0.2, 0.5, 1.0, 2.0)

# --- Moment priors ---
cat("\n============================================================\n")
cat("Diagnostics: Moment Priors\n")
cat("============================================================\n")

moment_diag <- data.frame(
  label  = paste0("M(d=", moment_d, ")"),
  d      = moment_d,
  kappa  = round(moment_kappa, 4)
)

for (i in seq_along(moment_d)) {
  d_i <- moment_d[i]
  k_i <- moment_kappa[i]
  pdf_fun <- function(x) prior_moment_density(x, d = d_i, kappa = k_i)
  moment_diag$median[i] <- round(prior_quantile(pdf_fun, 0.5, upper_start = k_i * PRIOR_SCALE), 3)
  moment_diag$q05[i]    <- round(prior_quantile(pdf_fun, 0.05, upper_start = k_i * PRIOR_SCALE), 3)
  moment_diag$q95[i]    <- round(prior_quantile(pdf_fun, 0.95, upper_start = k_i * PRIOR_SCALE), 3)
  for (e in eps_vals) {
    col <- paste0("Pr(d<", e, ")")
    moment_diag[[col]][i] <- round(
      integrate(pdf_fun, 0, e, rel.tol = 1e-8)$value, 4
    )
  }
}
print(moment_diag, row.names = FALSE)

# --- Inverse-moment priors ---
cat("\n============================================================\n")
cat(sprintf("Diagnostics: Inverse-Moment Priors (nu=%d)\n", NU_IMOMENT))
cat("============================================================\n")

imoment_diag <- data.frame(
  label = paste0("IM(d=", imoment_d, ")"),
  d     = imoment_d,
  tau   = round(imoment_tau, 4)
)

for (i in seq_along(imoment_d)) {
  d_i <- imoment_d[i]
  t_i <- imoment_tau[i]
  imoment_diag$median[i] <- round(imoment_quantile(0.5, d_i, t_i), 3)
  imoment_diag$q05[i]    <- round(imoment_quantile(0.05, d_i, t_i), 3)
  imoment_diag$q95[i]    <- round(imoment_quantile(0.95, d_i, t_i), 3)
  for (e in eps_vals) {
    col <- paste0("Pr(d<", e, ")")
    C_norm <- imoment_norm_const(d_i, t_i)
    imoment_diag[[col]][i] <- round(imoment_cdf(e, d_i, t_i, C_norm = C_norm), 4)
  }
}
print(imoment_diag, row.names = FALSE)


# ============================================================
# 4. Plots — all 8 priors overlaid
# ============================================================

m_cols  <- c("black", "#E41A1C", "#FF7F00", "#A65628")
im_cols <- c("#377EB8", "#4DAF4A", "#984EA3", "#66C2A5")

all_labels <- c(paste0("M(d=", moment_d, ")"), paste0("IM(d=", imoment_d, ")"))
all_cols   <- setNames(c(m_cols, im_cols), all_labels)
all_lty    <- setNames(c(rep("solid", 4), rep("dashed", 4)), all_labels)

# --- Build density data ---
delta_grid <- seq(0.001, 8, length.out = 1000)
delta_zoom <- seq(0.001, 2, length.out = 500)
delta_cdf  <- seq(0.001, 4, length.out = 400)

build_dens <- function(grid) {
  df_list <- list()

  for (i in seq_along(moment_d)) {
    df_list[[i]] <- data.frame(
      delta   = grid,
      density = prior_moment_density(grid, d = moment_d[i], kappa = moment_kappa[i]),
      label   = paste0("M(d=", moment_d[i], ")"),
      family  = "Moment"
    )
  }

  for (i in seq_along(imoment_d)) {
    C_norm <- imoment_norm_const(imoment_d[i], imoment_tau[i])
    df_list[[4 + i]] <- data.frame(
      delta   = grid,
      density = imoment_density(grid, d = imoment_d[i], tau = imoment_tau[i],
                                C_norm = C_norm),
      label   = paste0("IM(d=", imoment_d[i], ")"),
      family  = "Inverse-moment"
    )
  }

  df <- do.call(rbind, df_list)
  df$label <- factor(df$label, levels = all_labels)
  df
}

dens_full <- build_dens(delta_grid)
dens_zoom <- build_dens(delta_zoom)


# ---- Plot 1: Full range ----
p1 <- ggplot(dens_full, aes(x = delta, y = density, colour = label, linetype = label)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = all_cols, name = "Prior") +
  scale_linetype_manual(values = all_lty, name = "Prior") +
  labs(title = bquote("Prior density on δ scale (full range), IM:" ~ nu == .(NU_IMOMENT)),
       x = expression(delta), y = expression(pi(delta))) +
  theme_minimal(base_size = 12)
print(p1)

# ---- Plot 2: Zoomed to (0, 2) ----
p2 <- ggplot(dens_zoom, aes(x = delta, y = density, colour = label, linetype = label)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = all_cols, name = "Prior") +
  scale_linetype_manual(values = all_lty, name = "Prior") +
  labs(title = bquote("Prior density on δ scale — zoomed to (0, 2), IM:" ~ nu == .(NU_IMOMENT)),
       x = expression(delta), y = expression(pi(delta))) +
  theme_minimal(base_size = 12)
print(p2)

# ---- Plot 3: CDF ----
cdf_list <- list()

for (i in seq_along(moment_d)) {
  cdf_vals <- vapply(delta_cdf, function(e) {
    integrate(function(x) prior_moment_density(x, d = moment_d[i], kappa = moment_kappa[i]),
              0, e, rel.tol = 1e-8)$value
  }, numeric(1))
  cdf_list[[i]] <- data.frame(
    delta = delta_cdf, cdf = cdf_vals,
    label = paste0("M(d=", moment_d[i], ")"), family = "Moment"
  )
}

for (i in seq_along(imoment_d)) {
  C_norm <- imoment_norm_const(imoment_d[i], imoment_tau[i])
  cdf_vals <- vapply(delta_cdf, function(e) {
    imoment_cdf(e, imoment_d[i], imoment_tau[i], C_norm = C_norm)
  }, numeric(1))
  cdf_list[[4 + i]] <- data.frame(
    delta = delta_cdf, cdf = cdf_vals,
    label = paste0("IM(d=", imoment_d[i], ")"), family = "Inverse-moment"
  )
}

cdf_df <- do.call(rbind, cdf_list)
cdf_df$label <- factor(cdf_df$label, levels = all_labels)

p3 <- ggplot(cdf_df, aes(x = delta, y = cdf, colour = label, linetype = label)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  scale_colour_manual(values = all_cols, name = "Prior") +
  scale_linetype_manual(values = all_lty, name = "Prior") +
  labs(title = bquote("Prior CDF on δ scale, IM:" ~ nu == .(NU_IMOMENT)),
       x = expression(delta), y = expression(Pr(delta < epsilon))) +
  theme_minimal(base_size = 12)
print(p3)


# ---- Plot 4: Density on δ² scale ----
v_grid <- seq(0.005, 8, length.out = 500)

dens_sq_list <- list()

for (i in seq_along(moment_d)) {
  dv <- sqrt(v_grid)
  pi_d <- prior_moment_density(pmax(dv, 1e-12), d = moment_d[i], kappa = moment_kappa[i])
  dens_sq_list[[i]] <- data.frame(
    v = v_grid, density = pi_d / (2 * dv),
    label = paste0("M(d=", moment_d[i], ")"), family = "Moment"
  )
}

for (i in seq_along(imoment_d)) {
  C_norm <- imoment_norm_const(imoment_d[i], imoment_tau[i])
  dv <- sqrt(v_grid)
  pi_d <- imoment_density(pmax(dv, 1e-12), d = imoment_d[i], tau = imoment_tau[i],
                           C_norm = C_norm)
  dens_sq_list[[4 + i]] <- data.frame(
    v = v_grid, density = pi_d / (2 * dv),
    label = paste0("IM(d=", imoment_d[i], ")"), family = "Inverse-moment"
  )
}

dens_sq_df <- do.call(rbind, dens_sq_list)
dens_sq_df$label <- factor(dens_sq_df$label, levels = all_labels)

y_cap <- quantile(dens_sq_df$density[dens_sq_df$v > 0.1], 0.97, na.rm = TRUE)

p4 <- ggplot(dens_sq_df, aes(x = v, y = density, colour = label, linetype = label)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey50") +
  coord_cartesian(ylim = c(0, y_cap * 1.1)) +
  scale_colour_manual(values = all_cols, name = "Prior") +
  scale_linetype_manual(values = all_lty, name = "Prior") +
  labs(title = bquote("Induced prior density on " * delta^2 * " scale, IM:" ~ nu == .(NU_IMOMENT)),
       x = expression(delta^2), y = expression(pi(delta^2))) +
  theme_minimal(base_size = 12)
print(p4)

cat("\nDone. All plots displayed in RStudio.\n")
