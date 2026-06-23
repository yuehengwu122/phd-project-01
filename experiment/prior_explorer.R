# ============================================================
# Prior Explorer: Control Mass in (0, 1)
#
# Specify 4 priors via (d, target_mass) where
# target_mass = Pr(δ < 1). The script finds κ for each,
# prints a diagnostics table, and overlays all 4 on one plot.
#
# Usage:
#   source("moment_prior.R")  # adjust path as needed
#   source("prior_explorer.R")
#
# Then call:
#   explore_priors(priors)
# where priors is a data.frame with columns: label, d, target_mass
# ============================================================

library(ggplot2)

# ---- Core functions ----

#' Compute Pr(δ < eps) for given d and kappa (not using the cache)
prior_prob_direct <- function(d, kappa, eps) {
  if (eps <= 0) return(0)
  pdf <- function(x) prior_moment_density(x, d = d, kappa = kappa)
  integrate(pdf, lower = 0, upper = eps, rel.tol = 1e-8)$value
}

#' Find kappa such that Pr(δ < threshold | d, kappa) = target_mass
find_kappa_for_mass <- function(d, target_mass, threshold = 1,
                                kappa_range = c(0.01, 50)) {
  objective <- function(kappa) {
    prior_prob_direct(d, kappa, threshold) - target_mass
  }

  # Check that the interval brackets the root
  f_lo <- objective(kappa_range[1])
  f_hi <- objective(kappa_range[2])

  if (f_lo * f_hi > 0) {
    stop(sprintf(
      "Cannot bracket root for d=%d, target=%.3f.\n  f(%.2f)=%.4f, f(%.2f)=%.4f",
      d, target_mass, kappa_range[1], f_lo, kappa_range[2], f_hi
    ))
  }

  uniroot(objective, interval = kappa_range, tol = 1e-8)$root
}

#' Compute median of the prior given d and kappa
prior_median <- function(d, kappa) {
  pdf_fun <- function(x) prior_moment_density(x, d = d, kappa = kappa)
  prior_quantile(pdf_fun, p = 0.5, upper_start = kappa * PRIOR_SCALE)
}

#' Compute arbitrary quantile
prior_q <- function(d, kappa, p) {
  pdf_fun <- function(x) prior_moment_density(x, d = d, kappa = kappa)
  prior_quantile(pdf_fun, p = p, upper_start = kappa * PRIOR_SCALE)
}


# ---- Main function ----

#' Explore a set of priors
#'
#' @param priors data.frame with columns: label, d, target_mass
#'   label       — display name (e.g., "d=0, 27%")
#'   d           — moment order (integer >= 0)
#'   target_mass — target Pr(δ < 1)
#' @param threshold  the upper bound for the target mass (default 1)
#' @param colours    optional colour vector (length = nrow(priors))
explore_priors <- function(
    priors,
    threshold = 1,
    colours = NULL
) {

  n_priors <- nrow(priors)
  if (is.null(colours)) {
    colours <- c("black", "#E41A1C", "#377EB8", "#4DAF4A",
                 "#984EA3", "#FF7F00", "#A65628", "#F781BF")[1:n_priors]
  }

  # ---- Step 1: Find kappa for each prior ----
  cat("Finding kappa values...\n")
  priors$kappa <- NA_real_
  for (i in seq_len(n_priors)) {
    priors$kappa[i] <- find_kappa_for_mass(
      priors$d[i], priors$target_mass[i], threshold = threshold
    )
  }

  # ---- Step 2: Diagnostics table ----
  eps_vals <- c(0.2, 0.5, 1.0, 2.0)

  diag <- data.frame(
    label  = priors$label,
    d      = priors$d,
    kappa  = round(priors$kappa, 4),
    median = NA_real_,
    q05    = NA_real_,
    q95    = NA_real_
  )

  for (e in eps_vals) {
    diag[[paste0("Pr(d<", e, ")")]] <- NA_real_
  }

  for (i in seq_len(n_priors)) {
    d_i <- priors$d[i]
    k_i <- priors$kappa[i]

    diag$median[i] <- round(prior_median(d_i, k_i), 3)
    diag$q05[i]    <- round(prior_q(d_i, k_i, 0.05), 3)
    diag$q95[i]    <- round(prior_q(d_i, k_i, 0.95), 3)

    for (e in eps_vals) {
      col <- paste0("Pr(d<", e, ")")
      diag[[col]][i] <- round(prior_prob_direct(d_i, k_i, e), 4)
    }
  }

  cat("\n============================================================\n")
  cat("Diagnostics Table\n")
  cat("============================================================\n")
  print(diag, row.names = FALSE)
  cat("\n")

  # ---- Step 3: Build density and CDF data ----
  delta_grid_full <- seq(0.001, 25, length.out = 1000)
  delta_grid_zoom <- seq(0.001, 2, length.out = 500)
  delta_grid_cdf <- seq(0.001, 25, length.out = 500)
  delta_sq_grid <- seq(0.005, 25, length.out = 500)

  labels_factor <- factor(priors$label, levels = priors$label)

  # Density — full range
  dens_full <- do.call(rbind, lapply(seq_len(n_priors), function(i) {
    data.frame(
      delta   = delta_grid_full,
      density = prior_moment_density(delta_grid_full, d = priors$d[i],
                                     kappa = priors$kappa[i]),
      label   = factor(priors$label[i], levels = levels(labels_factor))
    )
  }))

  # Density — zoomed
  dens_zoom <- do.call(rbind, lapply(seq_len(n_priors), function(i) {
    data.frame(
      delta   = delta_grid_zoom,
      density = prior_moment_density(delta_grid_zoom, d = priors$d[i],
                                     kappa = priors$kappa[i]),
      label   = factor(priors$label[i], levels = levels(labels_factor))
    )
  }))

  # CDF
  cdf_df <- do.call(rbind, lapply(seq_len(n_priors), function(i) {
    cdf_vals <- vapply(delta_grid_cdf, function(e) {
      prior_prob_direct(priors$d[i], priors$kappa[i], e)
    }, numeric(1))
    data.frame(
      delta = delta_grid_cdf,
      cdf   = cdf_vals,
      label = factor(priors$label[i], levels = levels(labels_factor))
    )
  }))

  # Density on δ² scale
  dens_sq <- do.call(rbind, lapply(seq_len(n_priors), function(i) {
    dv <- sqrt(delta_sq_grid)
    pi_d <- prior_moment_density(pmax(dv, 1e-12), d = priors$d[i],
                                 kappa = priors$kappa[i])
    data.frame(
      delta_sq = delta_sq_grid,
      density  = pi_d / (2 * dv),
      label    = factor(priors$label[i], levels = levels(labels_factor))
    )
  }))

  col_scale <- scale_colour_manual(values = setNames(colours, levels(labels_factor)),
                                   name = "Prior")

  # ---- Plot 1: Density, full range ----
  p1 <- ggplot(dens_full, aes(x = delta, y = density, colour = label)) +
    geom_line(linewidth = 0.9) +
    col_scale +
    labs(title = "Prior density on δ scale (full range)",
         x = expression(delta), y = expression(pi(delta))) +
    theme_minimal(base_size = 12)
  print(p1)

  # ---- Plot 2: Density, zoomed to (0, 2) ----
  p2 <- ggplot(dens_zoom, aes(x = delta, y = density, colour = label)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    col_scale +
    labs(title = "Prior density on δ scale — zoomed to (0, 2)",
         x = expression(delta), y = expression(pi(delta))) +
    theme_minimal(base_size = 12)
  print(p2)

  # ---- Plot 3: CDF, zoomed to (0, 4) ----
  p3 <- ggplot(cdf_df, aes(x = delta, y = cdf, colour = label)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    col_scale +
    labs(title = "Prior CDF on δ scale",
         x = expression(delta), y = expression(Pr(delta < epsilon))) +
    theme_minimal(base_size = 12)
  print(p3)

  # ---- Plot 4: Density on δ² scale ----
  y_cap <- quantile(dens_sq$density[dens_sq$delta_sq > 0.1], 0.98, na.rm = TRUE)

  p4 <- ggplot(dens_sq, aes(x = delta_sq, y = density, colour = label)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    coord_cartesian(ylim = c(0, y_cap * 1.1)) +
    col_scale +
    labs(title = expression("Induced prior density on " * delta^2 * " scale"),
         x = expression(delta^2), y = expression(pi(delta^2))) +
    theme_minimal(base_size = 12)
  print(p4)

  # ---- Return results invisibly ----
  invisible(list(priors = priors, diagnostics = diag,
                 plots = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)))
}


# ============================================================
# Example usage (uncomment and modify):
# ============================================================

source("src/moment_prior.R")  # adjust path

priors <- data.frame(
  label       = c("d=0, 27%", "d=1, 18%", "d=0, 8%", "d=1, 5%"),
  d           = c(0,           1,           0,          1),
  target_mass = c(0.27,        0.18,        0.08,       0.05)
)

result <- explore_priors(priors)
