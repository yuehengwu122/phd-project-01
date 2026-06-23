# ============================================================
# Prior Mass Analysis: δ and δ² Scales
#
# Computes prior probabilities near the null and visualises
# prior densities / CDFs under the moment-prior family
# (d = 0, 1, 2, 3) on both the δ and δ² scales.
#
# Reference: Johnson & Rossell (2010), Section 3.4 and Fig. 2
#
# Outputs:
#   Table 1 — Pr(δ < ε) with J&R benchmark comparison
#   Table 2 — Ratio Pr(δ < ε | d) / Pr(δ < ε | d=0)
#   Plot 1 — Prior CDF on δ scale, zoomed to (0, 2)
#   Plot 2 — Prior density on δ scale, zoomed to (0, 1)
#   Plot 3 — Induced prior density on δ² scale
#   Plot 4 — Prior CDF on δ² scale, zoomed to (0, 4)
# ============================================================

source("src/_setup.R")
library(ggplot2)

# ---- Setup ----
d_orders <- 0:3
d_labels <- c("d=0 (local)", "d=1", "d=2", "d=3")
d_colours <- c("black", "#E41A1C", "#377EB8", "#4DAF4A")
names(d_colours) <- d_labels

kappas <- sapply(d_orders, get_kappa)
cat("Kappa values (median-matching scale):\n")
for (i in seq_along(d_orders)) {
  cat(sprintf("  d = %d : kappa = %.6f\n", d_orders[i], kappas[i]))
}
cat("\n")

# ============================================================
# TABLE 1: Pr(δ < ε) for discrete ε values
# ============================================================

eps_table <- c(0.1, 0.2, 0.5, 1.0)

prob_mat <- matrix(NA, nrow = 4, ncol = length(eps_table),
                   dimnames = list(d_labels, paste0("eps=", eps_table)))

for (i in seq_along(d_orders)) {
  for (j in seq_along(eps_table)) {
    prob_mat[i, j] <- prior_prob_delta(d_orders[i], eps_table[j])
  }
}

cat("============================================================\n")
cat("Table 1: Pr(δ < ε) under each prior\n")
cat("============================================================\n")
print(round(prob_mat, 6))
cat("\n")

# ---- J&R benchmark ----
cat("J&R benchmark (Section 3.4):\n")
cat("  They recommend Pr(|θ| < 0.2) ≈ 0.05 or 0.01 for the alternative prior.\n")
cat("  Our values at ε = 0.2 (interval-null threshold):\n")
for (i in seq_along(d_orders)) {
  cat(sprintf("    d=%d: Pr(δ < 0.2) = %.6f\n", d_orders[i], prob_mat[i, "eps=0.2"]))
}
cat("\n")

cat("  Our values at ε = 1.0 (minimum meaningful effect):\n")
for (i in seq_along(d_orders)) {
  p <- prob_mat[i, "eps=1"]
  cat(sprintf("    d=%d: Pr(δ < 1.0) = %.4f  (%.1f%%)\n", d_orders[i], p, 100 * p))
}
cat("\n")

# ---- Ratio table ----
ratio_mat <- sweep(prob_mat, 2, prob_mat[1, ], "/")

cat("============================================================\n")
cat("Table 2: Ratio Pr(δ < ε | d) / Pr(δ < ε | d=0)\n")
cat("============================================================\n")
print(round(ratio_mat, 4))
cat("\n")

# ---- Mass in (0.5, 1.0) ----
cat("Prior mass in (0.5, 1.0) — region where priors might still differ:\n")
for (i in seq_along(d_orders)) {
  diff <- prob_mat[i, "eps=1"] - prob_mat[i, "eps=0.5"]
  cat(sprintf("  d=%d: %.4f  (%.1f%%)\n", d_orders[i], diff, 100 * diff))
}
cat("\n")


# ============================================================
# Helper: induced density on δ² scale
#
#   If δ ~ π_d(δ), then v = δ² has density
#   π_d(√v) / (2√v)
# ============================================================

prior_density_delta_sq <- function(v, d, kappa) {
  v <- pmax(v, 1e-12)
  delta_vals <- sqrt(v)
  prior_moment_density(delta_vals, d = d, kappa = kappa) / (2 * delta_vals)
}


# ============================================================
# PLOT 1: Prior CDF on δ scale, zoomed to (0, 2)
# ============================================================

eps_grid <- seq(0.001, 2, length.out = 500)

cdf_df <- do.call(rbind, lapply(d_orders, function(d) {
  cdf_vals <- vapply(eps_grid, function(e) prior_prob_delta(d, e), numeric(1))
  data.frame(eps = eps_grid, cdf = cdf_vals,
             prior = factor(d_labels[d + 1], levels = d_labels))
}))

p1 <- ggplot(cdf_df, aes(x = eps, y = cdf, colour = prior)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.2, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 1.0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.23, y = 0.93, label = "ε = 0.2)",
           hjust = 0, size = 3, colour = "grey40") +
  annotate("text", x = 1.03, y = 0.93, label = "δ = 1",
           hjust = 0, size = 3, colour = "grey40") +
  scale_colour_manual(values = d_colours, name = "Prior order") +
  labs(title = "Prior CDF on δ scale — zoomed to (0, 2)",
       subtitle = "Pr(δ < ε) under moment priors d = 0, 1, 2, 3",
       x = expression(epsilon), y = expression(Pr(delta < epsilon))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.82, 0.85))

ggsave("plot1_prior_cdf_delta.png", p1, width = 7, height = 5, bg = "white")
cat("Saved: plot1_prior_cdf_delta.png\n")


# ============================================================
# PLOT 2: Prior density on δ scale, zoomed to (0, 1)
# ============================================================

delta_zoom <- seq(0.001, 1, length.out = 500)

dens_df <- do.call(rbind, lapply(d_orders, function(d) {
  kap <- get_kappa(d)
  data.frame(delta = delta_zoom,
             density = prior_moment_density(delta_zoom, d = d, kappa = kap),
             prior = factor(d_labels[d + 1], levels = d_labels))
}))

p2 <- ggplot(dens_df, aes(x = delta, y = density, colour = prior)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.2, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.22, y = max(dens_df$density) * 0.65,
           label = "ε = 0.2", hjust = 0, size = 3, colour = "grey40") +
  scale_colour_manual(values = d_colours, name = "Prior order") +
  labs(title = "Prior density on δ scale — zoomed to (0, 1)",
       subtitle = "Near-null region where non-local priors differ from local",
       x = expression(delta), y = expression(pi[d](delta))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.78, 0.35))

ggsave("plot2_prior_density_delta_zoom.png", p2, width = 7, height = 5, bg = "white")
cat("Saved: plot2_prior_density_delta_zoom.png\n")


# ============================================================
# PLOT 3: Induced prior density on δ² scale
# ============================================================

v_grid <- seq(0.005, 4, length.out = 500)

dens_sq_df <- do.call(rbind, lapply(d_orders, function(d) {
  kap <- get_kappa(d)
  data.frame(v = v_grid,
             density = prior_density_delta_sq(v_grid, d = d, kappa = kap),
             prior = factor(d_labels[d + 1], levels = d_labels))
}))

# Cap extreme densities near v=0 for d=0 (diverges as 1/√v)
y_cap <- quantile(dens_sq_df$density[dens_sq_df$v > 0.05], 0.99, na.rm = TRUE)

p3 <- ggplot(dens_sq_df, aes(x = v, y = density, colour = prior)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.04, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 1.0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.08, y = y_cap * 0.9,
           label = expression(epsilon^2 == 0.04), hjust = 0, size = 3, colour = "grey40") +
  annotate("text", x = 1.05, y = y_cap * 0.9,
           label = expression(delta^2 == 1), hjust = 0, size = 3, colour = "grey40") +
  scale_colour_manual(values = d_colours, name = "Prior order") +
  coord_cartesian(ylim = c(0, y_cap * 1.1)) +
  labs(title = expression("Induced prior density on " * delta^2 * " scale"),
       subtitle = expression(pi(v) == pi(sqrt(v)) / (2 * sqrt(v))),
       x = expression(delta^2), y = expression(pi(delta^2))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.78, 0.65))

ggsave("plot3_prior_density_delta_sq.png", p3, width = 7, height = 5, bg = "white")
cat("Saved: plot3_prior_density_delta_sq.png\n")


# ============================================================
# PLOT 4: Prior CDF on δ² scale, zoomed to (0, 4)
# ============================================================

v_cdf_grid <- seq(0.001, 4, length.out = 500)

cdf_sq_df <- do.call(rbind, lapply(d_orders, function(d) {
  # Pr(δ² < v) = Pr(δ < √v)
  cdf_vals <- vapply(v_cdf_grid, function(v) prior_prob_delta(d, sqrt(v)), numeric(1))
  data.frame(v = v_cdf_grid, cdf = cdf_vals,
             prior = factor(d_labels[d + 1], levels = d_labels))
}))

p4 <- ggplot(cdf_sq_df, aes(x = v, y = cdf, colour = prior)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.04, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.25, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 1.0, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.06, y = 0.93,
           label = expression(epsilon^2 == 0.04 ~ (epsilon == 0.2)),
           hjust = 0, size = 2.8, colour = "grey40") +
  annotate("text", x = 0.28, y = 0.80,
           label = expression(delta^2 == 0.25 ~ (delta == 0.5)),
           hjust = 0, size = 2.8, colour = "grey40") +
  annotate("text", x = 1.03, y = 0.93,
           label = expression(delta^2 == 1 ~ (delta == 1)),
           hjust = 0, size = 2.8, colour = "grey40") +
  scale_colour_manual(values = d_colours, name = "Prior order") +
  labs(title = expression("Prior CDF on " * delta^2 * " scale — zoomed to (0, 4)"),
       subtitle = expression(Pr(delta^2 < v) == Pr(delta < sqrt(v))),
       x = expression(delta^2), y = expression(Pr(delta^2 < v))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.82, 0.85))

ggsave("plot4_prior_cdf_delta_sq.png", p4, width = 7, height = 5, bg = "white")
cat("Saved: plot4_prior_cdf_delta_sq.png\n")


# ============================================================
# Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Summary for discussion\n")
cat("============================================================\n\n")

cat("1. Prior mass in interval-null region (δ < 0.2):\n")
for (i in seq_along(d_orders)) {
  p <- prob_mat[i, "eps=0.2"]
  cat(sprintf("   d=%d: %.6f  (%.4f%%)\n", d_orders[i], p, 100 * p))
}

cat("\n2. Prior mass below 'small effect' (δ < 1.0):\n")
for (i in seq_along(d_orders)) {
  p <- prob_mat[i, "eps=1"]
  cat(sprintf("   d=%d: %.4f  (%.1f%%)\n", d_orders[i], p, 100 * p))
}

cat("\n3. Prior mass in (0.5, 1.0) — beyond where priors visually separate:\n")
for (i in seq_along(d_orders)) {
  diff <- prob_mat[i, "eps=1"] - prob_mat[i, "eps=0.5"]
  cat(sprintf("   d=%d: %.4f  (%.1f%%)\n", d_orders[i], diff, 100 * diff))
}

cat("\nDone. Check png plots in working directory.\n")
