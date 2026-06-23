# ============================================================
# Plot Exploratory Simulation: BF vs N (5 methods, 6 priors)
#
# Methods: Bridge, Logspline, ECDF, Laplace, Jeffreys
# Priors:  M(d=0), M(d=1), M(d=2), M(d=3), IM(d=1), IM(d=2)
#
# Reads:
#   results/exploratory/bf_bridge.csv
#   results/exploratory/posterior_draws.csv
#
# Usage:
#   source("src/_setup.R")
#   source("analysis/exploratory/plot_exploratory_bf.R")
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

source("src/_setup.R")

# ============================================================
# SETTINGS
# ============================================================

EPS_FIXED  <- 0.2
data_dir   <- "results/exploratory"
report_dir <- "results/exploratory/report"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

N_vals     <- c(20, 50, 100, 200)
delta_vals <- c(0, 1, 2, 3)
eps_grid   <- seq(0.01, 1, length.out = 100)

DPI <- 200
W   <- 14
H   <- 5

save_png <- function(fname, plot, width = W, height = H) {
  path <- file.path(report_dir, fname)
  ggsave(path, plot = plot, width = width, height = height, dpi = DPI, bg = "white")
  message("  Saved: ", path)
}


# ============================================================
# INVERSE-MOMENT PRIOR FUNCTIONS (for prior probability)
# ============================================================

IM_NU <- 2

imoment_unnorm <- function(delta, d, tau, nu = IM_NU) {
  delta <- pmax(delta, 1e-300)
  exp(-(nu + 1) * log(delta) - (tau / delta^2)^d)
}

imoment_norm_const <- function(d, tau, nu = IM_NU) {
  integrate(function(x) imoment_unnorm(x, d, tau, nu),
            lower = 1e-10, upper = Inf,
            rel.tol = 1e-10, subdivisions = 1000)$value
}

imoment_cdf <- function(eps, d, tau, nu = IM_NU, C_norm = NULL) {
  if (eps <= 0) return(0)
  if (is.null(C_norm)) C_norm <- imoment_norm_const(d, tau, nu)
  integrate(function(x) imoment_unnorm(x, d, tau, nu) / C_norm,
            lower = 1e-10, upper = eps,
            rel.tol = 1e-10, subdivisions = 1000)$value
}

imoment_quantile <- function(p, d, tau, nu = IM_NU) {
  C_norm <- imoment_norm_const(d, tau, nu)
  upper <- 1
  while (imoment_cdf(upper, d, tau, nu, C_norm) < p) upper <- upper * 2
  uniroot(function(q) imoment_cdf(q, d, tau, nu, C_norm) - p,
          lower = 1e-10, upper = upper, tol = 1e-10)$root
}

find_tau_for_median <- function(d, target_median = 2, nu = IM_NU) {
  objective <- function(log_tau) {
    imoment_quantile(0.5, d, exp(log_tau), nu) - target_median
  }
  exp(uniroot(objective, interval = log(c(1e-6, 1e4)), tol = 1e-8)$root)
}

# Calibrate IM priors
message("Calibrating IM priors...")
im_tau <- c(
  "IM(d=1)" = find_tau_for_median(1),
  "IM(d=2)" = find_tau_for_median(2)
)
for (nm in names(im_tau)) {
  cat(sprintf("  %s: tau = %.4f\n", nm, im_tau[nm]))
}


# ============================================================
# LABELS & COLOURS
# ============================================================

prior_levels <- c("M(d=0)", "M(d=1)", "M(d=2)", "M(d=3)", "IM(d=1)", "IM(d=2)")
prior_colors <- c(
  "M(d=0)"  = "black",
  "M(d=1)"  = "#E41A1C",
  "M(d=2)"  = "#FF7F00",
  "M(d=3)"  = "#A65628",
  "IM(d=1)" = "#377EB8",
  "IM(d=2)" = "#4DAF4A"
)

add_labels <- function(df) {
  df |>
    mutate(
      delta_label = factor(paste0("δ = ", delta_true),
                           levels = paste0("δ = ", delta_vals)),
      N_label = factor(paste0("N = ", N), levels = paste0("N = ", N_vals)),
      prior = factor(prior, levels = prior_levels)
    )
}

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 8)
  )


# ============================================================
# LOAD DATA
# ============================================================

message("Loading bridge BF...")
bf_bridge <- read_csv(file.path(data_dir, "bf_bridge.csv"),
                      show_col_types = FALSE) |>
  add_labels()

message("Loading posterior draws...")
draws_all <- read_csv(file.path(data_dir, "posterior_draws.csv"),
                      show_col_types = FALSE)


# ============================================================
# COMPUTE PRIOR PROBABILITIES FOR ALL PRIORS
# ============================================================

message("Pre-computing prior probabilities...")

prior_table <- do.call(rbind, lapply(prior_levels, function(plab) {
  prior0 <- vapply(eps_grid, function(e) {
    if (grepl("^M", plab)) {
      d_ord <- as.integer(gsub("M\\(d=|\\)", "", plab))
      prior_prob_delta(d_ord, e)
    } else {
      d_ord <- as.integer(gsub("IM\\(d=|\\)", "", plab))
      imoment_cdf(e, d_ord, im_tau[plab])
    }
  }, numeric(1))
  data.frame(prior = plab, epsilon = eps_grid, prior0 = prior0,
             stringsAsFactors = FALSE)
}))


# ============================================================
# COMPUTE INTERVAL-NULL BFs FROM DRAWS
# ============================================================

message("Computing logspline & ECDF BFs from posterior draws...")
bf_interval <- draws_all |>
  mutate(prior = factor(prior, levels = prior_levels)) |>
  group_by(delta_true, N, prior) |>
  group_modify(function(df, keys) {
    x <- df$delta_draw[is.finite(df$delta_draw) & df$delta_draw >= 0]
    n_draws <- length(x)
    x_sort <- sort(x)

    p0 <- prior_table$prior0[prior_table$prior == as.character(keys$prior)]

    # ECDF
    k_vec <- findInterval(eps_grid, x_sort)
    post0_ecdf <- k_vec / n_draws

    # Logspline
    post0_ls <- vapply(eps_grid, function(e) {
      tryCatch(post_prob_delta_logspline(x, e),
               error = function(err) NA_real_)
    }, numeric(1))

    # Laplace & Jeffreys
    post0_laplace  <- (k_vec + 1) / (n_draws + 2)
    post0_jeffreys <- (k_vec + 0.5) / (n_draws + 1)

    data.frame(
      epsilon          = eps_grid,
      prior0           = p0,
      post0_logspline  = post0_ls,
      post0_ecdf       = post0_ecdf,
      post0_laplace    = post0_laplace,
      post0_jeffreys   = post0_jeffreys,
      log_bf10_logspline = log(p0 / post0_ls),
      log_bf10_ecdf      = log(p0 / post0_ecdf),
      log_bf10_laplace   = log(p0 / post0_laplace),
      log_bf10_jeffreys  = log(p0 / post0_jeffreys)
    )
  }) |>
  ungroup() |>
  add_labels()


# ============================================================
# EXTRACT BF AT FIXED EPSILON
# ============================================================

eps_sel <- eps_grid[which.min(abs(eps_grid - EPS_FIXED))]

bf_at_eps <- bf_interval |>
  filter(abs(epsilon - eps_sel) < 1e-9) |>
  select(delta_true, N, prior, delta_label, N_label,
         log_bf10_logspline, log_bf10_ecdf,
         log_bf10_laplace, log_bf10_jeffreys)

# ============================================================
# PLOT HELPER
# ============================================================

make_bf_vs_N <- function(df, bf_col, title_method, flag_inf = FALSE) {
  df <- df |> mutate(lbf = .data[[bf_col]])

  p <- ggplot(df, aes(x = N, y = lbf, colour = prior)) +
    facet_wrap(~ delta_label, nrow = 1, scales = "free_y") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_colour_manual(values = prior_colors, name = "Prior") +
    scale_x_continuous(breaks = N_vals) +
    labs(title = sprintf("log BF₁₀ vs N  [%s]", title_method),
         x = "N", y = "log BF₁₀") +
    base_theme

  if (flag_inf) {
    inf_pts <- df |> filter(is.infinite(lbf) & lbf > 0)
    if (nrow(inf_pts) > 0) {
      p <- p + geom_point(data = inf_pts,
                          shape = 4, size = 3.5, stroke = 1.2,
                          show.legend = FALSE)
    }
  }
  p
}


# ============================================================
# GENERATE PLOTS
# ============================================================

message("Generating plots...")

p_bridge <- make_bf_vs_N(
  bf_bridge, "log_bf10_bridge", "Bridge Sampling"
)
save_png("bf_vs_N_bridge.png", p_bridge)

p_logspline <- make_bf_vs_N(
  bf_at_eps, "log_bf10_logspline",
  sprintf("Logspline, ε = %.2f", EPS_FIXED)
)
save_png("bf_vs_N_logspline.png", p_logspline)

p_ecdf <- make_bf_vs_N(
  bf_at_eps, "log_bf10_ecdf",
  sprintf("ECDF, ε = %.2f", EPS_FIXED), flag_inf = TRUE
)
save_png("bf_vs_N_ecdf.png", p_ecdf)

p_laplace <- make_bf_vs_N(
  bf_at_eps, "log_bf10_laplace",
  sprintf("Laplace, ε = %.2f", EPS_FIXED)
)
save_png("bf_vs_N_laplace.png", p_laplace)

p_jeffreys <- make_bf_vs_N(
  bf_at_eps, "log_bf10_jeffreys",
  sprintf("Jeffreys, ε = %.2f", EPS_FIXED)
)
save_png("bf_vs_N_jeffreys.png", p_jeffreys)


# ---- Combined 5-method plot ----
p_all <- (p_bridge + p_logspline + p_ecdf + p_laplace + p_jeffreys) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

save_png("bf_vs_N_all_methods.png", p_all, W, H * 5)

message("\n=== Done ===")
