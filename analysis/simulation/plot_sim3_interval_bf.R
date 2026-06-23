#' ============================================================
#' Simulation Report: Interval-Null BF with Non-local Priors
#' ============================================================
#' Generates all PNG figures and an HTML report.
#' Run: source("report_sim3.R")
#' Output: analysis/simulation/report/

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# ============================================================
# SETTINGS
# ============================================================

EPS_FIXED <- 0.2 # chosen operational epsilon
cv_threshold <- 0.2 # CV criterion for adaptive epsilon

data_dir <- "data/simulated"
report_dir <- "analysis/simulation/report"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

delta_vals <- 0:3
N_vals <- c(20, 50, 100, 200)

# Plot dimensions
W <- 14
H <- 10 # panel grid plots
W2 <- 13
H2 <- 5 # single-row plots
DPI <- 200

save_png <- function(fname, plot, width = W2, height = H2) {
  path <- file.path(report_dir, fname)
  ggsave(path, plot = plot, width = width, height = height, dpi = DPI, bg = "white")
  message("  Saved: ", path)
  invisible(path)
}

# ============================================================
# LOAD DATA
# ============================================================

message("Loading data...")
bf_int <- lapply(delta_vals, function(d) {
  readr::read_csv(file.path(data_dir, sprintf("results_sim3_rstan_bf_interval_d%d.csv", d)), show_col_types = FALSE)
}) |>
  bind_rows()

bf_bridge <- lapply(delta_vals, function(d) {
  readr::read_csv(file.path(data_dir, sprintf("results_sim3_rstan_bf_bridge_d%d.csv", d)), show_col_types = FALSE)
}) |>
  bind_rows()

# ============================================================
# LABELS & COLORS
# ============================================================

order_labels <- c(d0 = "d=0 (local)", d1 = "d=1", d2 = "d=2", d3 = "d=3")
order_colors <- c("d=0 (local)" = "#E69F00", "d=1" = "#56B4E9", "d=2" = "#009E73", "d=3" = "#CC79A7")

add_labels <- function(df) {
  df |>
    mutate(
      delta_label = factor(paste0("delta = ", delta_true), levels = paste0("delta = ", 0:3)),
      N_label = factor(paste0("N = ", N), levels = paste0("N = ", N_vals)),
      order = factor(order, levels = names(order_labels), labels = order_labels)
    )
}

bf_int <- add_labels(bf_int)
bf_bridge <- add_labels(bf_bridge)

bridge_med <- bf_bridge |>
  group_by(delta_label, N_label, order) |>
  summarise(median_bridge = median(log_bf10_bridge, na.rm = TRUE), .groups = "drop")

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 8)
  )

# ============================================================
# FIGURE 1: BF vs epsilon convergence (logspline + raw)
# ============================================================

message("Figure 1: BF convergence curves...")

make_fig1 <- function(method_sel, show_ribbon = FALSE) {
  df <- bf_int |>
    filter(method == method_sel) |>
    group_by(delta_label, N_label, order, epsilon) |>
    summarise(
      med_raw = median(log_bf10, na.rm = TRUE),
      q25 = {
        fv <- log_bf10[is.finite(log_bf10)]
        if (length(fv) == 0) NA_real_ else quantile(fv, 0.25)
      },
      q75 = {
        fv <- log_bf10[is.finite(log_bf10)]
        if (length(fv) == 0) NA_real_ else quantile(fv, 0.75)
      },
      .groups = "drop"
    ) |>
    mutate(is_inf = is.infinite(med_raw) & med_raw > 0, med = ifelse(is_inf, NA_real_, med_raw))

  df_onset <- df |> filter(!is_inf) |> group_by(delta_label, N_label, order) |> slice_min(epsilon, n = 1) |> ungroup()

  df_plot <- df |> filter(!is_inf)

  p <- ggplot(df_plot, aes(x = epsilon, y = med, color = order, fill = order)) +
    facet_grid(delta_label ~ N_label, scales = "free_y") +
    geom_vline(xintercept = c(0.1, 0.2), linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_hline(
      data = bridge_med,
      aes(yintercept = median_bridge, color = order),
      linetype = "dotted",
      linewidth = 0.6,
      show.legend = FALSE
    )

  if (show_ribbon) {
    p <- p + geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA)
  }

  p +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_point(data = df_onset, aes(x = epsilon, y = med, color = order), shape = 17, size = 2.5, show.legend = FALSE) +
    scale_color_manual(values = order_colors, name = "Prior order") +
    scale_fill_manual(values = order_colors, name = "Prior order") +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(
      title = sprintf("Interval-null log BF10 vs epsilon  [%s]", method_sel),
      subtitle = "Dotted = bridge sampling. Dashed = eps 0.1, 0.2. Triangle = first finite BF (left = Inf).",
      x = "epsilon",
      y = "Median log BF10"
    ) +
    base_theme
}

save_png("fig1a_bf_vs_eps_logspline.png", make_fig1("logspline"), W, H)
save_png("fig1b_bf_vs_eps_raw.png", make_fig1("raw", show_ribbon = TRUE), W, H)

# ============================================================
# FIGURE 2: Effective convergence epsilon diagnostic
# ============================================================

message("Figure 2: Convergence epsilon diagnostic...")

eps_grid <- sort(unique(bf_int$epsilon))

df_eff <- bf_int |>
  filter(method == "raw") |>
  group_by(delta_label, N_label, order, replicate) |>
  summarise(
    eff_eps = {
      pe <- epsilon[post0 > 0 & is.finite(post0)]
      if (length(pe) == 0) max(eps_grid) else min(pe)
    },
    .groups = "drop"
  )

p2 <- ggplot(df_eff, aes(x = delta_label, y = eff_eps, fill = order, color = order)) +
  facet_wrap(~N_label, nrow = 1) +
  geom_violin(position = position_dodge(0.85), alpha = 0.4, linewidth = 0.5, scale = "width") +
  geom_boxplot(position = position_dodge(0.85), width = 0.15, outlier.size = 0.5, linewidth = 0.5, alpha = 0.8) +
  geom_hline(yintercept = c(0.1, 0.2), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = order_colors, name = "Prior order") +
  scale_color_manual(values = order_colors, name = "Prior order") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = "Effective convergence epsilon (smallest eps with P(delta < eps | y) > 0)",
    subtitle = "Raw method. Dashed = eps 0.1 and 0.2.",
    x = "True delta",
    y = "Effective convergence epsilon"
  ) +
  base_theme

save_png("fig2_convergence_eps.png", p2)

# ============================================================
# FIGURE 3: BF at fixed epsilon = 0.2 (logspline + raw)
# ============================================================

message("Figure 3: BF at fixed epsilon = 0.2...")

eps_sel <- eps_grid[which.min(abs(eps_grid - EPS_FIXED))]

df_fix <- bf_int |>
  filter(abs(epsilon - eps_sel) < 1e-9) |>
  group_by(delta_label, N_label, order, method) |>
  summarise(
    med = median(log_bf10, na.rm = TRUE),
    q25 = quantile(log_bf10[is.finite(log_bf10)], 0.25, na.rm = TRUE),
    q75 = quantile(log_bf10[is.finite(log_bf10)], 0.75, na.rm = TRUE),
    n_inf = sum(is.infinite(log_bf10) & log_bf10 > 0),
    n_total = n(),
    is_inf = is.infinite(median(log_bf10, na.rm = TRUE)),
    .groups = "drop"
  )

make_fig3 <- function(method_sel) {
  df <- df_fix |> filter(method == method_sel)
  finite_vals <- c(df$med[is.finite(df$med)], bridge_med$median_bridge[is.finite(bridge_med$median_bridge)])
  y_ceil <- max(c(finite_vals, 5), na.rm = TRUE) * 1.15

  df <- df |>
    mutate(
      med_plot = ifelse(is_inf, y_ceil, med),
      inf_label = ifelse(is_inf, sprintf("Inf\n(%d/%d fin)", n_total - n_inf, n_total), "")
    )

  ggplot(df, aes(x = delta_label, y = med_plot, fill = order)) +
    facet_wrap(~N_label, nrow = 1) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
    geom_errorbar(
      data = df |> filter(!is_inf),
      aes(ymin = q25, ymax = q75),
      position = position_dodge(0.8),
      width = 0.3,
      linewidth = 0.4
    ) +
    geom_text(
      data = df |> filter(is_inf),
      aes(label = inf_label),
      position = position_dodge(0.8),
      vjust = -0.1,
      size = 2,
      color = "grey30",
      lineheight = 0.8
    ) +
    geom_point(
      data = bridge_med,
      aes(x = delta_label, y = median_bridge, group = order),
      position = position_dodge(0.8),
      shape = 4,
      size = 2.5,
      stroke = 1.2,
      color = "black",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = order_colors, name = "Prior order") +
    coord_cartesian(ylim = c(NA, y_ceil * 1.2)) +
    labs(
      title = sprintf("Log BF10 at eps = %.2f  [%s]", EPS_FIXED, method_sel),
      subtitle = "Bars = median; errors = IQR (finite); x = bridge sampling; Inf labeled.",
      x = "True delta",
      y = "Median log BF10"
    ) +
    base_theme
}

save_png("fig3a_bf_eps02_logspline.png", make_fig3("logspline"))
save_png("fig3b_bf_eps02_raw.png", make_fig3("raw"))

# ============================================================
# FIGURE 4: Evidence accumulation — combined panels
# Bridge / Logspline / Raw, each with BF-vs-N (top) + BF-vs-delta (bottom)
# ============================================================

message("Figure 4: Evidence accumulation combined panels...")

# Helper: BF vs N
make_bf_vs_N <- function(df_summary, y_col, title_method) {
  ggplot(df_summary, aes(x = N, y = .data[[y_col]], color = order, fill = order)) +
    facet_wrap(~delta_label, nrow = 1) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = order_colors, name = "Prior order") +
    scale_fill_manual(values = order_colors, name = "Prior order") +
    scale_x_continuous(breaks = N_vals) +
    labs(title = sprintf("log BF10 vs N  [%s]", title_method), x = "N", y = "Median log BF10") +
    base_theme
}

# Helper: BF vs delta
make_bf_vs_delta <- function(df_summary, y_col, title_method) {
  ggplot(df_summary, aes(x = delta_true, y = .data[[y_col]], color = order, fill = order)) +
    facet_wrap(~N_label, nrow = 1) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = order_colors, name = "Prior order") +
    scale_fill_manual(values = order_colors, name = "Prior order") +
    scale_x_continuous(breaks = 0:3) +
    labs(title = sprintf("log BF10 vs delta  [%s]", title_method), x = "True delta", y = "Median log BF10") +
    base_theme
}

# Bridge
df_bs <- bf_bridge |>
  group_by(delta_label, N_label, delta_true, N, order) |>
  summarise(
    med = median(log_bf10_bridge, na.rm = TRUE),
    q25 = quantile(log_bf10_bridge, 0.25, na.rm = TRUE),
    q75 = quantile(log_bf10_bridge, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Interval at eps = 0.2
make_interval_summary <- function(method_sel) {
  bf_int |>
    filter(method == method_sel, abs(epsilon - eps_sel) < 1e-9) |>
    mutate(log_bf_finite = ifelse(is.finite(log_bf10), log_bf10, NA_real_)) |>
    group_by(delta_label, N_label, delta_true, N, order) |>
    summarise(
      med = median(log_bf_finite, na.rm = TRUE),
      q25 = quantile(log_bf_finite, 0.25, na.rm = TRUE),
      q75 = quantile(log_bf_finite, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

df_ls <- make_interval_summary("logspline")
df_rw <- make_interval_summary("raw")

# Combined: 2 rows x 1 col per method
library(patchwork)

for (info in list(
  list(df = df_bs, col = "med", lab = "Bridge Sampling", tag = "bridge"),
  list(df = df_ls, col = "med", lab = "Interval (logspline, eps=0.2)", tag = "logspline"),
  list(df = df_rw, col = "med", lab = "Interval (raw, eps=0.2)", tag = "raw")
)) {
  p_top <- make_bf_vs_N(info$df, info$col, info$lab)
  p_bot <- make_bf_vs_delta(info$df, info$col, info$lab)
  p_combined <- p_top / p_bot + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_png(sprintf("fig4_%s_combined.png", info$tag), p_combined, W2, H2 * 2)
}

# ============================================================
# Bridge sampling BF table (prior order differences)
# ============================================================

message("Bridge BF table...")

bridge_table <- bf_bridge |>
  group_by(delta_true, N, order) |>
  summarise(median_logBF = sprintf("%.3f", median(log_bf10_bridge, na.rm = TRUE)), .groups = "drop") |>
  pivot_wider(names_from = order, values_from = median_logBF) |>
  arrange(delta_true, N)

readr::write_csv(bridge_table, file.path(report_dir, "bridge_bf_table.csv"))

# ============================================================
# FIGURE 5: Method comparison scatter (logspline vs raw, eps=0.2)
# ============================================================

message("Figure 5: Method comparison scatter...")

df_scatter_base <- df_fix |>
  select(delta_label, N_label, order, method, med) |>
  mutate(med_finite = ifelse(is.finite(med), med, NA_real_)) |>
  select(-med) |>
  pivot_wider(names_from = method, values_from = med_finite) |>
  filter(!is.na(logspline) & !is.na(raw))

lim <- range(c(df_scatter_base$logspline, df_scatter_base$raw), na.rm = TRUE)
lim <- c(floor(lim[1]) - 0.5, ceiling(lim[2]) + 0.5)

p5 <- ggplot(df_scatter_base, aes(x = raw, y = logspline, color = delta_label, shape = order)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d(name = "True delta", option = "C", end = 0.85) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Prior order") +
  coord_fixed(xlim = lim, ylim = lim) +
  labs(
    title = "Logspline vs Raw at eps = 0.2 (finite medians only)",
    subtitle = "Dashed = identity. Points with Inf in either method excluded.",
    x = "Median log BF10 (raw)",
    y = "Median log BF10 (logspline)"
  ) +
  base_theme +
  theme(legend.position = "right")

save_png("fig5_method_comparison.png", p5, 7, 6)

# ============================================================
# FIGURE 6: CV-optimal epsilon (Plots 7a, 7b, 7c)
# ============================================================

message("Figure 6: CV-optimal epsilon analysis...")

df_cv <- bf_int |>
  filter(method == "raw") |>
  mutate(cv = ifelse(post0 > 0 & is.finite(post0) & is.finite(post0_mcse), post0_mcse / post0, Inf))

df_cv_opt <- df_cv |>
  filter(cv < cv_threshold) |>
  group_by(delta_label, N_label, delta_true, N, order, replicate) |>
  slice_min(epsilon, n = 1) |>
  ungroup() |>
  rename(eps_opt = epsilon, log_bf10_opt = log_bf10, cv_opt = cv)

# 6a: distribution of optimal epsilon
p6a <- ggplot(df_cv_opt, aes(x = delta_label, y = eps_opt, fill = order, color = order)) +
  facet_wrap(~N_label, nrow = 1) +
  geom_violin(position = position_dodge(0.85), alpha = 0.4, linewidth = 0.5, scale = "width") +
  geom_boxplot(position = position_dodge(0.85), width = 0.15, outlier.size = 0.5, linewidth = 0.5, alpha = 0.8) +
  geom_hline(yintercept = c(0.1, 0.2), linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = order_colors, name = "Prior order") +
  scale_color_manual(values = order_colors, name = "Prior order") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title = sprintf("CV-optimal epsilon (CV < %.0f%%)", cv_threshold * 100),
    subtitle = "Smallest eps where SE/p_hat < threshold.",
    x = "True delta",
    y = "CV-optimal epsilon"
  ) +
  base_theme

save_png("fig6a_cv_optimal_eps.png", p6a)

# 6b: BF at CV-optimal epsilon
df_cv_summary <- df_cv_opt |>
  group_by(delta_label, N_label, order) |>
  summarise(
    med_bf = median(log_bf10_opt, na.rm = TRUE),
    q25 = quantile(log_bf10_opt, 0.25, na.rm = TRUE),
    q75 = quantile(log_bf10_opt, 0.75, na.rm = TRUE),
    med_eps = median(eps_opt),
    n_valid = n(),
    .groups = "drop"
  )

p6b <- ggplot(df_cv_summary, aes(x = delta_label, y = med_bf, fill = order)) +
  facet_wrap(~N_label, nrow = 1) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  geom_errorbar(aes(ymin = q25, ymax = q75), position = position_dodge(0.8), width = 0.3, linewidth = 0.4) +
  geom_point(
    data = bridge_med,
    aes(x = delta_label, y = median_bridge, group = order),
    position = position_dodge(0.8),
    shape = 4,
    size = 2.5,
    stroke = 1.2,
    color = "black",
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = order_colors, name = "Prior order") +
  labs(
    title = sprintf("Interval BF at CV-optimal eps (CV < %.0f%%) [raw]", cv_threshold * 100),
    subtitle = "Bars = median; errors = IQR; x = bridge sampling.",
    x = "True delta",
    y = "Median log BF10"
  ) +
  base_theme

save_png("fig6b_bf_cv_optimal.png", p6b)

# 6c: scatter vs bridge (per replicate)
df_cv_scatter <- df_cv_opt |>
  left_join(
    bf_bridge |> select(delta_true, N, order, replicate, log_bf10_bridge),
    by = c("delta_true", "N", "order", "replicate")
  ) |>
  filter(is.finite(log_bf10_opt) & is.finite(log_bf10_bridge))

p6c <- ggplot(df_cv_scatter, aes(x = log_bf10_bridge, y = log_bf10_opt, color = delta_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.3, size = 1.2) +
  facet_wrap(~N_label, nrow = 1) +
  scale_color_viridis_d(name = "True delta", option = "C", end = 0.85) +
  labs(
    title = "Interval BF (CV-optimal eps) vs Bridge Sampling [per replicate]",
    subtitle = "Dashed = identity line.",
    x = "log BF10 (bridge)",
    y = "log BF10 (interval, CV-optimal eps)"
  ) +
  base_theme

save_png("fig6c_scatter_cv_vs_bridge.png", p6c)

# ============================================================
# GENERATE HTML REPORT
# ============================================================

message("Generating HTML report...")

bridge_table_html <- paste0(
  "<table border='1' cellpadding='4' cellspacing='0' style='border-collapse:collapse; font-size:12px;'>\n",
  "<tr><th>delta</th><th>N</th>",
  paste0("<th>", levels(bf_bridge$order), "</th>", collapse = ""),
  "</tr>\n",
  paste(
    apply(bridge_table, 1, function(row) {
      paste0("<tr>", paste0("<td>", row, "</td>", collapse = ""), "</tr>")
    }),
    collapse = "\n"
  ),
  "\n</table>"
)

html <- sprintf(
  '<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Simulation Report: Interval-Null BF</title>
<style>
  body { font-family: "Helvetica Neue", Arial, sans-serif; max-width: 1200px; margin: auto; padding: 20px; color: #333; }
  h1 { border-bottom: 2px solid #333; padding-bottom: 8px; }
  h2 { color: #555; margin-top: 40px; border-bottom: 1px solid #ccc; padding-bottom: 4px; }
  h3 { color: #777; }
  img { max-width: 100%%; border: 1px solid #ddd; margin: 10px 0; }
  .note { background: #f8f8f0; border-left: 3px solid #E69F00; padding: 10px 15px; margin: 10px 0; font-size: 14px; }
  .finding { background: #f0f8f0; border-left: 3px solid #009E73; padding: 10px 15px; margin: 10px 0; font-size: 14px; }
  table { margin: 10px 0; }
  td, th { text-align: right; padding: 4px 8px; }
  th { background: #f0f0f0; }
</style>
</head><body>

<h1>Simulation Report: Interval-Null BF with Non-local Priors</h1>
<p><strong>Date:</strong> %s &nbsp; | &nbsp; <strong>Simulation grid:</strong> N in {20, 50, 100, 200}, delta in {0, 1, 2, 3}, prior order d in {0, 1, 2, 3}, 50 replicates each.</p>

<div class="note">
<strong>Interval-Null BF:</strong> Instead of testing delta = 0, we test delta &lt; epsilon.
BF10(eps) = P(delta &lt; eps | prior) / P(delta &lt; eps | y).
Two posterior estimation methods: <em>logspline</em> (density-based) and <em>raw</em> (empirical proportion k/N).
Bridge sampling BF serves as the gold-standard benchmark.
</div>

<h2>1. BF Convergence with Epsilon</h2>
<p>How does log BF10 behave as epsilon increases? Triangles mark the first epsilon where the raw method yields a finite BF (left of triangle = Inf). Dotted horizontal lines = bridge sampling benchmark.</p>
<h3>Logspline method</h3>
<img src="fig1a_bf_vs_eps_logspline.png">
<h3>Raw method (with IQR ribbons)</h3>
<img src="fig1b_bf_vs_eps_raw.png">
<div class="finding">
<strong>Finding:</strong> Curves stabilize around eps = 0.2. At smaller eps, the raw method shows Inf (no posterior draws below eps). The logspline method is smoother but underestimates BF at large delta due to boundary extrapolation.
</div>

<h2>2. Effective Convergence Epsilon</h2>
<p>For the raw method: what is the smallest epsilon with any posterior mass below it? This diagnostic shows where the raw interval BF "turns on."</p>
<img src="fig2_convergence_eps.png">
<div class="finding">
<strong>Finding:</strong> Convergence eps increases with delta (posterior moves away from zero) and with prior order (non-local priors push mass away from zero). At delta = 0, convergence eps is near 0.01 for all priors and N.
</div>

<h2>3. BF at Fixed Epsilon = 0.2</h2>
<p>Bars = median log BF10 over 50 replicates at eps = 0.2. Error bars = IQR (finite replicates only). x = bridge sampling median. Bars labeled "Inf" indicate majority of replicates have zero posterior mass below eps.</p>
<h3>Logspline</h3>
<img src="fig3a_bf_eps02_logspline.png">
<h3>Raw</h3>
<img src="fig3b_bf_eps02_raw.png">
<div class="finding">
<strong>Finding:</strong> Logspline gives graded, finite evidence across all settings. Prior order effect is visible: d=0 (local) gives strongest BF10, d=3 weakest. Raw method saturates at Inf for delta >= 1 at small N, but gives finite values at N = 200 for delta = 1.
Under the null (delta = 0), both methods correctly return negative log BF (evidence against nonlinearity), and higher prior order gives more decisively negative values.
</div>

<h2>4. Evidence Accumulation Across Methods</h2>
<p>Each panel shows BF vs N (top row) and BF vs delta (bottom row). Comparing bridge sampling, logspline, and raw side by side.</p>
<h3>Bridge Sampling (benchmark)</h3>
<img src="fig4_bridge_combined.png">
<h3>Logspline (eps = 0.2)</h3>
<img src="fig4_logspline_combined.png">
<h3>Raw (eps = 0.2)</h3>
<img src="fig4_raw_combined.png">

<h3>Bridge BF by Prior Order (exact values)</h3>
<p>Bridge sampling BF is nearly prior-order-invariant. Differences are real but small (&lt; 0.6 on log scale):</p>
%s
<div class="finding">
<strong>Finding:</strong> Bridge sampling shows near-identical curves across prior orders (max difference &lt; 0.6 in log BF). This is because the priors are calibrated to the same median; they differ only near delta = 0, which barely affects the marginal likelihood integral. The interval BF retains sensitivity to prior order because it specifically evaluates evidence near the null.
</div>

<h2>5. Method Comparison: Logspline vs Raw</h2>
<p>Scatter of median log BF10 at eps = 0.2. Points with Inf in either method are excluded. Dashed = identity line.</p>
<img src="fig5_method_comparison.png">
<div class="finding">
<strong>Finding:</strong> Where both methods are finite, they agree well (points near the diagonal). The raw method is unbiased but has lower resolution; logspline provides smoother estimates but can be biased near the boundary.
</div>

<h2>6. CV-Optimal Epsilon Selection</h2>
<p>Principled epsilon: choose the smallest eps where SE(p_hat)/p_hat &lt; 20%%. This ensures the posterior probability estimate (denominator of BF) is precise enough to trust.</p>
<h3>6a. Distribution of CV-optimal epsilon</h3>
<img src="fig6a_cv_optimal_eps.png">
<h3>6b. BF at CV-optimal epsilon</h3>
<img src="fig6b_bf_cv_optimal.png">
<h3>6c. Per-replicate comparison vs bridge sampling</h3>
<img src="fig6c_scatter_cv_vs_bridge.png">
<div class="finding">
<strong>Finding:</strong> CV-optimal eps ranges from ~0.05 (delta=0, strong posterior mass near zero) to ~0.4 (delta=3, posterior far from zero). The interval BF at CV-optimal eps tracks bridge sampling well in moderate regimes but remains conservative for large effects.
</div>

<h2>7. Summary</h2>
<div class="finding">
<ol>
<li><strong>Bridge sampling</strong> is the gold standard: works everywhere, insensitive to prior order, but requires fitting a separate null model.</li>
<li><strong>Interval BF (logspline, eps = 0.2)</strong> is the practical workaround: requires only the alternative model, gives graded evidence, and retains sensitivity to prior order near the null.</li>
<li><strong>eps ~ 0.2</strong> is a defensible operational choice; the CV-based adaptive criterion provides a principled alternative.</li>
<li><strong>Non-local priors</strong> (d > 0) are more conservative under both H0 and H1: stronger evidence for the null when delta = 0, slightly weaker evidence for the alternative otherwise.</li>
<li><strong>Next steps:</strong> Run simulations for N = 500, 1000. Derive the L Hopital limit connecting interval BF to precise-null BF. Refit empirical applications with updated Stan model (v17, full-X projection).</li>
</ol>
</div>

</body></html>',
  Sys.Date(),
  bridge_table_html
)

writeLines(html, file.path(report_dir, "report_sim3.html"))
message("\n=== Report generated: ", file.path(report_dir, "report_sim3.html"), " ===")
