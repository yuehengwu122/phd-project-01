source("src/_setup.R")

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

source("src/moment_prior.R")
source("src/bayes_factor.R")

# ============================================================
# Paper figure: simulation Bayes factors vs sample size
#
# Rows:
# 1. Bridge sampling results.
# 2. d=0 bridge BF plus reweighted bridge BF for d=1,2,3.
# 3. d=0 interval-null BF (logspline, eps = 0.2) plus reweighted
#    non-local BF for d=1,2,3.
# ============================================================

direct_bridge_dirs <- c("data/simulated", "data/simulated_0623")
base_result_dirs <- c("data/simulated", "results/simulated_0627", "data/simulated_0623")
output_dir <- "results/simulated_0627"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

direct_bridge_tags <- c(
  "job01_N_small_delta0",
  "job04_N_small_delta1",
  "job07_N_small_delta2",
  "job08_N_small_delta3"
)

base_job_tags <- c(
  "job01_N_small_delta0",
  "job04_N_small_delta1",
  "job07_N_small_delta2",
  "job08_N_small_delta3"
)

eps_fixed <- 0.2
N_vals <- c(20, 50, 100, 200)
delta_vals <- c(0, 1, 2, 3)
order_levels <- sprintf("d%d", 0:3)
order_labels <- c(
  d0 = "d=0 (local)",
  d1 = "d=1",
  d2 = "d=2",
  d3 = "d=3"
)
order_colors <- c(
  "d=0 (local)" = "#E69F00",
  "d=1" = "#56B4E9",
  "d=2" = "#009E73",
  "d=3" = "#CC79A7"
)

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 8),
    plot.title = element_text(size = 12, face = "plain")
  )

load_csvs_from_dirs <- function(tags, prefix, dirs) {
  pieces <- lapply(tags, function(tag) {
    candidates <- file.path(dirs, sprintf("%s_%s.csv", prefix, tag))
    path <- candidates[file.exists(candidates)][1]
    if (is.na(path)) {
      warning("Missing file for tag: ", tag, call. = FALSE)
      return(tibble())
    }
    read_csv(path, show_col_types = FALSE)
  })

  bind_rows(pieces)
}

add_labels <- function(df) {
  df |>
    mutate(
      delta_label = factor(
        paste0("delta = ", delta_true),
        levels = paste0("delta = ", delta_vals)
      ),
      order = factor(order, levels = order_levels, labels = order_labels)
    )
}

build_base_draws_table <- function(draws, base_order = 0L) {
  draws |>
    filter(order == sprintf("d%d", base_order)) |>
    group_by(delta_true, N, replicate) |>
    summarise(delta_draws = list(delta_draw), .groups = "drop")
}

compute_log_reweight_adjustment_local <- function(
  delta_draws,
  target_order,
  base_order = 0L
) {
  draws <- as.numeric(delta_draws)
  draws <- draws[is.finite(draws) & draws >= 0]
  if (!length(draws)) {
    return(NA_real_)
  }

  log_target <- log(prior_moment_density(draws, target_order, get_kappa(target_order)))
  log_base <- log(prior_moment_density(draws, base_order, get_kappa(base_order)))
  log_weights <- log_target - log_base

  max_log_weight <- max(log_weights)
  max_log_weight + log(mean(exp(log_weights - max_log_weight)))
}

build_reweighted_bridge_table <- function(
  draws,
  bf_bridge,
  target_orders = 0:3,
  base_order = 0L
) {
  base_draws <- build_base_draws_table(draws, base_order = base_order)

  base_bridge <- bf_bridge |>
    filter(order == sprintf("d%d", base_order)) |>
    transmute(
      delta_true,
      N,
      replicate,
      log_bf10_bridge_base = log_bf10_bridge
    )

  base_tbl <- left_join(
    base_draws,
    base_bridge,
    by = c("delta_true", "N", "replicate")
  )

  bind_rows(lapply(target_orders, function(target_order) {
    base_tbl |>
      mutate(
        order = sprintf("d%d", target_order),
        log_reweight_adjustment = vapply(
          delta_draws,
          compute_log_reweight_adjustment_local,
          numeric(1),
          target_order = target_order,
          base_order = base_order
        ),
        log_bf10_reweighted = log_bf10_bridge_base + log_reweight_adjustment
      ) |>
      select(
        delta_true,
        N,
        order,
        replicate,
        log_bf10_bridge_base,
        log_reweight_adjustment,
        log_bf10_reweighted
      )
  }))
}

build_reweighted_logspline_interval_table <- function(
  draws,
  target_orders = 0:3,
  base_order = 0L,
  eps = eps_fixed
) {
  base_tbl <- build_base_draws_table(draws, base_order = base_order) |>
    mutate(
      log_bf10_base = vapply(
        delta_draws,
        function(x) {
          bf10_interval_null_delta(
            delta_draws = x,
            d_order = base_order,
            eps = eps,
            return_log = TRUE
          )$log_bf10[1]
        },
        numeric(1)
      )
    )

  bind_rows(lapply(target_orders, function(target_order) {
    base_tbl |>
      mutate(
        order = sprintf("d%d", target_order),
        log_reweight_adjustment = vapply(
          delta_draws,
          compute_log_reweight_adjustment_local,
          numeric(1),
          target_order = target_order,
          base_order = base_order
        ),
        log_bf10_reweighted = log_bf10_base + log_reweight_adjustment,
        epsilon = eps
      ) |>
      select(
        delta_true,
        N,
        order,
        replicate,
        epsilon,
        log_bf10_base,
        log_reweight_adjustment,
        log_bf10_reweighted
      )
  }))
}

summarise_bf <- function(df, bf_col) {
  df |>
    group_by(delta_true, N, order, delta_label) |>
    summarise(
      med = median(.data[[bf_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

make_column_limits <- function(df_summary) {
  df_summary |>
    group_by(delta_true, delta_label) |>
    summarise(
      ymin = min(c(0, med), na.rm = TRUE),
      ymax = max(c(0, med), na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      yrange = pmax(ymax - ymin, 0.50),
      ymin = ymin - 0.08 * yrange,
      ymax = ymax + 0.08 * yrange
    ) |>
    select(delta_label, ymin, ymax) |>
    tidyr::pivot_longer(
      cols = c(ymin, ymax),
      names_to = "bound",
      values_to = "y"
    ) |>
    mutate(N = ifelse(bound == "ymin", min(N_vals), max(N_vals)))
}

make_bf_row_plot <- function(
  df_summary,
  title_text,
  y_label = "Median log BF10",
  y_limits = NULL,
  show_x_label = TRUE
) {
  p <- ggplot(df_summary, aes(x = N, y = med, color = order)) +
    facet_wrap(~delta_label, nrow = 1, scales = "free_y") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = order_colors, name = "Prior order", drop = FALSE) +
    scale_x_continuous(breaks = N_vals) +
    labs(
      title = title_text,
      x = if (show_x_label) "Sample size N" else NULL,
      y = y_label
    ) +
    base_theme

  if (!is.null(y_limits)) {
    p <- p + geom_blank(
      data = y_limits,
      mapping = aes(x = N, y = y),
      inherit.aes = FALSE
    )
  }

  if (!show_x_label) {
    p <- p + theme(axis.title.x = element_blank())
  }

  p
}

message("Loading direct bridge files...")
bridge_direct <- load_csvs_from_dirs(
  direct_bridge_tags,
  "results_sim3_rstan_bf_bridge",
  dirs = direct_bridge_dirs
) |>
  filter(N %in% N_vals)

message("Loading base d=0 bridge files and posterior draws...")
bridge_base <- load_csvs_from_dirs(
  base_job_tags,
  "results_sim3_rstan_bf_bridge",
  dirs = base_result_dirs
) |>
  filter(N %in% N_vals)

draws_base <- load_csvs_from_dirs(
  base_job_tags,
  "results_sim3_rstan_posterior_draws",
  dirs = base_result_dirs
) |>
  filter(N %in% N_vals)

if (!nrow(bridge_base) || !nrow(draws_base)) {
  stop("Local bridge files or posterior draws are missing.")
}

message("Computing bridge-based reweighting...")
bf_bridge_reweighted <- build_reweighted_bridge_table(
  draws = draws_base,
  bf_bridge = bridge_base,
  target_orders = 0:3
)

message("Computing logspline-interval reweighting (eps = 0.2)...")
bf_logspline_reweighted <- build_reweighted_logspline_interval_table(
  draws = draws_base,
  target_orders = 0:3,
  eps = eps_fixed
)

bridge_row_raw <- bridge_direct |>
  transmute(
    delta_true,
    N,
    order,
    replicate,
    log_bf10_plot = log_bf10_bridge
  ) |>
  add_labels()

bridge_reweighted_row <- bf_bridge_reweighted |>
  transmute(
    delta_true,
    N,
    order,
    replicate,
    log_bf10_plot = log_bf10_reweighted
  ) |>
  add_labels()

logspline_reweighted_row <- bf_logspline_reweighted |>
  transmute(
    delta_true,
    N,
    order,
    replicate,
    log_bf10_plot = log_bf10_reweighted
  ) |>
  add_labels()

sum_bridge_row <- summarise_bf(bridge_row_raw, "log_bf10_plot")
sum_bridge_reweighted <- summarise_bf(bridge_reweighted_row, "log_bf10_plot")
sum_logspline_reweighted <- summarise_bf(logspline_reweighted_row, "log_bf10_plot")
column_limits <- make_column_limits(sum_bridge_row)

p_bridge <- make_bf_row_plot(
  df_summary = sum_bridge_row,
  title_text = "Bridge Sampling",
  y_limits = column_limits,
  show_x_label = FALSE
)

p_bridge_reweighted <- make_bf_row_plot(
  df_summary = sum_bridge_reweighted,
  title_text = "Reweighting from d=0 Posterior (Bridge BF Base)",
  y_limits = column_limits,
  show_x_label = FALSE
)

p_logspline_reweighted <- make_bf_row_plot(
  df_summary = sum_logspline_reweighted,
  title_text = sprintf(
    "Reweighting from d=0 Posterior (Interval-Null Base, eps = %.2f)",
    eps_fixed
  ),
  y_limits = column_limits,
  show_x_label = TRUE
)

combined_plot <- (
  p_bridge /
  p_bridge_reweighted /
  p_logspline_reweighted
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf_path <- file.path(output_dir, "sim3_reweighting_paper_figure.pdf")
png_path <- file.path(output_dir, "sim3_reweighting_paper_figure.png")

ggsave(pdf_path, combined_plot, width = 12, height = 10.5)
ggsave(png_path, combined_plot, width = 12, height = 10.5, dpi = 300)

message("Saved plot to:")
message("  ", pdf_path)
message("  ", png_path)
