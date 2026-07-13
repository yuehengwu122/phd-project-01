library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

source("src/moment_prior.R")
source("src/bayes_factor.R")

local_data_dir <- "results/simulated_0627"
legacy_data_dir <- "data/simulated"
output_dir <- "results/simulated_0627"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

local_job_tags <- c(
  "job01_N_small_delta0_local",
  "job04_N_small_delta1_local",
  "job07_N_small_delta2_local",
  "job08_N_small_delta3_local"
)

legacy_bridge_tags <- c(
  "job01_N_small_delta0",
  "job04_N_small_delta1"
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

load_csvs <- function(tags, prefix, dir) {
  paths <- file.path(dir, sprintf("%s_%s.csv", prefix, tags))
  existing <- paths[file.exists(paths)]
  if (!length(existing)) {
    return(tibble())
  }

  bind_rows(lapply(existing, read_csv, show_col_types = FALSE))
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

build_base_draws_table <- function(draws, base_order = 0L) {
  draws |>
    filter(order == sprintf("d%d", base_order)) |>
    group_by(delta_true, N, replicate) |>
    summarise(delta_draws = list(delta_draw), .groups = "drop")
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

summarise_bf <- function(df, bf_col) {
  df |>
    group_by(delta_true, N, order, delta_label) |>
    summarise(
      med = median(.data[[bf_col]], na.rm = TRUE),
      q25 = quantile(.data[[bf_col]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[bf_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

build_missing_annotation <- function(missing_delta_vals, label_text) {
  if (!length(missing_delta_vals)) {
    return(tibble())
  }

  tibble(
    delta_true = missing_delta_vals,
    delta_label = factor(
      paste0("delta = ", missing_delta_vals),
      levels = paste0("delta = ", delta_vals)
    ),
    x = 110,
    y = 0,
    label = label_text
  )
}

make_bf_row_plot <- function(
  df_summary,
  title_text,
  y_label,
  annotations = tibble()
) {
  p <- ggplot(df_summary, aes(x = N, y = med, color = order, fill = order)) +
    facet_wrap(~delta_label, nrow = 1, scales = "free_y") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = order_colors, name = "Prior order", drop = FALSE) +
    scale_fill_manual(values = order_colors, name = "Prior order", drop = FALSE) +
    scale_x_continuous(breaks = N_vals) +
    labs(
      title = title_text,
      x = "N",
      y = y_label
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(size = 8)
    )

  if (nrow(annotations) > 0) {
    p <- p + geom_text(
      data = annotations,
      mapping = aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      color = "grey35",
      size = 3
    )
  }

  p
}

message("Loading local-only bridge and posterior draws from simulated_0627...")
bridge_local <- load_csvs(
  local_job_tags,
  "results_sim3_rstan_bf_bridge",
  dir = local_data_dir
)
draws_local <- load_csvs(
  local_job_tags,
  "results_sim3_rstan_posterior_draws",
  dir = local_data_dir
)

if (!nrow(bridge_local) || !nrow(draws_local)) {
  stop(
    "Local-only bridge and posterior draw files were not found in ",
    local_data_dir
  )
}

message("Computing reweighted logspline interval-null BF (eps = 0.2)...")
bf_logspline_reweighted <- build_reweighted_logspline_interval_table(
  draws = draws_local,
  target_orders = 0:3,
  eps = eps_fixed
) |>
  add_labels()

sum_logspline_reweighted <- summarise_bf(
  bf_logspline_reweighted,
  "log_bf10_reweighted"
)

message("Computing reweighted bridge BF from latest local bridge results...")
bf_bridge_reweighted <- build_reweighted_bridge_table(
  draws = draws_local,
  bf_bridge = bridge_local,
  target_orders = 0:3
) |>
  add_labels()

sum_bridge_reweighted <- summarise_bf(
  bf_bridge_reweighted,
  "log_bf10_reweighted"
)

message("Loading legacy direct bridge results when available...")
bridge_legacy <- load_csvs(
  legacy_bridge_tags,
  "results_sim3_rstan_bf_bridge",
  dir = legacy_data_dir
) |>
  filter(N %in% N_vals) |>
  add_labels()

sum_bridge_legacy <- summarise_bf(bridge_legacy, "log_bf10_bridge")
missing_bridge_deltas <- setdiff(delta_vals, unique(bridge_legacy$delta_true))
bridge_annotations <- build_missing_annotation(
  missing_bridge_deltas,
  "legacy direct bridge\nnot available"
)

p_bridge_legacy <- make_bf_row_plot(
  df_summary = sum_bridge_legacy,
  title_text = "Legacy direct bridge sampling",
  y_label = "Median log BF10",
  annotations = bridge_annotations
)

p_bridge_reweighted <- make_bf_row_plot(
  df_summary = sum_bridge_reweighted,
  title_text = "Latest local bridge BF + reweighted nonlocal bridge BF",
  y_label = "Median log BF10"
)

p_logspline_reweighted <- make_bf_row_plot(
  df_summary = sum_logspline_reweighted,
  title_text = sprintf(
    "Reweighting from latest local results [logspline, eps = %.1f]",
    eps_fixed
  ),
  y_label = "Median log BF10"
)

combined_plot <- (
  p_bridge_legacy /
  p_bridge_reweighted /
  p_logspline_reweighted
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf_path <- file.path(output_dir, "sim3_reweighting_0627_overview.pdf")
png_path <- file.path(output_dir, "sim3_reweighting_0627_overview.png")

ggsave(pdf_path, combined_plot, width = 12, height = 11)
ggsave(png_path, combined_plot, width = 12, height = 11, dpi = 300)

message("Saved plot to:")
message("  ", pdf_path)
message("  ", png_path)
