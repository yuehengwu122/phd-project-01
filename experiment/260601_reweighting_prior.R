library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

source("src/moment_prior.R")
source("src/bayes_factor.R")

# Shared simulation settings for the BF-vs-N experiments.
data_dir <- "results/simulated"
job_tags <- c(
  "job01_N_small_delta0",
  "job04_N_small_delta1",
  "job07_N_small_delta2",
  "job08_N_small_delta3"
)
EPS_FIXED <- 0.2
BASE_PRIOR_NULL_MASS_GRID <- c(0.05, 0.10, 0.20)
N_vals <- c(20, 50, 100, 200)
delta_vals <- c(0, 1, 2, 3)
nlp_orders <- 0:3

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

y_limits_by_delta <- data.frame(
  delta_true = delta_vals,
  ymin = c(-3, -1, -1, -1),
  ymax = c(1, 10, 40, 62)
) |>
  mutate(
    delta_label = factor(
      paste0("delta = ", delta_true),
      levels = paste0("delta = ", delta_vals)
    )
  )

base_theme <- theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(size = 8)
  )

add_labels <- function(df) {
  df |>
    mutate(
      delta_label = factor(
        paste0("delta = ", delta_true),
        levels = paste0("delta = ", delta_vals)
      ),
      N_label = factor(
        paste0("N = ", N),
        levels = paste0("N = ", N_vals)
      ),
      order = factor(order, levels = names(order_labels), labels = order_labels)
    )
}

epsilon_from_base_prior_mass <- function(alpha, base_order = 0L) {
  stopifnot(length(alpha) >= 1, all(is.finite(alpha)), all(alpha > 0), all(alpha < 1))
  stopifnot(base_order %in% nlp_orders)

  kap <- get_kappa(base_order)
  pdf_fun <- function(x) prior_moment_density(x, d = base_order, kappa = kap)

  vapply(
    alpha,
    function(a) prior_quantile(pdf_fun, p = a, upper_start = PRIOR_MEDIAN_TARGET),
    numeric(1)
  )
}

load_bridge_results <- function(tags = job_tags, dir = data_dir) {
  lapply(tags, function(tag) {
    read_csv(
      file.path(dir, sprintf("results_sim3_rstan_bf_bridge_%s.csv", tag)),
      show_col_types = FALSE
    )
  }) |>
    bind_rows()
}

load_posterior_draws <- function(tags = job_tags, dir = data_dir) {
  lapply(tags, function(tag) {
    read_csv(
      file.path(dir, sprintf("results_sim3_rstan_posterior_draws_%s.csv", tag)),
      show_col_types = FALSE
    )
  }) |>
    bind_rows()
}

compute_log_reweight_adjustment <- function(
  delta_draws,
  target_order,
  base_order = 0L
) {
  stopifnot(target_order %in% nlp_orders, base_order %in% nlp_orders)

  draws <- delta_draws[is.finite(delta_draws) & delta_draws >= 0]
  if (!length(draws)) {
    return(NA_real_)
  }

  log_target <- log(prior_moment_density(draws, target_order, get_kappa(target_order)))
  log_base <- log(prior_moment_density(draws, base_order, get_kappa(base_order)))

  max_log_weight <- max(log_target - log_base)
  max_log_weight + log(mean(exp(log_target - log_base - max_log_weight)))
}

build_base_draws_table <- function(draws, base_order = 0L) {
  draws |>
    filter(order == sprintf("d%s", base_order)) |>
    group_by(delta_true, N, replicate) |>
    summarise(delta_draws = list(delta_draw), .groups = "drop")
}

compute_log_sdr_base_factor <- function(delta_draws, base_order = 0L) {
  prior0 <- prior_moment_density(0, d = base_order, kappa = get_kappa(base_order))
  post0 <- post_prob_delta_logspline(delta_draws, eps = 0)

  if (!is.finite(post0) || post0 <= 0) {
    return(NA_real_)
  }

  log(prior0) - log(post0)
}

build_reweighted_table <- function(
  base_tbl,
  base_log_bf_col,
  output_log_bf_col,
  target_orders = nlp_orders,
  base_order = 0L
) {
  bind_rows(lapply(target_orders, function(target_order) {
    base_tbl |>
      mutate(
        order = sprintf("d%s", target_order),
        log_reweight_adjustment = vapply(
          delta_draws,
          compute_log_reweight_adjustment,
          numeric(1),
          target_order = target_order,
          base_order = base_order
        ),
        !!output_log_bf_col := .data[[base_log_bf_col]] + log_reweight_adjustment
      )
  })) |>
    select(
      delta_true,
      N,
      order,
      replicate,
      all_of(base_log_bf_col),
      log_reweight_adjustment,
      all_of(output_log_bf_col)
    )
}

build_reweighted_bridge_table <- function(
  draws,
  bf_bridge,
  target_orders = nlp_orders,
  base_order = 0L
) {
  base_draws <- build_base_draws_table(draws, base_order = base_order)

  base_bridge <- bf_bridge |>
    filter(order == sprintf("d%s", base_order)) |>
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

  build_reweighted_table(
    base_tbl = base_tbl,
    base_log_bf_col = "log_bf10_bridge_base",
    output_log_bf_col = "log_bf10_reweighted",
    target_orders = target_orders,
    base_order = base_order
  )
}

build_reweighted_sdr_table <- function(
  draws,
  target_orders = nlp_orders,
  base_order = 0L
) {
  base_tbl <- build_base_draws_table(draws, base_order = base_order) |>
    mutate(
      log_bf10_sdr_base = vapply(
        delta_draws,
        compute_log_sdr_base_factor,
        numeric(1),
        base_order = base_order
      )
    )

  build_reweighted_table(
    base_tbl = base_tbl,
    base_log_bf_col = "log_bf10_sdr_base",
    output_log_bf_col = "log_bf10_sdr_reweighted",
    target_orders = target_orders,
    base_order = base_order
  )
}

compute_log_interval_base_factor <- function(
  delta_draws,
  eps = EPS_FIXED,
  base_order = 0L,
  method = "logspline"
) {
  bf10_interval_null_delta(
    delta_draws = delta_draws,
    d_order = base_order,
    eps = eps,
    method = method,
    return_log = TRUE
  )$log_bf10[1]
}

build_reweighted_interval_logspline_base_table <- function(
  draws,
  target_orders = nlp_orders,
  base_order = 0L,
  eps = EPS_FIXED
) {
  base_tbl <- build_base_draws_table(draws, base_order = base_order) |>
    mutate(
      epsilon = eps,
      log_bf10_interval_logspline_base = vapply(
        delta_draws,
        compute_log_interval_base_factor,
        numeric(1),
        eps = eps,
        base_order = base_order,
        method = "logspline"
      )
    )

  build_reweighted_table(
    base_tbl = base_tbl,
    base_log_bf_col = "log_bf10_interval_logspline_base",
    output_log_bf_col = "log_bf10_interval_logspline_reweighted",
    target_orders = target_orders,
    base_order = base_order
  ) |>
    mutate(epsilon = eps, .before = log_bf10_interval_logspline_base)
}

compute_weighted_interval_stats <- function(
  delta_draws,
  target_order,
  eps = EPS_FIXED,
  base_order = 0L,
  method = c("ecdf", "laplace", "jeffreys")
) {
  method <- match.arg(method)

  draws <- delta_draws[is.finite(delta_draws) & delta_draws >= 0]
  if (!length(draws)) {
    return(NA_real_)
  }

  weights <- prior_moment_density(draws, target_order, get_kappa(target_order)) /
    prior_moment_density(draws, base_order, get_kappa(base_order))

  total_weight <- sum(weights)

  if (!is.finite(total_weight) || total_weight <= 0) {
    return(c(
      p_hat = NA_real_,
      n_eff = NA_real_,
      post0 = NA_real_
    ))
  }

  weights_norm <- weights / total_weight
  p_hat <- sum(weights_norm[draws < eps])
  n_eff <- 1 / sum(weights_norm^2)

  post0 <- switch(
    method,
    ecdf = p_hat,
    laplace = {
      k_eff <- n_eff * p_hat
      (k_eff + 1) / (n_eff + 2)
    },
    jeffreys = {
      k_eff <- n_eff * p_hat
      (k_eff + 0.5) / (n_eff + 1)
    }
  )

  c(
    p_hat = p_hat,
    n_eff = n_eff,
    post0 = post0
  )
}

build_reweighted_interval_table <- function(
  draws,
  target_orders = nlp_orders,
  base_order = 0L,
  eps = EPS_FIXED
) {
  base_tbl <- build_base_draws_table(draws, base_order = base_order)
  method_map <- c(
    weighted_ecdf = "ecdf",
    weighted_laplace = "laplace",
    weighted_jeffreys = "jeffreys"
  )

  bind_rows(lapply(target_orders, function(target_order) {
    prior0 <- prior_prob_delta(target_order, eps)

    bind_rows(lapply(names(method_map), function(method_label) {
      stats_mat <- t(vapply(
        base_tbl$delta_draws,
        compute_weighted_interval_stats,
        numeric(3),
        target_order = target_order,
        eps = eps,
        base_order = base_order,
        method = method_map[[method_label]]
      ))

      base_tbl |>
        mutate(
          order = sprintf("d%s", target_order),
          epsilon = eps,
          prior0 = prior0,
          p_hat = stats_mat[, "p_hat"],
          n_eff = stats_mat[, "n_eff"],
          post0 = stats_mat[, "post0"],
          log_bf10 = log(prior0) - log(post0),
          method = method_label
        ) |>
        select(
          delta_true,
          N,
          order,
          replicate,
          epsilon,
          prior0,
          p_hat,
          n_eff,
          post0,
          log_bf10,
          method
        )
    }))
  }))
}

compare_reweighting_to_bridge <- function(reweighted, bf_bridge) {
  reweighted |>
    left_join(
      bf_bridge |>
        select(delta_true, N, order, replicate, log_bf10_bridge),
      by = c("delta_true", "N", "order", "replicate")
    ) |>
    mutate(
      error = log_bf10_reweighted - log_bf10_bridge,
      abs_error = abs(error)
    )
}

summarise_reweighting_error <- function(comparison_tbl) {
  comparison_tbl |>
    group_by(delta_true, N, order) |>
    summarise(
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      max_abs_error = max(abs_error, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_bf <- function(df, bf_col) {
  df |>
    mutate(lbf = ifelse(is.finite(.data[[bf_col]]), .data[[bf_col]], NA_real_)) |>
    group_by(delta_label, N_label, delta_true, N, order) |>
    summarise(
      med = median(lbf, na.rm = TRUE),
      q25 = quantile(lbf, 0.25, na.rm = TRUE),
      q75 = quantile(lbf, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_bf_with_inf <- function(df, bf_col) {
  df |>
    group_by(delta_label, N_label, delta_true, N, order) |>
    mutate(lbf = .data[[bf_col]]) |>
    summarise(
      med = median(ifelse(is.finite(lbf), lbf, NA_real_), na.rm = TRUE),
      q25 = quantile(ifelse(is.finite(lbf), lbf, NA_real_), 0.25, na.rm = TRUE),
      q75 = quantile(ifelse(is.finite(lbf), lbf, NA_real_), 0.75, na.rm = TRUE),
      maj_inf = sum(is.infinite(lbf) & lbf > 0) > dplyr::n() / 2,
      .groups = "drop"
    )
}

summarise_metric <- function(df, value_col) {
  df |>
    group_by(delta_label, N_label, delta_true, N, order) |>
    summarise(
      med = median(.data[[value_col]], na.rm = TRUE),
      q25 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

make_bf_vs_N_plot <- function(
  df_summary,
  title_method,
  use_shared_y_limits = TRUE,
  flag_inf = FALSE
) {
  p <- ggplot(df_summary, aes(x = N, y = med, color = order, fill = order)) +
    facet_wrap(~delta_label, nrow = 1, scales = "free_y") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = order_colors, name = "Prior order") +
    scale_fill_manual(values = order_colors, name = "Prior order") +
    scale_x_continuous(breaks = N_vals) +
    labs(
      title = sprintf("log BF10 vs N [%s]", title_method),
      x = "N",
      y = "Median log BF10"
    ) +
    base_theme

  if (use_shared_y_limits) {
    p <- p +
      geom_blank(
        data = y_limits_by_delta,
        aes(x = N_vals[1], y = ymin),
        inherit.aes = FALSE
      ) +
      geom_blank(
        data = y_limits_by_delta,
        aes(x = N_vals[1], y = ymax),
        inherit.aes = FALSE
      )
  }

  if (flag_inf && "maj_inf" %in% names(df_summary)) {
    p <- p + geom_point(
      data = dplyr::filter(df_summary, maj_inf),
      shape = 4,
      size = 3.5,
      stroke = 1.2,
      show.legend = FALSE
    )
  }

  p
}

make_metric_vs_N_plot <- function(
  df_summary,
  title_method,
  y_label,
  log10_y = FALSE
) {
  p <- ggplot(df_summary, aes(x = N, y = med, color = order, fill = order)) +
    facet_wrap(~delta_label, nrow = 1, scales = "free_y") +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = order_colors, name = "Prior order") +
    scale_fill_manual(values = order_colors, name = "Prior order") +
    scale_x_continuous(breaks = N_vals) +
    labs(
      title = title_method,
      x = "N",
      y = y_label
    ) +
    base_theme

  if (log10_y) {
    p <- p + scale_y_log10()
  }

  p
}

build_interval_method_plot_bundle <- function(
  draws,
  alpha,
  base_order = 0L
) {
  eps <- epsilon_from_base_prior_mass(alpha, base_order = base_order)

  bf_interval_raw <- build_reweighted_interval_table(
    draws = draws,
    target_orders = nlp_orders,
    base_order = base_order,
    eps = eps
  )
  bf_interval <- add_labels(bf_interval_raw)

  sum_ecdf <- bf_interval |>
    dplyr::filter(method == "weighted_ecdf") |>
    summarise_bf_with_inf("log_bf10")
  sum_laplace <- bf_interval |>
    dplyr::filter(method == "weighted_laplace") |>
    summarise_bf("log_bf10")
  sum_jeffreys <- bf_interval |>
    dplyr::filter(method == "weighted_jeffreys") |>
    summarise_bf("log_bf10")

  plot_obj <- (
    make_bf_vs_N_plot(
      sum_ecdf,
      sprintf("Weighted ECDF (alpha = %.2f, eps = %.3f)", alpha, eps),
      flag_inf = TRUE
    ) /
    make_bf_vs_N_plot(
      sum_laplace,
      sprintf("Weighted Laplace smoothing (alpha = %.2f, eps = %.3f)", alpha, eps)
    ) /
    make_bf_vs_N_plot(
      sum_jeffreys,
      sprintf("Weighted Jeffreys smoothing (alpha = %.2f, eps = %.3f)", alpha, eps)
    )
  ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  list(
    alpha = alpha,
    epsilon = eps,
    bf_interval_raw = bf_interval_raw,
    bf_interval = bf_interval,
    plot = plot_obj
  )
}

inspect_reweight_case <- function(
  delta_true_i,
  N_i,
  target_order,
  replicate_i,
  comparison_tbl = reweighting_comparison
) {
  comparison_tbl |>
    filter(
      delta_true == delta_true_i,
      N == N_i,
      order == sprintf("d%s", target_order),
      replicate == replicate_i
    ) |>
    select(
      delta_true,
      N,
      order,
      replicate,
      log_bf10_bridge,
      log_bf10_reweighted,
      error
    )
}

message("Loading bridge BF...")
bf_bridge_raw <- load_bridge_results()
bf_bridge <- add_labels(bf_bridge_raw)

message("Loading posterior draws...")
draws_all <- load_posterior_draws()

message("Computing reweighted BF from d=0 posterior draws...")
bf_reweighted_raw <- build_reweighted_bridge_table(draws_all, bf_bridge_raw)
bf_reweighted <- add_labels(bf_reweighted_raw)

message("Computing logspline-SDR base BF and reweighting...")
bf_sdr_reweighted_raw <- build_reweighted_sdr_table(draws_all)
bf_sdr_reweighted <- add_labels(bf_sdr_reweighted_raw)

message("Computing interval-null logspline base BF and reweighting...")
bf_interval_logspline_base_reweighted_raw <- build_reweighted_interval_logspline_base_table(
  draws_all,
  eps = EPS_FIXED
)
bf_interval_logspline_base_reweighted <- add_labels(
  bf_interval_logspline_base_reweighted_raw
)

bf_interval_reweighted_raw <- build_reweighted_interval_table(draws_all)
bf_interval_reweighted <- add_labels(bf_interval_reweighted_raw)

reweighting_comparison <- compare_reweighting_to_bridge(
  bf_reweighted_raw,
  bf_bridge_raw
)
reweighting_validation_summary <- summarise_reweighting_error(
  reweighting_comparison
)
reweighting_validation_overall <- reweighting_comparison |>
  summarise(
    mean_abs_error = mean(abs_error, na.rm = TRUE),
    max_abs_error = max(abs_error, na.rm = TRUE)
  )

sum_bridge <- summarise_bf(bf_bridge, "log_bf10_bridge")
sum_reweighted <- summarise_bf(bf_reweighted, "log_bf10_reweighted")
sum_sdr_reweighted <- summarise_bf(bf_sdr_reweighted, "log_bf10_sdr_reweighted")
sum_interval_logspline_base_reweighted <- summarise_bf(
  bf_interval_logspline_base_reweighted,
  "log_bf10_interval_logspline_reweighted"
)
sum_interval_reweighted_ecdf <- bf_interval_reweighted |>
  dplyr::filter(method == "weighted_ecdf") |>
  summarise_bf_with_inf("log_bf10")
sum_interval_reweighted_laplace <- bf_interval_reweighted |>
  dplyr::filter(method == "weighted_laplace") |>
  summarise_bf("log_bf10")
sum_interval_reweighted_jeffreys <- bf_interval_reweighted |>
  dplyr::filter(method == "weighted_jeffreys") |>
  summarise_bf("log_bf10")
diag_interval_weighted_ecdf <- bf_interval_reweighted |>
  dplyr::filter(method == "weighted_ecdf") |>
  mutate(
    p_hat_plot = pmax(p_hat, 1e-12)
  )
sum_interval_p_hat <- diag_interval_weighted_ecdf |>
  summarise_metric("p_hat_plot")
sum_interval_n_eff <- diag_interval_weighted_ecdf |>
  summarise_metric("n_eff")

p_bridge_row <- make_bf_vs_N_plot(sum_bridge, "Bridge Sampling")
p_reweighted_row <- make_bf_vs_N_plot(sum_reweighted, "Reweighting from d=0 posterior")
p_sdr_reweighted_row <- make_bf_vs_N_plot(
  sum_sdr_reweighted,
  "Reweighting with logspline SDR base"
)
p_interval_logspline_base_reweighted_row <- make_bf_vs_N_plot(
  sum_interval_logspline_base_reweighted,
  sprintf("Reweighting with logspline interval base (eps = %.2f)", EPS_FIXED)
)
p_sdr_reweighted_zoom_row <- make_bf_vs_N_plot(
  sum_sdr_reweighted,
  "Reweighting with logspline SDR base (zoomed y-scale)",
  use_shared_y_limits = FALSE
)

p_bridge_vs_reweighting <- (
  p_bridge_row /
  p_reweighted_row /
  p_sdr_reweighted_row /
  p_interval_logspline_base_reweighted_row /
  p_sdr_reweighted_zoom_row
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_reweighted_interval_methods <- (
  make_bf_vs_N_plot(
    sum_interval_reweighted_ecdf,
    sprintf("Weighted ECDF (eps = %.2f)", EPS_FIXED),
    flag_inf = TRUE
  ) /
  make_bf_vs_N_plot(
    sum_interval_reweighted_laplace,
    sprintf("Weighted Laplace smoothing (eps = %.2f)", EPS_FIXED)
  ) /
  make_bf_vs_N_plot(
    sum_interval_reweighted_jeffreys,
    sprintf("Weighted Jeffreys smoothing (eps = %.2f)", EPS_FIXED)
  )
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

p_reweighted_interval_diagnostics <- (
  make_metric_vs_N_plot(
    sum_interval_p_hat,
    sprintf("Weighted posterior interval mass p_hat (eps = %.2f, floored at 1e-12)", EPS_FIXED),
    "Median p_hat",
    log10_y = TRUE
  ) /
  make_metric_vs_N_plot(
    sum_interval_n_eff,
    "Importance-sampling effective sample size",
    "Median n_eff"
  )
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

interval_methods_alpha_010 <- build_interval_method_plot_bundle(
  draws = draws_all,
  alpha = 0.10
)
interval_methods_alpha_020 <- build_interval_method_plot_bundle(
  draws = draws_all,
  alpha = 0.20
)

p_reweighted_interval_methods_alpha_010 <- interval_methods_alpha_010$plot
p_reweighted_interval_methods_alpha_020 <- interval_methods_alpha_020$plot

if (interactive()) {
  print(p_bridge_vs_reweighting)
}

epsilon_calibration_table <- data.frame(
  alpha = BASE_PRIOR_NULL_MASS_GRID,
  epsilon = epsilon_from_base_prior_mass(BASE_PRIOR_NULL_MASS_GRID),
  stringsAsFactors = FALSE
)
