#' Simulation Experiment 2: Interval-Null Bayes Factor with Non-Local Priors
#'
#' Assessment of interval-null Bayes factor behavior under moment priors
#' (non-local priors) of different orders (d = 0, 1, 2, 3).
#'
#' @details
#' Simulation setup:
#' - Data domain: X in [-0.5, 0.5]
#' - Noise standard deviation: sigma = 0.1
#' - GP lengthscale: lambda = 0.3
#' - Sample sizes: N in {20, 50, 100, 200, 500, 1000}
#' - True effect sizes: delta in {0, 1, 2, 3}
#' - NLP orders: d in {0, 1, 2, 3}
#' - Epsilon grid: seq(0.01, 1, length.out = 100)
#' - Replicates per configuration: 50 independent datasets
#' - Output: Median log(BF10) across replicates for each (delta, N, d, epsilon)
#'
#' Expected results:
#' - For delta = 0 (H0 true): log(BF10) should decrease with N;
#'   higher d should show stronger evidence for H0
#' - For delta > 0 (H1 true): log(BF10) should increase with N;
#'   higher d should eventually show stronger evidence for H1 as N grows

source("_common.R")

library(parallel)
library(rstan)

# ============================================
# EXPERIMENT 2: PARAMETER SPECIFICATIONS
# ============================================

# Sample sizes to test
sample_sizes <- c(20, 50, 100, 200, 500, 1000)

# True effect sizes (delta) to vary
delta_values <- c(0, 1, 2, 3)

# NLP orders
nlp_orders <- 0:3

# GP lengthscale (fixed)
lambda_fixed <- 0.3

# Noise standard deviation (fixed)
sigma_fixed <- 0.1

# Domain range
x_min <- -0.5
x_max <- 0.5

# Epsilon grid for interval-null BF
eps_grid <- seq(0.01, 1, length.out = 100)

# Number of independent replicates per configuration
n_replicates <- 50

# Stan sampling settings
stan_iter <- 2000
stan_warmup <- 500
stan_chains <- 2
stan_seed <- 123

# Number of cores for parallelization
n_cores <- min(parallel::detectCores() - 1, 4)

# ============================================
# MOMENT PRIOR SPECIFICATIONS
# (Copied from 260223_nonlocal_prior.qmd)
# ============================================

df0 <- 6
s0 <- 2.79
m_target <- 2

g_base <- function(x, df = df0, s = s0) {
  2 * dt(x / s, df = df) / s
}

C_d_half_t <- function(d, df = df0, s = s0) {
  if (d >= df) {
    stop("Need d < df for the moment to exist (choose df > max d).")
  }
  s^d * (df^(d / 2)) * gamma((d + 1) / 2) * gamma((df - d) / 2) / (sqrt(pi) * gamma(df / 2))
}

g_scaled <- function(x, kappa, df = df0, s = s0) {
  (1 / kappa) * g_base(x / kappa, df = df, s = s)
}

prior_moment <- function(x, d, kappa, df = df0, s = s0) {
  if (any(x < 0)) {
    stop("x must be >= 0")
  }
  if (d == 0) {
    return(g_scaled(x, kappa = kappa, df = df, s = s))
  }
  Ck <- (kappa^d) * C_d_half_t(d, df = df, s = s)
  x^d * g_scaled(x, kappa = kappa, df = df, s = s) / Ck
}

prior_quantile <- function(pdf_fun, p, upper_start = NULL, rel.tol = 1e-8) {
  cdf <- function(q) {
    if (q <= 0) {
      return(0)
    }
    integrate(pdf_fun, lower = 0, upper = q, rel.tol = rel.tol)$value
  }
  upper <- if (is.null(upper_start)) 1 else upper_start
  while (cdf(upper) < p) {
    upper <- upper * 2
  }
  uniroot(function(q) cdf(q) - p, lower = 0, upper = upper, tol = 1e-8)$root
}

compute_kappa <- function(d) {
  pdf_k1 <- function(x) prior_moment(x, d = d, kappa = 1)
  med1 <- prior_quantile(pdf_k1, p = 0.5, upper_start = s0)
  m_target / med1
}

# Precompute kappa values
kappa <- setNames(vapply(nlp_orders, compute_kappa, numeric(1)), paste0("d", nlp_orders))

# ============================================
# HELPER FUNCTIONS
# ============================================

make_stan_data <- function(X, y, d_order) {
  X <- as.matrix(X)
  y <- as.vector(y)
  stopifnot(nrow(X) == length(y))

  list(
    N = nrow(X),
    P = ncol(X),
    X = X,
    y = y,
    d_order = as.integer(d_order),
    kappa_delta = as.numeric(kappa[paste0("d", d_order)]),
    df_delta = as.integer(df0),
    s_delta = as.numeric(s0),
    logC_base = if (d_order == 0) 0 else log(C_d_half_t(d_order, df = df0, s = s0))
  )
}

get_delta_draws_rstan <- function(fit) {
  delta <- rstan::extract(fit)$delta
  as.numeric(delta)
}

# Posterior interval probability via logspline
post_prob_delta_lt_epsilon_logspline <- function(x, eps, max_attempts = 5) {
  stopifnot(is.numeric(x), length(x) > 10, is.finite(eps), eps >= 0)
  x <- x[is.finite(x) & x >= 0]

  quantiles <- c(1, 0.999, 0.99, 0.9, 0.8)

  for (q in quantiles[1:max_attempts]) {
    cutoff <- as.numeric(stats::quantile(x, q, names = FALSE, type = 5))
    keep <- x <= cutoff
    x_trim <- x[keep]

    fit <- try(logspline::logspline(x_trim, lbound = 0, ubound = max(x_trim)), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      p_keep <- mean(keep)
      p_trim <- logspline::plogspline(eps, fit)
      return(p_keep * p_trim)
    }
  }

  warning("logspline failed for all tail-cuts. Returning NA.")
  NA_real_
}

# Prior interval probability
prior_prob_delta_lt_epsilon_moment <- function(d_order, eps, rel.tol = 1e-8) {
  stopifnot(is.finite(eps), eps >= 0)
  key <- paste0("d", d_order)

  if (eps == 0) {
    return(0)
  }

  pdf <- function(z) prior_moment(z, d = d_order, kappa = unname(kappa[key]), df = df0, s = s0)
  integrate(pdf, lower = 0, upper = eps, rel.tol = rel.tol)$value
}

# Compute interval-null BF10 for a single fit
compute_interval_bf <- function(delta_draws, d_order, eps_grid) {
  prior0 <- vapply(eps_grid, function(e) prior_prob_delta_lt_epsilon_moment(d_order, e), numeric(1))
  post0 <- vapply(eps_grid, function(e) post_prob_delta_lt_epsilon_logspline(delta_draws, e), numeric(1))

  log_bf10 <- log(prior0) - log(post0)

  data.frame(
    epsilon = eps_grid,
    prior0 = prior0,
    post0 = post0,
    log_bf10 = log_bf10,
    stringsAsFactors = FALSE
  )
}

# ============================================
# SINGLE REPLICATE FUNCTION (for parallelization)
# ============================================

run_single_replicate <- function(rep_id, delta, N, d_order, X_grid, model) {
  set.seed(rep_id + delta * 1000 + N * 10000 + d_order * 100000)

  # Generate data
  if (delta == 0) {
    y <- rnorm(N, mean = 0, sd = sigma_fixed)
  } else {
    y <- generate_gp_data(
      delta = delta,
      lambda = lambda_fixed,
      sigma = sigma_fixed,
      N = N,
      x = X_grid,
      n_samples = 1
    ) |>
      as.vector()
  }

  # Prepare Stan data
  stan_data <- make_stan_data(X_grid, y, d_order)

  # Fit model
  fit <- rstan::sampling(
    model,
    data = stan_data,
    iter = stan_iter,
    warmup = stan_warmup,
    chains = stan_chains,
    seed = stan_seed,
    refresh = 0,
    verbose = FALSE,
    open_progress = FALSE
  )

  # Extract draws and compute BF

  draws <- get_delta_draws_rstan(fit)
  bf_res <- compute_interval_bf(draws, d_order, eps_grid)

  bf_res$delta_true <- delta
  bf_res$N <- N
  bf_res$order <- paste0("d", d_order)
  bf_res$replicate <- rep_id

  bf_res
}

# ============================================
# MAIN SIMULATION LOOP
# ============================================

message("=== Starting Simulation Experiment 2 ===")
message(sprintf("Sample sizes: %s", paste(sample_sizes, collapse = ", ")))
message(sprintf("Delta values: %s", paste(delta_values, collapse = ", ")))
message(sprintf("NLP orders: %s", paste(nlp_orders, collapse = ", ")))
message(sprintf("Replicates per config: %d", n_replicates))
message(sprintf("Using %d cores for parallelization", n_cores))

# Compile Stan model once
message("Compiling Stan model...")
# Assumes gp_add_16 is available from _common.R

# Initialize results storage
all_results <- list()
result_counter <- 0

total_configs <- length(delta_values) * length(sample_sizes) * length(nlp_orders)
config_counter <- 0

for (delta in delta_values) {
  # Generate master dataset for this delta (N_max = max(sample_sizes))
  N_max <- max(sample_sizes)
  X_max <- seq(x_min, x_max, length = N_max) |> matrix(ncol = 1)

  for (N in sample_sizes) {
    X_grid <- X_max[1:N, , drop = FALSE]

    for (d_order in nlp_orders) {
      config_counter <- config_counter + 1
      message(sprintf(
        "\n[Config %d/%d] delta=%.1f, N=%d, d=%d - Running %d replicates...",
        config_counter,
        total_configs,
        delta,
        N,
        d_order,
        n_replicates
      ))

      # Run replicates in parallel
      rep_results <- mclapply(
        1:n_replicates,
        function(rep_id) {
          tryCatch(
            run_single_replicate(rep_id, delta, N, d_order, X_grid, gp_add_16),
            error = function(e) {
              warning(sprintf("Replicate %d failed: %s", rep_id, e$message))
              NULL
            }
          )
        },
        mc.cores = n_cores
      )

      # Filter out failed replicates
      rep_results <- Filter(Negate(is.null), rep_results)

      if (length(rep_results) > 0) {
        result_counter <- result_counter + 1
        all_results[[result_counter]] <- do.call(rbind, rep_results)
      }

      message(sprintf("  Completed %d/%d replicates successfully", length(rep_results), n_replicates))
    }
  }

  # Save intermediate results after each delta
  intermediate_df <- do.call(rbind, all_results)
  write.csv(
    intermediate_df,
    sprintf("data/simulated/results_sim2_delta%d_intermediate.csv", delta),
    row.names = FALSE
  )
  message(sprintf("Intermediate results saved for delta = %d", delta))
}

# ============================================
# COMBINE AND SUMMARIZE RESULTS
# ============================================

message("\n=== Combining results ===")

full_results <- do.call(rbind, all_results)

# Compute median log(BF10) across replicates
summary_results <- full_results |>
  dplyr::group_by(delta_true, N, order, epsilon) |>
  dplyr::summarise(
    median_log_bf10 = median(log_bf10, na.rm = TRUE),
    mean_log_bf10 = mean(log_bf10, na.rm = TRUE),
    sd_log_bf10 = sd(log_bf10, na.rm = TRUE),
    q25_log_bf10 = quantile(log_bf10, 0.25, na.rm = TRUE),
    q75_log_bf10 = quantile(log_bf10, 0.75, na.rm = TRUE),
    n_valid = sum(!is.na(log_bf10)),
    .groups = "drop"
  )

# Save results
write.csv(
  full_results,
  "data/simulated/results_sim2_full.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "data/simulated/results_sim2_summary.csv",
  row.names = FALSE
)

message("\n=== Experiment 2 COMPLETE ===")
message(sprintf("Total configurations: %d", total_configs))
message(sprintf("Total replicates attempted: %d", total_configs * n_replicates))
message(sprintf("Total rows in full results: %d", nrow(full_results)))
message("Results saved to:")
message("  - data/simulated/results_sim2_full.csv")
message("  - data/simulated/results_sim2_summary.csv")

# ============================================
# PLOT: Median logBF10 vs N at selected epsilon
# with 5%-95% error bars
# ============================================

library(ggplot2)
library(dplyr)
library(readr)

# Prefer full summary if available; otherwise use partial summary
summary_path_full <- "data/simulated/results_sim2_summary.csv"
summary_path_partial <- "data/simulated/results_sim2_partial_summary.csv"

summary_path <- if (file.exists(summary_path_full)) {
  summary_path_full
} else {
  summary_path_partial
}

sim_summary <- readr::read_csv(summary_path, show_col_types = FALSE)

# If 5%/95% columns are missing, compute them from full/partial raw results
if (!all(c("q05_log_bf10", "q95_log_bf10") %in% names(sim_summary))) {
  raw_path_full <- "data/simulated/results_sim2_full.csv"
  raw_path_partial <- "data/simulated/results_sim2_partial.csv"

  raw_path <- if (file.exists(raw_path_full)) {
    raw_path_full
  } else {
    raw_path_partial
  }

  sim_raw <- readr::read_csv(raw_path, show_col_types = FALSE)

  sim_summary <- sim_raw |>
    group_by(delta_true, N, order, epsilon) |>
    summarise(
      median_log_bf10 = median(log_bf10, na.rm = TRUE),
      mean_log_bf10 = mean(log_bf10, na.rm = TRUE),
      sd_log_bf10 = sd(log_bf10, na.rm = TRUE),
      q05_log_bf10 = quantile(log_bf10, 0.05, na.rm = TRUE),
      q25_log_bf10 = quantile(log_bf10, 0.25, na.rm = TRUE),
      q75_log_bf10 = quantile(log_bf10, 0.75, na.rm = TRUE),
      q95_log_bf10 = quantile(log_bf10, 0.95, na.rm = TRUE),
      n_valid = sum(!is.na(log_bf10)),
      .groups = "drop"
    )
}

# Select delta to plot (change if needed: 0, 1, 2, 3)
delta_to_plot <- 0

# Three representative epsilon values
eps_targets <- c(0.01, 0.05, 0.1)

# Match requested epsilons to nearest available values in the data
eps_available <- sort(unique(sim_summary$epsilon))
eps_selected <- eps_available[
  vapply(
    eps_targets,
    function(eps_i) which.min(abs(eps_available - eps_i)),
    integer(1)
  )
] |> unique()

plot_df <- sim_summary |>
  filter(delta_true == delta_to_plot, epsilon %in% eps_selected) |>
  mutate(
    order = factor(order, levels = c("d0", "d1", "d2", "d3")),
    N = factor(N, levels = c(20, 50, 100, 200, 500, 1000)),
    epsilon_lab = paste0("epsilon = ", format(round(epsilon, 3), nsmall = 2))
  )

p_median_logbf10 <- ggplot(
  plot_df,
  aes(x = N, y = median_log_bf10, color = order, group = order)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = q05_log_bf10, ymax = q95_log_bf10),
    width = 0.12,
    alpha = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~epsilon_lab, nrow = 1) +
  scale_color_manual(
    values = c("d0" = "black", "d1" = "red", "d2" = "blue", "d3" = "darkgreen")
  ) +
  labs(
    title = "Median log(BF10) across 50 replicates with 5%-95% intervals",
    subtitle = paste("delta =", delta_to_plot),
    x = "Sample size (N)",
    y = expression(median(log(BF[10]))),
    color = "NLP order"
  ) +
  theme_minimal()

p_median_logbf10