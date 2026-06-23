#' Simulation Experiment 3: Interval-Null BF + Bridge Sampling BF with NLP (RStan, Saturn-ready)
#'
#' Key features:
#' - Uses rstan for model fitting
#' - Fits both alternative (delta free) and null (delta = 0) models
#' - Computes interval-null BF via posterior/prior densities at delta=0
#' - Computes classical BF via bridge sampling
#' - Saves both BF methods to CSV (resume-friendly)
#' - Saves posterior draws incrementally to a separate CSV (optional thinning)
#' - Frees memory aggressively after each fit (rm + gc)
#'
#' Recommended Saturn Cloud machine:
#' - 4XLARGE (16 vCPU, 128 GiB RAM)
#'
#' Output files:
#' - data/simulated/results_sim3_rstan_bf_interval.csv      (interval-null BF)
#' - data/simulated/results_sim3_rstan_bf_bridge.csv        (bridge sampling BF)
#' - data/simulated/results_sim3_rstan_posterior_draws.csv
#' - data/simulated/results_sim3_rstan_bf_summary.csv
#' - data/simulated/results_sim3_rstan_completed_configs.csv

source("_common.R")

library(rstan)
library(bridgesampling)
library(dplyr)
library(readr)
library(logspline)

# ============================================
# SATURN/COMPUTE SETTINGS
# ============================================

# Detect available CPUs (Saturn often exposes this correctly)
cpu_available <- parallel::detectCores(logical = TRUE)

# Keep some headroom for OS and I/O
n_cores_use <- min(8L, max(2L, floor(cpu_available * 0.5)))
chains <- min(4L, n_cores_use)

rstan_options(auto_write = TRUE)
options(mc.cores = n_cores_use)

message(sprintf("Detected CPUs: %d | chains=%d | using cores=%d", cpu_available, chains, n_cores_use))
message("Running delta-specific script: delta=3")

# ============================================
# EXPERIMENT SETTINGS
# ============================================

sample_sizes <- c(20, 50, 100, 200)
delta_values <- c(3)
nlp_orders <- 0:3

eps_grid <- seq(0.01, 1, length.out = 100)
n_replicates <- 50

# Data-generating settings
lambda_fixed <- 0.3
sigma_fixed <- 0.1
x_min <- -0.5
x_max <- 0.5

# sampling settings
iter_warmup <- 500
iter_sampling <- 1500
adapt_delta <- 0.95
max_treedepth <- 12
seed_base <- 123

# Posterior draw storage settings
SAVE_POSTERIOR_DRAWS <- TRUE
POSTERIOR_THIN <- 5L # save every 5th draw; set 1L to save all draws

# ============================================
# PATHS
# ============================================

dir.create("data/simulated", recursive = TRUE, showWarnings = FALSE)

bf_interval_path <- "data/simulated/results_sim3_rstan_bf_interval_d3.csv"
bf_bridge_path <- "data/simulated/results_sim3_rstan_bf_bridge_d3.csv"
posterior_path <- "data/simulated/results_sim3_rstan_posterior_draws_d3.csv"
summary_path <- "data/simulated/results_sim3_rstan_bf_summary_d3.csv"
completed_path <- "data/simulated/results_sim3_rstan_completed_configs_d3.csv"

stan_file <- "stan/multiple_gp/fit_add_16.stan"
stan_file_null <- "stan/multiple_gp/fit_null.stan"

if (!file.exists(stan_file)) {
  stop("Stan file not found: ", stan_file)
}
if (!file.exists(stan_file_null)) {
  stop("Stan null model file not found: ", stan_file_null)
}

if (!exists("gp_simu")) {
  message("Compiling simulation data generator model (gp_simu)...")
  gp_simu <- compile_model("simu")
}

# ============================================
# MOMENT PRIOR HELPERS
# ============================================

df0 <- 6
s0 <- 2.79
m_target <- 2

g_base <- function(x, df = df0, s = s0) {
  2 * dt(x / s, df = df) / s
}

C_d_half_t <- function(d, df = df0, s = s0) {
  if (d >= df) {
    stop("Need d < df for finite moment")
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

kappa <- setNames(vapply(nlp_orders, compute_kappa, numeric(1)), paste0("d", nlp_orders))

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

# ============================================
# BF HELPERS (logspline + raw)
# ============================================

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

  NA_real_
}

post_prob_delta_lt_epsilon_raw <- function(x, eps) {
  stopifnot(is.numeric(x), length(x) > 0, is.finite(eps), eps >= 0)
  x <- x[is.finite(x) & x >= 0]

  p_hat <- mean(x < eps)
  mcse <- sqrt(p_hat * (1 - p_hat) / length(x))

  data.frame(post0 = p_hat, post0_mcse = mcse)
}

prior_prob_delta_lt_epsilon_moment <- function(d_order, eps, rel.tol = 1e-8) {
  stopifnot(is.finite(eps), eps >= 0)
  key <- paste0("d", d_order)
  if (!key %in% names(kappa)) {
    stop("Missing kappa for d_order")
  }
  if (eps == 0) {
    return(0)
  }

  pdf <- function(z) prior_moment(z, d = d_order, kappa = unname(kappa[key]), df = df0, s = s0)
  integrate(pdf, lower = 0, upper = eps, rel.tol = rel.tol)$value
}

post_prob_vec <- function(x, eps_vals, method = c("logspline", "raw"), max_attempts = 5) {
  method <- match.arg(method)

  if (method == "logspline") {
    return(
      data.frame(
        post0 = vapply(eps_vals, function(e) post_prob_delta_lt_epsilon_logspline(x, e, max_attempts), numeric(1)),
        post0_mcse = NA_real_
      )
    )
  }

  do.call(rbind, lapply(eps_vals, function(e) post_prob_delta_lt_epsilon_raw(x, e)))
}

bf10_interval_null_delta_moment <- function(
  delta_draws,
  d_order,
  eps,
  post0_method = c("logspline", "raw"),
  return_log = TRUE,
  max_attempts = 5
) {
  post0_method <- match.arg(post0_method)

  delta_draws <- as.numeric(delta_draws)
  delta_draws <- delta_draws[is.finite(delta_draws) & delta_draws >= 0]

  prior0 <- vapply(eps, function(e) prior_prob_delta_lt_epsilon_moment(d_order, e), numeric(1))
  post0_df <- post_prob_vec(delta_draws, eps_vals = eps, method = post0_method, max_attempts = max_attempts)
  post0 <- post0_df$post0

  bf10 <- prior0 / post0
  if (return_log) {
    bf10 <- log(bf10)
  }

  out <- data.frame(
    epsilon = eps,
    prior0 = prior0,
    post0 = post0,
    post0_mcse = post0_df$post0_mcse,
    stringsAsFactors = FALSE
  )

  if (return_log) {
    out$log_bf10 <- bf10
  } else {
    out$bf10 <- bf10
  }

  out
}

# ============================================
# RESUME SUPPORT
# ============================================

completed_configs <- data.frame(
  delta_true = numeric(),
  N = integer(),
  order = character(),
  replicate = integer(),
  stringsAsFactors = FALSE
)

if (file.exists(completed_path)) {
  completed_configs <- readr::read_csv(completed_path, show_col_types = FALSE)
}

append_csv <- function(df, path) {
  if (file.exists(path)) {
    readr::write_csv(df, path, append = TRUE)
  } else {
    readr::write_csv(df, path)
  }
}

# ============================================
# COMPILE RSTAN MODELS
# ============================================

message("Compiling alternative Stan model...")
stan_mod_alt <- rstan::stan_model(stan_file)

message("Compiling null Stan model...")
stan_mod_null <- rstan::stan_model(stan_file_null)

# ============================================
# MAIN LOOP
# ============================================

total_configs <- length(delta_values) * length(sample_sizes) * length(nlp_orders) * n_replicates
counter <- 0L

for (delta in delta_values) {
  N_max <- max(sample_sizes)
  X_master <- seq(x_min, x_max, length = N_max) |> matrix(ncol = 1)

  for (N in sample_sizes) {
    idx_sub <- round(seq(1, N_max, length.out = N)) |> as.integer()
    idx_sub[1] <- 1L
    idx_sub[length(idx_sub)] <- N_max
    X_sub <- X_master[idx_sub, , drop = FALSE]

    for (d_order in nlp_orders) {
      for (rep_id in seq_len(n_replicates)) {
        counter <- counter + 1L

        # resume check
        done <- nrow(
          completed_configs |>
            dplyr::filter(delta_true == delta, N == !!N, order == paste0("d", d_order), replicate == rep_id)
        ) >
          0
        if (done) {
          next
        }

        message(sprintf("[%d/%d] delta=%s N=%d d=%d rep=%d", counter, total_configs, delta, N, d_order, rep_id))

        # ---- generate one dataset
        set.seed(seed_base + rep_id + 1000 * d_order + 10000 * N + 100000 * delta)
        if (delta == 0) {
          y_sub <- rnorm(N, mean = 0, sd = sigma_fixed)
        } else {
          y_sub <- generate_gp_data(
            delta = delta,
            lambda = lambda_fixed,
            sigma = sigma_fixed,
            N = N,
            x = X_sub,
            n_samples = 1
          ) |>
            as.vector()
        }

        stan_data <- make_stan_data(X_sub, y_sub, d_order)

        # ---- fit alternative model with rstan
        fit_alt <- rstan::sampling(
          object = stan_mod_alt,
          data = stan_data,
          iter = iter_warmup + iter_sampling,
          warmup = iter_warmup,
          chains = chains,
          cores = n_cores_use,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
          seed = seed_base + rep_id,
          verbose = FALSE,
          refresh = 0
        )

        # ---- fit null model (delta = 0) with rstan
        fit_null <- rstan::sampling(
          object = stan_mod_null,
          data = stan_data,
          iter = iter_warmup + iter_sampling,
          warmup = iter_warmup,
          chains = chains,
          cores = n_cores_use,
          control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
          seed = seed_base + rep_id + 10000,
          verbose = FALSE,
          refresh = 0
        )

        # ---- extract posterior draws of delta from alternative model
        delta_draws <- rstan::extract(fit_alt, pars = "delta", inc_warmup = FALSE)$delta |>
          as.numeric()

        # ---- save posterior draws incrementally (separate dataframe)
        if (SAVE_POSTERIOR_DRAWS) {
          draw_index <- seq_along(delta_draws)
          keep <- ((draw_index - 1L) %% POSTERIOR_THIN) == 0L

          posterior_df <- data.frame(
            delta_true = delta,
            N = N,
            order = paste0("d", d_order),
            replicate = rep_id,
            draw_id = draw_index[keep],
            delta_draw = delta_draws[keep],
            stringsAsFactors = FALSE
          )

          append_csv(posterior_df, posterior_path)
        }

        # ---- BF via both methods
        bf_logspline <- bf10_interval_null_delta_moment(
          delta_draws = delta_draws,
          d_order = d_order,
          eps = eps_grid,
          post0_method = "logspline",
          return_log = TRUE
        ) |>
          dplyr::mutate(method = "logspline")

        bf_raw <- bf10_interval_null_delta_moment(
          delta_draws = delta_draws,
          d_order = d_order,
          eps = eps_grid,
          post0_method = "raw",
          return_log = TRUE
        ) |>
          dplyr::mutate(method = "raw")

        bf_df <- dplyr::bind_rows(bf_logspline, bf_raw) |>
          dplyr::mutate(
            delta_true = delta,
            N = N,
            order = paste0("d", d_order),
            replicate = rep_id
          )

        append_csv(bf_df, bf_interval_path)

        # ---- compute bridge sampling Bayes factor
        tryCatch(
          {
            bridge_alt <- bridgesampling::bridge_sampler(fit_alt)
            bridge_null <- bridgesampling::bridge_sampler(fit_null)

            log_bf10_bridge <- as.numeric(bridgesampling::bf(bridge_alt, bridge_null, log = TRUE)$bf)[1]
            bf10_bridge <- as.numeric(bridgesampling::bf(bridge_alt, bridge_null, log = FALSE)$bf)[1]

            bf_bridge_df <- data.frame(
              delta_true = delta,
              N = N,
              order = paste0("d", d_order),
              replicate = rep_id,
              log_bf10_bridge = log_bf10_bridge,
              bf10_bridge = bf10_bridge,
              bridge_alt_logml = as.numeric(bridge_alt$logml)[1],
              bridge_null_logml = as.numeric(bridge_null$logml)[1],
              stringsAsFactors = FALSE
            )

            append_csv(bf_bridge_df, bf_bridge_path)
          },
          error = function(e) {
            message(sprintf(
              "Bridge sampling failed for delta=%s N=%d d=%d rep=%d: %s",
              delta,
              N,
              d_order,
              rep_id,
              conditionMessage(e)
            ))
          }
        )

        # ---- save completion marker
        completed_row <- data.frame(
          delta_true = delta,
          N = N,
          order = paste0("d", d_order),
          replicate = rep_id,
          stringsAsFactors = FALSE
        )
        append_csv(completed_row, completed_path)

        # ---- memory cleanup (important for long Saturn jobs)
        rm(
          list = intersect(
            c(
              "stan_data",
              "y_sub",
              "fit_alt",
              "fit_null",
              "delta_draws",
              "bf_logspline",
              "bf_raw",
              "bf_df",
              "bridge_alt",
              "bridge_null",
              "bf_bridge_df",
              "posterior_df"
            ),
            ls()
          )
        )
        gc(FALSE)
      }
    }
  }
}

# ============================================
# SUMMARY TABLE
# ============================================

if (file.exists(bf_interval_path)) {
  bf_all <- readr::read_csv(bf_interval_path, show_col_types = FALSE)

  bf_summary <- bf_all |>
    dplyr::group_by(delta_true, N, order, method, epsilon) |>
    dplyr::summarise(
      median_log_bf10 = median(log_bf10, na.rm = TRUE),
      mean_log_bf10 = mean(log_bf10, na.rm = TRUE),
      sd_log_bf10 = sd(log_bf10, na.rm = TRUE),
      q05_log_bf10 = quantile(log_bf10, 0.05, na.rm = TRUE),
      q25_log_bf10 = quantile(log_bf10, 0.25, na.rm = TRUE),
      q75_log_bf10 = quantile(log_bf10, 0.75, na.rm = TRUE),
      q95_log_bf10 = quantile(log_bf10, 0.95, na.rm = TRUE),
      n_valid = sum(is.finite(log_bf10)),
      .groups = "drop"
    )

  readr::write_csv(bf_summary, summary_path)
}

# Return a useful object
list(
  bf_interval_file = bf_interval_path,
  bf_bridge_file = bf_bridge_path,
  posterior_file = posterior_path,
  summary_file = summary_path,
  completed_file = completed_path,
  chains = chains,
  cores_used = n_cores_use,
  cpu_available = cpu_available
)
