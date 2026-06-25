#' Job 04: N={20,50,100,200}, delta=1, local prior only with bridge BF

source("src/_setup.R")

library(bridgesampling)
library(dplyr)
library(readr)

target_delta   <- 1L
target_Ns      <- c(20L, 50L, 100L, 200L)
target_d_order <- 0L
order_label    <- "d0"
job_tag        <- "job04_N_small_delta1_local"

cpu_available <- parallel::detectCores(logical = TRUE)
n_jobs        <- 4L
reserve_cores <- 2L
cores_budget  <- max(2L, cpu_available - reserve_cores)
n_cores_use   <- max(1L, floor(cores_budget / n_jobs))
chains        <- min(2L, n_cores_use)

options(mc.cores = n_cores_use)

message(sprintf("[%s] CPUs: %d | chains=%d | cores/job=%d",
                job_tag, cpu_available, chains, n_cores_use))

sample_sizes   <- c(20L, 50L, 100L, 200L)
master_sample_size <- 1000L
n_replicates   <- 50L

lambda_fixed   <- 0.3
sigma_fixed    <- 0.1
x_min          <- -0.5
x_max          <- 0.5

iter_warmup    <- 500L
iter_sampling  <- 1500L
adapt_delta    <- 0.95
max_treedepth  <- 12L
seed_base      <- 123L

SAVE_POSTERIOR_DRAWS <- TRUE

dir.create("data/simulated", recursive = TRUE, showWarnings = FALSE)

bf_bridge_path   <- paste0("data/simulated/results_sim3_rstan_bf_bridge_", job_tag, ".csv")
posterior_path   <- paste0("data/simulated/results_sim3_rstan_posterior_draws_", job_tag, ".csv")
completed_path   <- paste0("data/simulated/results_sim3_rstan_completed_", job_tag, ".csv")

stan_file      <- "models/fit_full_model.stan"
stan_file_null <- "models/fit_linear_model.stan"

if (!file.exists(stan_file)) stop("Stan file not found: ", stan_file)
if (!file.exists(stan_file_null)) stop("Stan null model file not found: ", stan_file_null)

make_stan_data_alt <- function(X, y, d_order = target_d_order) {
  X <- as.matrix(X)
  y <- as.vector(y)
  stopifnot(nrow(X) == length(y))
  c(list(N = nrow(X), P = ncol(X), X = X, y = y), moment_prior_stan_data(d_order))
}

make_stan_data_null <- function(X, y) {
  X <- as.matrix(X)
  y <- as.vector(y)
  list(N = nrow(X), P = ncol(X), X = X, y = y)
}

append_csv <- function(df, path) {
  if (file.exists(path)) readr::write_csv(df, path, append = TRUE)
  else readr::write_csv(df, path)
}

completed_configs <- data.frame(
  delta_true = numeric(), N = integer(), order = character(), replicate = integer(),
  stringsAsFactors = FALSE
)
if (file.exists(completed_path)) {
  completed_configs <- readr::read_csv(completed_path, show_col_types = FALSE)
}

message("Compiling alternative Stan model...")
stan_mod_alt <- rstan::stan_model(stan_file)
message("Compiling null Stan model...")
stan_mod_null <- rstan::stan_model(stan_file_null)

total_configs <- length(target_Ns) * n_replicates
counter <- 0L

N_master <- master_sample_size
x_master_vec <- seq(x_min, x_max, length.out = N_master)
X_master <- matrix(x_master_vec, ncol = 1)

sorted_Ns <- sort(sample_sizes)
idx_list <- vector("list", length(sorted_Ns))
names(idx_list) <- as.character(sorted_Ns)
idx_list[[as.character(max(sorted_Ns))]] <-
  as.integer(round(seq(1L, N_master, length.out = max(sorted_Ns))))
for (k in rev(seq_len(length(sorted_Ns) - 1L))) {
  parent <- idx_list[[as.character(sorted_Ns[k + 1L])]]
  pick <- as.integer(round(seq(1L, length(parent), length.out = sorted_Ns[k])))
  idx_list[[as.character(sorted_Ns[k])]] <- parent[pick]
}

for (rep_id in seq_len(n_replicates)) {
  set.seed(seed_base + rep_id + 100000L * target_delta)
  y_master <- generate_gp_data(
    delta = target_delta,
    lambda = lambda_fixed,
    sigma = sigma_fixed,
    N = N_master,
    x = as.numeric(X_master),
    n_samples = 1L,
    seed = seed_base + rep_id + 100000L * target_delta
  ) |> as.vector()

  for (N in target_Ns) {
    idx_sub <- idx_list[[as.character(N)]]
    X_sub <- X_master[idx_sub, , drop = FALSE]
    y_sub <- y_master[idx_sub]
    counter <- counter + 1L

    done <- nrow(
      completed_configs |>
        dplyr::filter(
          delta_true == target_delta,
          N == !!N,
          order == order_label,
          replicate == rep_id
        )
    ) > 0
    if (done) next

    message(sprintf("[%s] [%d/%d] delta=%s N=%d d=%d rep=%d",
                    job_tag, counter, total_configs, target_delta, N, target_d_order, rep_id))

    fit_alt <- rstan::sampling(
      object = stan_mod_alt,
      data = make_stan_data_alt(X_sub, y_sub),
      iter = iter_warmup + iter_sampling,
      warmup = iter_warmup,
      chains = chains,
      cores = n_cores_use,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      seed = seed_base + rep_id,
      verbose = FALSE,
      refresh = 0
    )

    fit_null <- rstan::sampling(
      object = stan_mod_null,
      data = make_stan_data_null(X_sub, y_sub),
      iter = iter_warmup + iter_sampling,
      warmup = iter_warmup,
      chains = chains,
      cores = n_cores_use,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      seed = seed_base + rep_id + 10000L,
      verbose = FALSE,
      refresh = 0
    )

    delta_draws <- rstan::extract(fit_alt, pars = "delta", inc_warmup = FALSE)$delta |>
      as.numeric()

    if (SAVE_POSTERIOR_DRAWS) {
      append_csv(
        data.frame(
          delta_true = target_delta,
          N = N,
          order = order_label,
          replicate = rep_id,
          draw_id = seq_along(delta_draws),
          delta_draw = delta_draws,
          stringsAsFactors = FALSE
        ),
        posterior_path
      )
    }

    tryCatch({
      bridge_alt <- bridgesampling::bridge_sampler(fit_alt)
      bridge_null <- bridgesampling::bridge_sampler(fit_null)
      append_csv(
        data.frame(
          delta_true = target_delta,
          N = N,
          order = order_label,
          replicate = rep_id,
          log_bf10_bridge = as.numeric(bridgesampling::bf(bridge_alt, bridge_null, log = TRUE)$bf)[1],
          bf10_bridge = as.numeric(bridgesampling::bf(bridge_alt, bridge_null, log = FALSE)$bf)[1],
          bridge_alt_logml = as.numeric(bridge_alt$logml)[1],
          bridge_null_logml = as.numeric(bridge_null$logml)[1],
          stringsAsFactors = FALSE
        ),
        bf_bridge_path
      )
    }, error = function(e) {
      message(sprintf("Bridge sampling failed for delta=%s N=%d d=%d rep=%d: %s",
                      target_delta, N, target_d_order, rep_id, conditionMessage(e)))
    })

    append_csv(
      data.frame(
        delta_true = target_delta,
        N = N,
        order = order_label,
        replicate = rep_id,
        stringsAsFactors = FALSE
      ),
      completed_path
    )

    rm(list = intersect(c("fit_alt", "fit_null", "delta_draws", "bridge_alt", "bridge_null"), ls()))
    gc(FALSE)
  }

  rm(y_master)
  gc(FALSE)
}

list(
  bf_bridge_file = bf_bridge_path,
  posterior_file = posterior_path,
  completed_file = completed_path,
  chains = chains,
  cores_used = n_cores_use,
  cpu_available = cpu_available
)
