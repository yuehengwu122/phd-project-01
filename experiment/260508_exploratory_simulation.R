# ============================================================
# Exploratory Simulation: Moment vs Inverse-Moment Priors
#
# Design:
#   delta_true = 0, 1, 2, 3
#   N = 20, 50, 100, 200
#   6 priors: M(d=0), M(d=1), M(d=2), M(d=3), IM(d=1), IM(d=2)
#   1 replicate per (delta, N) cell
#
# Data generation:
#   Master grid of 1000 points on [-0.5, 0.5].
#   N=200 is a subset of the 1000, N=100 subset of N=200, etc.
#   Identical to the main simulation approach.
#
# Usage:
#   source("src/_setup.R")
#   source("analysis/exploratory_simulation.R")
# ============================================================

library(bridgesampling)
library(dplyr)
library(readr)

# ============================================================
# EXPERIMENT SETTINGS
# ============================================================

target_deltas    <- c(0, 1, 2, 3)
target_Ns        <- c(20L, 50L, 100L, 200L)
all_sample_sizes <- c(20, 50, 100, 200)    # used for nested index construction

# --- Moment priors (existing) ---
moment_orders <- 0:3

# --- Inverse-moment priors (new) ---
IM_NU    <- 2           # tail parameter (from exploration)
im_orders <- c(1, 2)    # d values for IM

# --- Shared simulation parameters ---
lambda_fixed  <- 0.3
sigma_fixed   <- 0.1
x_min         <- -0.5
x_max         <- 0.5

seed_base     <- 2026

# --- MCMC settings ---
iter_warmup   <- 500
iter_sampling <- 1500
adapt_delta   <- 0.95
max_treedepth <- 12
chains        <- 2L

SAVE_POSTERIOR_DRAWS <- TRUE
POSTERIOR_THIN       <- 5L

# ============================================================
# PATHS
# ============================================================

dir.create("results/exploratory", recursive = TRUE, showWarnings = FALSE)

bf_bridge_path   <- "results/exploratory/bf_bridge.csv"
posterior_path   <- "results/exploratory/posterior_draws.csv"

stan_file_moment  <- "models/fit_full_model.stan"
stan_file_imoment <- "models/fit_full_model_imoment.stan"
stan_file_null    <- "models/fit_linear_model.stan"

if (!file.exists(stan_file_moment))  stop("Stan file not found: ", stan_file_moment)
if (!file.exists(stan_file_imoment)) stop("Stan file not found: ", stan_file_imoment)
if (!file.exists(stan_file_null))    stop("Stan null file not found: ", stan_file_null)


# ============================================================
# INVERSE-MOMENT PRIOR FUNCTIONS
# ============================================================

#' Unnormalised inverse-moment density on delta > 0
imoment_unnorm <- function(delta, d, tau, nu = IM_NU) {
  delta <- pmax(delta, 1e-300)
  log_dens <- -(nu + 1) * log(delta) - (tau / delta^2)^d
  exp(log_dens)
}

#' Normalising constant
imoment_norm_const <- function(d, tau, nu = IM_NU) {
  integrate(function(x) imoment_unnorm(x, d, tau, nu),
            lower = 1e-10, upper = Inf,
            rel.tol = 1e-10, subdivisions = 1000)$value
}

#' Normalised density
imoment_density <- function(delta, d, tau, nu = IM_NU, C_norm = NULL) {
  if (is.null(C_norm)) C_norm <- imoment_norm_const(d, tau, nu)
  imoment_unnorm(delta, d, tau, nu) / C_norm
}

#' CDF: Pr(delta < eps)
imoment_cdf <- function(eps, d, tau, nu = IM_NU, C_norm = NULL) {
  if (eps <= 0) return(0)
  if (is.null(C_norm)) C_norm <- imoment_norm_const(d, tau, nu)
  integrate(function(x) imoment_unnorm(x, d, tau, nu) / C_norm,
            lower = 1e-10, upper = eps,
            rel.tol = 1e-10, subdivisions = 1000)$value
}

#' Quantile by numerical inversion
imoment_quantile <- function(p, d, tau, nu = IM_NU) {
  C_norm <- imoment_norm_const(d, tau, nu)
  upper <- 1
  while (imoment_cdf(upper, d, tau, nu, C_norm) < p) upper <- upper * 2
  uniroot(function(q) imoment_cdf(q, d, tau, nu, C_norm) - p,
          lower = 1e-10, upper = upper, tol = 1e-10)$root
}

#' Find tau for median = target
find_tau_for_median <- function(d, target_median = 2, nu = IM_NU,
                                tau_range = c(1e-6, 1e4)) {
  objective <- function(log_tau) {
    tau <- exp(log_tau)
    imoment_quantile(0.5, d, tau, nu) - target_median
  }
  exp(uniroot(objective, interval = log(tau_range), tol = 1e-8)$root)
}

#' Prepare Stan data for IM prior
imoment_stan_data <- function(d, tau, nu = IM_NU) {
  C_norm <- imoment_norm_const(d, tau, nu)
  list(
    d_order       = as.integer(d),
    tau_delta     = tau,
    nu_delta      = nu,
    logC_imoment  = log(C_norm)
  )
}


# ============================================================
# CALIBRATE IM PRIORS
# ============================================================

cat("Calibrating inverse-moment priors (nu =", IM_NU, ")...\n")
im_tau_cache <- setNames(numeric(length(im_orders)), paste0("im", im_orders))
for (d in im_orders) {
  tau_d <- find_tau_for_median(d, target_median = 2)
  im_tau_cache[paste0("im", d)] <- tau_d
  med_check <- imoment_quantile(0.5, d, tau_d)
  cat(sprintf("  IM(d=%d): tau = %.4f, median = %.4f\n", d, tau_d, med_check))
}


# ============================================================
# HELPERS
# ============================================================

make_stan_data_alt <- function(X, y, d_order) {
  X <- as.matrix(X); y <- as.vector(y)
  c(list(N = nrow(X), P = ncol(X), X = X, y = y), moment_prior_stan_data(d_order))
}

make_stan_data_imoment <- function(X, y, d, tau) {
  X <- as.matrix(X); y <- as.vector(y)
  c(list(N = nrow(X), P = ncol(X), X = X, y = y), imoment_stan_data(d, tau))
}

make_stan_data_null <- function(X, y) {
  X <- as.matrix(X); y <- as.vector(y)
  list(N = nrow(X), P = ncol(X), X = X, y = y)
}

append_csv <- function(df, path) {
  if (file.exists(path)) readr::write_csv(df, path, append = TRUE)
  else                   readr::write_csv(df, path)
}


# ============================================================
# COMPILE MODELS
# ============================================================

message("Compiling moment-prior Stan model...")
stan_mod_moment  <- rstan::stan_model(stan_file_moment)
message("Compiling inverse-moment-prior Stan model...")
stan_mod_imoment <- rstan::stan_model(stan_file_imoment)
message("Compiling null Stan model...")
stan_mod_null    <- rstan::stan_model(stan_file_null)


# ============================================================
# BUILD PRIOR CONFIG TABLE
# ============================================================

prior_configs <- dplyr::bind_rows(
  # Moment priors
  lapply(moment_orders, function(d) {
    data.frame(family = "M", d = d, label = paste0("M(d=", d, ")"),
               stringsAsFactors = FALSE)
  }),
  # Inverse-moment priors
  lapply(im_orders, function(d) {
    data.frame(family = "IM", d = d, label = paste0("IM(d=", d, ")"),
               stringsAsFactors = FALSE)
  })
)

cat("\nPrior configurations:\n")
print(prior_configs, row.names = FALSE)
cat("\n")


# ============================================================
# MAIN LOOP
# ============================================================

N_master <- 1000L
x_master_vec <- seq(x_min, x_max, length.out = N_master)
X_master <- matrix(x_master_vec, ncol = 1)

# Build nested index list (largest N first, then subsample)
sorted_Ns <- sort(all_sample_sizes)
idx_list <- vector("list", length(sorted_Ns))
names(idx_list) <- as.character(sorted_Ns)
idx_list[[as.character(max(sorted_Ns))]] <-
  as.integer(round(seq(1L, N_master, length.out = max(sorted_Ns))))
for (k in rev(seq_len(length(sorted_Ns) - 1L))) {
  parent <- idx_list[[as.character(sorted_Ns[k + 1L])]]
  pick <- as.integer(round(seq(1L, length(parent), length.out = sorted_Ns[k])))
  idx_list[[as.character(sorted_Ns[k])]] <- parent[pick]
}

total_fits <- length(target_deltas) * length(target_Ns) * nrow(prior_configs)
counter <- 0L

for (delta in target_deltas) {

  # ---- Generate master data ----
  set.seed(seed_base + 100000L * delta)
  if (delta == 0) {
    y_master <- rnorm(N_master, mean = 0, sd = sigma_fixed)
  } else {
    y_master <- generate_gp_data(
      delta     = delta,
      lambda    = lambda_fixed,
      sigma     = sigma_fixed,
      N         = N_master,
      x         = as.numeric(X_master),
      n_samples = 1
    ) |> as.vector()
  }

  for (N in target_Ns) {
    idx_sub <- idx_list[[as.character(N)]]
    X_sub <- X_master[idx_sub, , drop = FALSE]
    y_sub <- y_master[idx_sub]

    # ---- Fit null model (once per (delta, N)) ----
    message(sprintf("\n=== delta=%d, N=%d: fitting null model ===", delta, N))
    fit_null <- rstan::sampling(
      object = stan_mod_null,
      data = make_stan_data_null(X_sub, y_sub),
      iter = iter_warmup + iter_sampling,
      warmup = iter_warmup,
      chains = chains,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      seed = seed_base + 10000L,
      verbose = FALSE,
      refresh = 0,
      open_progress = FALSE
    )

    bridge_null <- tryCatch(
      bridgesampling::bridge_sampler(fit_null),
      error = function(e) { message("  Bridge null failed: ", e$message); NULL }
    )

    # ---- Loop over priors ----
    for (pc in seq_len(nrow(prior_configs))) {
      fam   <- prior_configs$family[pc]
      d_ord <- prior_configs$d[pc]
      lab   <- prior_configs$label[pc]

      counter <- counter + 1L
      message(sprintf("[%d/%d] delta=%d N=%d %s", counter, total_fits, delta, N, lab))

      # ---- Fit alternative model ----
      if (fam == "M") {
        stan_data <- make_stan_data_alt(X_sub, y_sub, d_ord)
        stan_mod  <- stan_mod_moment
      } else {
        tau_d <- im_tau_cache[paste0("im", d_ord)]
        stan_data <- make_stan_data_imoment(X_sub, y_sub, d_ord, tau_d)
        stan_mod  <- stan_mod_imoment
      }

      fit_alt <- rstan::sampling(
        object = stan_mod,
        data = stan_data,
        iter = iter_warmup + iter_sampling,
        warmup = iter_warmup,
        chains = chains,
        control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
        seed = seed_base,
        verbose = FALSE,
        refresh = 0,
        open_progress = FALSE
      )

      delta_draws <- rstan::extract(fit_alt, pars = "delta")$delta |> as.numeric()

      # ---- Save posterior draws ----
      if (SAVE_POSTERIOR_DRAWS) {
        draw_index <- seq_along(delta_draws)
        keep <- ((draw_index - 1L) %% POSTERIOR_THIN) == 0L
        append_csv(
          data.frame(
            delta_true = delta, N = N, prior = lab,
            draw_id = draw_index[keep], delta_draw = delta_draws[keep]
          ),
          posterior_path
        )
      }

      # ---- Bridge sampling BF ----
      tryCatch({
        bridge_alt <- bridgesampling::bridge_sampler(fit_alt)
        if (!is.null(bridge_null)) {
          bf_result <- bridgesampling::bf(bridge_alt, bridge_null, log = TRUE)
          append_csv(
            data.frame(
              delta_true = delta, N = N, prior = lab,
              log_bf10_bridge   = as.numeric(bf_result$bf)[1],
              bridge_alt_logml  = as.numeric(bridge_alt$logml)[1],
              bridge_null_logml = as.numeric(bridge_null$logml)[1]
            ),
            bf_bridge_path
          )
        }
      }, error = function(e) {
        message(sprintf("  Bridge sampling failed for %s: %s", lab, e$message))
      })

      rm(fit_alt, delta_draws); gc(FALSE)
    } # end prior

    rm(fit_null, bridge_null); gc(FALSE)
  } # end N

  rm(y_master); gc(FALSE)
} # end delta

cat("\n============================================================\n")
cat("Simulation complete.\n")
cat("Results saved to:\n")
cat("  ", bf_bridge_path, "\n")
cat("  ", posterior_path, "\n")
cat("============================================================\n")
cat("Interval-null BFs can be computed from posterior draws later.\n")
