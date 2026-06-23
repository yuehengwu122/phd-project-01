# ============================================================
# Experiment: Diabetes — compare post0 methods for nonlinear SDR
#
# For fit_d0 (local prior, d_order = 0) on the diabetes dataset,
# compute log BF10 for each predictor's nonlinear effect (delta)
# using four posterior-at-null estimation methods at eps = 0.2:
#
#   1. logspline   — post_prob_delta_logspline()
#   2. raw (ecdf)  — post_prob_delta_ecdf()
#   3. laplace     — beta(k+1, n-k+1) mode: (k+1)/(n+2)
#   4. jeffreys    — beta(k+0.5, n-k+0.5) mode: (k+0.5)/(n+1)
#
# All four share the same prior0 = prior_prob_delta(d_order=0, eps).
# Output: one row per predictor, one column per method (log BF10).
# ============================================================

source("src/_setup.R")

EPS <- 0.2

# Load fit ---------------------------------------------------------------

fit_d0 <- readRDS("results/fit_03_diabetes_d0.rds")
fit_d0 <- fit_ig

post       <- rstan::extract(fit_d0$fit)
delta_mat  <- post$delta            # [S x P] matrix of posterior draws
d_order    <- fit_d0$d_order        # 0
pred_names <- fit_d0$predictor_names
P          <- fit_d0$stan_data$P

# Shared prior probability -----------------------------------------------

prior0 <- prior_prob_delta(d_order, EPS)

# compute_group_bf (laplace + jeffreys) ----------------------------------

compute_group_bf <- function(x, eps, prior0) {
  x      <- x[is.finite(x) & x >= 0]
  n      <- length(x)
  x_sort <- sort(x)
  k      <- findInterval(eps, x_sort)   # count of draws <= eps
  data.frame(
    k               = k,
    n               = n,
    prior0          = prior0,
    post0_laplace   = (k + 1)   / (n + 2),
    post0_jeffreys  = (k + 0.5) / (n + 1)
  )
}

# Build comparison table -------------------------------------------------

rows <- lapply(seq_len(P), function(j) {
  x <- delta_mat[, j]

  post0_ls  <- post_prob_delta_logspline(x, EPS)
  post0_raw <- post_prob_delta_ecdf(x, EPS)$post0
  gb        <- compute_group_bf(x, EPS, prior0)

  data.frame(
    predictor          = pred_names[j],
    log_bf10_logspline = round(log(prior0 / post0_ls),              4),
    log_bf10_raw       = round(log(prior0 / post0_raw),             4),
    log_bf10_laplace   = round(log(prior0 / gb$post0_laplace),      4),
    log_bf10_jeffreys  = round(log(prior0 / gb$post0_jeffreys),     4)
  )
})

tbl <- do.call(rbind, rows)

# Print ------------------------------------------------------------------

cat(sprintf("\neps = %.1f,  prior0 = %.6f,  d_order = %d\n\n", EPS, prior0, d_order))
print(tbl, row.names = FALSE)
