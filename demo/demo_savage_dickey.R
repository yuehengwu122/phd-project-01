# ============================================================
# Demo: Savage-Dickey Ratio Computation for GP Model Selection
# ============================================================
#
#   1. Using the helper function compute_savage_dickey_ratios()
#   2. Manual computation step-by-step
#
# SDR = prior_density(0) / posterior_density(0)
#   - SDR > 1: evidence FOR the null (no effect)
#   - SDR < 1: evidence AGAINST the null (effect exists)
#
# Two priors are compared:
#   - p1: delta ~ half-t(4, 2.7)    [original prior]
#   - p2: delta2 ~ half-t(1.5, 6)   [squared parameterization]
# ============================================================

source("_common.R")

# ============================================================
# Part 1: Simulate Data with Known GP Effect
# ============================================================
N <- 50
X <- seq(-0.5, 0.5, length = N) |> matrix(ncol = 1)

# Generate data with true delta = 2 (nonlinear effect exists)
y <- generate_gp_data(
  delta = 2,
  lambda = 0.3,
  sigma = 0.1,
  N = 50,
  x = X,
  n_samples = 1
) |>
  as.vector()

# ============================================================
# Part 2: Fit GP Models with Two Different Priors
# ============================================================

# Fit with p1 prior (delta ~ half-t(4, 2.7))
fit_gp_p1 <- fit_model(
  gp_add_p1_15,
  X = X,
  y = y,
  chains = 2,
  iter = 2000,
  warmup = 500
)

# Fit with p2 prior (delta2 ~ half-t(1.5, 6))
fit_gp_p2 <- fit_model(
  gp_add_p2_15,
  X = X,
  y = y,
  chains = 2,
  iter = 2000,
  warmup = 500,
  delta_prior = "p2"
)

# View summaries (includes SDR in logSDR10 column)
summary_gp_fit(fit_gp_p1)
summary_gp_fit(fit_gp_p2)

# ============================================================
# Part 3: Compute SDR Using Helper Function
# ============================================================
sdr_p1 <- compute_savage_dickey_ratios(fit_gp_p1, delta_prior = "p1")
sdr_p2 <- compute_savage_dickey_ratios(fit_gp_p2, delta_prior = "p2")

# ============================================================
# Part 4: Manual SDR Computation (Step-by-Step)
# ============================================================

# --- For p1 prior: delta ~ half-t(4, 2.7) ---
post_p1 <- rstan::extract(fit_gp_p1$fit)
delta_p1 <- post_p1$delta[, 1] # posterior draws of delta

# Prior density at 0: half-t(df=4, scale=2.7)
# Formula: 2 * dt(x/scale, df) / scale
prior_delta_at_0 <- 2 * dt(0 / 2.7, df = 4) / 2.7

# Posterior density at 0: estimated via logspline (lbound=0 for half distribution)
post_delta_at_0 <- logspline::dlogspline(
  0,
  logspline::logspline(delta_p1, lbound = 0)
)

# SDR = prior(0) / posterior(0)
sdr_nonlinear_p1 <- prior_delta_at_0 / post_delta_at_0


# --- For p2 prior: delta2 ~ half-t(1.5, 6) ---
post_p2 <- rstan::extract(fit_gp_p2$fit)
delta2_p2 <- post_p2$delta2[, 1] # posterior draws of delta^2

# Prior density at 0: half-t(df=1.5, scale=6)
prior_delta2_at_0 <- 2 * dt(0 / 6, df = 1.5) / 6

# Posterior density at 0
post_delta2_at_0 <- logspline::dlogspline(
  0,
  logspline::logspline(delta2_p2, lbound = 0)
)

# SDR
sdr_nonlinear_p2 <- prior_delta2_at_0 / post_delta2_at_0

# ============================================================
# Part 5: Application to Real Data (Diabetes Dataset)
# ============================================================
fit_a3_dep <- readRDS("analysis/application/fits/fit_a3_dep.rds")
fit_a3_dep_p2 <- readRDS("analysis/application/fits/fit_a3_dep_p2.rds")

# Compute SDR for all predictors
compute_savage_dickey_ratios(fit_a3_dep)
compute_savage_dickey_ratios(fit_a3_dep_p2, delta_prior = "p2")

# --- Manual SDR for s5 predictor (8th column) ---

# p1 prior
post_a3_p1 <- rstan::extract(fit_a3_dep$fit)
delta_s5_p1 <- post_a3_p1$delta[, 8]
post_delta_at_0 <- logspline::dlogspline(
  0,
  logspline::logspline(delta_s5_p1, lbound = 0, max = max(delta_s5_p1))
)
sdr_s5_p1 <- prior_delta_at_0 / post_delta_at_0

# p2 prior
post_a3_p2 <- rstan::extract(fit_a3_dep_p2$fit)
delta_s5_p2 <- post_a3_p2$delta2[, 8]
post_delta2_at_0 <- logspline::dlogspline(
  0,
  logspline::logspline(delta_s5_p2, lbound = 0, max = max(delta_s5_p2))
)
sdr_s5_p2 <- prior_delta2_at_0 / post_delta2_at_0

# Note: If logspline warns about "maxknots reduced", can trim heavy tail:
delta_s5_p2_trimmed <- delta_s5_p2[delta_s5_p2 <= quantile(delta_s5_p2, 0.90)]
post_delta2_at_0 <- logspline::dlogspline(
  0,
  logspline::logspline(delta_s5_p2_trimmed, lbound = 0)
)
sdr_s5_p2_trimmed <- prior_delta2_at_0 / post_delta2_at_0
