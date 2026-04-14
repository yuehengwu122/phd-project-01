# ============================================================
# Application 4 (BF): Boston Housing
# Leave-one-nonlinear-predictor-out via bridge sampling
#
# For each predictor j in 1..P, fit a model that includes the
# full linear part for all predictors but drops the GP component
# of predictor j. Compare to the full model (all GP components)
# via bridge sampling:  BF_j = m(full) / m(drop_j).
#
# A LARGE BF_j means predictor j's nonlinear component is
# supported by the data. We only do the nonlinear sweep here
# because (a) each fit is expensive and (b) the linear BFs are
# already well handled by the Savage-Dickey approach in
# compute_bf_table().
# ============================================================

source("src/_setup.R")

# --- Data ---
library(MASS)
data(Boston)

X <- as.matrix(Boston[, -c(4, 9, 14)])   # drop chas, rad, medv
y <- Boston$medv

P <- ncol(X)
pred_names <- colnames(X)
cat("N =", nrow(X), ", P =", P, "\n")
cat("Predictors:", paste(pred_names, collapse = ", "), "\n\n")

# --- Common MCMC controls ---
mcmc <- list(
  chains = 2, iter = 5500, warmup = 500,
  adapt_delta = 0.95, max_treedepth = 12
)

# --- Choose moment-prior order for the BF sweep ---
# Use d_order = 2 (non-local) so the BF sweep is consistent with
# the non-local analysis in the main application. Change if desired.
D_ORDER <- 2L

# ============================================================
# Step 1: Full reference model (all nonlinear components active)
# ============================================================
cat("---- Fitting FULL reference model ----\n")
fit_full <- fit_model_select(
  X, y,
  include_linear    = rep(1L, P),
  include_nonlinear = rep(1L, P),
  d_order = D_ORDER,
  chains = mcmc$chains, iter = mcmc$iter, warmup = mcmc$warmup,
  adapt_delta = mcmc$adapt_delta, max_treedepth = mcmc$max_treedepth
)
saveRDS(fit_full, "results/fit_04_boston_full_select.rds")

bridge_full <- bridgesampling::bridge_sampler(fit_full$fit, silent = TRUE)
cat("logml(full) =", bridge_full$logml, "\n\n")

# ============================================================
# Step 2: Loop over predictors, dropping the GP component of j
# ============================================================
results <- data.frame(
  predictor   = pred_names,
  logml_drop  = NA_real_,
  log_BF10    = NA_real_,   # full vs drop_j
  BF10        = NA_real_,
  stringsAsFactors = FALSE
)

fits_drop <- vector("list", P)
names(fits_drop) <- pred_names

for (j in seq_len(P)) {
  cat(sprintf("---- Fitting DROP-NL[%d] = %s ----\n", j, pred_names[j]))

  mask_nl <- rep(1L, P); mask_nl[j] <- 0L

  fit_j <- fit_model_select(
    X, y,
    include_linear    = rep(1L, P),
    include_nonlinear = mask_nl,
    d_order = D_ORDER,
    chains = mcmc$chains, iter = mcmc$iter, warmup = mcmc$warmup,
    adapt_delta = mcmc$adapt_delta, max_treedepth = mcmc$max_treedepth
  )
  fits_drop[[j]] <- fit_j

  bridge_j <- bridgesampling::bridge_sampler(fit_j$fit, silent = TRUE)
  log_bf <- bridgesampling::bf(bridge_full, bridge_j, log = TRUE)$bf

  results$logml_drop[j] <- bridge_j$logml
  results$log_BF10[j]   <- log_bf
  results$BF10[j]       <- exp(log_bf)

  cat(sprintf("   logml(drop_%s) = %.3f   log BF10 = %.3f\n\n",
              pred_names[j], bridge_j$logml, log_bf))
}

# --- Save everything ---
saveRDS(fits_drop, "results/fit_04_boston_drop_nl_list.rds")
saveRDS(results,   "results/bf_04_boston_nl_loo.rds")

cat("\n================ Leave-one-nonlinear-out BF table ================\n")
print(results, row.names = FALSE)
