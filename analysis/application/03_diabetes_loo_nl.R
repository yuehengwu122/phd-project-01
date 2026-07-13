# ============================================================
# Application 3 (BF): Diabetes
# Leave-one-nonlinear-predictor-out via bridge sampling
# ============================================================

source("src/_setup.R")

# --- Data ---
library(care)
data(efron2004, package = "care")

x_mat <- as.matrix(efron2004$x)
class(x_mat) <- "matrix"
X <- x_mat[, -2]         # drop binary "sex"
y <- efron2004$y[, 1]

P <- ncol(X)
pred_names <- colnames(X)
cat("N =", nrow(X), ", P =", P, "\n")
cat("Predictors:", paste(pred_names, collapse = ", "), "\n\n")

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# --- Common MCMC controls ---
mcmc <- list(
  chains = 2, iter = 5500, warmup = 500,
  adapt_delta = 0.95, max_treedepth = 12
)

D_ORDER <- 2L

# ============================================================
# Step 1: Full reference model
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
saveRDS(fit_full, "results/fit_03_diabetes_full_select.rds")

bridge_full <- bridgesampling::bridge_sampler(fit_full$fit, silent = TRUE)
cat("logml(full) =", bridge_full$logml, "\n\n")

# ============================================================
# Step 2: Drop-one-nonlinear sweep
# ============================================================
results <- data.frame(
  predictor   = pred_names,
  logml_drop  = NA_real_,
  log_BF10    = NA_real_,
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
  saveRDS(
    fit_j,
    sprintf("results/fit_03_diabetes_drop_nl_%02d_%s.rds", j, pred_names[j])
  )

  bridge_j <- bridgesampling::bridge_sampler(fit_j$fit, silent = TRUE)
  log_bf <- bridgesampling::bf(bridge_full, bridge_j, log = TRUE)$bf

  results$logml_drop[j] <- bridge_j$logml
  results$log_BF10[j]   <- log_bf
  results$BF10[j]       <- exp(log_bf)

  saveRDS(fits_drop, "results/fit_03_diabetes_drop_nl_list.rds")
  saveRDS(results,   "results/bf_03_diabetes_nl_loo.rds")

  cat(sprintf("   logml(drop_%s) = %.3f   log BF10 = %.3f\n\n",
              pred_names[j], bridge_j$logml, log_bf))
}

saveRDS(fits_drop, "results/fit_03_diabetes_drop_nl_list.rds")
saveRDS(results,   "results/bf_03_diabetes_nl_loo.rds")

cat("\n================ Leave-one-nonlinear-out BF table ================\n")
print(results, row.names = FALSE)
