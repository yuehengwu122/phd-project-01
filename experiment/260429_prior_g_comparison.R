# ============================================================
# Experiment: Compare lognormal-g vs IG-g prior on theta
#
# Fits the same 2D simulated dataset under four models:
#   Model A:    fit_full_model.stan         (log_g ~ N(log N, 0.5))
#   Model B0.5: Stan_code_for_full_model_IG (g ~ IG(0.5, N/2), Zellner-Siow)
#   Model B1.0: Stan_code_for_full_model_IG (g ~ IG(1.0, N/2), Moderate)
#   Model B2.0: Stan_code_for_full_model_IG (g ~ IG(2.0, N/2), Weakly informative)
#
# For each model and each predictor j, compute:
#   prior0[j]       = marginal prior density at theta_j = 0
#   post0[j]        = posterior density at theta_j = 0 (logspline)
#   log_SD_ratio[j] = log(prior0 / post0)  (Savage-Dickey log BF10)
#
# X1 has a true linear effect; X2 has a true nonlinear effect.
# We expect the SD ratio to favour H1 for X1 and H0 for X2.
# ============================================================

source("src/_setup.R")

SEED    <- 42L
IG_VALS <- c(0.5, 1.0, 2.0)   # ig_a values to compare for Model B

# ── 1. Simulate data ──────────────────────────────────────────────────────────

set.seed(SEED)

N_sim <- 100L
P_sim <- 2L

X_sim <- matrix(runif(N_sim * P_sim, -2, 2), N_sim, P_sim)
colnames(X_sim) <- c("X1_linear", "X2_nonlinear")

f1 <- function(x) 0 * x
f2 <- function(x) sin(2 * x)

y_sim <- f1(X_sim[, 1]) + f2(X_sim[, 2]) + rnorm(N_sim, 0, 0.3)

X_sim_scaled <- apply(X_sim, 2, rescale_to_unit)
colnames(X_sim_scaled) <- colnames(X_sim)

# ── 2. Build shared Stan data ─────────────────────────────────────────────────

d_order <- 0L
mp      <- moment_prior_stan_data(d_order)

stan_data_base <- c(
  list(N = N_sim, P = P_sim, X = X_sim_scaled, y = y_sim),
  mp
)

ig_b <- 0.5 * N_sim   # shared across all IG variants (Zellner-Siow scale)

# ── 3. Compile Stan models (once each) ───────────────────────────────────────

message("Compiling Model A ...")
mod_A <- rstan::stan_model(file = "models/fit_full_model.stan")

message("Compiling Model B ...")
mod_B <- rstan::stan_model(file = "models/fit_full_model_IG.stan ")

# ── 4. Helper: fit one model ──────────────────────────────────────────────────

run_sampling <- function(mod, data, seed = SEED) {
  rstan::sampling(
    object        = mod,
    data          = data,
    iter          = 2000L,
    warmup        = 500L,
    chains        = 2L,
    control       = list(adapt_delta = 0.95, max_treedepth = 12),
    seed          = seed,
    refresh       = 500L,
    verbose       = FALSE,
    open_progress = FALSE
  )
}

# ── 5. Fit all models ─────────────────────────────────────────────────────────

message("\n=== Fitting Model A: lognormal-g ===\n")
fit_A <- run_sampling(mod_A, stan_data_base)

fits_B <- setNames(
  lapply(IG_VALS, function(a) {
    message(sprintf("\n=== Fitting Model B: IG-g (ig_a = %.1f) ===\n", a))
    run_sampling(mod_B, c(stan_data_base, list(ig_a = a, ig_b = ig_b)))
  }),
  paste0("B_ig", IG_VALS)
)

# ── 6. Convergence diagnostics ────────────────────────────────────────────────

print_convergence <- function(fit, label) {
  diag       <- extract_diagnostics(fit)
  summary_df <- rstan::summary(fit)$summary
  pars <- c("theta[1]", "theta[2]", "delta[1]", "delta[2]", "sigma", "log_g_theta")
  pars <- pars[pars %in% rownames(summary_df)]

  cat(strrep("-", 55), "\n")
  cat(label, "\n")
  cat(strrep("-", 55), "\n")
  elapsed    <- rstan::get_elapsed_time(fit)
  total_sec  <- max(rowSums(elapsed))
  time_str   <- if (total_sec > 3600)
    sprintf("%.2f hours", total_sec / 3600)
  else if (total_sec > 60)
    sprintf("%.1f minutes", total_sec / 60)
  else
    sprintf("%.1f seconds", total_sec)

  cat(sprintf(
    "  Fit time:         %s\n  Divergences:      %d\n  Max treedepth:    %d\n  Mean accept stat: %.3f\n",
    time_str, diag$n_divergent, diag$n_max_treedepth, diag$mean_accept_stat
  ))
  cat(sprintf(
    "  Rhat max:         %.4f\n  Min ESS ratio:    %.3f\n",
    max(diag$rhat, na.rm = TRUE),
    min(diag$n_eff_ratio, na.rm = TRUE)
  ))
  cat("  Per-parameter (key params):\n")
  sub    <- summary_df[pars, c("mean", "sd", "Rhat", "n_eff"), drop = FALSE]
  sub_df <- as.data.frame(round(sub, 3))
  sub_df$n_eff <- round(sub_df$n_eff)
  print(sub_df)
  cat("\n")
}

cat("\n=== Convergence Diagnostics ===\n\n")
print_convergence(fit_A, "Model A: lognormal-g (fit_full_model.stan)")
for (nm in names(fits_B)) {
  a <- as.numeric(sub("B_ig", "", nm))
  print_convergence(fits_B[[nm]], sprintf("Model B: IG-g (ig_a = %.1f, ig_b = %.1f)", a, ig_b))
}

# ── 7. Compute prior0 and post0 for theta ────────────────────────────────────

# Prior at theta = 0 -----------------------------------------------------------

prior0_A <- prior0_theta_analytical_lognormal(
  X        = X_sim_scaled,
  mu_log_g = log(N_sim),
  sd_log_g = 0.5
)

prior0_B <- lapply(IG_VALS, function(a) {
  prior0_theta_analytical(X = X_sim_scaled, ig_a = a, ig_b = ig_b)
})
names(prior0_B) <- paste0("B_ig", IG_VALS)

# Posterior at theta = 0 (logspline) -------------------------------------------

post0_for_fit <- function(fit) {
  draws <- rstan::extract(fit)$theta
  vapply(seq_len(P_sim), function(j) post0_theta_logspline(draws[, j]), numeric(1))
}

post0_A <- post0_for_fit(fit_A)
post0_B <- lapply(fits_B, post0_for_fit)

# ── 8. Assemble and print results ─────────────────────────────────────────────

pred_names <- colnames(X_sim)

# Build one wide data.frame: predictor | Model A cols | Model B cols (x3)
results <- data.frame(predictor = pred_names)

add_model_cols <- function(df, prior0, post0, suffix) {
  df[[paste0("prior0_",       suffix)]] <- round(prior0,                    6)
  df[[paste0("post0_",        suffix)]] <- round(post0,                     6)
  df[[paste0("log_SD_ratio_", suffix)]] <- round(log(prior0 / post0),       4)
  df
}

results <- add_model_cols(results, prior0_A, post0_A, "A")
for (nm in names(fits_B)) {
  suffix <- sub("B_ig", "B", nm)   # "B0.5", "B1.0", "B2.0"
  results <- add_model_cols(results, prior0_B[[nm]], post0_B[[nm]], suffix)
}

cat("\n")
cat("Savage-Dickey log BF10 for theta (linear effect)\n")
cat("  Model A:    lognormal-g  (mu=log(N), sd=0.5)\n")
for (a in IG_VALS)
  cat(sprintf("  Model B%.1f:  IG-g         (ig_a=%.1f, ig_b=%.1f)\n", a, a, ig_b))
cat("\n")
print(results, row.names = FALSE)

cat("\nPositive log_SD_ratio -> evidence for linear effect (H1: theta != 0)\n")
cat("Negative log_SD_ratio -> evidence for null          (H0: theta  = 0)\n")
