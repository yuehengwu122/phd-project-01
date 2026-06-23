#' Simulation Experiment 1: Bayes Factor Behavior in Single-Predictor GP Model
#'
#' Assessment of Bayes factor behavior for detecting nonlinear (GP) effects
#' via a single continuous predictor. Evaluates two Bayes factor approaches:
#' bridge sampling (BF_BS) and Savage-Dickey ratio (BF_SDR).
#'
#' @details
#' Simulation setup:
#' - Data domain: range(X) = 1
#' - Noise standard deviation: sigma = 0.1
#' - Sample sizes: N in {30, 100, 200}
#' - Effect sizes: delta in {0, 0.5, 1, 1.5, 2, 2.5, 3}
#' - GP lengthscale: lambda = 1/3 (medium scale for range = 1)
#' - Prior on delta: 
#'   * Prior 1: Horseshoe-type prior (continuous support)
#'   * Prior 2: Nonlocal prior with zero density at null (delta = 0)
#' - Replicates per configuration: 50 independent datasets
#' - Output: BF_10^BS (bridge sampling) and BF_10^SDR (Savage-Dickey)
#'
#' Expected results (from paper):
#' - BF_BS increases monotonically with delta, steeper growth at larger N
#' - Prior 2 shows mild effects: smaller for small signals (delta <= 1),
#'   slightly larger for strong signals (delta > 1)
#' - BF_SDR tracks BF_BS for delta <= 1.5, then levels off
#' - Leveling reflects difficulty estimating pointwise posterior density
#'   when posterior mass is far from zero

source("_common.R")

# ============================================
# EXPERIMENT 1: PARAMETER SPECIFICATIONS
# ============================================

# Sample sizes to test
sample_sizes <- c(30, 100, 200)

# True effect sizes (delta) to vary
delta_values <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)

# GP lengthscale (medium scale for range(X) = 1)
lambda_fixed <- 1/3

# Noise standard deviation (fixed)
sigma_fixed <- 0.1

# Domain range
domain_range <- 1

# Number of independent replicates per configuration
n_replicates <- 50

# ============================================
# PRIOR SPECIFICATIONS
# ============================================

# Prior 1: Horseshoe-type prior (continuous support at zero)
# Density at zero: p(delta=0|Prior1)
prior_1_density_at_zero <- (1 / 2.7) * dt(0 / 2.7, df = 4) * 2

# Prior 2: Nonlocal prior (zero density at null)
# Parametrized on delta^2 with zero probability at delta^2 = 0
# Density at zero: p(delta^2=0|Prior2)
prior_2_density_at_zero <- (1 / 7) * dt(0 / 7, df = 1.15) * 2

# ============================================
# INITIALIZE RESULTS STORAGE
# ============================================

results <- data.frame(
  N = numeric(),                 # Sample size
  prior = numeric(),             # Prior version (1 or 2)
  delta = numeric(),             # True effect size
  lambda = numeric(),            # GP lengthscale
  bf_bs = numeric(),             # Bayes factor (bridge sampling)
  bf_sdr = numeric(),            # Bayes factor (Savage-Dickey ratio)
  delta_posterior_median = numeric()  # Posterior median of delta estimate
)

# ============================================
# COMPILE STAN MODELS
# ============================================



# ============================================
# MAIN SIMULATION LOOP
# ============================================

total_configs <- length(sample_sizes) * length(delta_values)
config_counter <- 0

for (N in sample_sizes) {
  # Create covariate grid for this sample size
  x_grid <- seq(-domain_range/2, domain_range/2, length = N)
  
  for (delta in delta_values) {
    config_counter <- config_counter + 1
    
    # ========== GENERATE DATA REPLICATES ==========
    
    y_matrix <- NULL
    if (delta == 0) {
      # Null case: pure noise
      message(sprintf(
        "[Config %d/%d] N=%d, delta=%.2f (null case), lambda=%.3f",
        config_counter, total_configs, N, delta, lambda_fixed
      ))
    } else {
      # Alternative case: use simulation model to generate data
      message(sprintf(
        "[Config %d/%d] N=%d, delta=%.2f, lambda=%.3f - Generating %d datasets",
        config_counter, total_configs, N, delta, lambda_fixed, n_replicates
      ))
      
      y_matrix <- generate_gp_data(
        delta = delta,
        lambda = lambda_fixed,
        sigma = sigma_fixed,
        N = N,
        x = x_grid,
        n_samples = n_replicates
      )
    }
    
    # ========== FIT MODELS TO EACH REPLICATE ==========
    
    for (rep in 1:n_replicates) {
      # Generate response variable
      if (delta == 0) {
        # Under null: pure noise realization
        set.seed(rep)
        y <- rnorm(N, mean = 0, sd = sigma_fixed)
      } else {
        # Under alternative: simulated GP realization
        y <- as.vector(y_matrix[rep, ])
      }
      
      # Prepare data list for Stan
      stan_data <- list(
        N = N,
        y = y,
        X = x_grid  |> matrix(ncol = 1),
        P = 1
      )
      
      # Fit null model (linear only) - compute once per replicate
      fit_null <- fit_model(gp_null, stan_data, 
        chains = 2,
        iter = 5500,
        warmup = 500)
      
      # ========== PRIOR 1: HORSESHOE ==========
      fit_gp_p1 <- fit_model(
        gp_add_p1_11,
        stan_data,
        chains = 2,
        iter = 5500,
        warmup = 500
      )
      
      # Extract posterior samples for Prior 1
      post_p1 <- rstan::extract(fit_gp_p1$fit)
      
      # Compute BF_BS (bridge sampling)
      bf_bs_p1 <- compute_bayes_factor_bridge(fit_gp_p1$fit, fit_null$fit)
      
      # Compute BF_SDR (Savage-Dickey ratio)
      bf_sdr_p1 <- compute_savage_dickey_ratio(
        post_p1$delta,
        prior_1_density_at_zero
      )
      
      # Store Prior 1 results
      results <- rbind(results, data.frame(
        N = N,
        prior = 1,
        delta = delta,
        lambda = lambda_fixed,
        bf_bs = bf_bs_p1,
        bf_sdr = bf_sdr_p1,
        delta_posterior_median = median(post_p1$delta)
      ))
      
      # ========== PRIOR 2: NONLOCAL ==========
      fit_gp_p2 <- fit_model(
        gp_single_p2,
        stan_data,
        chains = 2,
        iter = 5500,
        warmup = 500
      )
      
      # Extract posterior samples for Prior 2
      post_p2 <- rstan::extract(fit_gp_p2$fit)
      
      # Compute BF_BS (bridge sampling)
      bf_bs_p2 <- compute_bayes_factor_bridge(fit_gp_p2$fit, fit_null$fit)
      
      # Compute BF_SDR (Savage-Dickey ratio)
      # For Prior 2: delta^2 is the sampled parameter
      # Compute density at zero for delta^2
      bf_sdr_p2 <- compute_savage_dickey_ratio(
        post_p2$delta2,
        prior_2_density_at_zero
      )
      
      # Store Prior 2 results
      results <- rbind(results, data.frame(
        N = N,
        prior = 2,
        delta = delta,
        lambda = lambda_fixed,
        bf_bs = bf_bs_p2,
        bf_sdr = bf_sdr_p2,
        delta_posterior_median = median(sqrt(post_p2$delta2))
      ))
      
      # Save results periodically to disk
      write.csv(
        results,
        "results/simulation/01_single_gp_comparison.csv",
        row.names = FALSE
      )
      
      if (rep %% 10 == 0) {
        message(sprintf("  Completed %d/%d replicates", rep, n_replicates))
      }
    }
  }
}

message("\n=== Experiment 1 COMPLETE ===")
message(sprintf("Total configurations: %d", total_configs))
message(sprintf("Total replicates: %d", total_configs * n_replicates))
message(sprintf("Total model fits: %d", nrow(results)))
message(sprintf("Results saved to: results/simulation/01_single_gp_comparison.csv"))