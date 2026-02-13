#' Compute Bayes Factor via Bridge Sampling
#'
#' Computes the Bayes Factor (BF10) comparing a focal model to a reference 
#' (null) model using bridge sampling approximation of marginal likelihoods.
#'
#' @param fit1 Fitted Stan model object (focal/alternative model)
#' @param fit0 Fitted Stan model object (reference/null model)
#'
#' @return Bayes Factor BF10 = P(data | M1) / P(data | M0).
#' Values > 1 favor the focal model; values < 1 favor the null model.
#'
#' @details
#' Bridge sampling estimates the marginal likelihood for each model by
#' using samples from the posterior combined with importance weights.
#' This method is more numerically stable than thermodynamic integration
#' and works well in moderate to high dimensions.
#'
#' @references
#' Gronau, Q. F., Singmann, H., & Wagenmakers, E. J. (2020).
#' "bridgesampling: An R package for estimating normalizing constants."
#' Journal of Statistical Software, 92(10), 1-29.
#'
#' @examples
#' \dontrun{
#' # Compare GP model to linear null model
#' bf_gp_vs_null <- compute_bayes_factor(fit_gp, fit_null)
#' log_bf <- log(bf_gp_vs_null)
#' }
compute_bayes_factor <- function(result1, result0) {
  fit1 <- result1$fit
  fit0 <- result0$fit
  # Compute marginal likelihoods via bridge sampling
  bridge1 <- bridgesampling::bridge_sampler(fit1)
  bridge0 <- bridgesampling::bridge_sampler(fit0)
  
  # Return Bayes Factor (ratio of marginal likelihoods)
  bridgesampling::bf(bridge1, bridge0)$bf
}


# ============================================================
# Savage-Dickey Ratio Functions
# ============================================================

#' Estimate Posterior Density at Zero using Logspline
#' 
#' Fits a logspline density estimator to posterior samples and estimates 
#' the density at zero. Attempts multiple tail cutoffs if the initial fit fails.
#' Used for computing the Savage-Dickey ratio Bayes factor.
#'
#' @param x Posterior samples (numeric vector)
#' @param lbound Lower bound for logspline support (default: 0)
#' @param ubound Upper bound for logspline support (numeric, optional)
#' @param max_attempts Maximum number of tail cutoff attempts (default: 7)
#'
#' @return Estimated posterior density at zero. Returns NA if all attempts fail.
#'
#' @details
#' The logspline density estimator can support bounded intervals via lbound and ubound.
#' When fitting fails, the function progressively trims the tail of the sample
#' by using increasingly stringent quantile cutoffs (1.0, 0.999, 0.99, ..., 0.6).
#'
#' @references
#' Stone, C. J., Hansen, M. H., Kooperberg, C., & Truong, Y. K. (1997).
#' "Polynomial splines and their tensor products in extended linear modeling."
#' Annals of Statistics, 25(4), 1371-1470.
#'
#' @examples
#' \dontrun{
#' # Simulate posterior samples
#' post_samples <- rnorm(1000, mean = 2, sd = 0.5)
#' 
#' # Estimate density at zero
#' density_at_zero <- estimate_post_density_at_zero(
#'   post_samples, 
#'   lbound = 0
#' )
#' }
estimate_post_density_at_zero <- function(x, lbound = 0, ubound = NULL, 
                                          max_attempts = 7) {
  # Candidate quantile cutoffs for tail trimming
  quantiles <- c(1, 0.999, 0.99, 0.9, 0.8, 0.7, 0.6)
  
  for (q in quantiles[1:max_attempts]) {
    # Trim right tail at specified quantile
    x_trimmed <- x[x <= quantile(x, q)]
    
    # Build logspline call with optional upper bound
    logspline_args <- list(x_trimmed, lbound = lbound)
    if (!is.null(ubound)) {
      logspline_args$ubound <- ubound
    }
    
    # Attempt logspline fit
    try_fit <- try(do.call(logspline::logspline, logspline_args), 
                   silent = TRUE)
    
    # If successful, return density estimate at zero
    if (!inherits(try_fit, "try-error")) {
      return(logspline::dlogspline(0, try_fit))
    }
  }
  
  # If all attempts fail, issue warning and return NA
  warning("logspline failed for all tail-cuts. Returning NA.")
  return(NA)
}

#' Compute Prior Density at Zero for Theta (Linear Coefficients)
#'
#' Computes the prior density at zero for Cauchy prior on theta.
#'
#' @param X Predictor matrix (N x P)
#' @return Vector of prior densities at 0 for each predictor
compute_prior_density_theta <- function(X) {
  N <- nrow(X)
  
  # Center predictors (matching Stan's transformed data)
  X_c <- scale(X, center = TRUE, scale = FALSE)
  
  # Compute JZS scale factor for each predictor
  scales_theta <- apply(X_c, 2, function(x_col) {
    sqrt(N) / sqrt(sum(x_col^2))
  })
  
  # Cauchy density at 0
  dcauchy(0, location = 0, scale = scales_theta)
}



#' Compute Savage-Dickey Ratios for All Parameters
#'
#' Computes the Savage-Dickey Bayes Factors for linear (theta) and 
#' nonlinear (delta) effects. Prior density for delta depends on the specified prior.
#'
#' @param result Result object from run_analysis()
#' @param delta_prior Prior specification for delta ("p1" or "p2"). Default is "p1".
#'   - "p1": half-t(df=4, s=2.7) prior, uses "delta" parameter
#'   - "p2": half-t(df=1.5, s=6) prior, uses "delta2" parameter
#' @return Data frame with SDR for linear (theta) and nonlinear (delta) effects
#'
#' @details
#' The function extracts X from result$data and computes:
#' - SDR_linear: Savage-Dickey ratio for theta (linear coefficients)
#' - SDR_nonlinear: Savage-Dickey ratio for delta (nonlinear effects)
#' 
#' For p1: Prior density = (1/2.7) * dt(0/2.7, df = 4) * 2
#' For p2: Prior density = (1/6) * dt(0/6, df = 1.5) * 2
compute_savage_dickey_ratios <- function(result, delta_prior = "p1") {
  
  fit <- result$fit
  post <- rstan::extract(fit)
  X <- result$data$X
  P <- ncol(X)
  
  # --- Linear effect (theta) ---
  prior_d0_theta <- compute_prior_density_theta(X)
  
  theta_draws <- post$theta
  if (is.matrix(theta_draws)) {
    post_d0_theta <- sapply(1:P, function(j) {
      estimate_post_density_at_zero(
        theta_draws[, j],
        lbound = min(0, min(theta_draws[, j])),
        ubound = max(0, max(theta_draws[, j]))
      )
    })
  } else {
    # Single predictor case
    post_d0_theta <- estimate_post_density_at_zero(
      theta_draws,
      lbound = min(0, min(theta_draws)),
      ubound = max(0, max(theta_draws))
    )
  }
  
  sdr_linear <- prior_d0_theta / post_d0_theta
  
  # --- Nonlinear effect (delta) ---
  # Set prior and parameter based on delta_prior argument
  if (delta_prior == "p1") {
    # Prior: half-t(df=4, s=2.7)
    prior_d0_delta <- (1/2.7) * dt(0/2.7, df = 4) * 2
    delta_draws <- post$delta
  } else if (delta_prior == "p2") {
    # Prior: half-t(df=1.5, s=6)
    prior_d0_delta <- (1/6) * dt(0/6, df = 1.5) * 2
    delta_draws <- post$delta2
  } else {
    stop("delta_prior must be either 'p1' or 'p2'")
  }
  if (is.matrix(delta_draws)) {
    post_d0_delta <- sapply(1:P, function(j) {
      estimate_post_density_at_zero(
        delta_draws[, j],
        lbound = 0,
        ubound = max(delta_draws[, j])
      )
    })
  } else {
    # Single predictor case
    post_d0_delta <- estimate_post_density_at_zero(
      delta_draws,
      lbound = 0,
      ubound = max(delta_draws)
    )
  }
  
  sdr_nonlinear <- prior_d0_delta / post_d0_delta
  
  
  data.frame(
    # predictor = predictor_names,
    SDR_linear = format(sdr_linear, digits = 3, scientific = TRUE),
    SDR_nonlinear = format(sdr_nonlinear, digits = 3, scientific = TRUE)
  )
}


#' # Compute Savage-Dickey ratio
#' sdr <- compute_savage_dickey_ratio(posterior_delta, prior_at_zero)
#' log_sdr <- log(sdr)
#' }
compute_savage_dickey_ratio <- function(posterior_samples, 
                                        prior_density_at_zero) {
  # Estimate posterior density at zero using logspline
  post_density_at_zero <- estimate_post_density_at_zero(posterior_samples)
  
  # Compute and return Savage-Dickey ratio
  prior_density_at_zero / post_density_at_zero
}