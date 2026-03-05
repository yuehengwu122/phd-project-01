#' Extract and Summarize Posterior
#'
#' Extracts posterior samples and computes summary statistics
#' (median, mean) for specified parameters.
#'
#' @param fit Stan sampling object (output from fit_model)
#' @param param_names Vector of parameter names to extract (default: c("delta", "lambda"))
#'
#' @return List containing median and mean for each parameter
#'
#' @examples
#' \dontrun{
#' post_summary <- extract_posterior_summary(fit,
#'                                           param_names = c("delta", "lambda"))
#' }
extract_posterior_summary <- function(fit, param_names = c("delta", "lambda")) {
  # Extract all posterior samples
  post <- rstan::extract(fit)

  # Compute summaries for each requested parameter
  result <- list()
  for (param in param_names) {
    if (param %in% names(post)) {
      result[[paste0(param, "_median")]] <- median(post[[param]])
      result[[paste0(param, "_mean")]] <- mean(post[[param]])
    }
  }

  return(result)
}

#' Create Single Result Row for Simulation Study
#'
#' Computes Bayes Factor and Savage-Dickey ratio for a single simulated dataset
#' and returns results as a data frame row.
#'
#' @param N Sample size
#' @param prior Prior version (1 or 2)
#' @param delta True delta parameter value
#' @param lambda True lambda parameter value
#' @param fit1 Fitted GP model (Stan sampling object)
#' @param fit0 Fitted null model (Stan sampling object)
#' @param prior_density Prior density evaluated at zero
#'
#' @return Data frame with one row containing: N, prior, delta, lambda,
#'         delta_median, lambda_median, bf10, sdr10
#'
#' @details
#' Extracts parameter estimates and computes both bridge-sampling Bayes Factor
#' and Savage-Dickey ratio for model comparison.
#'
#' @examples
#' \dontrun{
#' result_row <- create_result_row(
#'   N = 100, prior = 1, delta = 2.5, lambda = 1/3,
#'   fit1 = fit_gp, fit0 = fit_null,
#'   prior_density = 0.1487
#' )
#' }
create_result_row <- function(
  N,
  prior,
  delta,
  lambda,
  fit1,
  fit0,
  prior_density
) {
  # Extract posterior samples
  post <- rstan::extract(fit1)

  # Extract delta parameter (differs by prior version)
  if (prior == 1) {
    # Prior 1: delta is direct parameter
    param_delta <- post$delta
  } else {
    # Prior 2: delta2 is sampled, take square root
    param_delta <- sqrt(post$delta2)
  }

  # Compute Bayes Factor via bridge sampling
  bf10 <- compute_bayes_factor(fit1, fit0)

  # Compute Savage-Dickey Bayes Factor
  sdr10 <- compute_savage_dickey_ratio(param_delta, prior_density)

  # Return results as data frame row
  data.frame(
    N = N,
    prior = prior,
    delta = delta,
    lambda = lambda,
    delta_median = median(param_delta),
    lambda_median = median(post$lambda),
    bf10 = bf10,
    sdr10 = round(sdr10, 3)
  )
}
