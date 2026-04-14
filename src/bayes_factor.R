# ============================================================
# Bayes Factor Computation
#
# Two approaches:
# 1. Interval-null BF via Savage-Dickey density/interval ratios
# 2. Bridge sampling BF via marginal likelihood estimation
# ============================================================


# --- Posterior density/interval estimation (logspline) ---

#' Estimate posterior P(delta < eps | y) via logspline
#'
#' Uses tail-trimming for robustness when logspline fails on the
#' full distribution (heavy tails).
#'
#' @param x Posterior draws of delta (non-negative)
#' @param eps Epsilon threshold (>= 0)
#' @param max_attempts Number of tail-trimming levels to try
#' @return Scalar probability estimate, or NA on failure
post_prob_delta_logspline <- function(x, eps, max_attempts = 5) {
  x <- x[is.finite(x) & x >= 0]
  stopifnot(length(x) > 10, is.finite(eps), eps >= 0)

  quantiles <- c(1, 0.999, 0.99, 0.9, 0.8)

  for (q in quantiles[seq_len(max_attempts)]) {
    cutoff <- as.numeric(stats::quantile(x, q, names = FALSE, type = 5))
    keep <- x <= cutoff
    x_trim <- x[keep]

    fit <- try(
      logspline::logspline(x_trim, lbound = 0, ubound = max(x_trim)),
      silent = TRUE
    )
    if (!inherits(fit, "try-error")) {
      p_keep <- mean(keep)
      if (eps == 0) {
        return(p_keep * 2 * logspline::dlogspline(0, fit))
      }
      return(p_keep * logspline::plogspline(eps, fit))
    }
  }

  warning("logspline failed for all tail-cuts. Returning NA.")
  NA_real_
}


#' Estimate posterior P(delta < eps | y) via empirical CDF
#'
#' Fallback method when logspline is unreliable.
#' Returns both the point estimate and its Monte Carlo standard error.
#'
#' @param x Posterior draws of delta
#' @param eps Epsilon threshold
#' @return Named list with `post0` (estimate) and `mcse` (MC standard error)
post_prob_delta_ecdf <- function(x, eps) {
  x <- x[is.finite(x) & x >= 0]
  p_hat <- mean(x < eps)
  mcse <- sqrt(p_hat * (1 - p_hat) / length(x))
  list(post0 = p_hat, mcse = mcse)
}


#' Vectorized posterior probability computation
#'
#' @param x Posterior draws of delta
#' @param eps_vals Vector of epsilon thresholds
#' @param method "logspline" or "ecdf"
#' @return data.frame with columns `post0` and (for ecdf) `mcse`
post_prob_delta_vec <- function(
  x, eps_vals,
  method = c("logspline", "ecdf")
) {
  method <- match.arg(method)

  if (method == "logspline") {
    post0 <- vapply(eps_vals, function(e) post_prob_delta_logspline(x, e), numeric(1))
    return(data.frame(post0 = post0, mcse = NA_real_))
  }

  res <- lapply(eps_vals, function(e) post_prob_delta_ecdf(x, e))
  data.frame(
    post0 = vapply(res, `[[`, numeric(1), "post0"),
    mcse  = vapply(res, `[[`, numeric(1), "mcse")
  )
}


# --- Interval-null BF for delta (nonlinear effect) ---

#' Compute interval-null BF10 for delta
#'
#' BF10 = P(delta < eps) / P(delta < eps | y)
#'
#' where the prior is the moment prior of order d_order,
#' and the posterior is estimated via logspline.
#'
#' @param delta_draws Posterior draws of delta (vector, non-negative)
#' @param d_order Moment prior order (0, 1, 2, 3)
#' @param eps Vector of epsilon thresholds
#' @param method Posterior estimation method ("logspline" or "ecdf")
#' @param return_log If TRUE, return log(BF10) (default: TRUE)
#' @return data.frame with columns: epsilon, prior0, post0, log_bf10 (or bf10).
#'   When method = "ecdf", also includes post0_mcse.
bf10_interval_null_delta <- function(
  delta_draws,
  d_order,
  eps,
  method = c("logspline", "ecdf"),
  return_log = TRUE
) {
  method <- match.arg(method)

  delta_draws <- as.numeric(delta_draws)
  delta_draws <- delta_draws[is.finite(delta_draws) & delta_draws >= 0]

  prior0 <- prior_prob_delta(d_order, eps)
  post0_df <- post_prob_delta_vec(delta_draws, eps, method = method)

  bf10 <- prior0 / post0_df$post0
  if (return_log) bf10 <- log(bf10)

  out <- data.frame(epsilon = eps, prior0 = prior0, post0 = post0_df$post0)
  if (method == "ecdf") out$post0_mcse <- post0_df$mcse
  if (return_log) out$log_bf10 <- bf10 else out$bf10 <- bf10
  out
}


# --- Savage-Dickey ratio for theta (linear effect) ---

#' Compute analytical prior density at theta = 0
#'
#' Under the mixture-of-g prior:
#'   theta | g ~ N(0, g * (X'X)^{-1})
#'   log(g) ~ N(log(N), 0.5^2)
#'
#' @param X Scaled predictor matrix (as passed to Stan)
#' @param mu_log_g Prior mean of log(g) (default: log(N))
#' @param sd_log_g Prior SD of log(g) (default: 0.5)
#' @return Vector of prior densities at zero, one per predictor
prior0_theta_analytical <- function(
  X,
  mu_log_g = log(nrow(X)),
  sd_log_g = 0.5
) {
  X_c <- scale(X, center = TRUE, scale = FALSE)
  XtX_inv <- solve(crossprod(X_c))
  v <- diag(XtX_inv)
  Eg_inv_sqrt_g <- exp(-0.5 * mu_log_g + 0.125 * sd_log_g^2)
  (1 / sqrt(2 * pi * v)) * Eg_inv_sqrt_g
}


#' Estimate posterior density at theta = 0 via logspline
#'
#' @param theta_draws Posterior draws of theta (vector for one predictor)
#' @return Scalar density estimate at zero
post0_theta_logspline <- function(theta_draws) {
  fit <- logspline::logspline(theta_draws)
  logspline::dlogspline(0, fit)
}


# --- Main entry points for application analysis ---

#' Compute BF10 table for all predictors (linear + nonlinear)
#'
#' For each predictor j:
#'   - Linear BF10: Savage-Dickey ratio for theta_j = 0
#'   - Nonlinear BF10: interval-null BF for delta_j < eps
#'
#' @param result A gp_fit object from fit_model()
#' @param eps Epsilon threshold for nonlinear BF (scalar, default: 0 = point null)
#' @param p Optional vector of predictor indices to include
#' @return data.frame with columns: predictor, log_bf10_linear, log_bf10_nonlinear
compute_bf_table <- function(result, eps = 0, p = NULL) {
  stopifnot(inherits(result, "gp_fit"))
  stopifnot(result$model_type == "full")

  post <- rstan::extract(result$fit)
  X <- result$stan_data$X
  P <- result$stan_data$P
  d_order <- result$d_order
  pred_names <- result$predictor_names

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  # --- Linear BF10: Savage-Dickey at theta = 0 ---
  prior0_theta <- prior0_theta_analytical(X)
  theta_draws <- post$theta

  post0_theta <- vapply(seq_len(P), function(j) {
    post0_theta_logspline(theta_draws[, j])
  }, numeric(1))

  log_bf10_linear <- log(prior0_theta / post0_theta)

  # --- Nonlinear BF10: interval-null for delta ---
  delta_draws <- post$delta

  log_bf10_nonlinear <- vapply(seq_len(P), function(j) {
    draws_j <- if (is.matrix(delta_draws)) delta_draws[, j] else delta_draws
    bf <- bf10_interval_null_delta(draws_j, d_order, eps = eps, return_log = TRUE)
    bf$log_bf10[1]
  }, numeric(1))

  data.frame(
    predictor = pred_names[p_idx],
    log_bf10_linear = round(log_bf10_linear[p_idx], 4),
    log_bf10_nonlinear = round(log_bf10_nonlinear[p_idx], 4)
  )
}


#' Compute BF10 via bridge sampling (full model vs linear model)
#'
#' @param result_full gp_fit object from fit_model(..., model = "full")
#' @param result_linear gp_fit object from fit_model(..., model = "linear")
#' @param log If TRUE, return log(BF10)
#' @return Scalar BF10 (or log BF10)
compute_bf_bridge <- function(result_full, result_linear, log = FALSE) {
  bridge_full <- bridgesampling::bridge_sampler(result_full$fit)
  bridge_lin <- bridgesampling::bridge_sampler(result_linear$fit)
  bridgesampling::bf(bridge_full, bridge_lin, log = log)$bf
}
