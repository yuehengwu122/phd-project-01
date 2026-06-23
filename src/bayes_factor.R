# ============================================================
# Bayes Factor Computation
#
# Working decisions reflected here:
# 1. Linear effects use point-null Savage-Dickey ratios at theta = 0.
#    Posterior density at zero is estimated with dlogspline().
# 2. Nonlinear effects use interval-null Bayes factors for delta < eps,
#    with eps fixed by default at 0.2.
#    Posterior interval mass is estimated with plogspline().
# 3. Direct fitting under each non-local prior is retained, but we also
#    support reweighting a local/base-prior BF to a target non-local prior.
# ============================================================

BF_EPSILON_DEFAULT <- 0.2


# --- Nonlinear posterior estimation via logspline ---

#' Estimate posterior P(delta < eps | y) via logspline
#'
#' Uses tail-trimming for robustness when logspline fails on the
#' full distribution (heavy tails).
#'
#' @param x Posterior draws of delta (non-negative)
#' @param eps Interval-null threshold (>= 0)
#' @param max_attempts Number of tail-trimming levels to try
#' @return Scalar probability estimate, or NA on failure
post_prob_delta_logspline <- function(x, eps = BF_EPSILON_DEFAULT, max_attempts = 5) {
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


# --- Nonlinear interval-null BF ---

#' Compute interval-null BF10 for delta
#'
#' BF10 = P(delta < eps) / P(delta < eps | y)
#'
#' The prior is the moment prior of order d_order. The posterior interval
#' probability is estimated by logspline.
#'
#' @param delta_draws Posterior draws of delta (vector, non-negative)
#' @param d_order Moment prior order (0, 1, 2, 3)
#' @param eps Interval-null threshold. Defaults to 0.2.
#' @param return_log If TRUE, return log(BF10)
#' @return data.frame with columns epsilon, prior0, post0, and log_bf10 or bf10
bf10_interval_null_delta <- function(
  delta_draws,
  d_order,
  eps = BF_EPSILON_DEFAULT,
  return_log = TRUE
) {
  delta_draws <- as.numeric(delta_draws)
  delta_draws <- delta_draws[is.finite(delta_draws) & delta_draws >= 0]

  prior0 <- prior_prob_delta(d_order, eps)
  post0 <- post_prob_delta_logspline(delta_draws, eps = eps)

  bf10 <- prior0 / post0
  if (return_log) bf10 <- log(bf10)

  out <- data.frame(epsilon = eps, prior0 = prior0, post0 = post0)
  if (return_log) out$log_bf10 <- bf10 else out$bf10 <- bf10
  out
}


# --- Reweighting from a base/local prior ---

#' Compute the log reweighting adjustment from a base prior to a target prior
#'
#' For posterior draws sampled under the base prior,
#'
#'   log adjustment = log E_base[ pi_target(delta) / pi_base(delta) | y ]
#'
#' @param delta_draws Posterior draws of delta under the base fit
#' @param target_order Target moment-prior order
#' @param base_order Base moment-prior order
#' @return Scalar log adjustment
compute_log_reweight_adjustment <- function(
  delta_draws,
  target_order,
  base_order = 0L
) {
  stopifnot(target_order %in% 0:3, base_order %in% 0:3)

  draws <- as.numeric(delta_draws)
  draws <- draws[is.finite(draws) & draws >= 0]
  if (!length(draws)) {
    return(NA_real_)
  }

  log_target <- log(prior_moment_density(draws, target_order, get_kappa(target_order)))
  log_base <- log(prior_moment_density(draws, base_order, get_kappa(base_order)))
  log_weights <- log_target - log_base

  max_log_weight <- max(log_weights)
  max_log_weight + log(mean(exp(log_weights - max_log_weight)))
}


#' Reweight a base log-BF to a target moment-prior order
#'
#' @param base_log_bf Log BF computed under the base prior
#' @param delta_draws Posterior draws of delta from the base fit
#' @param target_order Target moment-prior order
#' @param base_order Base moment-prior order
#' @return data.frame with the base log-BF, log adjustment, and reweighted log-BF
reweight_log_bf <- function(
  base_log_bf,
  delta_draws,
  target_order,
  base_order = 0L
) {
  log_adjustment <- compute_log_reweight_adjustment(
    delta_draws = delta_draws,
    target_order = target_order,
    base_order = base_order
  )

  data.frame(
    base_order = base_order,
    target_order = target_order,
    log_bf10_base = base_log_bf,
    log_reweight_adjustment = log_adjustment,
    log_bf10_reweighted = base_log_bf + log_adjustment
  )
}


#' Compute a reweighted interval-null BF for one predictor from a base fit
#'
#' @param result_base gp_fit object, typically fitted with d_order = 0
#' @param target_order Target moment-prior order
#' @param predictor Predictor index
#' @param eps Interval-null threshold. Defaults to 0.2.
#' @return data.frame with base and reweighted log-BFs
compute_reweighted_nonlinear_bf <- function(
  result_base,
  target_order,
  predictor = 1L,
  eps = BF_EPSILON_DEFAULT
) {
  stopifnot(inherits(result_base, "gp_fit"))

  post <- rstan::extract(result_base$fit)
  delta_draws_all <- post$delta
  predictor <- as.integer(predictor)

  delta_draws <- if (is.matrix(delta_draws_all)) delta_draws_all[, predictor] else delta_draws_all
  base_order <- result_base$d_order
  base_log_bf <- bf10_interval_null_delta(
    delta_draws = delta_draws,
    d_order = base_order,
    eps = eps,
    return_log = TRUE
  )$log_bf10[1]

  out <- reweight_log_bf(
    base_log_bf = base_log_bf,
    delta_draws = delta_draws,
    target_order = target_order,
    base_order = base_order
  )
  out$predictor <- predictor
  out$epsilon <- eps
  out
}


# --- Linear point-null BF ---

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
#' Linear effects use point-null Savage-Dickey ratios at theta = 0.
#' Nonlinear effects use interval-null BFs for delta < eps with logspline.
#'
#' @param result A gp_fit object from fit_model()
#' @param eps Interval-null threshold for nonlinear BF. Defaults to 0.2.
#' @param p Optional vector of predictor indices to include
#' @return data.frame with columns predictor, log_bf10_linear, log_bf10_nonlinear
compute_bf_table <- function(result, eps = BF_EPSILON_DEFAULT, p = NULL) {
  stopifnot(inherits(result, "gp_fit"))
  stopifnot(result$model_type == "full")

  post <- rstan::extract(result$fit)
  X <- result$stan_data$X
  P <- result$stan_data$P
  d_order <- result$d_order
  pred_names <- result$predictor_names

  if (is.null(p)) p_idx <- seq_len(P) else p_idx <- as.integer(p)

  prior0_theta <- prior0_theta_analytical(X)
  theta_draws <- post$theta

  post0_theta <- vapply(seq_len(P), function(j) {
    post0_theta_logspline(theta_draws[, j])
  }, numeric(1))

  log_bf10_linear <- log(prior0_theta / post0_theta)

  delta_draws <- post$delta
  log_bf10_nonlinear <- vapply(seq_len(P), function(j) {
    draws_j <- if (is.matrix(delta_draws)) delta_draws[, j] else delta_draws
    bf <- bf10_interval_null_delta(
      delta_draws = draws_j,
      d_order = d_order,
      eps = eps,
      return_log = TRUE
    )
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
#' This is the direct-fit marginal-likelihood BF for the currently fitted prior.
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
