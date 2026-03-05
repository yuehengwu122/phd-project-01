#' Compute BF10 by bridge sampling (inputs: two fit results; output: scalar BF10).
compute_bayes_factor_bridge <- function(result1, result0, log = FALSE) {
  fit1 <- result1$fit
  fit0 <- result0$fit
  # Compute marginal likelihoods via bridge sampling
  bridge1 <- bridgesampling::bridge_sampler(fit1)
  bridge0 <- bridgesampling::bridge_sampler(fit0)

  # Return Bayes Factor (ratio of marginal likelihoods)
  bridgesampling::bf(bridge1, bridge0, log = log)$bf
}


#' Compute prior density at zero for theta (input: X; output: vector of densities).
prior0_theta_analytical <- function(
  X,
  N = nrow(X),
  mu_log_g = log(N),
  sd_log_g = 0.5
) {
  # input rescaled X in [-0.5, 0.5]
  X_c <- scale(X, center = TRUE, scale = FALSE)

  XtX_inv <- solve(crossprod(X_c)) # (X'X)^{-1}
  v <- diag(XtX_inv)

  Eg_inv_sqrt_g <- exp(-0.5 * mu_log_g + 0.125 * sd_log_g^2) # E[g^{-1/2}]
  prior_d0 <- (1 / sqrt(2 * pi * v)) * Eg_inv_sqrt_g

  prior_d0
}


#' Logspline estimate posterior density at eps for delta (input: draws; output: scalar density).
post0_delta_logspline <- function(x, eps = 0, max_attempts = 5) {
  quantiles <- c(1, 0.999, 0.99, 0.9, 0.8)

  for (q in quantiles[1:max_attempts]) {
    cutoff <- as.numeric(stats::quantile(x, q, names = FALSE, type = 5))
    keep <- x <= cutoff
    x_trim <- x[keep]

    logspline_args <- list(x_trim, lbound = 0, ubound = max(x_trim))

    fit <- try(do.call(logspline::logspline, logspline_args), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      # density under conditional (trimmed) distribution
      # renormalize back to the original distribution
      p_keep <- mean(keep) # ~ q, but explicit

      if (eps == 0) {
        f0_trim <- 2 * logspline::dlogspline(0, fit)
        return(p_keep * f0_trim)
      } else {
        f0_trim <- logspline::plogspline(eps, fit)
        return(p_keep * f0_trim)
      }
    }
  }
  warning("logspline failed for all tail-cuts. Returning NA.")
  NA_real_
}


#' Compute prior mass at zero or interval (inputs: prior_id, eps/d; output: vector).
prior_prob_delta_lt_epsilon <- function(
  prior_id = c("p1", "p2", "p3", "p4"),
  eps = NULL,
  d = NULL,
  mad_x = NULL
) {
  prior_id <- match.arg(prior_id)

  spec <- switch(
    prior_id,
    p1 = list(power = 1, df = 4, scale = 2.7),
    p2 = list(power = 2, df = 1.5, scale = 6),
    p3 = list(power = 3, df = 1, scale = 12),
    p4 = list(power = 4, df = 0.7, scale = 25)
  )

  prior_from_eps <- function(eps_val) {
    if (!is.finite(eps_val) || eps_val < 0) {
      return(NA_real_)
    }
    if (eps_val == 0) {
      return((1 / spec$scale) * dt(0 / spec$scale, df = spec$df) * 2)
    }
    2 * pt(eps_val / spec$scale, df = spec$df) - 1
  }

  if (!is.null(eps)) {
    if (!is.numeric(eps)) {
      stop("`eps` must be numeric.")
    }
    if (length(eps) < 1L) {
      stop("`eps` must have length >= 1.")
    }
    if (anyNA(eps)) {
      stop("`eps` must not contain NA.")
    }
    if (any(!is.finite(eps))) {
      stop("`eps` must be finite.")
    }
    if (any(eps < 0)) stop("`eps` must be nonnegative.")
  }

  if (!is.null(d)) {
    if (!is.numeric(d)) {
      stop("`d` must be numeric when provided.")
    }
    if (length(d) < 1L) {
      stop("`d` must have length >= 1.")
    }
    if (anyNA(d)) {
      stop("`d` must not contain NA.")
    }
    if (any(!is.finite(d))) {
      stop("`d` must be finite.")
    }
    if (any(d <= 0)) {
      stop("`d` must be positive.")
    }
    if (is.null(mad_x) || !is.finite(mad_x) || mad_x <= 0) {
      stop("`mad_x` must be provided and > 0 when `d` is used.")
    }
  }

  if (is.null(eps) && (is.null(d) || all(d == 0))) {
    eps <- 0
  }
  eps_vals <- c(
    if (is.null(eps)) numeric(0) else eps,
    if (is.null(d) || all(d == 0)) numeric(0) else mad_x / d
  )
  vapply(eps_vals, prior_from_eps, numeric(1))
}


#' Compute posterior mass at zero or interval (inputs: draws, eps/d; output: vector).
post_prob_delta_lt_epsilon <- function(x, eps = NULL, d = NULL) {
  d0_from_eps <- function(eps_val) {
    if (!is.finite(eps_val) || eps_val < 0) {
      return(NA_real_)
    }
    post0_delta_logspline(x, eps = eps_val)
  }

  if (!is.null(eps)) {
    if (!is.numeric(eps)) {
      stop("`eps` must be numeric.")
    }
    if (length(eps) < 1L) {
      stop("`eps` must have length >= 1.")
    }
    if (anyNA(eps)) {
      stop("`eps` must not contain NA.")
    }
    if (any(!is.finite(eps))) {
      stop("`eps` must be finite.")
    }
    if (any(eps < 0)) stop("`eps` must be nonnegative.")
  }

  if (!is.null(d)) {
    if (!is.numeric(d)) {
      stop("`d` must be numeric when provided.")
    }
    if (length(d) < 1L) {
      stop("`d` must have length >= 1.")
    }
    if (anyNA(d)) {
      stop("`d` must not contain NA.")
    }
    if (any(!is.finite(d))) {
      stop("`d` must be finite.")
    }
    if (any(d <= 0)) stop("`d` must be positive.")
  }

  if (is.null(eps) && (is.null(d) || all(d == 0))) {
    eps <- 0
  }
  eps_vals <- c(
    if (is.null(eps)) numeric(0) else eps,
    if (is.null(d) || all(d == 0)) numeric(0) else mad(x) / d
  )
  vapply(eps_vals, d0_from_eps, numeric(1))
}


#' Compute BF10 for delta with point/interval null (inputs: draws, eps/d; output: data.frame).
bf10_interval_null_delta_1d <- function(
  delta_draws,
  eps = NULL,
  d = NULL,
  prior_id = c("p1", "p2", "p3", "p4"),
  return_log = TRUE
) {
  prior_id <- match.arg(prior_id)
  if (!prior_id %in% c("p1", "p2", "p3", "p4")) {
    stop("`prior_id` must be one of 'p1', 'p2', 'p3', 'p4'.")
  }

  post_prob <- post_prob_delta_lt_epsilon(delta_draws, eps = eps, d = d)
  prior_prob <- prior_prob_delta_lt_epsilon(
    prior_id,
    eps = eps,
    d = d,
    mad_x = mad(delta_draws)
  )

  if (is.null(eps) && (is.null(d) || all(d == 0))) {
    eps <- 0
  }
  eps_vals <- c(
    if (is.null(eps)) numeric(0) else eps,
    if (is.null(d) || all(d == 0)) numeric(0) else mad(delta_draws) / d
  )
  tag <- c(
    if (is.null(eps)) character(0) else paste0("eps=", eps),
    if (is.null(d) || all(d == 0)) character(0) else paste0("mad/", d)
  )

  bf10 <- prior_prob / post_prob
  if (return_log) {
    bf10 <- log(bf10)
  }

  out <- data.frame(
    epsilon = eps_vals,
    prior0 = prior_prob,
    post0 = post_prob,
    bf10 = bf10
  ) |>
    round(4)
  out$tag <- tag
  if (return_log) {
    out$log_bf10 <- out$bf10
    out$bf10 <- NULL
  }
  out
}


#' Compute linear + nonlinear BF10s from a fit result (inputs: result, eps/d; output: data.frame).
compute_sdr_from_fit <- function(
  result,
  eps = NULL,
  d = NULL,
  prior_id = c("p1", "p2", "p3", "p4"),
  p = NULL,
  predictor_names = NULL
) {
  if (!is.list(result) || is.null(result$fit) || is.null(result$data$X)) {
    stop(
      "`result` must be a fit_model() result containing `$fit` and `$data$X`."
    )
  }
  if (!is.null(eps)) {
    eps <- as.numeric(eps)
    if (!is.numeric(eps)) {
      stop("`eps` must be numeric.")
    }
    if (length(eps) < 1L) {
      stop("`eps` must have length >= 1.")
    }
    if (anyNA(eps)) {
      stop("`eps` must not contain NA.")
    }
    if (any(!is.finite(eps))) {
      stop("`eps` must be finite.")
    }
    if (any(eps < 0)) stop("`eps` must be nonnegative.")
  }

  if (!is.null(d)) {
    if (!is.numeric(d)) {
      stop("`d` must be numeric when provided.")
    }
    if (length(d) < 1L) {
      stop("`d` must have length >= 1.")
    }
    if (anyNA(d)) {
      stop("`d` must not contain NA.")
    }
    if (any(!is.finite(d))) {
      stop("`d` must be finite.")
    }
    if (any(d <= 0)) stop("`d` must be positive.")
  }
  prior_id <- match.arg(prior_id)
  if (!prior_id %in% c("p1", "p2", "p3", "p4")) {
    stop("`prior_id` must be one of 'p1', 'p2', 'p3', 'p4'.")
  }

  X <- result$data$X
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  P <- ncol(X)

  if (is.null(p)) {
    p_idx <- seq_len(P)
  } else {
    if (!is.numeric(p) || anyNA(p)) {
      stop("`p` must be an integer-like vector with no NA.")
    }
    p_idx <- as.integer(p)
    if (any(p_idx < 1L | p_idx > P)) stop("`p` indices out of range.")
  }

  if (is.null(predictor_names)) {
    predictor_names <- colnames(X)
    if (is.null(predictor_names)) predictor_names <- paste0("X", seq_len(P))
  }
  if (!is.character(predictor_names) || length(predictor_names) != P) {
    stop("`predictor_names` must be a character vector of length P.")
  }

  # extract posterior draws matrix
  post <- rstan::extract(result$fit)

  delta_draws <- switch(
    prior_id,
    p1 = post$delta,
    p2 = post$delta2,
    p3 = post$delta3,
    p4 = post$delta4
  )

  theta_draws <- post$theta

  # compute prior0, post0, and sdr for theta
  prior0_theta <- prior0_theta_analytical(X)

  post0_theta <- sapply(1:P, function(j) {
    fit <- logspline::logspline(
      theta_draws[, j],
      lbound = min(min(theta_draws[, j]), 0),
      ubound = max(max(theta_draws[, j]), 0)
    )
    logspline::dlogspline(0, fit)
  })

  log_sdr10_linear <- log(prior0_theta / post0_theta)

  if (is.null(eps) && (is.null(d) || all(d == 0))) {
    eps <- 0
  }
  log_sdr10_nonlinear <- sapply(1:P, function(j) {
    out <- bf10_interval_null_delta_1d(
      delta_draws[, j],
      eps = eps,
      d = d,
      prior_id = prior_id,
      return_log = TRUE
    )
    out$log_bf10[1]
  })

  data.frame(
    Predictor = predictor_names,
    SDR_linear = format(log_sdr10_linear, digits = 3, scientific = TRUE),
    SDR_nonlinear = format(log_sdr10_nonlinear, digits = 3, scientific = TRUE)
  )
}
