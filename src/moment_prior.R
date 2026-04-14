# ============================================================
# Moment-Prior Family for delta (Nonlinear Magnitude)
#
# Implements the non-local moment prior:
#   pi_d(delta) = delta^d / C_d * g_kappa(delta)
#
# where g_kappa is a scaled half-t base density and kappa is
# chosen so that the median equals PRIOR_MEDIAN_TARGET.
#
# Reference: Johnson & Rossell (2010)
# ============================================================

# --- Base prior family parameters ---
PRIOR_DF <- 6L
PRIOR_SCALE <- 2.79
PRIOR_MEDIAN_TARGET <- 2


#' Half-t base density (unnormalized on positive reals)
#' @param x Non-negative values
#' @param df Degrees of freedom
#' @param s Scale parameter
#' @return Density values
g_base <- function(x, df = PRIOR_DF, s = PRIOR_SCALE) {
  2 * dt(x / s, df = df) / s
}


#' Normalizing constant C_d for the moment prior
#'
#' C_d = E[delta^d] under the base half-t(df, s) distribution.
#' Requires d < df for finiteness.
#'
#' @param d Moment order (0, 1, 2, or 3)
#' @param df Degrees of freedom
#' @param s Scale parameter
#' @return Scalar normalizing constant
C_d_half_t <- function(d, df = PRIOR_DF, s = PRIOR_SCALE) {
  if (d == 0) return(1)
  if (d >= df) stop("Need d < df for finite moment")
  s^d * (df^(d / 2)) * gamma((d + 1) / 2) * gamma((df - d) / 2) /
    (sqrt(pi) * gamma(df / 2))
}


#' Scaled base density g_kappa(x) = (1/kappa) * g(x/kappa)
g_scaled <- function(x, kappa, df = PRIOR_DF, s = PRIOR_SCALE) {
  (1 / kappa) * g_base(x / kappa, df = df, s = s)
}


#' Moment prior density pi_d(x)
#'
#' @param x Non-negative values
#' @param d Moment order (0, 1, 2, 3)
#' @param kappa Median-matching scale
#' @param df Degrees of freedom
#' @param s Scale parameter
#' @return Density values
prior_moment_density <- function(x, d, kappa, df = PRIOR_DF, s = PRIOR_SCALE) {
  stopifnot(all(x >= 0))
  if (d == 0) return(g_scaled(x, kappa = kappa, df = df, s = s))
  Ck <- (kappa^d) * C_d_half_t(d, df = df, s = s)
  x^d * g_scaled(x, kappa = kappa, df = df, s = s) / Ck
}


#' Compute quantile of a density by numerical inversion
#' @param pdf_fun Density function (single argument)
#' @param p Probability level
#' @param upper_start Initial upper bound for search
#' @return Quantile value
prior_quantile <- function(pdf_fun, p, upper_start = NULL, rel.tol = 1e-8) {
  cdf <- function(q) {
    if (q <= 0) return(0)
    integrate(pdf_fun, lower = 0, upper = q, rel.tol = rel.tol)$value
  }
  upper <- if (is.null(upper_start)) 1 else upper_start
  while (cdf(upper) < p) upper <- upper * 2
  uniroot(function(q) cdf(q) - p, lower = 0, upper = upper, tol = 1e-8)$root
}


#' Compute median-matching kappa for moment order d
#'
#' Finds kappa such that the median of pi_d equals PRIOR_MEDIAN_TARGET.
#'
#' @param d Moment order
#' @return Scalar kappa value
compute_kappa <- function(d) {
  pdf_k1 <- function(x) prior_moment_density(x, d = d, kappa = 1)
  med1 <- prior_quantile(pdf_k1, p = 0.5, upper_start = PRIOR_SCALE)
  PRIOR_MEDIAN_TARGET / med1
}


# --- Precompute kappa for d = 0, 1, 2, 3 ---
.kappa_cache <- setNames(
  vapply(0:3, compute_kappa, numeric(1)),
  paste0("d", 0:3)
)


#' Get precomputed kappa value
#' @param d_order Moment order (0, 1, 2, 3)
#' @return Kappa value
get_kappa <- function(d_order) {
  key <- paste0("d", d_order)
  if (!key %in% names(.kappa_cache)) stop("Invalid d_order: ", d_order)
  unname(.kappa_cache[key])
}


#' Compute prior interval probability P(delta < eps) under moment prior
#'
#' @param d_order Moment order (0, 1, 2, 3)
#' @param eps Epsilon threshold (scalar or vector)
#' @param rel.tol Relative tolerance for integration
#' @return Vector of prior probabilities
prior_prob_delta <- function(d_order, eps, rel.tol = 1e-8) {
  kap <- get_kappa(d_order)

  compute_one <- function(e) {
    if (!is.finite(e) || e < 0) return(NA_real_)
    if (e == 0) {
      if (d_order > 0) return(0)
      return(prior_moment_density(0, d = 0, kappa = kap))
    }
    pdf <- function(z) prior_moment_density(z, d = d_order, kappa = kap)
    integrate(pdf, lower = 0, upper = e, rel.tol = rel.tol)$value
  }

  vapply(eps, compute_one, numeric(1))
}


#' Prepare moment-prior Stan data fields
#'
#' Computes the d_order, kappa_delta, df_delta, s_delta, logC_base
#' fields needed by fit_full_model.stan.
#'
#' @param d_order Moment order (0, 1, 2, 3)
#' @return Named list of prior data fields
moment_prior_stan_data <- function(d_order) {
  list(
    d_order = as.integer(d_order),
    kappa_delta = get_kappa(d_order),
    df_delta = as.integer(PRIOR_DF),
    s_delta = PRIOR_SCALE,
    logC_base = if (d_order == 0) 0 else log(C_d_half_t(d_order))
  )
}
