# 2026-02-17 — Modify `estimate_post_density_at_zero()`

## Before
``` {r}
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
```

- Once you discard draws, you are no longer fitting the posterior density f(x); you are fitting the conditional density given the trimming event.


## Problem

Let the true target density be f(x) on $[0,\infty)$. Suppose you trim at a cutoff $c$ (e.g., the empirical q-quantile) and keep only $x \le c$. The distribution you are fitting is
$$
f_{\text{trim}}(x) \;=\; f(x \mid x \le c) \;=\; \frac{f(x)}{\Pr(x \le c)} \quad \text{for } x \le c.
$$

If $0$ is inside the kept region (it is), then
$$
f_{\text{trim}}(0) \;=\; \frac{f(0)}{\Pr(x \le c)}.
$$

So if you return `dlogspline(0, fit_trim)`, you are estimating $f_{\text{trim}}(0)$, which is inflated relative to $f(0)$ by the factor $1/\Pr(x \le c)$.
- If you trim at the 0.9 quantile, $\Pr(x\le c)\approx 0.9$, so your estimate at 0 is inflated by $\approx 1/0.9 \approx 1.11$ (≈11% too large).
- If you trim at 0.7, inflation is $\approx 1/0.7 \approx 1.43$.
- If you trim at 0.99, inflation is $\approx 1.01$ (tiny).

## Fix

If you trim at quantile level q, then $\Pr(x \le c)\approx q$. To recover the original-scale density at 0,
$$
f(0) \approx q \, f_{\text{trim}}(0).
$$


## After
```{r}
estimate_post_density_at_zero <- function(x, lbound = 0, ubound = NULL,
                                          max_attempts = 7) {

  quantiles <- c(1, 0.999, 0.99, 0.9, 0.8)

  for (q in quantiles[1:max_attempts]) {

    cutoff <- as.numeric(stats::quantile(x, q, names = FALSE, type = 7))
    keep   <- x <= cutoff
    x_trim <- x[keep]

    logspline_args <- list(x_trim, lbound = lbound)
    if (!is.null(ubound)) logspline_args$ubound <- ubound

    fit <- try(do.call(logspline::logspline, logspline_args), silent = TRUE)
    if (!inherits(fit, "try-error")) {

      # density under conditional (trimmed) distribution
      f0_trim <- logspline::dlogspline(0, fit)

      # renormalize back to the original distribution
      p_keep <- mean(keep)  # ~ q, but explicit
      return(p_keep * f0_trim)
    }
  }

  warning("logspline failed for all tail-cuts. Returning NA.")
  NA_real_
}
```

