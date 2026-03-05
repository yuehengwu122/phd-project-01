---
title: "Test on Linear Effects"
date: "2026-02-17"
format: 
  html:
    toc: true
    toc-depth: 3
    theme: cosmo
    embed-resources: true
---

## Prior

$$
\theta \mid g \sim \mathcal N\!\left(0,\; g (X_c^\top X_c)^{-1}\right), \\
\log g \sim \mathcal N(\log N,\, 0.5^2)
$$

So conditional on $g$, the **marginal prior** for component $\theta_j$ is $$
\theta_j \mid g \sim \mathcal N\!\left(0,\; g\, v_j\right),\qquad v_j = \left[(X_c^\top X_c)^{-1}\right]_{jj}.
$$

Crucially, $v_j$ comes from the full $(X_c^\top X_c)^{-1}$, which reflects predictor correlations (multicollinearity inflates $v_j$).

## Prior density at 0 for $\theta_j$

For SDR we need the marginal prior ordinate $\pi(\theta_j=0)$, i.e. after integrating out g.

Because $$
\pi(\theta_j=0\mid g)=\frac{1}{\sqrt{2\pi g v_j}},
$$

we have $$
\pi(\theta_j=0)=\mathbb E_g\!\left[\frac{1}{\sqrt{2\pi g v_j}}\right]
=\frac{1}{\sqrt{2\pi v_j}}\; \mathbb E_g[g^{-1/2}].
$$

With $\log g \sim \mathcal N(m,s^2)$, $$
\mathbb E[g^{-1/2}] = \exp\!\left(-\tfrac12 m + \tfrac18 s^2\right).
$$

In our prior, $m=\log N$, $s=0.5$, so $$
\mathbb E[g^{-1/2}] = \exp(0.03125)\, N^{-1/2}.
$$

Therefore $$
\boxed{
\pi(\theta_j=0)=\frac{\exp(0.03125)}{\sqrt{2\pi\, v_j\, N}}
}
\qquad\text{where } v_j=\left[(X_c^\top X_c)^{-1}\right]_{jj}.
$$

## R Implementation

``` r
prior0_theta_analytical <- function(X, N = nrow(X), 
                         mu_log_g = log(N), sd_log_g = 0.5) {
  # input rescaled X in [-0.5, 0.5]
  X_c <- scale(X, center = TRUE, scale = FALSE)

  XtX_inv <- solve(crossprod(X_c))  # (X'X)^{-1}
  v <- diag(XtX_inv)

  Eg_inv_sqrt_g <- exp(-0.5 * mu_log_g + 0.125 * sd_log_g^2)  # E[g^{-1/2}]
  prior_d0 <- (1 / sqrt(2 * pi * v)) * Eg_inv_sqrt_g

  prior_d0
}
```

￼

## Experiments

**1. Monte Carlo sanity check of the analytic ordinate**\
Sample $g\sim \log\mathcal N(m,s^2)$, then for each g evaluate `dnorm`(0,0,$\sqrt{g v_j}$), average over draws; compare to the closed form above. This verifies the prior-ordinate code independently.

``` r
X <- fit_a3_dep$data$X
X_c <- scale(X, center = TRUE, scale = FALSE)
v <- diag(solve(t(X_c) %*% X_c))

N <- dim(X)[1]
P <- length(v)

N_g <- 1e3
log_g_theta_draws <- rnorm(N_g, mean = log(N), sd = 0.5)

prior0_simu <- prior0_simu_se <- numeric(P)
prior0_analy <- prior0_theta_analytical(X)

for (i in 1:P) {
    vals <- dnorm(0, mean=0, sd=sqrt( exp(log_g_theta_draws) * v[i] ))
    prior0_simu[i] <- mean(vals)
    prior0_simu_se[i] <- sd(vals) / sqrt(N_g)
} 

round(data.frame(prior0_analy, prior0_simu, prior0_simu_se), 4)
```

| Variable | Analytic Prior Ordinate | Simulated Prior Ordinate | Simulation SE |
|---------------|---------------------|---------------------|---------------|
| age      | 0.0818                  | 0.0808                   | 0.0007        |
| bmi      | 0.0617                  | 0.0609                   | 0.0005        |
| bp       | 0.0675                  | 0.0666                   | 0.0005        |
| s1       | 0.0091                  | 0.0090                   | 0.0001        |
| s2       | 0.0100                  | 0.0098                   | 0.0001        |
| s3       | 0.0176                  | 0.0174                   | 0.0001        |
| s4       | 0.0252                  | 0.0249                   | 0.0002        |
| s5       | 0.0238                  | 0.0234                   | 0.0002        |
| s6       | 0.0589                  | 0.0582                   | 0.0005        |

The Monte Carlo estimator targets the same integral as the analytic formula.

**2. Check whether your posterior ordinate is on the same parameter**\
You must estimate $\pi(\theta_j=0\mid y)$ using draws of $\theta_j$, not $\beta_{1j}$ (since $\beta_{1j}=\sigma\theta_j$ changes the null scaling). Your Stan defines $\theta$ explicitly, so keep SDR on $\theta$.

``` r
theta_draws <- rstan::extract(fit_a3_dep$fit)$theta

post0 <- sapply(1:P, function(j) {
    fit <- logspline::logspline(theta_draws[, j], 
                                lbound = min(min(theta_draws[, j]), 0), 
                                ubound = max(max(theta_draws[, j]), 0))
    logspline::dlogspline(0, fit)
})

log_sdr_10 <- log(prior0_analy / post0)

format(data.frame(prior0_analy, post0, log_sdr_10, v), digit = 3, scientific = TRUE)
```

| Variable | Analytic Prior Ordinate | Posterior Ordinate | log(SDR10) | v        |
|--------------|-------------------|--------------|--------------|--------------|
| age      | 8.18e-02                | 1.60e+00           | -2.98e+00  | 5.73e-02 |
| bmi      | 6.17e-02                | 1.32e-08           | 1.54e+01   | 1.01e-01 |
| bp       | 6.75e-02                | 7.14e-04           | 4.55e+00   | 8.41e-02 |
| s1       | 9.08e-03                | 3.98e-02           | -1.48e+00  | 4.65e+00 |
| s2       | 9.96e-03                | 4.85e-02           | -1.58e+00  | 3.87e+00 |
| s3       | 1.76e-02                | 1.43e-01           | -2.10e+00  | 1.24e+00 |
| s4       | 2.52e-02                | 3.67e-01           | -2.68e+00  | 6.05e-01 |
| s5       | 2.38e-02                | 3.38e-02           | -3.51e-01  | 6.80e-01 |
| s6       | 5.89e-02                | 8.37e-01           | -2.65e+00  | 1.10e-01 |

**3. Watch** $v_j$ as a “collinearity amplifier”\
If some predictors are correlated, $v_j$ inflates, making $\pi(\theta_j=0)$ smaller (more diffuse prior on $\theta_j$). That will directly affect `SDR10`. Printing $v_j$ alongside your SDR table often immediately explains “why this variable’s BF is weaker/stronger than expected.”

-   $v_j$ is the prior variance “inflation factor” for $\theta_j$ under the multivariate g-prior.
-   Large v_j typically means predictor j is close to a linear combination of others (collinearity), so the prior for \theta\_j becomes more diffuse and the data need to work harder to pin down that coefficient.
-   `corr(s1, s2) = 0.90` → classic collinearity inflation: huge (v1,v2) = (4.65, 3.87)