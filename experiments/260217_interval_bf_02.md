---
title: "Interval-Null Bayes Factor"
author: "Your Name"
date: "2026-02-17"
format: 
  html:
    toc: true
    toc-depth: 3
    code-fold: true
    theme: cosmo
    embed-resources: true
execute:
  echo: true
  warning: false
  message: false
---

-   [ ] TODO: write functions for the new method and try on a small dataset
-   [ ] TODO: modify the simulation code and run
-   [ ] TODO: derive the mathematical justification of SDR via limit

## Interval-Null Bayes Factor

Instead of testing

$$
H_0: \delta = 0
$$

approximating it by

$$
H_{0,\varepsilon}: \delta < \varepsilon
$$

and computing $$
BF_{(0,\varepsilon),1}
=
\frac{P(\delta < \varepsilon \mid y)}
{P(\delta < \varepsilon)}
$$

for different choices of ε.

### Advantages

#### 1. Avoid extrapolation bias in SDR

In SDR, you estimate

$$
\frac{p(\delta=0 \mid y)}{p(\delta=0)}
$$

But:

-   You never have posterior draws at exactly 0
-   You extrapolate density to the boundary
-   `logspline` at an exterior point may be biased

Interval probability avoids this entirely:

-   No extrapolation
-   Just count posterior draws below ε
-   Prior probability is analytic

So it is:

-   unbiased
-   Monte Carlo error easily computed (binomial variance)

#### 2. Works naturally for non-local priors

SDR is undefine for non-local priors $p(\delta=0) = 0$. But the interval BF works. And as $\varepsilon \to 0$, it converges to the true precise-null BF, but not the SDR (both numerator and denominator go to zero). So this gives you:

-   A principled way to compare local vs non-local priors
-   Without relying on unstable density-at-zero estimates

#### 3. Mathematical justification of SDR via limit

The interval-based BF is not a hack — it is theoretically justified.

To justify SDR as:

$$
BF_{01}
=
\lim_{\varepsilon \to 0}
\frac{P(\delta < \varepsilon \mid y)}
{P(\delta < \varepsilon)}
$$

This provides a clean mathematical motivation of SDR — especially important since you test an exterior point (δ ≥ 0). Please refer to Berger & Delampady (1987) [^1] and Wetzels et al. (2010) [^2]

[^1]: Berger, J.O., & Delampady, M. (1987). Testing Precise Hypotheses. Statistical Science, 2, 317-335.

[^2]: Wetzels, R., Grasman, R. P. P. P., & Wagenmakers, E.-J. (2010). An encompassing prior generalization of the Savage–Dickey density ratio. Computational Statistics & Data Analysis, 54(9), 2094–2102. https://doi.org/10.1016/j.csda.2010.03.016.

## Related Papers

### Berger & Delampady (1987): *Testing Precise Hypotheses* [^3]

[^3]: Berger, J.O., & Delampady, M. (1987). Testing Precise Hypotheses. Statistical Science, 2, 317-335.

-   They argued that testing a precise null is best understood as the limit of testing small intervals around the null.

    $$
      BF_{01}
      =
      \lim_{\varepsilon \to 0}
      \frac{P(|\theta| < \varepsilon \mid y)}
      {P(|\theta| < \varepsilon)}
      $$

-   They also propose practical ε choices: ε proportional to standard error (e.g. ε = SE/2, SE/4, SE/6). In Bayesian context, use posterior SD as proxy for uncertainty.

-   This paper gives theoretical grounding for interval-null approach and a principle way to choose $\epsilon$.

### Wetzels, Grasman, Wagenmakers (2010): *Encompassing prior generalization of SDR* [^4]

[^4]: Wetzels, R., Grasman, R. P. P. P., & Wagenmakers, E.-J. (2010). An encompassing prior generalization of the Savage–Dickey density ratio. Computational Statistics & Data Analysis, 54(9), 2094–2102. https://doi.org/10.1016/j.csda.2010.03.016.

-   Exact equality constrains are obtained as limits.

-   SDR is a special case of the encompassing prior approach.

-   Exact equality testing is delicate due to:

    -   Borel–Kolmogorov paradox
    -   Dependence on parameterization

-   Inequality Bayes factors: $$
    BF = \frac{\text{posterior proportion satisfying constraint}}
    {\text{prior proportion satisfying constraint}}
    $$

-   They show that $$
      BF
      =
      \frac{p(\theta_0 \mid y)}
      {p(\theta_0)}
      $$

    only if numerator and denominator go to zero at same rate.

-   For non-local priors, both densities are zero at the null. So SR is not valid. But interval approach remains valid.

-   This paper gives mathematical link between interval BF and SDR, a way to motivate my method rigorously, and a warning about non-local priors.

### Ly, Verhagen, & Wagenmakers (2016): *Harold Jeffreys’s default Bayes factor hypothesis tests: Explanation, extension, and application in psychology*. [^5]

[^5]: Ly, A., Verhagen, J., & Wagenmakers, E.-J. (2016). Harold Jeffreys’s default Bayes factor hypothesis tests: Explanation, extension, and application in psychology. Psychological Methods, 21(4), 463–478. https://doi.org/10.1037/met0000041

-   This paper explains interval null as “peri-null” testing.
-   Treat precise null as limit of shrinking interval
-   Provides computationally stable implementation
-   Argues interval null is philosophically coherent

## Implementation

1.  The choice of $\epsilon$

    (A) SD-based epsilons $$
         \varepsilon =
         \frac{\text{post.sd}(\delta)}{2},
         \frac{\text{post.sd}(\delta)}{4},
         \frac{\text{post.sd}(\delta)}{6}
         $$

    (B) Fixed small epsilons $$
         \varepsilon =
         0.05,\ 0.1,\ 0.2,\ 0.3,\ 0.4
         $$

2.  Compare across priors

    Compute interval BFs for:

    -   δ \~ t+
    -   δ² \~ t+
    -   δ³ \~ t+
    -   δ⁴ \~ t+

3.  Compute Monte Carlo error

    -   Posterior probability: $\hat p = k/N$
        -   $k$: the number of posterior draws that satisfy the interval constraint ($\delta^{(s)} < \epsilon$).
        -   $N$: the total number of posterior draws.
    -   Variance: $\hat p(1-\hat p)/N$\
        So you can give standard errors for BF estimates, which is a very strong methodological contribution.

4.  Check limit behavior

    Plot: $$
     BF_{(0,\varepsilon),1}
     \quad \text{vs} \quad \varepsilon
     $$

    As ε → 0:

    -   For local prior → converges to SDR
    -   For non-local prior → converges to true precise-null BF

    This is extremely publishable.

5.  Add mathematical section: “Precise Null as Limit of Interval Null”

    Derive: $$
     BF_{01}
     =
     \lim_{\varepsilon \to 0}
     \frac{P(\delta < \varepsilon \mid y)}
     {P(\delta < \varepsilon)}
     $$

**Suggestions**:

-   Do NOT remove SDR

-   Use bridge sampling as benchmark.

-   If results align:

-   You can argue SDR is safe in moderate regimes.

-   Interval BF supports your conclusions.

-   Non-local prior behaves better.

------------------------------------------------------------------------