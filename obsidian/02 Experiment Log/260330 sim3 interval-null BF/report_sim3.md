 

# Simulation Report: Interval-Null BF with Non-local Priors

**Date:** 2026-03-30   |   **Simulation grid:** N in {20, 50, 100, 200}, delta in {0, 1, 2, 3}, prior order d in {0, 1, 2, 3}, 50 replicates each.

**Interval-Null BF:** Instead of testing delta = 0, we test delta < epsilon. BF10(eps) = P(delta < eps | prior) / P(delta < eps | y). Two posterior estimation methods: _logspline_ (density-based) and _raw_ (empirical proportion k/N). Bridge sampling BF serves as the gold-standard benchmark.

## 1. BF Convergence with Epsilon

How does log BF10 behave as epsilon increases? Triangles mark the first epsilon where the raw method yields a finite BF (left of triangle = Inf). Dotted horizontal lines = bridge sampling benchmark.

### Logspline method

![](HTML%20import/Attachments/fig1a_bf_vs_eps_logspline.png)

### Raw method (with IQR ribbons)

![](HTML%20import/Attachments/fig1b_bf_vs_eps_raw.png)

**Finding:** Curves stabilize around eps = 0.2. At smaller eps, the raw method shows Inf (no posterior draws below eps). The logspline method is smoother but underestimates BF at large delta due to boundary extrapolation.

## 2. Effective Convergence Epsilon

For the raw method: what is the smallest epsilon with any posterior mass below it? This diagnostic shows where the raw interval BF "turns on."

![](HTML%20import/Attachments/fig2_convergence_eps.png)

**Finding:** Convergence eps increases with delta (posterior moves away from zero) and with prior order (non-local priors push mass away from zero). At delta = 0, convergence eps is near 0.01 for all priors and N.

## 3. BF at Fixed Epsilon = 0.2

Bars = median log BF10 over 50 replicates at eps = 0.2. Error bars = IQR (finite replicates only). x = bridge sampling median. Bars labeled "Inf" indicate majority of replicates have zero posterior mass below eps.

### Logspline

![](02 Experiment Log/260330 sim3 interval-null BF/Attachments/fig3a_bf_eps02_logspline.png)

### Raw

![](HTML%20import/Attachments/fig3b_bf_eps02_raw.png)

**Finding:** Logspline gives graded, finite evidence across all settings. Prior order effect is visible: d=0 (local) gives strongest BF10, d=3 weakest. Raw method saturates at Inf for delta >= 1 at small N, but gives finite values at N = 200 for delta = 1. Under the null (delta = 0), both methods correctly return negative log BF (evidence against nonlinearity), and higher prior order gives more decisively negative values.

## 4. Evidence Accumulation Across Methods

Each panel shows BF vs N (top row) and BF vs delta (bottom row). Comparing bridge sampling, logspline, and raw side by side.

### Bridge Sampling (benchmark)

![](HTML%20import/Attachments/fig4_bridge_combined.png)

### Logspline (eps = 0.2)

![](HTML%20import/Attachments/fig4_logspline_combined.png)

### Raw (eps = 0.2)

![](HTML%20import/Attachments/fig4_raw_combined.png)

### Bridge BF by Prior Order (exact values)

Bridge sampling BF is nearly prior-order-invariant. Differences are real but small (< 0.6 on log scale):

|delta|N|d=0 (local)|d=1|d=2|d=3|
|---|---|---|---|---|---|
|0|20|-0.797|-0.898|-0.814|-0.934|
|0|50|-1.087|-1.147|-1.234|-1.247|
|0|100|-1.290|-1.508|-1.632|-1.695|
|0|200|-1.427|-1.874|-2.007|-1.999|
|1|20|0.038|0.096|0.105|0.099|
|1|50|1.111|1.209|1.219|1.217|
|1|100|-0.365|-0.428|-0.472|-0.502|
|1|200|7.146|7.238|7.229|7.215|
|2|20|0.751|0.857|0.884|0.881|
|2|50|7.753|7.940|7.997|8.001|
|2|100|9.481|9.673|9.724|9.729|
|2|200|27.655|27.883|27.944|27.956|
|3|20|1.515|1.643|1.681|1.683|
|3|50|16.199|16.350|16.399|16.397|
|3|100|27.228|27.460|27.530|27.538|
|3|200|54.096|54.352|54.429|54.438|

**Finding:** Bridge sampling shows near-identical curves across prior orders (max difference < 0.6 in log BF). This is because the priors are calibrated to the same median; they differ only near delta = 0, which barely affects the marginal likelihood integral. The interval BF retains sensitivity to prior order because it specifically evaluates evidence near the null.

## 5. Method Comparison: Logspline vs Raw

Scatter of median log BF10 at eps = 0.2. Points with Inf in either method are excluded. Dashed = identity line.

![](HTML%20import/Attachments/fig5_method_comparison.png)

**Finding:** Where both methods are finite, they agree well (points near the diagonal). The raw method is unbiased but has lower resolution; logspline provides smoother estimates but can be biased near the boundary.

## 6. CV-Optimal Epsilon Selection

Principled epsilon: choose the smallest eps where SE(p_hat)/p_hat < 20%. This ensures the posterior probability estimate (denominator of BF) is precise enough to trust.

### 6a. Distribution of CV-optimal epsilon

![](HTML%20import/Attachments/fig6a_cv_optimal_eps.png)

### 6b. BF at CV-optimal epsilon

![](HTML%20import/Attachments/fig6b_bf_cv_optimal.png)

### 6c. Per-replicate comparison vs bridge sampling

![](HTML%20import/Attachments/fig6c_scatter_cv_vs_bridge.png)

**Finding:** CV-optimal eps ranges from ~0.05 (delta=0, strong posterior mass near zero) to ~0.4 (delta=3, posterior far from zero). The interval BF at CV-optimal eps tracks bridge sampling well in moderate regimes but remains conservative for large effects.

## 7. Summary

1. **Bridge sampling** is the gold standard: works everywhere, insensitive to prior order, but requires fitting a separate null model.
2. **Interval BF (logspline, eps = 0.2)** is the practical workaround: requires only the alternative model, gives graded evidence, and retains sensitivity to prior order near the null.
3. **eps ~ 0.2** is a defensible operational choice; the CV-based adaptive criterion provides a principled alternative.
4. **Non-local priors** (d > 0) are more conservative under both H0 and H1: stronger evidence for the null when delta = 0, slightly weaker evidence for the alternative otherwise.
5. **Next steps:** Run simulations for N = 500, 1000. Derive the L Hopital limit connecting interval BF to precise-null BF. Refit empirical applications with updated Stan model (v17, full-X projection).