---
date: 2026-06-01
allDay: true
title: Reweighting Method for Nonlocal Priors
tags:
  - non-local-prior
files:
  - experiment/260601_reweighting_prior.R
---



> [!SUMMARY] Highlight
> It solves two things:
> - We only need to fit the model once under the local/base prior and can then obtain the corresponding BF under the non-local prior.
> - We no longer need to directly estimate the non-local posterior behavior near 0, which is helpful because of the very little mass there.
> 
> However, one important problem remains:
> - If the base posterior itself is already far from the null, the estimation is still unstable.


# Computing Moment-Prior Bayes Factors via Reweighting

## Setup

The moment prior of order $d$ on the nonlinear magnitude $\delta$ is

$$\pi^M(\delta) = \frac{\delta^d \; \pi^B(\delta)}{C}, \quad C = \int \delta^d \; \pi^B(\delta) \, d\delta,$$

where $\pi^B$ is a local base prior (half-$t$, positive at zero).

## Core identity (supervisor's derivation)

Using Bayes' theorem to substitute $p(y|\delta)\,\pi^B(\delta) = p^B(y|M_1)\,\pi^B(\delta|y)$:

$$p^M(y \mid M_1) = \int p(y|\delta)\,\pi^M(\delta)\,d\delta = p^B(y \mid M_1) \cdot \mathbb{E}^B_{\delta|y}\!\left[\frac{\pi^M(\delta)}{\pi^B(\delta)}\right].$$

Dividing both sides by $p(y \mid M_0)$:

$$BF^M_{10} = BF^B_{10} \cdot \mathbb{E}^B_{\delta|y}\!\left[\frac{\pi^M(\delta)}{\pi^B(\delta)}\right].$$

When the base prior used for fitting equals the base inside $\pi^M$, the ratio simplifies to $\delta^d / C$, recovering the formula on the handwritten note:

$$BF^M_{10} = BF^B_{10} \cdot \frac{\mathbb{E}^B_{\delta|y}[\delta^d]}{\mathbb{E}^B[\delta^d]}.$$

## Our implementation: fixed base, general reweighting

In our construction, $\pi^M$ of order $d$ uses a median-calibrated base $g_{\kappa(d)}$, where $\kappa$ depends on $d$. Rather than refitting the base model for each $d$, we fit **once** under a fixed base $\pi^{B_0} = g_{\kappa(0)}$ (the $d{=}0$ local prior) and apply the general form:

$$BF^{\pi_d}_{10} = BF^{B_0}_{10} \cdot \frac{1}{S}\sum_{s=1}^{S} \frac{\pi^M_d(\delta^{(s)})}{\pi^{B_0}(\delta^{(s)})},$$

where $\{\delta^{(s)}\}$ are posterior draws from the $d{=}0$ base fit. Both $\pi^M_d$ and $\pi^{B_0}$ have known analytic densities, so the weight at each draw is a simple density ratio evaluation.

## Verification

Across simulated datasets ($\delta_{\mathrm{true}} \in \{0, 1, 3\}$, $N \in \{20, 50, 100, 200\}$, $d \in \{1, 2, 3\}$), the reweighted $BF$ agrees with bridge-sampling estimates to within $\sim$0.03 on the log-BF scale.


## Current conclusion

- Reweighting itself seems correct: using bridge sampling as the reference, the reweighted BF is almost identical to the directly computed BF.
- The main practical gain is that I only need to fit the **local/base prior** model once, and can then obtain the corresponding BF under the non-local prior by reweighting.
- Another gain is that I no longer need to directly estimate the non-local posterior behavior near 0. This avoids the worst instability caused by the very small near-zero mass of the non-local prior itself.
- However, this still does **not** solve the deeper problem that when the posterior is already far from the null region, the null-neighborhood quantity is difficult to estimate accurately. So the reweighting solves the prior-order conversion problem, but not the general extrapolation / near-null estimation problem.

Empirically:

- Using bridge sampling as the benchmark, the reweighted BF and the directly computed BF are essentially the same, so the reweighting step itself looks correct.
- For the current working version, the most satisfactory route is to fit the **local/base prior** model and estimate its SDR with **logspline**. This still has the usual extrapolation concern, but among the methods explored so far it behaves best.
- For the **non-local prior**, the current plan is therefore to use the reweighting identity above to obtain its BF from the base fit, instead of directly estimating the non-local posterior near 0.
![[Pasted image 20260616111919.png]]

- I also explored interval-BF variants. The raw interval method is well defined, but becomes unstable when there is essentially no posterior mass in the null interval.
- For Laplace and Jeffreys smoothing, I have **not** found a satisfactory reweighted analogue: the raw weighted interval probability is clear, but once smoothing is introduced, it is unclear how to reweight it in a way that both preserves the intended regularization and recovers the corresponding direct-fit results.
![[Pasted image 20260616111930.png]]

## Working decision

- Local prior: fit directly and use logspline to compute the SDR.
- Non-local prior: compute the BF by reweighting from the local/base fit.
- Treat the interval-BF route, including the Laplace/Jeffreys smoothing variants, as exploratory only for now.
