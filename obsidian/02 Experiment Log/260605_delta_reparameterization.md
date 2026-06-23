---
date: 2026-06-05
allDay: true
title: Redefining `delta` after projection
tags:
  - non-local-prior
  - model-parameterization
  - figure1
files:
  - experiment/260605_delta_reparameterization_figures.R
  - experiment/260605_delta_reparameterization_fit_check.R
  - models/fit_full_model_projected_scaled.stan
---

# 260605_delta_reparameterization

## Question

I originally interpreted `delta` as the magnitude of the nonlinear effect **after** the orthogonalization step. But the current model does not define it that way. The question is whether `delta` should instead be redefined so that it directly controls the magnitude of the projected nonlinear component.

## Old definition of `delta`

In the current fitted model, `delta` is attached to the GP kernel **before** projection. The order is:

1. build the raw GP kernel using `delta` and `lambda`;
2. project the GP component onto the orthogonal complement of the linear space;
3. use the projected GP as the nonlinear effect.

So the old `delta` is a **pre-projection scale parameter**, not the magnitude of the final nonlinear effect.

## Why this became a problem

Under the old definition, fixing `delta` does **not** fix the size of the final nonlinear effect. Once the raw GP is projected, the visible nonlinear magnitude depends strongly on `lambda`.

This was first visible in the raw-vs-orthogonalized comparison:

📁 `experiment/260605_delta_reparameterization_figures.R`

![[Pasted image 20260605170756.png]]

For large `lambda`, the raw GP becomes very smooth and much of its variation lies close to the constant/linear subspace. After projection, that part is removed, so the remaining nonlinear effect can be much smaller than the input `delta` suggests.

The decomposition view makes the same point more explicitly:

📁 `experiment/260605_delta_reparameterization_figures.R`

![[Pasted image 20260605170851.png|607]]

Each raw draw can be split into:

- a projected linear part;
- an orthogonal residual.

For larger `lambda`, a substantial fraction of the raw GP variation is absorbed by the projected linear part.

## Calibration check under the old definition

To check whether this was only a data-generation artifact or also a fitting-interpretation problem, I ran a small calibration experiment:

- generate data under the old orthogonalized setup;
- fit the current model back;
- compare the projected nonlinear magnitude with the fitted posterior `delta`.

Source: `experiment/260605_delta_reparameterization_fit_check.R`
![[Pasted image 20260605170914.png|572]]
For fixed input `delta`, the expected projected RMS decreases sharply as `lambda` increases. So the same nominal `delta` does not correspond to the same final nonlinear magnitude.

Source: `experiment/260605_delta_reparameterization_fit_check.R`

![[Pasted image 20260605170944.png|500]]
This means the old fitted `delta` is also difficult to interpret directly as an effect size. It is upstream of the final projected nonlinear effect.

## New definition of `delta`

I then introduced a reparameterized model:

📁 `models/fit_full_model_projected_scaled.stan`

The idea is:

1. build a unit-amplitude kernel from `lambda`;
2. project that kernel onto the nonlinear subspace;
3. standardize the projected kernel so that its average marginal variance is 1;
4. let `delta` scale this **already projected and standardized** kernel.

Under this definition, `delta` is intended to represent the magnitude of the **final projected nonlinear effect**.

## What changed after redefining `delta`

The first check was whether the projected RMS is now stable across `lambda`:

📁 `experiment/260605_delta_reparameterization_figures.R`
![[Pasted image 20260605171143.png|777]]
This is the key result. Under the new definition, the projected nonlinear RMS is flat in `lambda` for fixed `delta`. That is exactly the interpretation I originally wanted.

The draw-level comparison also became cleaner:

Source: `experiment/260605_delta_reparameterization_figures.R`
![[Pasted image 20260605171023.png|1176]]
Here:

- `lambda` mainly changes shape/smoothness;
- `delta` mainly changes magnitude.

This separation was not true under the old definition.

## Figure 1 replacement candidate

I then made a Figure-1-like plot under the new definition, using:

- `delta = 0, 1, 2, 3`
- `lambda = 0.15, 0.25, 0.4`
- symmetric `y` scale `[-10, 10]`

📁 `experiment/260605_delta_reparameterization_figures.R`
![[Pasted image 20260605171159.png|836]]
At this stage, this is the strongest candidate for a Figure 1 replacement, because it matches the interpretation I want:

- `delta = 0, 1, 2, 3` correspond to null / small / medium / large nonlinear effect;
- `lambda` controls the smoothness of that nonlinear effect;
- the plotted functions are already in the nonlinear subspace used by the fitted model.

## Current takeaway

The old `delta` was not wrong mathematically, but it was not the parameter I thought I was interpreting. It controlled the GP **before** projection, while the scientific quantity of interest is the nonlinear effect **after** projection.

The reparameterized `delta` is therefore preferable if I want:

- a clean interpretation of posterior `delta`;
- a Figure 1 that can be used for prior elicitation;
- a clearer separation between nonlinear magnitude (`delta`) and smoothness (`lambda`).
