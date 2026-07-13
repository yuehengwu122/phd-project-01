---
date: 2026-05-28
allDay: true
title: ARD alternative
tags:
  - application
file: analysis/application/03_diabetes_ard_gpboost.R, analysis/application/04_boston_housing_ard_gpboost.R
---
## Goal
- Add a standard ARD-based GP baseline for application datasets.
- Compare ARD patterns against the current additive GP method.

## Datasets
- App3: Diabetes
- App4: Boston Housing

## Method
- Package: `gpboost`
- Covariance: `gaussian_ard`
- Fixed effect: intercept only
- Output focus: predictor-wise estimated range / inverse range

## Preprocessing
- Use the same predictor sets as the current application scripts.
- Exclude binary / non-continuous predictors already dropped in the main analyses.
- Standardize each predictor by z-score before fitting.

## Outputs
### Diabetes
ARD ranking by inverse range:

| **predictor** | **inv_range_ard** |
| ------------- | ----------------- |
| s5            | 0.262             |
| bmi           | 0.173             |
| age           | 0.139             |
| bp            | 0.107             |
| s6            | 0.083             |
| s3            | 0.055             |
| s1            | 0.050             |
| s4            | 0.000             |
| s2            | 0.000             |
- Shortest range / largest inverse range: `s5`, `bmi`, `age`, `bp`, `s6`
- Longest range / smallest inverse range: `s4`, `s2`

Our method:
- strong linear evidence: `bmi`, `bp`, `s5`
- nonlinear evidence mostly weak; only `s6` is close to neutral, others are mostly negative

### Boston Housing
ARD ranking by inverse range:

| **predictor** | **inv_range_ard** |
| ------------- | ----------------- |
| tax           | 1.054             |
| dis           | 0.792             |
| lstat         | 0.566             |
| nox           | 0.450             |
| black         | 0.312             |
| rm            | 0.243             |
| age           | 0.172             |
| indus         | 0.133             |
| ptratio       | 0.098             |
| crim          | 0.093             |
| zn            | 0.003             |
- Shortest range / largest inverse range: `tax`, `dis`, `lstat`, `nox`, `black`, `rm`
- Very long range: `zn`

Our method:
- strong nonlinear evidence: `crim`, `indus`, `nox`, `rm`, `dis`, `tax`, `lstat`
- weak nonlinear evidence: `zn`, `age`, `ptratio`, `black`
- strong linear evidence: `crim`, `zn`, `nox`, `rm`, `dis`, `ptratio`, `black`, `lstat`

## Initial Observations
Diabetes:
- ARD does not provide clear nonlinear separation.
- `s5`, `bmi`, and `bp` are ranked highly by ARD, but in our method they are mainly supported through strong linear evidence rather than strong nonlinear evidence.
- This suggests that short ARD lengthscales can reflect overall predictor influence without distinguishing whether that influence is linear or nonlinear.

Boston Housing:
- ARD shows clearer qualitative agreement with the nonlinear results.
- Predictors such as `tax`, `dis`, `lstat`, `nox`, and `rm` have relatively short ARD ranges and also show strong nonlinear evidence under our method.
- The agreement is not exact. `black` is a useful counterexample: it has a relatively short ARD range, but weak nonlinear evidence and stronger support for a linear effect in our method.

## Comparison To Current Method
- ARD provides one summary per predictor through the GP lengthscale.
- The current method provides predictor-specific assessment of linear and nonlinear effects separately.
- The comparison is therefore qualitative / partial: check whether short ARD ranges tend to align with predictors showing stronger nonlinear evidence.
- Exact one-to-one agreement is not expected, because the two quantities target different aspects of predictor contribution.

## Next Steps
- Check stability of ARD rankings across datasets.
- Decide how to present ARD vs. Bayes-factor-based evidence in the paper.

## In Paper
> *As a GP-based baseline, we compare our method to automatic relevance determination (ARD) implemented through an anisotropic Gaussian process. As discussed in the introduction, ARD does not formally separate linear and nonlinear effects, and its lengthscale parameters should not be interpreted as definitive measures of variable relevance. We therefore use ARD only as a coarse GP-based sensitivity summary, rather than as a direct analogue of our predictor-specific evidence measures.*

> *We do not expect exact one-to-one agreement between ARD ranges and our nonlinear Bayes factors, because the two quantities target different aspects of predictor contribution. ARD summarizes overall covariance sensitivity in a GP, whereas our method explicitly tests linear and nonlinear effects separately.*

### Dataset-specific wording
> *In the diabetes data, ARD highlighted predictors that our method also identified as important, but it did not distinguish whether their contribution was primarily linear or nonlinear. In particular, predictors such as `bmi`, `bp`, and `s5` showed strong linear evidence in our analysis, whereas ARD summarized them only through short GP lengthscales.*

> *In the Boston housing data, ARD showed clearer qualitative agreement with our nonlinear results: predictors with strong nonlinear evidence under our method, such as `nox`, `rm`, `dis`, `tax`, and `lstat`, also tended to have shorter ARD ranges. However, the agreement was not exact. For example, `black` received a relatively short ARD range despite weak nonlinear evidence in our method and stronger support for a linear effect. This illustrates that ARD can capture overall sensitivity of the GP surface, but it does not decompose that sensitivity into linear versus nonlinear contributions.*
