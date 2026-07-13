
## Paper Refine Standards
1. **Academic Flow & Coherence**: Ensure the modified paragraph reads naturally and seamlessly. The newly introduced idea must be logically integrated, properly contextualized, and smoothly transitioned.
2. **Eliminate Redundancy & Contradictions**: Scan the surrounding sentences to ensure the new addition does not create repetitive arguments, awkward phrasing, or logical inconsistencies with the rest of the paragraph.
3. **Translate Concept into Prose**: If I pitch a rough or conceptual idea, translate it into precise, rigorous academic prose that fits the sophisticated tone of the paper.
4. **Output Format**: Provide the fully revised, complete paragraph directly so I can copy-paste it. Underneath, briefly point out the strategic adjustments you made to preserve the logical flow.


## Todo Round1 June29-July2
- [x] 2.2 Lack of an introduction of delta. See the narrative of of 2.2 Nonlinear Component (below), there is no introduction of delta. Or let's say, there is a gap between the covariance of GP and the kernel. In the formula at the beginning of section 2, we define the covariance as sigma2\*delta2\*K(lambda). Therefore, it's important to mention it here and emphasis the the "standardization" of delta here. 
- [x] 3.1.2 Follow the point of 2.2, because delta is standardized, we can propose a more "objective" prior on it, but do not use the word "objective" 
- [x] 3.1.2 For nonlocal priors, I want to add a point to make readers now we propose it just as an alternative priors (not in a strong feeling, because in my simulation study, the distinguishability is not VERY large)
- [x] 3.1.2 Be careful to make the four equivalent concepts clear: half-t prior, horseshoe prior, local prior, base prior.
- [x] 3.2.2 We have made the decision: use `logspline` to estimate interval-null BF. Therefore, the raw method is not needed anymore. We can remove it from this part. 
- [x] 3.2.2 The choice of epsilon: we still do not find a proper way to explain why we choose eps=0.2 (maybe add a prior-predictive / function-draw argument，instead of only eyeballing observed plot)
- [x] Add a new section after 3.2.2, the new reweighing method for nonlocal priors 
- [x] Add Simulation section
- [x] Add and explain the results for 5.1
- [x] complete abstract and discussion
5.3 &5.4
- [x] add results
- [x] add results of alternative method (ADR)
- [x] also compare with the results mentioned in previous work
- [ ] add bridge sampling results (not available yet)


## Round2
- [x] `rasmussen2006gaussian` → should be **`williams2006gaussian`** (Rasmussen & Williams 2006, _Gaussian Processes for Machine Learning_, MIT Press).
- [ ] Add prior and posterior samples in Appendix

- [x] **Full compile-health pass** — I can check the whole document (not just the intro) for any remaining undefined refs, orphaned labels, figure files that won't resolve (e.g. one figure points to `../results/fig_01_neuro_trend.pdf`, an outside-tree path worth confirming), and math/environment issues.
- [x] **Bibliography audit across the whole paper** — verify every one of the 59 cited entries against its source (authors/year/venue), flag any remaining mismatches, and check the 33 defined-but-uncited entries so you can decide what to prune.
	- Overall the bibliography is in good shape. No missing citations, no undefined keys, and author metadata is clean across all 59. 
- [x] **Section-by-section literature check** — extend what we did for the intro to the Discussion (which cites the scalability and extension literature) and the Empirical Application sections, so every claim is correctly attributed.
- [x] **Narrative / argument review** — read a full section (Modeling framework, or the Simulation/Application sections) for clarity, logical flow, and whether the claims are supported by what's actually shown.
	- Issues worth fixing
	- [x] **θ dimension error (L124)** — `[\theta_1,...,\theta_j]^T` uses the summation index `j` as the vector's terminal index. Since **X** is n×(p+1), θ is a (p+1)-vector; should be `[\theta_0,\theta_1,\dots,\theta_p]^T`.
	- [x] **Intercept handling never stated** — X carries an intercept column and the g-prior covariance g(XᵀX)⁻¹ shrinks it, which is unusual; and the linear tests presumably skip the intercept. One clarifying sentence in §2.1 would close this, and it also reconciles the `β_0 + β_1 x` notation used later in the mother's IQ model (L987).
	- [x] **Unsupported cross-reference (L1036)** — "Following the analysis setup used in Section 2, we exclude the binary predictor `sex`" — §2 never describes dropping binary predictors. Drop the cross-ref or point it where the convention is actually set.
	- [x] **Appendix title renders as "A Appendix A" (L1204)** — `\appendix` auto-numbers, so give it a real title like "Limit of the interval-null Bayes factor."
	- [x] **"Multivariate normal" prior imprecise (L594)** — the SDR uses the _marginal_ g-prior (a scaled t(4) after integrating g out), not the conditional MVN. The computation is fine; the wording can mislead.


## Round3 go through the paper



## Supplements
#### section 2.2 
The narrative of of 2.2 Nonlinear Component:
- GP prior on f
	- mean and covariance
	- SE kernel and the role of lambda
	- the role of delta, control the covariance
	- Figure 1 shows the interaction between delta and lambda
- Issue: confounding between L and NL
- propose orthogonalized GP



#### Thoughts
- [ ] App: BF10 or logBF10?
AA00FNM4IN



