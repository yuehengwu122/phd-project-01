# Project 2 构想:HSGP 近似下的 Bayes Factor 检验

日期:2026-05-28

## 一句话 pitch

> exact-GP 的 linear-vs-nonlinear Bayes factor 检验不可扩展;HSGP 能解决速度,但近似会不会破坏 BF 的结论?我们刻画它,给出实用准则,并在真实数据上验证。

这样它既是一篇 application paper,又带一个站得住的方法论 claim。

## 背景与动机

- Project 1 的核心:additive GP + Bayes factor,用于检验 linear / nonlinear effect(bridge sampling + Savage-Dickey)。
- 问题:exact GP fit 非常慢,缺乏实用价值。
- 想法:用 approximated GP(HSGP, Hilbert Space GP)替代,写一篇偏 application 的 paper。
- 前期探索:`experiments/260104_hsgp_experiment.R` + `stan/multiple_gp/fit_add_hsgp_01.stan`,用 **Matérn kernel** 结果很好;**SE kernel**(project 1 当前所用)fit 不出来。

## 方向判断

**方向合理。** HSGP 是把 exact GP 慢的问题变得"真的能跑在真实数据上"的自然切入点。

**要避开的陷阱:** 如果只是"把 exact GP 换成 HSGP 然后跑几个数据集",HSGP 不是自己提出的,容易被 reviewer 当成 incremental。

**真正的贡献点:** 在 HSGP 近似下,Bayes factor / 模型选择结论还成立吗?
- marginal likelihood 对近似的超参(基函数个数 M、边界 L、kernel)敏感;
- 而 BF 又恰恰依赖 marginal likelihood;
- "近似对 BF 结论的影响 + 什么时候会崩"这个 angle,把工程性工作抬成研究。

## SE vs Matérn(可写进 paper 的现象,而非 bug)

- SE kernel 的 spectral density 指数衰减 → 有限基函数逼近时高频项几乎为零,对 length-scale 极度敏感,容易塌缩 / 不可识别。
- Matérn 的 spectral tail 多项式衰减 → 有限基函数逼近稳健得多。
- 这也是文献里推荐 HSGP 配 Matérn 的常见理由,可作为论证里的支撑点。

## 空白确认(2026-05-28 web 搜索结论)

- 搜 "HSGP + Bayes factor / marginal likelihood / 检验":**没有任何论文把两者结合**。
- Riutort-Mayol (2023) 明确说近似精度取决于基函数个数 M,但只在**预测 / 拟合精度**层面讨论,**完全没碰 marginal likelihood / BF**。
- marginal likelihood 对模型有效维度(≈ 基函数个数)出了名地敏感。
- → 空白清晰且真实:**没人量化过 HSGP 的近似 knob(M、boundary、kernel)如何影响 Bayes factor 结论的稳定性。**

## 推荐阅读(按相关性排序)

### 1. HSGP 本身(打地基)
- **Solin & Särkkä (2020)**, *Hilbert space methods for reduced-rank Gaussian process regression*, Statistics and Computing 30, 419–446. arXiv:1401.5508 — HSGP 理论原始文献(sine basis、boundary L、spectral density 来源)。
- **Riutort-Mayol, Bürkner, Andersen, Solin, Vehtari (2023)**, *Practical Hilbert space approximate Bayesian Gaussian processes for probabilistic programming*, Statistics and Computing 33(1). arXiv:2004.11408 — **先读这篇**;选 M / boundary factor c / length-scale 可识别区间 + 诊断准则。

### 2. Bayes factor 的计算与可靠性(贡献点所在)
- **Gronau et al. (2017)**, *A tutorial on bridge sampling*, J. Mathematical Psychology(对应 `bridgesampling` R 包)。
- Savage-Dickey:Verdinelli & Wasserman (1995);Wagenmakers et al. (2010, Cognitive Psychology) 教程。
- 待搜:"marginal likelihood sensitivity basis function / rank approximation" —— 确认有没有人碰过(很可能是空白)。

### 3. 方法论邻居(linear vs nonlinear 检验、变量选择)
- **"Bayesian Testing of Linear Versus Nonlinear Effects Using Gaussian Process Priors"**, The American Statistician (2022) — 几乎是 project 1 的方法论母体,用的是 **exact GP**,正好坐实"需要近似版"的动机缺口。**Related work 必引,需确认确切作者 / 年份。**
- Piironen & Vehtari, projection predictive(`projpred`)— 和 BF 路线竞争 / 互补的主流做法,需对比讨论。参见 `experiments/260113_brms_projpred_gp_selection.R`。
- Bürkner, `brms` GP 文档 / 相关论文。

### 4. Kernel 选择与 spectral density
- 不必单独找新论文;Solin & Särkkä 和 Riutort-Mayol 两篇里关于 spectral density 衰减的部分即可支撑 SE-脆弱 / Matérn-稳健的论证。

## 建议读法

1. 先读 Riutort-Mayol (2023),把近似的 knobs 搞清楚。
2. 回头看 Solin & Särkkä 补理论。
3. 带着"这些 knob 怎么影响 BF"的问题去读 bridge sampling 那几篇。
4. 重点记:有没有人量化过 rank / basis 数对 marginal likelihood 的影响 —— 决定贡献是"新发现"还是"已知问题的应用"。

## 链接

- Solin & Särkkä (2020): https://arxiv.org/abs/1401.5508
- Riutort-Mayol et al. (2023) arXiv: https://arxiv.org/abs/2004.11408
- Riutort-Mayol et al. Springer: https://link.springer.com/article/10.1007/s11222-022-10167-2
- A tutorial on bridge sampling: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5699790/
- Bayesian Testing of Linear vs Nonlinear Effects (GP priors): https://www.tandfonline.com/doi/full/10.1080/00031305.2022.2028675
