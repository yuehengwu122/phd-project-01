# Project

This repository contains code and notes for a PhD paper on additive Gaussian process regression, with Bayes factors used to assess linear and nonlinear effects.

The latest paper draft is `draft.pdf`.

## Repo Structure

This project uses three GitHub repositories:

- **phd-project-01** (this repo, private): full workspace — experiments, obsidian notes, data, results, code, models
- **phd-code-01** (lean, shareable): only `src/` and `models/*.stan`; synced from this repo via `scripts/sync_core.sh`
- **phd-paper-01** (paper): LaTeX source synced with Overleaf; lives at `phd-paper-01/` inside this repo but is a separate git repo (gitignored here)

To sync core code to phd-code-01 (first clone it to `~/Documents/phd-code-01` once):
```bash
bash scripts/sync_core.sh
```

## Core Entry Points

- Latest paper draft: `draft.pdf`
- Project overview and high-level progress map: `obsidian/00 Project Home/🏠 Project Overview.md`
- Setup entry point: `src/_setup.R`
- Main fitting code: `src/model_fitting.R`
- Bayes factor code: `src/bayes_factor.R`
- Main Stan models: `models/`
  - `models/fit_full_model.stan`: main model
  - `models/fit_linear_model.stan`: model without GP term, used for bridge sampling comparison
  - `models/fit_full_model_select.stan`: selection version used to test linear/nonlinear effects of specific predictors by bridge sampling
  - `models/fit_interaction_model.stan` and `models/fit_interaction_only_model.stan`: only for application 2
  - `models/simulation.stan`: simulation data generation
  - `models/fit_full_model_imoment.stan`: inverse-moment prior variant under exploration

## Current Repo State

- The repository does not yet cleanly separate stable interfaces from exploratory scripts.
- `src/` should be treated as the main reusable code layer.
- `analysis/` contains analysis scripts and is not yet fully standardized.
- `260511_exploratory/` is a temporary folder previously prepared for supervisor communication. Do not modify, reorganize, or use it as the basis for refactoring unless the user explicitly asks.
- `analysis/application/` should not be reorganized or refactored unless the user explicitly asks.

## Collaboration Rules

- Do not rewrite or restructure personal planning notes, task-management notes, or reflective notes unless the user explicitly asks for that exact change.
- Before changing notes or project-organization files, first read the relevant source files and summarize your understanding.
- If project context is still insufficient, ask instead of inventing structure.
- Prefer discussing and clarifying the repo structure before making broad organizational edits.
- Preserve the user's wording and intent in notes unless the user asks for rewriting.

## Working Style for This Repo

- When discussing the paper or code, ground claims in the repository files rather than generic statistical templates.
- Prefer small, local, reversible edits over broad rewrites.
- When proposing structural improvements, separate:
  - facts about the current repo
  - suggestions for future cleanup
- Treat model-comparison logic, prior specification, and bridge sampling as active research decisions, not settled implementation details.

## Default First Read

If context is needed, read these files first before making substantive suggestions:

1. `AGENTS.md`
2. `README.md`
3. `src/_setup.R`
4. `src/model_fitting.R`
5. `src/bayes_factor.R`
6. the relevant file under `models/`
7. the relevant script under `analysis/`

## What To Avoid

- Do not assume the current obstacle is only technical; some tasks may be waiting on supervisor input or long-running computation.
- Do not treat missing structure in notes as permission to replace them with a generic productivity system.
- Do not move or rename directories as a cleanup step unless the user explicitly requests it.
