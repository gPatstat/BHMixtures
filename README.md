# BHMixtures

R code accompanying the paper:

> G. Patanè, S. Greven, A. Menafoglio. *Random mixtures in Bayes Hilbert spaces*. [arXiv:2609.03523](https://arxiv.org/abs/2609.03523)

The paper introduces a framework for defining and studying the identifiability of random mixtures of densities in the Bayes Hilbert space B²(I), and proposes a penalised maximum-likelihood estimation approach (via an MCEM algorithm) for unmixing such mixtures, with an application to hyperspectral unmixing (AVIRIS Indian Pines dataset).

## Repository contents

**Simulation of the mixtures**
- `BHM_simulator.R` — simulator of Bayes Hilbert mixtures (Beta-distributed vertices, ilr-normal proportions).
- `LM_simulator.R` — simulator of linear (probabilistic) mixtures, used as a model-misspecification scenario.
- `clr2density.R` — inverse clr transformation (clr coefficients → density).
- `second_derivative_fd.R` — second-derivative operator on an evenly spaced grid (smoothing penalty).
- `coef_function.R` — projection of densities onto the compositional B-spline basis (clr coefficient representation).

**Estimation (MCEM algorithm, Section 5 of the paper)**

The paper discusses two options for sampling the conditional distribution of the proportions p in the E-step (Section 5.3): importance sampling and Hamiltonian Monte Carlo (HMC). The repository provides both, but they are not used interchangeably in the actual study scripts:

- `Sampling_p.R` — importance-sampling routine for p | data.
- `update_functions.R` — closed-form M-step updates for H, Σ_ε, μ_p, Σ_p built on top of `Sampling_p.R` (importance sampling). Used only in the demonstrative script `examples_simulations.R`, not in the Study A/B/C scripts.
- `stan_intro.R` — Stan/`cmdstanr` model used to draw HMC samples of p in the E-step.
- `update_functions_new.R` — closed-form M-step updates for H, Σ_ε, μ_p, Σ_p built on top of the HMC sampler from `stan_intro.R` (via the helper `safe_stan_draws`, which also discards divergent transitions and non-finite draws). This is the version actually sourced by all `study*_BHM.R` / `study*_LM.R` scripts (Studies A, B, C).
- `without_pen.R` — a further variant of the update functions without the smoothing/trace penalisation, not called by the study scripts.
- `EM_updating.R` — main MCEM loop, with a stochastic-batch subsampling option (Algorithm 1 of the paper).

In short: **the simulations reported in the paper (Studies A, B, C) run on `update_functions_new.R` + `stan_intro.R`, i.e. the HMC-based E-step**; the importance-sampling version (`update_functions.R` + `Sampling_p.R`) is only illustrated in `examples_simulations.R`.

**Simulation studies (Section 6 and Appendices C, F of the paper)**
- `studyA_BHM.R`, `studyA_LM.R` — Study A: effect of increasing the residual variance σ².
- `studyB_BHM.R`, `studyB_LM.R` — Study B: effect of moving the mean proportion μ_p towards one vertex.
- `studyC_BHM.R`, `studyC_LM.R` — Study C (Appendix F): effect of the number of vertices m.
- `vertex_error.R` — L² estimation error between estimated and true vertices, with optimal permutation realignment.
- `examples_simulations.R` — illustrative examples of the simulator and estimation routines.

**Result aggregation and plots**
- `data_boxes.R`, `data_boxes_lm.R` — aggregation of results across simulation replicates.
- `simC_appendix_boxplot.R` — aggregation of the Study C results.
- `plot_boxes.R` — boxplots of the estimation error (Figures 3 and 12 of the paper).

## Requirements

R code. Packages used across the scripts: `FDboost`, `compositions`, `clue`, `ggplot2`; the HMC variant of the E-step additionally requires `cmdstanr` (and a CmdStan installation).

## Note

Several scripts contain absolute paths (`setwd`) referring to the authors' original development environment and need to be adapted before running. The AVIRIS Indian Pines case study (Section 7 of the paper) is not included in this repository as a separate script.

## Citation

If you use this code, please cite:

```
Patanè, G., Greven, S., Menafoglio, A. (2026). Random mixtures in Bayes Hilbert spaces.
arXiv:2609.03523. https://arxiv.org/abs/2609.03523
```
