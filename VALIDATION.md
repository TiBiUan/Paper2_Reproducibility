# Computational Validation

This document describes the validation of the MATLAB reproducibility workflow accompanying the manuscript **“When are population connectome harmonics adequate representations of individual brains?”**

## Scope and computational environment

Following final code cleanup and preparation of the reproducibility package, the complete Paper 2 workflow was rerun using **MATLAB R2025b** on macOS, **without parallel computation**.

The validation reused the archived `Inputs/PopSplits.mat`, preserving the family-aware outer partitions and Monte Carlo population samples used in the reference analyses. The workflow was executed with the archived analysis-ready inputs.

Validation covered:

- Robustness analyses, including realization-level trajectories, sample-size contrasts, power-law fits, and convergence estimates;
- Representativeness metrics and statistical analyses;
- Practical Utility analyses, including functional measures and Cognition prediction;
- Cross-pillar association analyses;
- Comparisons of population-connectome construction choices;
- Generation of all main and supplementary figures.

## Archive organization

The accompanying data archive distinguishes the following directories:

```text
Paper2_Data/
├── Inputs/       # Analysis-ready inputs and precomputed population splits
├── Original/     # Results from the original manuscript analysis runs
└── Validation/   # Results from the final validation of the released workflow
```

`Original/` preserves the historical outputs of the manuscript analyses. `Validation/` contains the results produced by rerunning the cleaned, released workflow and is therefore the appropriate reference when checking a fresh execution of that workflow.

New executions of `RunPaper2.m` write to the repository's local `results/` directory, leaving the archived `Original/` and `Validation/` outputs unchanged.

## Validation outcome

The full validation run completed successfully and reproduced the numerical results, statistical analyses, and scientific conclusions reported in the manuscript. The principal analysis stages and the main and supplementary figures were checked against the corresponding reference results.

The family-aware permutation analysis comparing population- and individual-level Cognition prediction used an explicit random seed in the finalized workflow, enabling deterministic repetition of its permutation draws under the validated environment. This should not be confused with permutation-based significance testing of each individual prediction model's R², which was not part of Paper 2.

**Reproducibility qualification:** The validation established agreement of the reported results and conclusions; this document does not claim bit-for-bit equality of every stored array or specify an unverified universal numerical tolerance. Minor floating-point differences can depend on software and computational environment.

## Visualization corrections identified during validation

Validation identified an aggregation issue affecting the visualization of realization-level data in selected heat maps. The affected plotting routines were corrected to show the intended aggregation across realizations.

These corrections changed the appearance of relevant heat maps in **Figures 2 and 6**, but did not change the underlying numerical results, statistical analyses, or scientific conclusions.

Minor plotting issues in **Supplementary Figure 1** were also corrected without changing the underlying analyses.

Accordingly, readers comparing historical `Original/` figure outputs with `Validation/` should use the finalized manuscript and validation figures as the references for the corrected visualizations. Differences in those particular plots do not, by themselves, indicate a discrepancy in the numerical analyses.

## Checking a new execution

1. Obtain the GitHub code release and the corresponding data archive, following `README.md` and `DATA.md`.
2. Configure the repository and data paths in `RunPaper2.m`.
3. Run the desired analysis stages using the archived inputs and population splits.
4. Inspect newly generated outputs in the local `results/` directory.
5. Compare the corresponding metrics, statistical summaries, and figures against `Paper2_Data/Validation/`, using the finalized manuscript to interpret expected plotting corrections.

Because the workflow includes computationally intensive repeated-sampling analyses, readers may inspect archived reference outputs without rerunning every stage. A fresh execution should be compared with the validation outputs for the same analysis and configuration.

## Interpretation

The purpose of the validation is to demonstrate that the cleaned, publicly released analysis workflow reproduces the findings of the manuscript. The `Original/` and `Validation/` directories are retained separately to provide a transparent record of the original computations and the final reproducibility check.
