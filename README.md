# Reproducibility Package
## When are population connectome harmonics adequate representations of individual brains?

This repository accompanies the manuscript

> **"When are population connectome harmonics adequate representations of individual brains?"**

and contains the MATLAB code required to reproduce the numerical analyses reported in the manuscript.

The study evaluates population-derived connectome harmonics along three complementary dimensions:

- **Robustness** — reproducibility across independently sampled population representations.
- **Representativeness** — correspondence between population harmonics and subject-specific harmonic representations.
- **Practical Utility** — consequences of replacing individual harmonics with population harmonics in downstream functional analyses.

The repository additionally contains the analyses linking these three properties and assessing their dependence on population-connectome construction.

The accompanying reproducibility dataset is distributed separately and contains the analysis-ready inputs, the population splits used in the study, the original manuscript analysis outputs, and the outputs obtained during the final computational validation of the repository.

---

# Repository at a glance

- **Language:** MATLAB
- **Main analysis script:** `workflows/Paper2/analyses/RunPaper2.m`
- **Data:** distributed separately with the accompanying reproducibility dataset
- **License:** MIT
- **Archived code and data:** [Zenodo DOI: 10.5281/zenodo.23428343](https://doi.org/10.5281/zenodo.23428343)

---

# Purpose

This repository combines the research version of the **Graph Signal Processing Toolbox (GSPTB)** with the complete analysis workflow developed for the accompanying study.

The analysis is organized into five main stages:

1. **Robustness**
2. **Representativeness**
3. **Practical Utility**
4. **Cross-pillar analyses**
5. **Connectome construction analyses**

Each stage is implemented as a separate analysis script under `workflows/Paper2/analyses/`, while `RunPaper2.m` provides the entry point for executing the complete workflow.

The reusable graph signal processing methods are implemented within GSPTB, whereas study-specific configuration, analysis and figure-generation code is maintained under `workflows/Paper2/`.

---

# Reproducibility and validation

Following finalization of the analysis code, the complete computational workflow was rerun from the analysis-ready inputs. All numerical results, statistical analyses, and main and supplementary figures reported in the manuscript were reproduced.

During this validation, minor visualization issues affecting the aggregation of realization-level heat maps were identified and corrected. These corrections affected the visualization of the corresponding quantities but did not alter the underlying numerical results, statistical analyses, or scientific conclusions.

The accompanying reproducibility dataset therefore distinguishes between:

- **Original outputs** — results generated during the original manuscript analyses.
- **Validation outputs** — results generated during the final clean validation of the released workflow.

New executions of the pipeline write their outputs to the local `results/` directory and do not modify either reference set.

---

# Repository structure

```text
Paper2_Reproducibility/
├── external/
├── src/
├── workflows/
│   └── Paper2/
│       ├── analyses/
│       │   ├── RunPaper2.m
│       │   ├── Run_Robustness.m
│       │   ├── Run_Representativeness.m
│       │   ├── Run_PracticalUtility_Part1.m
│       │   ├── Run_PracticalUtility_Part2.m
│       │   ├── Run_PracticalUtility_Part3.m
│       │   ├── Run_CrossPillarAnalysis.m
│       │   └── Run_ConstructionSchemeAnalysis.m
│       └── figures/
├── results/
├── README.md
├── DATA.md
├── VALIDATION.md
├── LICENSE
└── CITATION.cff
```

The `results/` directory stores newly generated outputs. The master script creates it if it does not already exist.

---

# Reproducibility data

The accompanying data package is organized separately from the code repository:

```text
Paper2_Data/
├── Inputs/
│   ├── Data.mat
│   ├── PopSplits.mat
│   └── X_REST1.mat
├── Original/
└── Validation/
```

`Inputs/` contains the analysis-ready quantities required to execute the workflow. `PopSplits.mat` contains the population sampling used for the final validated analysis, allowing the released workflow to reproduce the same sampling realization rather than generating a new one.

`Original/` contains the reference outputs from the original manuscript analyses, whereas `Validation/` contains the corresponding outputs obtained from the final clean rerun of the released workflow.

See `DATA.md` for a detailed description of the data files, provenance, access conditions, and installation instructions.

---

# Requirements

The reproducibility workflow requires MATLAB together with the toolboxes used by the analyses.

The released code was validated using:

- **MATLAB:** R2025b
- **Operating system:** macOS
- **Execution:** serial (no parallel computation used)

Other MATLAB toolbox requirements have not yet been independently verified. The presence of third-party functions in `external/` does not by itself establish which MATLAB toolboxes are required.

Third-party MATLAB functions required by the workflow are included under `external/` and retain their original licenses.

The complete workflow can be computationally demanding, particularly for analyses involving repeated population sampling and subject-specific harmonic comparisons. Runtime and memory requirements therefore depend on the analyses being executed.

---

# Installation

## 1. Download the code

Clone or download the `Paper2_Reproducibility` repository.

## 2. Download the reproducibility data

Download the accompanying `Paper2_Data` package and retain its directory structure.

The code and data directories do not need to be located inside one another.

## 3. Configure local paths

Open `workflows/Paper2/analyses/RunPaper2.m` and specify:

```matlab
RepositoryDir = '/path/to/Paper2_Reproducibility';
DataDir       = '/path/to/Paper2_Data';
```

No other machine-specific paths should be required.

---

# Running the complete workflow

The complete analysis can be executed from `RunPaper2.m` after setting `RepositoryDir` and `DataDir` at the beginning of the script.

The workflow sequentially executes Robustness, Representativeness, Practical Utility, cross-pillar analyses, and connectome construction analyses.

The released `PopSplits.mat` file is loaded by default so that the analyses use the same family-aware population sampling as the validated workflow. The cognition-prediction procedure uses a separate family-aware cross-validation scheme, as described in the manuscript.

Newly generated outputs are written to `Paper2_Reproducibility/results/`, which is created automatically if it does not already exist. The reference results contained in `Paper2_Data/Original/` and `Paper2_Data/Validation/` are never modified.

---

# Running individual analyses

| Analysis | Script |
|---|---|
| Robustness | `Run_Robustness.m` |
| Representativeness | `Run_Representativeness.m` |
| Practical Utility, Part 1 | `Run_PracticalUtility_Part1.m` |
| Practical Utility, Part 2 | `Run_PracticalUtility_Part2.m` |
| Practical Utility, Part 3 | `Run_PracticalUtility_Part3.m` |
| Cross-pillar analyses | `Run_CrossPillarAnalysis.m` |
| Connectome construction analyses | `Run_ConstructionSchemeAnalysis.m` |

Individual analysis scripts assume that the common initialization section of `RunPaper2.m` has first been executed, providing the configuration, analysis-ready data, population splits, and relevant directory paths.

---

# Code and data availability

The complete reproducibility package is archived on Zenodo:

**https://doi.org/10.5281/zenodo.23428343**

The archive contains two complementary files:

- `Paper2_Reproducibility.zip` — MATLAB source code, analysis workflows, dependencies, and documentation.
- `Paper2_Data.zip` — Analysis-ready inputs, predefined population splits, original manuscript results, and computational validation outputs.

The source code is also available on [GitHub](https://github.com/TiBiUan/Paper2_Reproducibility).

See [DATA.md](DATA.md) for details of the data organization and installation instructions.
