# Reproducibility Data

The reproducibility data accompanying this manuscript are available in the Zenodo archive:

**https://doi.org/10.5281/zenodo.23428343**

Download `Paper2_Data.zip` and extract it while preserving the original directory structure. The archive also contains `Paper2_Reproducibility.zip`, a snapshot of the corresponding MATLAB code.

The archive contains the analysis-ready inputs required to execute the reproducibility workflow, together with reference outputs from the original manuscript analyses and the final computational validation of the released code.

---

# Directory structure

```text
Paper2_Data/
├── Inputs/
│   ├── Data.mat
│   ├── PopSplits.mat
│   └── X_REST1.mat
├── Original/
└── Validation/
```

The data directory does not need to be located inside the `Paper2_Reproducibility` repository. Its location is specified through `DataDir` at the beginning of `RunPaper2.m`.

---

# Analysis-ready inputs

## `Data.mat`

`Data.mat` contains a MATLAB structure named `Data`, holding analysis-ready structural-connectivity, graph-spectral, demographic/relationship, behavioral, and head-motion information for 875 participants.

```matlab
load(fullfile(DataDir,'Inputs','Data.mat'),'Data');
```

The fields visible in the released structure are:

| Field | Stored size/type | Description |
|---|---|---|
| `Parents` | `875 × 2 double` | Parent/family identifiers used in family-related information |
| `SubjectIndices` | `1 × 875 double` | Participant identifiers |
| `RegionDistances` | `216 × 216 double` | Pairwise regional distances for the analyzed cortical regions |
| `HemiID` | `216 × 1 double` | Hemisphere assignments for analyzed cortical regions |
| `FactorScores` | `1 × 1 struct` | Behavioral factor scores |
| `Adj` | `1 × 875 cell` | Subject-specific structural adjacency matrices |
| `DegMat` | `1 × 875 cell` | Subject-specific graph degree matrices |
| `U` | `1 × 875 cell` | Subject-specific connectome harmonic bases |
| `Lambda` | `1 × 875 cell` | Corresponding graph Laplacian eigenvalues |
| `NumSubjects` | scalar, `875` | Number of participants |
| `NumRegions` | scalar, `274` | Number of parcels recorded in the input structure |
| `FamilyKernel` | `875 × 875 double` | Family-relatedness information |
| `FD_TimeCourses` | `1 × 9 cell` | Framewise-displacement time courses |
| `MeanFD` | `875 × 9 double` | Mean framewise-displacement summaries |
| `MeanFD_NonScrubbedOnly` | `875 × 9 double` | Mean framewise displacement for non-scrubbed frames |
| `PercentageScrubbed` | `875 × 9 double` | Percentage of frames flagged for scrubbing |

**Parcel-count distinction.** The input structure records 274 parcels, while the analyses described in the manuscript use 216 cortical regions. The `216 × 216` regional-distance matrix and `216 × 1` hemisphere identifiers reflect the cortical analysis space. The exact parcel-selection operation is defined by the analysis scripts; the stored `NumRegions` value should not be interpreted as the number of regions used in every downstream analysis.

## `PopSplits.mat`

`PopSplits.mat` contains a MATLAB structure named `PopSplits` that preserves the population-sampling scheme used in the study.

```matlab
load(fullfile(DataDir,'Inputs','PopSplits.mat'),'PopSplits');
```

Its top-level fields are:

| Field | Value/size | Description |
|---|---|---|
| `Outer` | `1 × 10 struct` | Ten outer retained/left-out partitions |
| `SampleSizes` | `1 × 22 double` | Population sample sizes from 5 to 700 |
| `NumOuterSplits` | `10` | Number of outer partitions |
| `NumMonteCarlo` | `20` | Monte Carlo realizations per outer partition |
| `LeftOutFraction` | `0.2` | Target fraction of participants held out |

Each `PopSplits.Outer(k)` has fields `MonteCarlo`, `Train`, `Test`, `NumTrain`, and `NumTest`. For example, the first outer split contains 700 training participants and 175 test participants, as well as a `1 × 20` Monte Carlo structure array. Exact retained/left-out counts may vary across splits because families are kept together.

The archived split definitions allow readers to reproduce the same sampling trajectories and comparisons. Regenerating random partitions may yield different numerical results. The separate cognition-prediction cross-validation is described in the manuscript and should not be confused with these repeated population-sampling splits.

## `X_REST1.mat`

`X_REST1.mat` contains `X_REST1`, a `1 × 875` cell array of regional resting-state activity time courses. Each participant's cell contains a `1 × 2` cell array, with each entry holding a `1194 × 274` matrix of type `single` (time points × parcels).

The two matrices correspond to the two acquisitions/directions represented in this input structure. The manuscript's functional analyses focus on REST1-LR; readers should follow the pipeline's acquisition-selection and cortical-parcel-selection steps rather than treating all 274 input parcels or both matrices as the analyzed signal.

This file is loaded by stages requiring subject-level functional time courses. The archive starts from processed regional signals rather than raw HCP imaging acquisitions.

---

# Reference outputs

## `Original/`

`Original/` contains the saved outputs generated during the original manuscript analyses.

These files provide a record of the numerical results underlying the manuscript prior to the final cleanup and validation of the released reproducibility workflow.

## `Validation/`

`Validation/` contains the outputs generated during the final clean rerun of the released analysis workflow.

During this validation, all manuscript analyses were rerun and the numerical results, statistical analyses, and main and supplementary figures were reproduced. Minor visualization issues identified during this process were corrected without altering the underlying numerical results or scientific conclusions.

The `Validation/` outputs therefore provide the reference state against which a new execution of the released workflow can be compared.

Neither `Original/` nor `Validation/` is modified when running the reproducibility workflow.

---

# Installing the data

1. Download or clone the `Paper2_Reproducibility` repository.
2. Download the accompanying `Paper2_Data` archive from Zenodo and extract it without modifying the directory structure.
3. At the beginning of `RunPaper2.m`, specify:

```matlab
RepositoryDir = '/path/to/Paper2_Reproducibility';
DataDir       = '/path/to/Paper2_Data';
```

The code and data directories may be stored independently and do not need to share a common parent directory.

No files from `Paper2_Data` need to be copied into the code repository.

---

# Data provenance

The analysis-ready inputs distributed with this package are derived from the Human Connectome Project Young Adult dataset and correspond to the processed quantities used by the accompanying study.

The reproducibility workflow begins from these analysis-ready inputs. It does not reproduce the complete preprocessing of the original Human Connectome Project acquisitions.

Further details regarding participant selection, image preprocessing, connectome construction, and functional-data processing are provided in the accompanying manuscript and its supplementary information.

---

# Newly generated results

Executing `RunPaper2.m` creates new analysis outputs under:

```text
Paper2_Reproducibility/results/
```

rather than writing into the distributed data archive.

This separation ensures that the supplied inputs and reference outputs remain unchanged while allowing newly generated results to be inspected or compared against the validated reference results.
