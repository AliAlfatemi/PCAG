# PCAG: Multi-omics Patient Subgrouping

MATLAB implementation of patient subgrouping through PCA-based representations,
patient-similarity graphs, and subspace integration on a Grassmann manifold.
The workflow integrates gene expression, DNA methylation, and microRNA views.

**Paper:** [Patient subgrouping with distinct survival rates via integration of multiomics data on a Grassmann manifold](https://doi.org/10.1186/s12911-022-01938-y)

**Authors:** Ali Alfatemi, Hong Peng, Wentao Rong, Bin Zhang, Hongmin Cai

**Publication:** BMC Medical Informatics and Decision Making, volume 22, article 190 (2022).

## Workflow

1. Align patients across the three omics matrices.
2. Standardize each view and construct its PCA-based representation.
3. Build a patient-similarity graph for each view.
4. Integrate the graph subspaces and perform spectral clustering.

## Requirements and input

Use MATLAB with the Statistics and Machine Learning Toolbox for `pca`, `zscore`,
and `kmeans`. The entry script expects a project-specific `BRCA` loader returning
`BRCAGE`, `BRCAMETH`, `BRCAMIR`, and `SUR`. Supply that loader and the aligned
cohort matrices before running the experiment.

```bash
git clone https://github.com/AliAlfatemi/PCAG.git
```

From the repository directory in MATLAB:

```matlab
addpath(pwd)
GrassMann_Fold_Total_Cleaning
```

The script configures the representation dimension, graph-neighbor count,
heat-kernel parameter, number of clusters, and subspace regularization.
Adapt these settings to the cohort and experiment protocol.

## Code reference

| File | Purpose |
|---|---|
| [GrassMann_Fold_Total_Cleaning.m](GrassMann_Fold_Total_Cleaning.m) | Experiment entry point |
| [pca_whit.m](pca_whit.m) | Standardization and PCA-based reconstruction |
| [constructW.m](constructW.m) | Similarity-graph construction |
| [sc_ml.m](sc_ml.m) | Multilayer spectral clustering |
| [displayClusters.m](displayClusters.m) | Cluster visualization |

## Citation and licensing

See [CITATION.cff](CITATION.cff) for citation metadata and [LICENSE](LICENSE)
for the repository license. Retain the third-party notices embedded in the
MATLAB source files.
