# Tau-PET Heterogeneity Links Structural, Functional, and Spatially Dependent Pathological Changes to Discordant Clinical Perception and Divergent Anti-Amyloid Therapy Responses


This repository contains three integrated analysis pipelines for
discovering and characterizing tau-PET subtypes in Alzheimer's disease.

The framework includes:

1.  **Subtype Discovery**\
    Self-supervised deep learning to identify spatial tau-PET subtypes.

2.  **Brain Connectivity (BC) Modeling**\
    Subtype-specific functional connectivity characterization using
    Bayesian modeling.

3.  **Treatment Response Analysis**\
    Evaluation of differential solanezumab treatment effects across
    discovered subtypes in the A4 Study cohort.



## Methodological Overview

<p align="center">
  <img src="docs/Arch_5.png" alt="Analysis Pipeline" width="800">
</p>

**Figure 1.** Multi-stage framework for tau PET subtype discovery and validation:\
(A) Self-supervised deep learning framework for tau pattern identification. PET-tau input data are encoded into 10-dimensional latent representations and clustered using Gaussian Mixture Model.,\
(B) Multidimensional subtype characterization pipeline. Identified subtypes undergo comprehensive validation through clinical, genetic, cognitive, and functional connectivity profiling.


------------------------------------------------------------------------

## Repository Structure

    TauHeterogeneity/
    ├── README.md
    ├── docs/
    ├── src/
    │   ├── 01_subtype_discovery/
    │   ├── 02_bc_model/
    │   └── 03_treatment_response_analysis/




## Data Access

-   **ADNI:** https://adni.loni.usc.edu\
-   **A4 Study:** https://www.a4studydata.org

Data access requires approved data use agreements.


## System Requirements

### Python

-   Python 3.11.4
-   PyTorch
-   CUDA 12.4 (GPU recommended)

### R

-   R ≥ 4.0
-   Required packages: coda, MASS, pROC, corrplot, ggplot2

### License
This repository is licensed under the MIT License.

Note: ADNI and A4 study datasets are not distributed with this repository.
Users must obtain data access approval from the respective data providers.


