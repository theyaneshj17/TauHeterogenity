# Bayesian Connectivity Modeling of fMRI Subtypes

This repository contains the analysis pipeline for subtype-specific
Bayesian Connectivity (BC) modeling of functional connectivity data.

The workflow consists of three main stages:

1.  Latent dimension selection (model selection via train/test split)
2.  Correlation-based model evaluation
3.  Final model fitting using the full dataset

------------------------------------------------------------------------

# Repository Structure

## Repository Structure

```
bc_model/
│
├── code/
│   └── BC model implementation and supporting utility functions
│
├── 01_ModelSelection_Split.R
├── 02_CorrelationEvaluation.R
├── 03_FinalModel_FullData.R
│
├── README.md
```


------------------------------------------------------------------------

# Data Description

**X.rds**\
A list of subject-level functional connectivity matrices.

**SubjectMapping.csv**\
Contains subject IDs and subtype cluster assignments derived from the prior subtype discovery procedure.

Each subject is assigned to one of K clusters (0, 1, ..., K−1), where K represents the optimal number of subtypes determined during clustering.

Subsequent connectivity modeling is conducted separately within each identified subtype.



------------------------------------------------------------------------

# Required R Packages

Install required packages before running the analysis:

```r
install.packages(c(
  "coda",
  "MASS",
  "pROC",
  "corrplot",
  "ggplot2",
  "reshape2",
  "grid",
  "gridExtra"
))
```
------------------------------------------------------------------------

# Analysis Pipeline

## Step 1 --- Latent Dimension Selection

Script: 01_ModelSelection_Split.R

This script: - Splits data into 80% training / 20% testing - Fits BC
models for multiple latent dimensions (e.g., K = 2--8) - Saves model
objects and held-out test data - Repeats across experiments if specified

Used to determine the optimal latent dimension (K).

------------------------------------------------------------------------

## Step 2 --- Correlation-Based Model Evaluation

Script: 02_CorrelationEvaluation.R

This script: - Loads saved BC model outputs - Computes correlation
between estimated connectivity (training data) and average held-out
connectivity (testing data) - Aggregates results across experiments -
Produces summary PDF and CSV outputs

The optimal latent dimension is selected based on predictive
consistency.

------------------------------------------------------------------------

## Step 3 --- Final Model Fitting (Full Dataset)

Script: 03_FinalModel_FullData.R

This script: - Uses the optimal latent dimension selected in Steps
1--2 - Fits BC models separately for each subtype cluster - Uses the
full dataset (no train/test split) - Generates cluster-specific
connectivity matrices and difference maps

These outputs correspond to the final results reported in the
manuscript.

------------------------------------------------------------------------
# Reproducibility Notes

-   All paths are relative to the project root directory.
-   Random seeds are set within scripts for reproducibility.
-   Latent dimensionality is selected via cross-validation before final
    model fitting.
-   Final connectivity estimates are derived from models fit on the full
    dataset.

------------------------------------------------------------------------

# Methodological Summary

1.  Subtypes were identified in a prior clustering analysis.
2.  Subjects were assigned to subtype clusters.
3.  BC models were fit independently within each cluster.
4.  Latent dimensionality (K) was selected using train/test validation.
5.  The final BC model was refit on the full dataset using the optimal
    K.

------------------------------------------------------------------------

# Output

Outputs include: - .rds model objects - Test data files - Correlation
summary tables - Connectivity heatmaps - Cluster difference maps


---

---

# Acknowledgment and Citation

The Bayesian Connectivity (BC) model implementation used in this repository is adapted from the framework proposed in:

Wang, S., Wang, Y., Xu, F. H., Shen, L., & Zhao, Y. (2025).  
*Establishing group-level brain structural connectivity incorporating anatomical knowledge under latent space modeling.*  
Medical Image Analysis, 99, 103309.  
https://doi.org/10.1016/j.media.2024.103309

We acknowledge the original authors for the theoretical development and foundational implementation of the model.  
For methodological details and the original framework, please refer to the cited publication.
