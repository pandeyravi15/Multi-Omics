Integrative Multi-Omics Analysis
================
2025-03-23

Author: Ravi Pandey, Jackson Laboratory

# Introduction

Alzheimer’s disease is a complex, multifactorial pathology with high
heterogeneity in biological alternations.One of the biggest challenges
in AD is to identify pathways and markers of disease progression, which
can be easily accessible, in asymptomatic at-risk individuals.
Multi-omics data could provide complementary information, which might be
helpful to reveal the underlying biology of the associations. Moreover,
integrating the information from single-omics studies provides an
opportunity for a thorough exploration of endophenotype networks and
biological interactions related to disease.

<img src="figures/intro_fig2.png" width="1000px" align="center" alt="Overview" >

## Challenges in multi-omics analysis

- Heterogeneity, Sparsity and outliers
- Omics datasets can differ vastly in size (number of features)
- More features than data (p \>\> n)
- Class imbalance and overfitting
- Computation and storage cost
- Additionally, biological datasets are complex, noisy, with potential
  errors due to measurement mistakes or unique biological deviations.

## Methods/Tools

<img src="figures/method_fig1.png" width="1000px" align="center" alt="tools" >

One can choose tools based on its ability to address biological question
of interests and approaches.

### The biological questions are broadly categorized into 3 different case studies:

- Disease subtyping and classification based on multi-omics profiles.

- Prediction of biomarkers for various applications including
  diagnostics and driver genes for diseases.

- Deriving insights into disease biology.

<img src="figures/method_fig3.jpg" width="1000px" align="center" alt="tools" >

In this lesson, we are going to integrate data from multiple omics
platforms (transcriptomics, proteomics, and metabolomics) in an unbiased
fashion, considering interaction between modalities using **multi-omics
factor analysis (MOFA)**.

## Overview of MOFA (multi-omics factor analysis)

[MOFA](https://biofam.github.io/MOFA2/) is a factor analysis model that
provides a general framework for the integration of multi-omic data sets
in an unsupervised fashion. Intuitively, MOFA can be viewed as a
versatile and statistically rigorous generalization of principal
component analysis to multi-omics data. Given several data matrices with
measurements of multiple -omics data types on the same or on overlapping
sets of samples, MOFA infers an interpretable low-dimensional
representation in terms of a few latent factors. These learnt factors
represent the driving sources of variation across data modalities, thus
facilitating the identification of cellular states or disease subgroups.
[MOFA2](https://biofam.github.io/MOFA2/NEWS.html)
<img src="figures/mofa_overview_fig1.png" width="1000px" align="center" alt="Overview of MOFA" >

## Multi-Omics Data

For this lesson, we are going to use multi-omics data from [LOAD2 mice
cohort](https://www.synapse.org/Synapse:syn51534997). LOAD2 mice cohort
consist of mouse models expressing humanized Abeta and two genetic risk
factors (APOE4 and Trem2\*R47H) at multiple ages for both sexes. The
data are being released as part of [MODEL-AD (Model Organism Development
& Evaluation for Late-Onset Alzheimer’s
Disease)](https://www.model-ad.org) consortium. We conducted brain RNA
sequencing, TMT-based brain proteomics, and targeted metabolomics in
both brain and plasma samples from over a hundred mice.

<img src="figures/intro_fig1.png" width="1000px" align="center" alt="workflow" >

**Let’s start …**

### Load Libraries

``` r
library("MOFA2")
library("MOFAdata")
library("AnnotationDbi")
library("MultiAssayExperiment")
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(cowplot)
library(reticulate)
library(gt)
#library(synapser)
```

### Synapse Download

You can download the data from Synapse data repository. API clients
provide a way to use Synapse programmatically. Installation instructions
are available at [Synapse API Documentation
Site](https://help.synapse.org/docs/Installing-Synapse-API-Clients.1985249668.html).
