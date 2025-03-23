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
