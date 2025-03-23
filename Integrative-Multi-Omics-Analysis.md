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

The Synapse command line client is implemented in Python and comes with
the Synapse Python package. To install the Synapse command line client,
make sure that you have Python and pip installed. For more information,
see the [Python](https://www.python.org/downloads/) and
[pip](https://pip.pypa.io/en/stable/installing/) installation
instructions.

``` bash
pip install synapseclient
synapse login -u SYNAPSEUSER -p SYNAPSE_AUTH_TOKEN
synapse -h

# example:download single file
synapse get syn34114003 --version 2

# download multiple files together
synapse get-download-list
```

``` r
# Proteomics
synapse get syn33605372 --version 2

df.protein <- read.csv("data/TMT_normAbundances_Batch_regressed_Jax.IU.Pitt_LOAD2.csv") %>% dplyr::select(-"specimenID") %>% column_to_rownames(.,var="individualID")

# Brain metabolomics
synapse get syn34114003

df.brainmet <- read.csv("data/Q500-brain_5833_Jax.IU.Pitt_LOAD2.csv")[,c(1:2,6:248)] %>% select(-"specimenID") %>% column_to_rownames(.,var="individualID")

# Plasma metabolomics 
synapse get syn34114002

df.plasmamet <- read.csv("data/Q500-plasma_5833_Jax.IU.Pitt_LOAD2.csv")[,c(1:2,6:459)] 

# Transcriptomics 
synapse get syn26195571

df.rna <- read.delim2("data/rnaseq_rsem.merged.gene_counts_Jax.IU.Pitt_LOAD2.tsv",check.names = F) %>% dplyr::select(-"transcript_id(s)") %>% column_to_rownames("gene_id")
```

We need to re-process the data to prepare it for input into MOFA. After
re-processing, you can store data in a list format like below and can be
inputted to MOFA.

``` r
AD_data <- list("Brain.Trans"= as.matrix(df.rna.norm),"PlasmaMeta"= as.matrix(df.plasmamet),"BrainMeta"= as.matrix(df.brainmet),"Proteins"= as.matrix(df.protein))
saveRDS(AD_data,file="data/AD_OmicsData.rds")
```

We’ll skip this step for now, but you can complete it in your free time.
We’ll be using the processed data and metadata for integrative analysis.

### Load data

``` r
AD_data <- readRDS("data/AD_OmicsData.rds")
lapply(AD_data, dim) 
```

    ## $Brain.Trans
    ## [1] 19878   106
    ## 
    ## $PlasmaMeta
    ## [1] 454 105
    ## 
    ## $BrainMeta
    ## [1] 243 106
    ## 
    ## $Proteins
    ## [1] 8602  106

Let’s check the data

``` r
AD_data$Brain.Trans[1:5,1:8]
```

    ##          50130    50316    50323    50389    50394    50395    51003    51007
    ## Gnai3 3.382053 3.407262 3.399493 3.414759 3.423258 3.410657 3.425562 3.416852
    ## Cdc45 2.840115 2.862579 2.850313 2.861727 2.846612 2.870444 2.878602 2.899075
    ## H19   2.786492 2.729832 2.746410 2.773358 2.809781 2.700492 2.702488 2.738903
    ## Scml2 2.891104 2.856333 2.815852 2.894488 2.885750 2.883277 2.875253 2.814438
    ## Apoh  2.729897 2.733669 2.703876 2.685259 2.729252 2.709516 2.699244 2.747691

``` r
AD_data$Proteins[1:5,1:8]
```

    ##                      50130       50316        50323        50389        50394
    ## Dync1h1.Q9JHU4  0.03593564  0.06332429 -0.028716182  0.069717972  0.000371381
    ## Sptan1.B9EKJ1  -0.01621585 -0.04678289  0.003484264 -0.025884450  0.042879525
    ## Sptan1.E9Q447   0.04036904 -0.06265565  0.114494945  0.152829910 -0.065677761
    ## Sptan1.B7ZWK3  -0.02807768 -0.08684510 -0.028227685 -0.003363745  0.030999364
    ## Ank2.Q8C8R3.2  -0.02666336 -0.03879798  0.009258663 -0.033007445  0.044784192
    ##                       50395       51003        51007
    ## Dync1h1.Q9JHU4 -0.006129719  0.03589722  0.001746378
    ## Sptan1.B9EKJ1   0.017133969 -0.06152933 -0.020950508
    ## Sptan1.E9Q447   0.212985616  0.16678477 -0.004633678
    ## Sptan1.B7ZWK3   0.027687324 -0.10505185 -0.037864239
    ## Ank2.Q8C8R3.2   0.035806980 -0.04287130  0.059415673

``` r
AD_data$PlasmaMeta[1:5,1:8]
```

    ##                      50130      50316      50323      50389      50394
    ## C0             4.773192040  4.7227381  5.1494215  4.9496897  5.0672718
    ## C2             3.782306844  4.0666836  4.2304026  4.0752836  3.8881590
    ## C3            -0.693605707 -0.3551559 -0.7073895 -0.5356703 -0.9493844
    ## C3.DC..C4.OH. -3.075408776 -2.3855597 -2.5899819 -2.3312477 -2.4073752
    ## C4             0.001051582  0.4977249  0.5294306  0.7924650  1.1846088
    ##                    50395      51003      51007
    ## C0             5.1212723  4.1466576  4.6358948
    ## C2             4.3667201  2.9239377  3.3212564
    ## C3             0.2624915 -1.2431486 -1.5094841
    ## C3.DC..C4.OH. -1.7204870 -3.2793603 -3.3123107
    ## C4             1.5563730  0.7464116 -0.4079952

``` r
AD_data$BrainMeta[1:5,1:8]
```

    ##                   50130     50316     50323     50389     50394     50395
    ## C0             3.628861  3.491175  3.335184  3.777141  3.492472  3.434089
    ## C2             2.019508  2.462707  1.593381  1.962925  2.414993  2.282543
    ## C3            -3.767542 -2.194707 -4.131880 -4.894299 -3.344102 -2.280190
    ## C3.DC..C4.OH. -3.898077 -3.511717 -4.602416 -3.952525 -4.250866 -3.678960
    ## C4            -3.859741 -2.535561 -4.089066 -4.503175 -3.315805 -2.932323
    ##                   51003     51007
    ## C0             2.908451  3.137472
    ## C2             1.742815  1.927878
    ## C3            -3.079559 -3.237546
    ## C3.DC..C4.OH. -3.874287 -4.102768
    ## C4            -2.811956 -3.256912

### Load metadata

``` r
AD_covariates <- readRDS("data/metadata.rds")
dim(AD_covariates)
```

    ## [1] 106   4

Let’s check the metadata

``` r
head(AD_covariates)
```

    ##          Sex Age Genotype Diet
    ## 50130 Female   4    LOAD1   CD
    ## 50316 Female  18    LOAD1   CD
    ## 50323   Male  12    LOAD1   CD
    ## 50389 Female  18    LOAD2   CD
    ## 50394   Male  18    LOAD2   CD
    ## 50395   Male  18    LOAD1   CD

``` r
dplyr::count(AD_covariates, Sex, Genotype, Age,Diet)  %>% gt() 
```

<div id="frvjdrwqfl" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#frvjdrwqfl table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#frvjdrwqfl thead, #frvjdrwqfl tbody, #frvjdrwqfl tfoot, #frvjdrwqfl tr, #frvjdrwqfl td, #frvjdrwqfl th {
  border-style: none;
}
&#10;#frvjdrwqfl p {
  margin: 0;
  padding: 0;
}
&#10;#frvjdrwqfl .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#frvjdrwqfl .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#frvjdrwqfl .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#frvjdrwqfl .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#frvjdrwqfl .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#frvjdrwqfl .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#frvjdrwqfl .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#frvjdrwqfl .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#frvjdrwqfl .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#frvjdrwqfl .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#frvjdrwqfl .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#frvjdrwqfl .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#frvjdrwqfl .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#frvjdrwqfl .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#frvjdrwqfl .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#frvjdrwqfl .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#frvjdrwqfl .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#frvjdrwqfl .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#frvjdrwqfl .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#frvjdrwqfl .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#frvjdrwqfl .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#frvjdrwqfl .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#frvjdrwqfl .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#frvjdrwqfl .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#frvjdrwqfl .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#frvjdrwqfl .gt_left {
  text-align: left;
}
&#10;#frvjdrwqfl .gt_center {
  text-align: center;
}
&#10;#frvjdrwqfl .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#frvjdrwqfl .gt_font_normal {
  font-weight: normal;
}
&#10;#frvjdrwqfl .gt_font_bold {
  font-weight: bold;
}
&#10;#frvjdrwqfl .gt_font_italic {
  font-style: italic;
}
&#10;#frvjdrwqfl .gt_super {
  font-size: 65%;
}
&#10;#frvjdrwqfl .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#frvjdrwqfl .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#frvjdrwqfl .gt_indent_1 {
  text-indent: 5px;
}
&#10;#frvjdrwqfl .gt_indent_2 {
  text-indent: 10px;
}
&#10;#frvjdrwqfl .gt_indent_3 {
  text-indent: 15px;
}
&#10;#frvjdrwqfl .gt_indent_4 {
  text-indent: 20px;
}
&#10;#frvjdrwqfl .gt_indent_5 {
  text-indent: 25px;
}
&#10;#frvjdrwqfl .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#frvjdrwqfl div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Sex">Sex</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Genotype">Genotype</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Age">Age</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="Diet">Diet</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="n">n</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">B6</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">4</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">12</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">HFD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">4</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">12</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Female</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">HFD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">B6</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">4</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">4</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">12</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD1</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">HFD</td>
<td headers="n" class="gt_row gt_right">5</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">4</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">7</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">12</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">CD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
    <tr><td headers="Sex" class="gt_row gt_center">Male</td>
<td headers="Genotype" class="gt_row gt_center">LOAD2</td>
<td headers="Age" class="gt_row gt_center">18</td>
<td headers="Diet" class="gt_row gt_center">HFD</td>
<td headers="n" class="gt_row gt_right">6</td></tr>
  </tbody>
  &#10;  
</table>
</div>

This table tells you all age group, sex, genotype, diet and how many
samples are in each group.

### Create MOFA object

``` r
# Create MultiAssayExperiment object 
mae_AD <- MultiAssayExperiment(
  experiments = AD_data,
  colData = AD_covariates)

# Create MOFA object ##
model <- create_mofa(mae_AD)
model
```

### Build and Train Object

This step can take upto 15-20 minutes. We’ll skip this step for now, but
you can complete it later.

``` r
## Define options ##
##data options
data_opts <- get_default_data_options(model)
data_opts

#Define model options
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 10
model_opts

#Define training options
train_opts <- get_default_training_options(model)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 1111
train_opts

## Prepare MOFA object ##
model <- prepare_mofa(model,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

## Train the model ##
model <- run_mofa(model,use_basilisk = TRUE)

## Save the model ##
saveRDS(model,"data/MOFA2_Object_10factor.rds")
#micromamba activate /Users/pandera/Library/Caches/org.R-project.R/R/basilisk/1.18.0/0
```

### Load trained model

We will load the trained model and explore the MOFA object for
downstream analyses.

``` r
MOFAobject <- readRDS("data/MOFA2_Object_15factor.rds")
MOFAobject
```

    ## Trained MOFA with the following characteristics: 
    ##  Number of views: 4 
    ##  Views names: Brain.Trans PlasmaMeta BrainMeta Proteins 
    ##  Number of features (per view): 19878 454 243 8602 
    ##  Number of groups: 1 
    ##  Groups names: group1 
    ##  Number of samples (per group): 106 
    ##  Number of factors: 15

### Overview of the data

The function `plot_data_overview` can be used to obtain an overview of
the input data. It shows how many views (rows) and how many groups
(columns) exist, what are their corresponding dimensionalities and how
many missing information they have (grey bars).

``` r
plot_data_overview(MOFAobject)
```

<img src="Integrative-Multi-Omics-Analysis_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

### Add metadata to the model

The metadata is stored as a data.frame object in
`model@samples_metadata`, and it requires at least the column sample.
The number of rows must match the total number of samples in the model.

``` r
head(AD_covariates)
```

    ##          Sex Age Genotype Diet
    ## 50130 Female   4    LOAD1   CD
    ## 50316 Female  18    LOAD1   CD
    ## 50323   Male  12    LOAD1   CD
    ## 50389 Female  18    LOAD2   CD
    ## 50394   Male  18    LOAD2   CD
    ## 50395   Male  18    LOAD1   CD

``` r
AD_covariates$sample <- rownames(AD_covariates)
samples_metadata(MOFAobject) <- AD_covariates
head(MOFAobject@samples_metadata, n=5)
```

    ##          Sex Age Genotype Diet sample  group
    ## 50130 Female   4    LOAD1   CD  50130 group1
    ## 50316 Female  18    LOAD1   CD  50316 group1
    ## 50323   Male  12    LOAD1   CD  50323 group1
    ## 50389 Female  18    LOAD2   CD  50389 group1
    ## 50394   Male  18    LOAD2   CD  50394 group1

### Variance decomposition

The first step in the MOFA analysis is to quantify the amount of
variance explained (𝑅^2) by each factor in each data modality. The
variance explained estimates are stored in the hdf5 file and loaded in
`MOFAobject@cache[["variance_explained"]]`:

##### Total variance explained per view and group

``` r
MOFAobject@cache$variance_explained$r2_total[[1]]
```

    ## Brain.Trans  PlasmaMeta   BrainMeta    Proteins 
    ##    42.76213    66.47333    32.57267    43.03760

Let’s plot variance explained estimates:

``` r
plot_variance_explained(MOFAobject , x="view", y="factor")
```

<img src="Integrative-Multi-Omics-Analysis_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

#### Renaming dimensions

The user can rename the dimensions of the model as follows:

``` r
views_names(MOFAobject)
```

    ## [1] "Brain.Trans" "PlasmaMeta"  "BrainMeta"   "Proteins"

``` r
views_names(MOFAobject) <- c("Transcriptomics", "Plasma Metabolomics","Brain Metabolomics","Proteomics")
views_names(MOFAobject)
```

    ## [1] "Transcriptomics"     "Plasma Metabolomics" "Brain Metabolomics" 
    ## [4] "Proteomics"

``` r
factors_names(MOFAobject)
```

    ##  [1] "Factor1"  "Factor2"  "Factor3"  "Factor4"  "Factor5"  "Factor6" 
    ##  [7] "Factor7"  "Factor8"  "Factor9"  "Factor10" "Factor11" "Factor12"
    ## [13] "Factor13" "Factor14" "Factor15"

``` r
factors_names(MOFAobject) <- paste0("LF",c(1:15))
factors_names(MOFAobject)
```

    ##  [1] "LF1"  "LF2"  "LF3"  "LF4"  "LF5"  "LF6"  "LF7"  "LF8"  "LF9"  "LF10"
    ## [11] "LF11" "LF12" "LF13" "LF14" "LF15"

``` r
plot_variance_explained(MOFAobject)
```

<img src="Integrative-Multi-Omics-Analysis_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

##### Variance explained for every factor in per view and group

``` r
MOFAobject@cache$variance_explained$r2_per_factor[[1]]
```

    ##      Transcriptomics Plasma Metabolomics Brain Metabolomics   Proteomics
    ## LF1     1.107528e-04        45.023971719        0.001599441 1.056614e-04
    ## LF2     2.388946e-04         0.103648409        2.401866961 2.124522e+01
    ## LF3     1.620314e+01         0.002631586        2.610556945 5.607056e-03
    ## LF4     1.993100e-01         7.611620670        8.371110022 2.543883e+00
    ## LF5     3.985905e+00         1.319170924        7.694674229 3.404121e+00
    ## LF6     8.252737e-05         0.001617342        0.223916448 1.102318e+01
    ## LF7     7.127392e-01         2.802452460        5.609232985 1.299765e+00
    ## LF8     1.943470e+00         4.649598197        2.356880180 1.943161e-01
    ## LF9     4.464782e-01         6.106694043        0.768684349 1.788879e+00
    ## LF10    8.252979e+00         0.026071594        0.772557042 1.652078e-02
    ## LF11    4.867525e+00         0.196413618        0.001227749 2.043101e-04
    ## LF12    3.136182e+00         0.192051373        0.090483299 5.067614e-03
    ## LF13    2.615062e+00         0.197868079        0.557388585 1.853962e-02
    ## LF14    8.827759e-04         0.276999158        0.260847252 2.569058e+00
    ## LF15    1.454285e+00         0.012322165        0.001368526 8.209902e-05

``` r
#Total variance explained per view
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
```

<img src="Integrative-Multi-Omics-Analysis_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

##### Correlation between factors

We can also visualize correlation between factors

``` r
#Correlation between factors
plot_factor_cor(MOFAobject)
```

<img src="Integrative-Multi-Omics-Analysis_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />
