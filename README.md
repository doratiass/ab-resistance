# Estimating Antibiotic Resistance Following Antibiotic Treatment in Outpatients

This repository contains the R code and data preparation scripts used for the analysis presented in the manuscript:

> **Estimating antibiotic resistance following antibiotic treatment in outpatients: a retrospective study**\
> *Chowers Michal, Atias Dor, Gottesman Batsheva, Low Marcelo, Obolski Uri*

The study evaluates the impact of different antibiotic treatments on the future risk of antibiotic resistance using advanced epidemiological methods. These methods include propensity score (PS) matching, risk difference (RD) estimation via g-methods (standardization), and inverse probability of censoring weighting (IPCW). The code is organized into modular functions and objects to enable efficient, reproducible, and repeated analyses.

------------------------------------------------------------------------

## Overview

This repository implements a full analysis pipeline designed to emulate a clinical trial comparing the effects of different antibiotic treatments on future antibiotic resistance. Key features include:

-   **Advanced Epidemiological Methods:**
    -   **Propensity Score Matching:** Patients receiving different antibiotic treatments are matched based on a rich set of pre-treatment covariates (e.g., age, sex, comorbidities, demographic sector, and season) to reduce confounding.
    -   **Risk Difference Estimation via G-Methods:** The g-formula is applied to standardized risk estimates, allowing the calculation of risk differences (RD) and number needed to harm (NNH) between exposure groups.
    -   **Inverse Probability of Censoring Weighting (IPCW):** Used to adjust for selection bias due to patients not obtaining a follow-up culture within the defined time frame.
-   **Efficient Repeated Jobs:**
    -   The code is structured around custom R objects and functions (e.g., `summarise_resist`, `match_ab`, `rd_calc`, `rd_fig`) that encapsulate common tasks such as data cleaning, matching, estimation, and plotting. This modular approach promotes reproducibility and efficient re-use in repeated analyses for the different exposures and outcomes.

------------------------------------------------------------------------

## Understanding the Repository

### Prerequisites

Ensure you have R (version ≥ 4.0) installed. Install the required packages by running the following commands in R:

`install.packages(c("tidyverse", "MatchIt", "WeightIt", "sandwich", "marginaleffects", "janitor", "ggbreak", "patchwork"))`

**Whats in the Analysis**

The analysis pipeline is divided into several steps. It is recommended to follow the scripts in this order:

1.  **Data Creation & Cleaning:**

`Rscripts/01 - bigfiles.R` is used to clean the EMR data in a format more suitable for coding, later in `Rscripts/02 - df_create.R` we load and preprocess the raw microbiology, antibiotic consumption, covariate, hospitalization, and demographic data to apply inclusion and exclusion criteria and create our cohort. The results are saved in a cohort object with designated functions to inspect the cohort creation process easily later.

2.  **Propensity Score Matching:**

In `Rscripts/03 - matching.R` we perform PS matching between different antibiotic exposure groups. The matching function (match_ab) is designed to efficiently create matched cohorts using the pre-processed data.

3.  **Risk Difference Estimation:**

Taking our matched cohorts, in `Rscripts/04 - analysis.R` we estimate risk differences (RD) using g-methods. The function rd_calc calculates the RD (and NNH) between matched groups, adjusting for residual confounding and censoring.

4.  **Plot Generation:**

Then, using the plotting functions (e.g., `rd_fig` in `Rscripts/04 - analysis.R`) we generate publication-ready figures that illustrate the risk differences and matching diagnostics.

5.  **Utility Functions:**

The functions in `Rscripts/05 - funcs.R` (and any additional functions) support variable labeling, transformation, and visualization, ensuring consistency and efficiency across analyses.

------------------------------------------------------------------------

## Methods

**Propensity Score Matching**

Patients were exactly matched on key baseline culture results and further balanced using PS matching. This process aimed to emulate a randomized trial by balancing measured covariates between patients receiving different antibiotic treatments.

**Risk Difference Estimation via G-Methods**

After matching, a multivariable logistic regression model (weighted by IPCW) was used to standardize risk estimates. The difference in predicted probabilities between scenarios (e.g., all patients receiving the control antibiotic versus the treatment antibiotic) was computed to obtain the RD and its confidence intervals via bootstrap methods.

------------------------------------------------------------------------

**Citation**

If you use this repository in your research, please cite the following article:

```         
Chowers, M., Atias, D., Gottesman, B., Low, M., & Obolski, U. (2025). *Estimating antibiotic resistance following antibiotic treatment in outpatients: a retrospective study.* [Journal Name, Volume(Issue), Pages]. DOI: [doi]
```

**License**

This project is licensed under the [MIT License](LICENSE).

**Contact**

For questions or further information, please contact:

**Dor Atias**

Email: [dor.atiass\@gmail.com](mailto:dor.atiass@gmail.com)
