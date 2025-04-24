# Multinomial Logit with Elastic Net Regularization

This repository provides an end-to-end pipeline to estimate a **Multinomial Logit (MNL)** discrete choice model with **Elastic Net Regularization** applied to the log likelihood. The model is applied to survey data regarding willingness to pay for a sustainable management plan for [**Dogger Bank**](https://en.wikipedia.org/wiki/Dogger_Bank). It analyzes interactions between marginal and demographic covariates, BIC-based model selection, and cross-validated tuning of the regularization (elastic net) parameter.

------------------------------------------------------------------------

## Features

-   Uses long-format discrete choice experiment data.
-   Constructs interaction terms between user-specific covariates and user-variant covariates.
-   Applies both BIC-based and cross-validation-based tuning of the L1 + L2 regularization parameter.
-   Fits the customized MNL model.
-   Provides coefficient summaries, statistical significance and shrinkage status for covariates.
-   Uses parallelization for scalability.

------------------------------------------------------------------------

## Workflow

1.  **Data Preparation**
    -   Input data in wide format is converted to long format.
    -   Columns with excessive missing data are dropped.
    -   Dummy variables imputed.
    -   Interaction terms computed
2.  **Model Estimation**
    -   The best MNL model is estimated using user selected covariates.
    -   Penalized log likelihood estimation is performed with L1 and L2 penalties.
    -   Parameter tuning is done using:
        -   BIC
        -   5-fold cross-validated out-of-sample log-likelihood.
4.  **Result Interpretation**
    -   Final coefficient estimates are reported with standard errors, t-scores, and significance markers.
    -   Shrunk (near-zero) coefficients are flagged.

------------------------------------------------------------------------

## Dependencies

-   `glmnet`
-   `maxLik`
-   `dplyr`
-   `tidyr`
-   `Rfast`
-   `matrixStats`
-   `bgw`
-   `glmnet`
-   `future.apply`
-   `future`

------------------------------------------------------------------------

## Usage

Please set working directory and file paths as needed. To run, use:

```         
source("scripts/mnl_execution_L1L2.R")
```
