# Aiding the Analysis of Systematic Preference Heterogeneity Using Scalable and Regularized Multinomial Logit Models

This repository reproduces an end-to-end pipeline for estimating a multinomial logit (MNL) model with 
elastic net regularization applied to the log-likelihood, with: 

(i) discovery of systematic preference heterogeneity

(ii) screening tool for practitioners to find relevant interactions

(iii) regularization parameter ($\lambda$) tuning using 5-fold out-of-sample log-likelihood

(iii) semi-synthetic coefficient recovery experiments using Gumbel-noise based choice generation

(iv) choice of optimizer between BFGS and BHHH 

(v) scalability using parallel execution

## Expected data format

The pipeline expects the input data in **wide format**, 
with one row per respondent and repeated choice-task information stored across columns. 
The scripts internally reshape the data into **long format** before estimation.

## Real Data (Dogger Bank)

-  Penalized MNL estimates for high-dimensional utility space.
-  $\lambda$ selection diagnostics: CV out-of-sample log-likelihood.
-  Coefficient table for screening.

## Semi-synthetic data

-  Controlled induction of interaction effects into the Dogger Bank data.
-  Recovery analysis of those induced coefficients.
-  Distribution of $\lambda$ across iterations.
-  Detection results under a threshold grid.
-  Receiver Operating Characteristic (ROC) curves for different SNR (signal-to-noise ratios) 

## Using own dataset

- The input must be in wide format, with one row per respondent-choice task.
- Use the script `scripts/run_model_general.R` for the purpose. Some examples are already provided along with data.
- Running the script needs the user to provide a `data_schema` to identify relevant columns. 
The scripts do not create `ASC` columns. 
- If column names contain any special characters other than underscores `_`, please remove them,  
as underscore `_` is used to denote interacting features. A line of code already in the script can help achieve it.
- The best way to name columns will be `<variablename><choicenumber>`, for example, `spec101` or `beach501`. 
- When providing choice specific variables in `data_schema`, just the `<variablename>` should be provided. 
- No need to provide a full `lambda_grid` but at times, a lower `lambda_min` or a higher `lambda_max` may be 
needed depending on the data.
- A threshold validation is performed after the optimal $\lambda$ has been selected. 
`threshold_rule = "best_ll"` selects the threshold with the highest held-out log-likelihood after pruning 
, `threshold_rule = "sparsest_within_tolerance"` selects the sparsest threshold whose 
held-out log-likelihood is within `threshold_ll_tolerance` of the best held-out log-likelihood and 
`threshold_rule = "sparsest_within_tolerance"` with `threshold_ll_tolerance = 1` 
selects the threshold retaining the fewest coefficients, provided its held-out log-likelihood 
is no more than `1` unit worse than the best threshold.
- Before running the full model, please ensure that all required columns exist, choice indicators are valid, 
`n_alt` is correct, `nrep` is correct, and the number of generated features is plausible.


## Quick start guide

1. Download or clone this repository.

2. Open the project in R or RStudio.

3. Set the working directory to the repository root.

4. Make sure the working directory is set to the repository root and 
   that the required input files are available in the expected folders.

5. Run one of the following scripts:

- For real-data:

  `source("scripts/run_model_penalized_mnl.R")`

- For semi-synthetic simulation study:

  `source("scripts/MNL_interactions_functional_v3.R")`
  
6. The Receiver Operating Characteristic (ROC) curves in the paper can be created using `scripts/roc_plot.R`.

The scripts will create output files in the `data/` and `plots/` folders with appropriate 
folder and file naming based on parameters used.

## Repository structure

- `scripts/`  

  Main scripts to run penalized MNL for real data and semi-synthetic simulations.

- `functions/`  

  Helper functions used by the main scripts.

- `data/`  

  Stores input data and generated output files and folders.

- `plots/`  

  Stores generated figures.
  
## Notes

Installing the R packages, if not already installed, may take some time. 
It is recommended to run the scripts on a machine with multiple cores (e.g. cluster).