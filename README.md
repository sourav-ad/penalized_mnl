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
- Respondent ID must be named `id` and choice task ID must be named `line`.
- The alternatives must be indexed with `y`, such as `y1`, `y2`,... and exactly one of these should equal 1 per row.
- If the base variable is to be named, for example, `variable`, then the columns must be `variable1`, `variable2`, ...
- Alternative-specific constants (ASCs) must also follow the same pattern.


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