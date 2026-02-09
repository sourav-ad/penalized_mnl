## Aiding the Analysis of Systematic Preference Heterogeneity Using Scalable and Regularized Multinomial Logit Models

This repository reproduces an end-to-end pipeline for estimating a multinomial logit (MNL) model with 
elastic net regularization applied to the log-likelihood, with: 

(i) interaction construction for systematic 
preference heterogeneity

(ii) $\lambda$ tuning via BIC and K-fold out-of-sample log-likelihood

(iii) semi-synthetic recovery experiments via Gumbel-noise choice generation

(iv) parallel execution for scalability. 

# Real Data (Dogger Bank)

-  Penalized MNL estimates for a high-dimensional utility specification.
-  $\lambda$ selection diagnostics: BIC($\lambda$) and CV out-of-sample log-likelihood($\lambda$).
-  Coefficient table with shrinkage flags / threshold screening results.

# Semi-synthetic data

-  Controlled induction of interaction effects + recovery evaluation
-  Distribution of $\lambda$ across iterations.
-  Detection results under a threshold grid.



# Usage

Please set working directory and file paths as needed. To run the script for real data and analysis, use:

```         
source("scripts/mnl_execution_L1L2.R") 

```

and for semisynthetic data analysis, use:

```         
source("scripts/MNL_interactions_functional_v3.R") 

```
