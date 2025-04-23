#Libraries

required_packages <- c("maxLik", "matrixStats", "tidyr", "dplyr", "glmnet", "bgw", 
                       "Rfast")

install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  if(length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)  # Install missing packages
  }
  
  # Load all packages
  lapply(packages, require, character.only = TRUE)
}

install_if_missing(required_packages)

#Source files
source("functions/utility_functions.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")
source("functions/MNL_functions_execution.R")

#Data
data <- read.csv("data/doggerbank_full_973_wide.csv")

#Implement functions

output <- data_wide_to_long(data, n_alt = 3)
df_demo <- output$df_demo
df_long <- output$df_long 

constant_vars <- c('cost', 'spec10', 'spec25', 'prot25', 'prot50', 'invasive')

changing_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
                   'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
                   'job5', 'job6', 'job7', 'job8')


final_df_scaled <- create_interaction_features(df_long, constant_vars, changing_vars)

# X <- as.matrix(final_df_scaled)
# y <- as.numeric(df_long$chosen)

#great for comparing our results !!
#selected_features <- run_elastic_net(X, y, alpha = 0.5, n = top_n)
# so I think the ultimate goal would be to replace glment with our own L1/L2 routine and we are that close !!

num_covariates <- 15
selected_features <- colnames(final_df_scaled)[1:num_covariates]


alt_matrices <- create_alt_matrices(df_demo, selected_features)

#utility function matrices
alt1 <- alt_matrices$alt1
alt2 <- alt_matrices$alt2
alt3 <- alt_matrices$alt3

#Adjust grid as needed
lambda_grid <- seq(0.001, 0.01, 0.001)

results <- lasso_lambda_bic(
  lambda_grid = lambda_grid,
  alt_matrices = alt_matrices,
  df_long = df_long,
  n = num_covariates
)

lambda_results <- results$lambda_results

# Plot BIC values
plot(lambda_results$lambda, lambda_results$BIC, type = "b", col = "blue",
     xlab = "L1 lambda)", ylab = "BIC",
     main = "Lasso parameter with BIC")

# legend("topleft", legend = c("BIC", "Non zero coeff"),
#        col = c("blue", "red"), lty = 1, bty = "n")

results_cv <- tune_lambda_cv(
  df_demo, 
  selected_features, 
  lambda_grid, 
  n_alt = 3, 
  n = num_covariates, 
  n_folds = 5)

lambda_results_cv <- results_cv$lambda_results

#Plot mean LL

plot(lambda_results_cv$lambda, lambda_results_cv$mean_LL,
     type = "b", xlab = "Lambda",
     ylab = "Mean Out-of-Sample Log-Likelihood",
     main = "CV Log-Likelihood vs Lambda")



# To do:
# 1. Add L2 (1-alpha)/2 * sum(coeff2) + alpha * sum(abs(coeff)) (done)
# 2. glmnet elastic net as comparison only and not for choosing number of covariates (too slow)
# 3. Just throw everything at regularizing, all the covariates (too slow)
# 4. compact script, all functions here to a source file (done)
# 5. Plots outside the functions (done)


#Show the final model
#Using CV tuned parameter

start.values <- rep(0, length(selected_features))

final_model <- maxLik(
  function(coeff) MNL(coeff, alt1, alt2, alt3, lambda = results_cv$best_lambda, alpha = 0.5),
  start = start.values,
  method = "BHHH",
  iterlim = 200,
  print.level = 0,
  finalHessian = TRUE
)

coefficients_table <- summary_table_mnl(final_model, selected_features)
