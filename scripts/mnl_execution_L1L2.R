#Libraries

required_packages <- c("maxLik", "matrixStats", "tidyr", "dplyr", "glmnet", "bgw", 
                       "Rfast", "future.apply", "future", "ggplot2")

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
source("functions/utility_functions_generalized.R") #has function create_alt_matrices2()
source("functions/mnl_function.R")
source("functions/pre_process.R")
source("functions/MNL_functions_execution.R")

#Thresholding
#BHHH + L1/L2 cannot introduce exact sparsity 
SCREENING_THRESHOLD <- 0.01

#Data
data <- read.csv("data/doggerbank_full_973_wide.csv")
#Alternative
#alternate specific constant for choices 2 and 3
data$ASC21 <- 0
data$ASC22 <- 1
data$ASC23 <- 0
data$ASC31 <- 0
data$ASC32 <- 0
data$ASC33 <- 1

#scaling (can be implemented if needed)
# data$income <- scale(data$income)
# data$age <- scale(data$age)
# data$cost <- scale(data$cost)

#Implement functions

output <- data_wide_to_long(data, n_alt = 3)
df_demo <- output$df_demo
df_long <- output$df_long 

choice_vars <- c('ASC2', 'ASC3', 'cost', 'spec10', 'spec25', 'prot25', 'prot50', 'invasive')

demographic_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
                   'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
                   'job5', 'job6', 'job7', 'job8')


final_df_scaled <- create_interaction_features(df_long, choice_vars, demographic_vars)
#Uncomment to consider ONLY interactions
#final_df_scaled <- final_df_scaled[, grepl("_", colnames(final_df_scaled))]

# X <- as.matrix(final_df_scaled)
# y <- as.numeric(df_long$chosen)

#great for comparing our results !!
#selected_features <- run_elastic_net(X, y, alpha = 0.5, n = top_n)
# so I think the ultimate goal would be to replace glment with our own L1/L2 routine and we are that close !!

#ncol(final_df_scaled) gives the maximum number of possible attributes
#marginal + interactions
num_covariates <- 25
selected_features <- colnames(final_df_scaled)[1:num_covariates]

#to use all the covariates
# selected_features <- colnames(final_df_scaled)
# num_covariates <- length(selected_features)


alt_matrices <- create_alt_matrices2(df_demo, 
                                     selected_features = selected_features, 
                                     demographic_vars = demographic_vars, 
                                     n_alt = 3
                                    )

#utility function matrices
alt1 <- alt_matrices$alt1
#print(nrow(alt1))
alt2 <- alt_matrices$alt2
alt3 <- alt_matrices$alt3

n_alt <- 3
alt_list <- lapply(1:n_alt, function(j) alt_matrices[[j]])
choice_list <- lapply(1:n_alt, function(j) df_demo[[paste0("choice", j)]])

plan(multisession)
options(future.rng.onMisuse = "ignore")

#Adjust grid as needed
#lambda_grid <- seq(0.005, 0.03, 0.005)
lambda_grid <- exp(seq(log(6e-4), log(3e-2), length.out = 10))
N <- nrow(df_long)

##Sequential
# results <- lasso_lambda_bic(
#   lambda_grid = lambda_grid,
#   alt_matrices = alt_matrices,
#   df_long = df_long,
#   n = num_covariates,
#   threshold = 1e-4,
#   N
# )

##Parallel
results <- lasso_lambda_bic_parallel(
  lambda_grid = lambda_grid,
  alt_list,
  choice_list,
  n = num_covariates,
  threshold = SCREENING_THRESHOLD,
  N
)

lambda_results <- results$lambda_results

# Plot BIC values
plot_bic <- ggplot(lambda_results, aes(x = lambda, y = BIC)) +
  geom_line(size = 1.5, colour = "grey") +
  geom_point(size = 3) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 26)
  ) +
  labs(
    x = "Lambda",
    y = "BIC"
    #title = "BIC vs Lambda (Elastic Net Regularization)"
  )

# Save (comment/uncomment as needed)
#ggsave("bic_vs_lambda.pdf", plot_bic, width = 7, height = 5, dpi = 300)
ggsave("plots/bic_vs_lambda.png", plot_bic, width = 10, height = 7, dpi = 400)



# legend("topleft", legend = c("BIC", "Non zero coeff"),
#        col = c("blue", "red"), lty = 1, bty = "n")

# #Sequential
# results_cv <- tune_lambda_cv(
#   df_demo, 
#   selected_features, 
#   lambda_grid, 
#   n_alt = 3, 
#   n = num_covariates, 
#   n_folds = 5
#   )

#Parallel
results_cv <- tune_lambda_cv_parallel(
  df_demo, 
  selected_features, 
  lambda_grid,
  demographic_vars, 
  n_alt = 3, 
  n = num_covariates, 
  n_folds = 5
)

lambda_results_cv <- results_cv$lambda_results

#Plot mean LL

plot_cv <-ggplot(lambda_results_cv, aes(x = lambda, y = mean_LL)) +
  geom_line(size = 1.5, colour = "grey") +
  geom_point(size = 3) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 26)
  ) +
  labs(
    x = "Lambda",
    y = "Mean Out-of-Sample LL"
    #title = "Cross-Validated Log-Likelihood vs Lambda"
  )

# Save (comment/uncomment as needed)
#ggsave("plots/cv_ll_vs_lambda.pdf", plot_cv, width = 7, height = 5, dpi = 300)
ggsave("plots/cv_ll_vs_lambda.png", plot_cv, width = 10, height = 7, dpi = 400)



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
  function(coeff) MNL(coeff, alt_list, choice_list, lambda = results_cv$best_lambda, 
                      alpha = 0.5,
                      final_eval = FALSE, nrep = 6, intercept_index = 1),
  start = start.values,
  method = "BHHH",
  iterlim = 200,
  print.level = 0,
  finalHessian = TRUE
)

coefficients_table <- summary_table_mnl(final_model, selected_features)
