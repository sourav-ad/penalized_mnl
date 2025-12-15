##INDUCE AND TEST IF INTERACTIONS CAN BE RECOVERED
#Libraries

required_packages <- c("maxLik", "matrixStats", "tidyr", "dplyr", "glmnet", "bgw", 
                       "Rfast", "future.apply", "future", "evd", "ggplot2")

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

data <- read.csv("data/doggerbank_full_973_wide.csv")

#Taking subset of persons using their id 
n_persons <- 100   # or 500,600,700,800,900
selected_ids <- sample(unique(data$id), n_persons)

data <- data[data$id %in% selected_ids, ]
#Alternative
data$ASC21 <- 0
data$ASC22 <- 1
data$ASC23 <- 0
data$ASC31 <- 0
data$ASC32 <- 0
data$ASC33 <- 1
#scaling
data$income <- data$income/10000

### extract the experimental design / choice alternatives from the dataset + 
alt1 = cbind(0, 0, data$spec101, data$spec251, data$prot251, data$prot501, data$invasive1, data$cost1, 
             data$ASC21 * data$male, data$ASC31 * data$age, data$cost1 * data$income
             , data$prot251 * data$child, data$prot501 * data$child, data$invasive1 * data$child)

alt2 = cbind(1, 0, data$spec102, data$spec252, data$prot252, data$prot502, data$invasive2, data$cost2, 
             data$ASC22 * data$male, data$ASC32 * data$age, data$cost2 * data$income
             , data$prot252 * data$child, data$prot502 * data$child, data$invasive2 * data$child)

alt3 = cbind(0, 1, data$spec103, data$spec253, data$prot253, data$prot503, data$invasive3, data$cost3, 
             data$ASC23 * data$male, data$ASC33 * data$age, data$cost3 * data$income
             , data$prot253 * data$child, data$prot503 * data$child, data$invasive3 * data$child)


resp.vars = cbind(data$age, data$male, data$income) 

### Use approx. betas from Table 5a in Doggerbank paper (may be changed afterwards)
true.betas.main = c(-0.4, -0.4, 0.2, 0.3, 1, 1.3, -0.9, -0.04)
### Assume some interactions
true.betas.ia   = c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
true.betas = c(true.betas.main, true.betas.ia)

#Does not change
choice_vars <- c('ASC2', 'ASC3', 'cost', 'prot25', 'prot50', 'invasive')
demographic_vars <- c('male', 'age', 'income', 'child')

num_covariates <- 24

n_alt <- 3
lambda_grid <- seq(0.001, 0.01, 0.001)

#Iterative function

plan(multisession)
options(future.rng.onMisuse = "ignore")

run_once <- function(rep_id) {
  
  set.seed(1231 + rep_id) ; gumbel_alt1 <- rgumbel(nrow(alt1), loc=0, scale=1)
  set.seed(2231 + rep_id) ; gumbel_alt2 <- rgumbel(nrow(alt1), loc=0, scale=1)
  set.seed(3231 + rep_id) ; gumbel_alt3 <- rgumbel(nrow(alt1), loc=0, scale=1)
  
  ### Compute utilities
  util1 = alt1 %*% true.betas + gumbel_alt1
  util2 = alt2 %*% true.betas + gumbel_alt2
  util3 = alt3 %*% true.betas + gumbel_alt3
  utils = cbind(util1, util2, util3)
  synth.choices = apply(utils, 1, which.max)
  
  data$syn_choice <- synth.choices
  data$y1 <- ifelse(data$syn_choice == 1, 1, 0)
  data$y2 <- ifelse(data$syn_choice == 2, 1, 0)
  data$y3 <- ifelse(data$syn_choice == 3, 1, 0)
  
  output <- data_wide_to_long(data, n_alt = 3)
  df_demo <- output$df_demo
  df_long <- output$df_long 
  
  final_df_scaled <- create_interaction_features(df_long, choice_vars, demographic_vars)
  final_df_scaled <- final_df_scaled[, grepl("_", colnames(final_df_scaled))]  
  
  selected_features <- colnames(final_df_scaled)[1:num_covariates]
  
  alt_matrices <- create_alt_matrices2(df_demo, 
                                       selected_features = selected_features, 
                                       demographic_vars = demographic_vars, 
                                       n_alt = 3)
  
  alt1 <- alt_matrices$alt1
  alt2 <- alt_matrices$alt2
  alt3 <- alt_matrices$alt3
  
  alt_list <- lapply(1:n_alt, function(j) alt_matrices[[j]])
  choice_list <- lapply(1:n_alt, function(j) df_demo[[paste0("choice", j)]])
  
  N <- nrow(df_long)
  
  results_cv <- tune_lambda_cv_parallel(
    df_demo, 
    selected_features, 
    lambda_grid,
    demographic_vars, 
    n_alt = 3, 
    n = num_covariates, 
    n_folds = 5
  )
  
  start.values <- rep(0, length(selected_features))
  final_model <- maxLik(
    function(coeff) MNL(coeff, alt_list, choice_list,
                        lambda = results_cv$best_lambda,
                        alpha = 0.5, final_eval = FALSE,
                        nrep = 6, intercept_index = 1),
    start = start.values,
    method = "BHHH",
    iterlim = 200,
    print.level = 0,
    finalHessian = TRUE
  )
  
  #return(coef(final_model))
  coefs <- coef(final_model)
  names(coefs) <- selected_features
  return(coefs)
}

# Outer sequential repetition
iterations <- 5
results_list <- vector("list", iterations)
for (i in 1:iterations) {
  results_list[[i]] <- run_once(rep_id = i)
}

#rownames(results) <- selected_features

# results is a matrix: rows = coefficients, cols = repetition
results <- do.call(cbind, results_list)
results_df <- as.data.frame(t(results))

# Boxplot of coefficients
#boxplot(results_df, las = 2, main = "Coefficient Recovery Across 100 Runs")

interaction_features <- c("ASC2_male", "ASC3_age", "cost_income",
                          "prot25_child", "prot50_child", "invasive_child")

results_df_interactions <- results_df[, interaction_features, drop = FALSE]

# boxplot(results_df_interactions, las = 2,
#         main = "Recovery of True Interaction Coefficients")

plot_data <- tidyr::pivot_longer(
  results_df_interactions,
  cols = everything(),
  names_to = "Feature",
  values_to = "Coefficient"
)

ggplot(plot_data, aes(x = Feature, y = Coefficient, fill = Feature)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(title = "Recovery of True Interaction Coefficients",
       x = "Interaction Feature",
       y = "Coefficient Estimate")
