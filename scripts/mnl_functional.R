#Libraries

set.seed(123)

required_packages <- c("maxLik", "matrixStats", "tidyr", "dplyr", "glmnet", "bgw", "Rfast")

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
#Set appropriate working directory
source("functions/utility_functions.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")


choice_model <- function(n_alt, vars_vary, vars_constant, n){
  #Data (change path as needed)
  data <- read.csv("data/doggerbank_full_973_wide.csv")
  #person id
  data$person <- rep(1:(dim(data)[1]/6), each=6)
  
  #Remove columns 'choice1' to 'choice 6' and create new columns instead
  #We modify the data to have an idea about each respondent's choice among the 3 options provided 
  data <- data[, !names(data) %in% c('choice1', 'choice2', 'choice3', 'choice4', 'choice5', 'choice6')]
  data$choice <- 0
  data <- data %>%
    rename('q227' = 'q22_7', 'q229' = 'q22_9')
  
  for (i in 1:n_alt) {
    data$choice[data[[paste0("y", i)]] == 1] <- i
  }
  
  for (i in 1:n_alt) {
    data[[paste0("choice", i)]] <- ifelse(data$choice == i, 1, 0)
  }
  
  df_demo <- data[, c('id', 'line',
                      'invasive1', 'cost1', 'spec101', 'spec251', 'prot251', 'prot501',
                      'invasive2','cost2', 'spec102', 'spec252', 'prot252', 'prot502',
                      'invasive3', 'cost3', 'spec103', 'spec253', 'prot253', 'prot503', 'q227', 'q229',
                      'edu', 'male', 'job', 'age', 'choice', 'y1', 'y2', 'y3',
                      'choice1', 'choice2', 'choice3', 'job1', 'job2', 
                      'job3', 'job4', 'job5', 'job6', 'job7', 'job8', 'q1', 'q2', 'q6', 'q7', 'q10')]
  
  # Convert wide data format to long
  df_long <- df_demo %>%
    pivot_longer(
      cols = starts_with('cost')|starts_with('spec10')|starts_with('spec25')|starts_with('prot25')|starts_with('prot50')|starts_with('invasive')|starts_with('protest'),
      names_to = c(".value", "choice_option"),
      names_pattern = "(.*)([1-3])"
    ) %>%
    mutate(
      choice_option = as.integer(choice_option),  # Convert choice_option to integer
      chosen = as.integer(choice_option == choice)  # Create binary indicator for chosen option
    ) %>%
    arrange(id, line, choice_option)  # Arrange the data for clarity
  
  
  df_long <- df_long %>%
    mutate(
      chosen = ifelse(choice_option == choice, 1, 0)  # 1 if choice_option matches choice, else 0
    )
  
  df_interactions <- df_long[, vars_vary]
  df_interactions_with <- df_long[, vars_constant]
  
  interaction_df <- data.frame(matrix(nrow = nrow(df_interactions), ncol = 0))
  
  # Create interactions **ONLY between df_interactions and df_interactions_with**
  for (col1 in colnames(df_interactions)) {
    for (col2 in colnames(df_interactions_with)) {
      interaction_name <- paste0(col1, "_", col2)
      
      # Create interaction term as the product
      interaction_df[[interaction_name]] <- df_interactions[[col1]] * df_interactions_with[[col2]]
    }
  }
  final_df <- cbind(df_interactions, interaction_df)
  
  #Scale and remove NA if any, before Elastic Net
  final_df_scaled <- final_df %>%
    mutate(across(where(is.numeric), scale))
  
  final_df_scaled <- final_df_scaled[, colSums(is.na(final_df_scaled)) == 0]
  
  X <- as.matrix(final_df_scaled)
  y <- as.numeric(df_long$chosen)
  
  elastic_net_model <- cv.glmnet(
    x = X,
    y = y,
    alpha = 0.5,            # Elastic Net (0.5 = balance between Lasso and Ridge)
    family = "binomial",    # For binary classification
    nfolds = 5,             # 5-fold cross-validation
    maxit = 5000,          # High number of iterations for convergence
    type.measure = "class"  # Classification accuracy
  )
  
  best_lambda <- elastic_net_model$lambda.min
  final_model <- glmnet(
    x = X,
    y = y,
    alpha = 0.5,
    family = "binomial",
    lambda = best_lambda
  )
  
  coefficients <- coef(final_model)
  coefficients_df <- as.data.frame(as.matrix(coefficients))
  colnames(coefficients_df)[1] <- "coefficient"
  coefficients_df$feature <- rownames(coefficients_df)
  
  sorted_coefficients <- coefficients_df %>%
    filter(feature != "(Intercept)") %>%
    arrange(desc(abs(coefficient)))
  
  selected_features <- sorted_coefficients$feature[1:n]
  
  alt_matrices <- create_alt_matrices(df_demo, selected_features, n_alt)
  
  #create_alt_matrices() to be generalized to any number of utility functions first
  
  #utility function matrices (applicable when 3 choices are presented)
  alt1 <- alt_matrices$alt1
  alt2 <- alt_matrices$alt2
  alt3 <- alt_matrices$alt3
  
  #alt_list <- mget(paste0("alt", 1:n_alt))
  alt_list <- list(alt1, alt2, alt3)
  
  nset <- nrow(df_demo)
  start.values <- rep(0, n)
  choice_list <- lapply(1:n_alt, function(i) df_demo[[paste0("choice", i)]])
  
  # estimate model (without inverting Hessian)
  res = maxBFGS(
    #MNL,
    function(coeff) MNL(coeff, alt_list, choice_list),
    grad=NULL,
    hess=NULL,
    start=start.values,
    fixed=NULL,
    print.level=1,
    iterlim=200,
    constraints=NULL,
    tol=1e-25, reltol=1e-25,
    finalHessian=FALSE,
    parscale=rep(1, length=length(start))
  )
  
  # estimate model (with inverting Hessian)
  #new start values are the values obtained above as starting values
  start.values = coef(res)
  
  res = maxLik(
    #MNL,
    function(coeff) MNL(coeff, alt_list, choice_list),  
    grad=NULL, 
    hess=NULL, 
    start=start.values, 
    fixed=NULL, 
    print.level=1, 
    method="BHHH", 
    iterlim=2,
    #iterlim=0, 
    constraints=NULL, 
    tol=1e-04, reltol=1e-04,
    finalHessian=TRUE
  )
  
  #Give names to the betas, top n candidates
  names(res$estimate) = selected_features
  cat(paste(",",res$estimate))
  MNL.res = res
  summary(MNL.res)
  
  #Extracting the significant interactions
  coeff_table <- summary(MNL.res)$estimate
  significant_interactions <- coeff_table[coeff_table[, "Pr(> t)"] < 0.001, ]
  #List of significant (***) covariates AND interactions
  print(significant_interactions)
  
  #Bayesian Information Criteria
  N <- nrow(df_long)
  BIC_value <- -2 * logLik(res) + length(coef(res)) * log(N)
  cat("BIC value for elastic-net + MNL: ", BIC_value, "\n")
  
}


choice_model(
    n_alt = 3, 
    vars_vary = c('cost', 'spec10', 'spec25', 'prot25', 'prot50', 'invasive'),
    vars_constant = c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2', 
                      'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4', 'job5', 
                      'job6', 'job7', 'job8'),
    n = 20
    
)

















