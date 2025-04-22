#Functional MNL with Lasso regularized LL

#Libraries

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
source("functions/utility_functions.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")

#Data
data <- read.csv("data/doggerbank_full_973_wide.csv")


data_wide_to_long <- function(data, n_alt = 3){
  data <- data[, !names(data) %in% c('choice1', 'choice2', 'choice3', 'choice4', 'choice5', 'choice6')]
  data$choice <- 0
  
  for (i in 1:n_alt) {
    data$choice[data[[paste0("y", i)]] == 1] <- i
  }
  
  for (i in 1:n_alt) {
    data[[paste0("choice", i)]] <- ifelse(data$choice == i, 1, 0)
  }
  
  data <- data %>%
    rename('q227' = 'q22_7', 'q229' = 'q22_9') %>%
    select(-which(colMeans(is.na(.)) > 0.5))
  
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
  
  return(list(df_demo = df_demo, df_long = df_long))
}



create_interaction_features <- function(df_long, constant_vars, changing_vars){
  df_interactions <- df_long[, constant_vars]
  df_interactions_with <- df_long[, changing_vars]
  interaction_df <- data.frame(matrix(nrow = nrow(df_interactions), ncol = 0))
  
  for (col1 in colnames(df_interactions)) {
    for (col2 in colnames(df_interactions_with)) {
      interaction_name <- paste0(col1, "_", col2)
      
      # Create interaction term as the product
      interaction_df[[interaction_name]] <- df_interactions[[col1]] * df_interactions_with[[col2]]
    }
  }
  
  final_df <- cbind(df_interactions, interaction_df)
  final_df_scaled <- final_df %>%
    mutate(across(where(is.numeric), scale))
  final_df_scaled <- final_df_scaled[, colSums(is.na(final_df_scaled)) == 0]
  
  return(final_df_scaled)
}

# X <- as.matrix(final_df_scaled)
# y <- as.numeric(df_long$chosen)

run_elastic_net <- function(X, y, alpha = 0.5, n = 15){
  set.seed(123)
  elastic_net_model <- cv.glmnet(
    x = X,
    y = y,
    alpha = 0.5,            # Elastic Net (0.5 = balance between Lasso and Ridge)
    family = "binomial",    # For binary classification
    nfolds = 5,             # 5-fold cross-validation
    maxit = 5000,           # High number of iterations for convergence
    type.measure = "class"  # Classification accuracy
  )
  elastic_net_parameter <- elastic_net_model$lambda.min
  
  final_model <- glmnet(
    x = X,
    y = y,
    alpha = 0.5,
    family = "binomial",
    lambda = elastic_net_parameter
  )
  
  coefficients <- coef(final_model)
  coefficients_df <- as.data.frame(as.matrix(coefficients))
  colnames(coefficients_df)[1] <- "coefficient"
  coefficients_df$feature <- rownames(coefficients_df)
  
  sorted_coefficients <- coefficients_df %>%
    filter(feature != "(Intercept)") %>%
    arrange(desc(abs(coefficient)))
  
  selected_features <- sorted_coefficients$feature[1:n]
  return(selected_features)
}


lasso_lambda_bic <- function(lambda_grid, alt1, alt2, alt3, df_long, n, 
                             threshold = 1e-3) {
  best_lambda <- NULL
  best_BIC <- Inf
  best_res <- NULL
  lambda_results <- data.frame(lambda = lambda_grid, BIC = NA)
  N <- nrow(df_long)
  
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    start.values <- rep(0, n)
    
    res <- maxBFGS(
      function(coeff) MNL(coeff, alt1, alt2, alt3, lambda, final_eval = FALSE),
      grad = NULL,
      hess = NULL,
      start = start.values,
      fixed = NULL,
      print.level = 0,
      iterlim = 200,
      constraints = NULL,
      tol = 1e-25,
      reltol = 1e-25,
      finalHessian = FALSE,
      parscale = rep(1, length(start.values))
    )
    
    invisible(MNL(res$estimate, alt1, alt2, alt3, lambda, final_eval = TRUE))
    start.values <- coef(res)
    
    res <- maxLik(
      function(coeff) MNL(coeff, alt1, alt2, alt3, lambda, final_eval = TRUE),
      grad = NULL,
      hess = NULL,
      start = start.values,
      fixed = NULL,
      print.level = 0,
      method = "BHHH",
      iterlim = 2,
      constraints = NULL,
      tol = 1e-04,
      reltol = 1e-04,
      finalHessian = TRUE
    )
    
    LL_unpenalized <- sum(MNL(res$estimate, alt1, alt2, alt3, lambda = 0, final_eval = FALSE))
    active_coeffs <- coef(res)[abs(coef(res)) >= threshold]
    k <- length(active_coeffs)
    BIC_lasso <- -2 * LL_unpenalized + k * log(N)
    
    lambda_results$BIC[i] <- BIC_lasso
    
    if (BIC_lasso < best_BIC) {
      best_lambda <- lambda
      best_BIC <- BIC_lasso
      best_res <- res
    }
  }
  
  # Plot the BIC curve
  plot(lambda_results$lambda, lambda_results$BIC, type = "b",
       xlab = "Lambda (L1 Penalty)", ylab = "BIC",
       main = "Model Selection using BIC")
  
  LL_unpenalized_best <- sum(MNL(best_res$estimate, alt1, alt2, alt3, lambda = 0, final_eval = FALSE))
  
  return(list(
    best_lambda = best_lambda,
    best_BIC = best_BIC,
    best_model = best_res,
    best_LL = LL_unpenalized_best,
    lambda_results = lambda_results
  ))
}


tune_lambda_cv <- function(df_demo, selected_features, lambda_grid, n_alt = 3, n = 15, n_folds = 5) {
  #Create folds (respondent-wise split)
  set.seed(123)
  id_list <- unique(df_demo$id)
  folds <- cut(seq_along(id_list), breaks = n_folds, labels = FALSE)
  id_folds <- split(id_list, folds)
  
  lambda_results <- data.frame(lambda = lambda_grid, mean_LL = NA)
  best_lambda <- NULL
  best_LL <- -Inf
  
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    fold_lls <- numeric(n_folds)
    
    for (fold in 1:n_folds) {
      test_ids <- id_folds[[fold]]
      train_ids <- setdiff(id_list, test_ids)
      
      train_df <- df_demo[df_demo$id %in% train_ids, ]
      test_df <- df_demo[df_demo$id %in% test_ids, ]
      
      alt_train <- create_alt_matrices(train_df, selected_features)
      alt_test  <- create_alt_matrices(test_df, selected_features)
      
      alt_list_train <- lapply(1:n_alt, function(j) alt_train[[j]])
      alt_list_test  <- lapply(1:n_alt, function(j) alt_test[[j]])
      
      choice_list_train <- lapply(1:n_alt, function(j) train_df[[paste0("choice", j)]])
      choice_list_test  <- lapply(1:n_alt, function(j) test_df[[paste0("choice", j)]])
      
      start.values <- rep(0, n)
      
      res <- maxBFGS(
        function(coeff) MNL_cv(coeff, alt_list_train, choice_list_train, lambda),
        start = start.values,
        print.level = 0,
        iterlim = 200,
        finalHessian = FALSE
      )
      
      # Evaluate unpenalized LL on test data
      ll_out_sample <- MNL_cv(res$estimate, alt_list_test, choice_list_test, lambda = 0)
      fold_lls[fold] <- sum(ll_out_sample)
    }
    
    mean_LL <- mean(fold_lls)
    lambda_results$mean_LL[i] <- mean_LL
    
    if (mean_LL > best_LL) {
      best_LL <- mean_LL
      best_lambda <- lambda
    }
  }
  
  # Plot
  plot(lambda_results$lambda, lambda_results$mean_LL,
       type = "b", xlab = "Lambda",
       ylab = "Mean Out-of-Sample Log-Likelihood",
       main = "CV Log-Likelihood vs Lambda")
  
  cat("\n===== Lambda tuning summary (CV) =====\n")
  print(lambda_results)
  cat("\nBest lambda based on mean out-of-sample LL:", best_lambda, "\n")
  
  return(list(best_lambda = best_lambda, lambda_results = lambda_results))
}

#Implement functions

output <- data_wide_to_long(data, n_alt = 3)
df_demo <- output$df_demo
df_long <- output$df_long 

constant_vars <- c('cost', 'spec10', 'spec25', 'prot25', 'prot50', 'invasive')

changing_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
                   'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
                   'job5', 'job6', 'job7', 'job8')


final_df_scaled <- create_interaction_features(df_long, constant_vars, changing_vars)

X <- as.matrix(final_df_scaled)
y <- as.numeric(df_long$chosen)

selected_features <- run_elastic_net(X, y, alpha = 0.5, n = 10)


alt_matrices <- create_alt_matrices(df_demo, selected_features)

#utility function matrices
alt1 <- alt_matrices$alt1
alt2 <- alt_matrices$alt2
alt3 <- alt_matrices$alt3

lambda_grid <- seq(0.001, 0.01, 0.001)

results <- lasso_lambda_bic(
  lambda_grid = lambda_grid,
  alt1 = alt1,
  alt2 = alt2,
  alt3 = alt3,
  df_long = df_long,
  n = 10
)


results_cv <- tune_lambda_cv(
  df_demo, 
  selected_features, 
  lambda_grid, 
  n_alt = 3, 
  n = 10, 
  n_folds = 5)