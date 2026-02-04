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

#Function
force_feature_interaction <- function(n_persons = 400,
                                      lambda_grid = exp(seq(log(1e-4), log(5e-1), length.out = 20)),
                                      iterations = 5) {
  
  # data
  data <- read.csv("data/doggerbank_full_973_wide.csv")
  selected_ids <- sample(unique(data$id), n_persons)
  data <- data[data$id %in% selected_ids, ]
  
  # alternatives
  data$ASC21 <- 0; data$ASC22 <- 1; data$ASC23 <- 0
  data$ASC31 <- 0; data$ASC32 <- 0; data$ASC33 <- 1
  #scaling

  data$income <- data$income/10000
  data$age <- data$age/100
  data$cost2 <- data$cost2/100
  data$cost3 <- data$cost3/100
  
  alt1 = cbind(0, 0, data$spec101, data$spec251, data$prot251, data$prot501, data$invasive1, data$cost1,
               data$ASC21 * data$male, data$ASC31 * data$age, data$cost1 * data$income,
               data$prot251 * data$child, data$prot501 * data$child, data$invasive1 * data$child)
  
  alt2 = cbind(1, 0, data$spec102, data$spec252, data$prot252, data$prot502, data$invasive2, data$cost2,
               data$ASC22 * data$male, data$ASC32 * data$age, data$cost2 * data$income,
               data$prot252 * data$child, data$prot502 * data$child, data$invasive2 * data$child)
  
  alt3 = cbind(0, 1, data$spec103, data$spec253, data$prot253, data$prot503, data$invasive3, data$cost3,
               data$ASC23 * data$male, data$ASC33 * data$age, data$cost3 * data$income,
               data$prot253 * data$child, data$prot503 * data$child, data$invasive3 * data$child)
  
  true.betas.main <- c(-0.4, -0.4, 0.2, 0.3, 1, 1.3, -0.9, -0.04)
  true.betas.ia   <- c(-0.4, -0.2, 0, 0.5, 0, -0.4) #inducing 4 signals, 2 noise
  #true.betas.ia   <- c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
  true.betas = c(true.betas.main, true.betas.ia)
  
  choice_vars <- c('ASC2', 'ASC3', 'cost', 'prot25', 'prot50', 'invasive')
  demographic_vars <- c('male', 'age', 'income', 'child')
  n_alt <- 3
  #true_betas_full <- NULL
  
  # ---- BUILD FEATURE UNIVERSE (ONCE) ----
  output_tmp <- data_wide_to_long(data, n_alt = 3)
  df_long_tmp <- output_tmp$df_long
  
  final_df_scaled_tmp <- create_interaction_features(
    df_long_tmp,
    choice_vars,
    demographic_vars
  )
  
  # keep only interaction terms
  final_df_scaled_tmp <- final_df_scaled_tmp[, grepl("_", colnames(final_df_scaled_tmp))]
  
  all_features <- colnames(final_df_scaled_tmp)
  
  # ---- DEFINE GROUND TRUTH (ONCE) ----
  true_betas_full <- setNames(
    numeric(length(all_features)),
    all_features
  )
  
  true_betas_full[c(
    "ASC2_male",
    "ASC3_age",
    "cost_income",
    "prot25_child",
    "prot50_child",
    "invasive_child"
  )] <- c(-0.4, -0.2, 0, 0.5, 0, -0.4)
  
  plan(multisession, workers = 50)
  options(future.rng.onMisuse = "ignore")
  
  run_once <- function(rep_id) {
    
    set.seed(1231 + rep_id); gumbel_alt1 <- rgumbel(nrow(alt1), loc=0, scale=1)
    set.seed(2231 + rep_id); gumbel_alt2 <- rgumbel(nrow(alt2), loc=0, scale=1)
    set.seed(3231 + rep_id); gumbel_alt3 <- rgumbel(nrow(alt3), loc=0, scale=1)
    
    util1 = alt1 %*% true.betas + gumbel_alt1
    util2 = alt2 %*% true.betas + gumbel_alt2
    util3 = alt3 %*% true.betas + gumbel_alt3
    synth.choices = apply(cbind(util1, util2, util3), 1, which.max)
    
    data$syn_choice <- synth.choices
    data$y1 <- ifelse(data$syn_choice == 1, 1, 0)
    data$y2 <- ifelse(data$syn_choice == 2, 1, 0)
    data$y3 <- ifelse(data$syn_choice == 3, 1, 0)
    
    output <- data_wide_to_long(data, n_alt = 3)
    df_demo <- output$df_demo
    df_long <- output$df_long
    
    final_df_scaled <- create_interaction_features(df_long, choice_vars, demographic_vars)
    #keeps only interactions
    final_df_scaled <- final_df_scaled[, grepl("_", colnames(final_df_scaled))]
    
    required <- c("ASC2_male","ASC3_age","cost_income",
                  "prot25_child","prot50_child","invasive_child")
    missing <- setdiff(required, colnames(final_df_scaled))
    
    if (length(missing)) stop("Missing interaction columns: ", paste(missing, collapse=", "))
    
    #selected_features <- colnames(final_df_scaled)[1:num_covariates]
    selected_features <- colnames(final_df_scaled)
    
    alt_matrices <- create_alt_matrices2(df_demo,
                                         selected_features = selected_features,
                                         demographic_vars = demographic_vars,
                                         n_alt = 3)
    alt_list <- lapply(1:n_alt, function(j) alt_matrices[[j]])
    choice_list <- lapply(1:n_alt, function(j) df_demo[[paste0("choice", j)]])
    
    results_cv <- tune_lambda_cv_parallel(df_demo,
                                          selected_features,
                                          lambda_grid,
                                          demographic_vars,
                                          n_alt = 3,
                                          #n = num_covariates,
                                          n = length(selected_features),
                                          n_folds = 5,
                                          alpha = 1, # Lasso
                                          #alpha = 0, #ridge
                                          #alpha = 0.5, #elastic net
                                          )
    
    start.values <- rep(0, length(selected_features))
    final_model <- maxLik(
      function(coeff) MNL(coeff, alt_list, choice_list,
                          lambda = results_cv$best_lambda,
                          alpha = 1, # Lasso
                          #alpha = 0, #ridge
                          #alpha = 0.5, #elastic net
                          final_eval = FALSE,
                          nrep = 6, 
                          #intercept_index = 1
                          intercept_index = NULL # let all coefficients be penalized
                          ),
      start = start.values,
      method = "BHHH",
      iterlim = 200,
      print.level = 0,
      finalHessian = TRUE
    )

    ##p-value based threshold detection for confusion matrix
    beta_hat <- coef(final_model)
    vcov_mat <- vcov(final_model)
    
    se_hat <- sqrt(diag(vcov_mat))
    z_stat <- beta_hat / se_hat
    p_val  <- 2 * (1 - pnorm(abs(z_stat)))
    ##
    
    coefs <- as.numeric(coef(final_model))
    names(coefs) <- selected_features
    coefs_full <- setNames(numeric(length(all_features)), all_features)
    coefs_full[selected_features] <- coefs
    #return(list(coefs_full = coefs_full, best_lambda = results_cv$best_lambda))
    return(list(
      coefs_full  = coefs_full,
      p_val       = p_val,
      best_lambda = results_cv$best_lambda
    ))
    
  }
  
  results_list <- vector("list", iterations)
  best_lambdas <- numeric(iterations)
  
  for (i in 1:iterations) {
    res <- run_once(rep_id = i)
    #results_list[[i]] <- run_once(rep_id = i)
    results_list[[i]] <- res$coefs_full
    best_lambdas[i] <- res$best_lambda
  }
  
  #Save best lambda per iteration
  lambda_df <- data.frame(Iteration = 1:iterations, BestLambda = best_lambdas)
  lambda_filename <- file.path("data",
                               paste0(n_persons, "n_", iterations, "_iter_lambda_values.csv"))
  write.csv(lambda_df, lambda_filename, row.names = FALSE)
  
  results <- do.call(cbind, results_list)
  results_df <- as.data.frame(t(results))
  
  
  alpha_sig <- 0.01   # or 0.001
  #detected_mat <- p_val < alpha_sig
  detected_mat[iter, ] <- p_val < alpha_sig
  #Everything works upto this
  #TP TN FP FN
  active_features   <- names(true_betas_full[true_betas_full != 0])
  inactive_features <- names(true_betas_full[true_betas_full == 0])
  
  #stopping criteria, failsafe
  stopifnot(setequal(names(true_betas_full), colnames(results_df)))
  
  # tp_rate <- colMeans(results_df[, active_features, drop = FALSE] != 0, na.rm = TRUE)
  # fp_rate <- colMeans(results_df[, inactive_features, drop = FALSE] != 0, na.rm = TRUE)
  
  tp_rate <- colMeans(detected_mat[, active_features, drop = FALSE], na.rm = TRUE)
  fp_rate <- colMeans(detected_mat[, inactive_features, drop = FALSE], na.rm = TRUE)
  
  fn_rate <- 1 - tp_rate
  tn_rate <- 1 - fp_rate
  
  # TP <- sum(results_df[, active_features, drop = FALSE] != 0, na.rm = TRUE)
  # FN <- sum(results_df[, active_features, drop = FALSE] == 0, na.rm = TRUE)
  # 
  # FP <- sum(results_df[, inactive_features, drop = FALSE] != 0, na.rm = TRUE)
  # TN <- sum(results_df[, inactive_features, drop = FALSE] == 0, na.rm = TRUE)
  
  TP <- sum(detected_mat[, active_features, drop = FALSE], na.rm = TRUE)
  FN <- sum(!detected_mat[, active_features, drop = FALSE], na.rm = TRUE)

  FP <- sum(detected_mat[, inactive_features, drop = FALSE], na.rm = TRUE)
  TN <- sum(!detected_mat[, inactive_features, drop = FALSE], na.rm = TRUE)
  
  confusion_matrix <- data.frame(
    TrueStatus = c("Induced", "Induced", "Not Induced", "Not Induced"),
    Predicted  = c("Detected", "Not Detected", "Detected", "Not Detected"),
    Count      = c(TP, FN, FP, TN)
  )
  
  write.csv(
    confusion_matrix,
    file.path("data", paste0(n_persons, "n_", iterations, "_confusion_matrix.csv")),
    row.names = FALSE
  )
  
  #confusion rates
  confusion_rates <- data.frame(
    TPR = TP / (TP + FN),
    FPR = FP / (FP + TN),
    TNR = TN / (TN + FP),
    FNR = FN / (FN + TP)
  )
  
  write.csv(
    confusion_rates,
    file.path("data", paste0(n_persons, "n_", iterations, "_confusion_rates.csv")),
    row.names = FALSE
  )
  
  
  #cat("\n=== Detection Summary ===\n")
  
  #cat("True Positives (active recovered):\n"); print(tp_rate)
  #cat("\nFalse Negatives (active missed):\n"); print(fn_rate)
  #cat("\nFalse Positives (inactive wrongly picked):\n"); print(head(fp_rate, 10))  # print first few
  #cat("\nTrue Negatives (inactive correctly zero):\n"); print(head(tn_rate, 10))
  
  all_feats <- names(true_betas_full)
  
  detection_summary <- data.frame(
    Feature    = all_feats,
    TrueBeta   = as.numeric(true_betas_full[all_feats]),
    TrueStatus = ifelse(all_feats %in% active_features, "Active", "Inactive"),
    DetectRate = colMeans(detected_mat[, all_feats, drop = FALSE], na.rm = TRUE)
    #DetectRate = colMeans(results_df[, all_feats, drop = FALSE] != 0, na.rm = TRUE)
  )
  
  # Add interpretable rate labels:
  # Active: DetectRate = TPR per interaction
  # Inactive: DetectRate = FPR per interaction
  
  summary_filename <- file.path("data",
                                paste0(n_persons, "n_", iterations, "_iter_detections.csv"))
  write.csv(detection_summary, summary_filename, row.names = FALSE)
  
  results_active <- results_df[, active_features, drop = FALSE]
  results_inactive <- results_df[, inactive_features, drop = FALSE]
  
  interaction_features <- c("ASC2_male", "ASC3_age", "cost_income",
                            "prot25_child", "prot50_child", "invasive_child")
  results_df_interactions <- results_df[, interaction_features, drop = FALSE]
  #Save the values as csv
  csv_filename <- file.path("data", paste0(n_persons, "n_", iterations, "_iter_results.csv"))
  write.csv(results_df_interactions, csv_filename, row.names = FALSE)
  
  true_values <- data.frame(
    Feature = c("ASC2_male", "ASC3_age", "cost_income",
                "prot25_child", "prot50_child", "invasive_child"),
    Forced = c(-0.4, -0.2, 0, 0.5, 0, -0.4)
    #Forced = c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
  )
  
  true_map <- setNames(true_values$Forced, true_values$Feature)
  
  error_df <- results_df_interactions
  for (f in names(true_map)) {
    error_df[[f]] <- error_df[[f]] - true_map[f]
  }
  
  plot_data <- tidyr::pivot_longer(error_df,
                                   cols = everything(),
                                   names_to = "Feature",
                                   values_to = "Error")
  
  p <- ggplot(plot_data, aes(x = Feature, y = Error, fill = Feature)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8) + #not the elastic net alpha
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(
      title = "Recovery of Induced Interaction Effects",
      subtitle = "Boxplots show estimation error (estimated − true)",
      x = "Interaction Feature",
      y = "Estimation Error"
    )
  
  filename <- file.path("plots", paste0(n_persons, "n_", iterations, "iter.png"))
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  
  print(p)
}


for (i in c(50)) {
  force_feature_interaction(
    n_persons = i,
    iterations = 3
  )
}
