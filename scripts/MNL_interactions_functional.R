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
                           #num_covariates = 24,
                           lambda_grid = seq(0.001, 0.01, 0.001),
                           iterations = 5) {
  
  # data
  data <- read.csv("data/doggerbank_full_973_wide.csv")
  selected_ids <- sample(unique(data$id), n_persons)
  data <- data[data$id %in% selected_ids, ]
  
  # alternatives
  data$ASC21 <- 0; data$ASC22 <- 1; data$ASC23 <- 0
  data$ASC31 <- 0; data$ASC32 <- 0; data$ASC33 <- 1
  data$income <- data$income/10000
  
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
  true.betas.ia   <- c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
  true.betas = c(true.betas.main, true.betas.ia)
  
  choice_vars <- c('ASC2', 'ASC3', 'cost', 'prot25', 'prot50', 'invasive')
  demographic_vars <- c('male', 'age', 'income', 'child')
  n_alt <- 3
  #lambda_grid <- lambda_grid
  #num_covariates <- num_covariates
  #data <- data
  #all_features <- NULL
  true_betas_full <- NULL
  
  plan(multisession, workers = 45)
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
    
    #For finding false positives and true negatives
    all_features <<- colnames(final_df_scaled)
    true_betas_full <<- setNames(numeric(length(all_features)), all_features)
    true_betas_full[c("ASC2_male","ASC3_age","cost_income",
                      "prot25_child","prot50_child","invasive_child")] <-
                                            c(-0.4,-0.2,0.3,
                                              0.5,0.3,-0.4)
    
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
                                          n_folds = 5)
    #Plot BIC
    # if ("BIC" %in% names(results_cv$lambda_results)) {
    #   bic_data <- results_cv$lambda_results
    #   p_bic <- ggplot(bic_data, aes(x = lambda, y = BIC)) +
    #     geom_line() + geom_point() +
    #     theme_minimal(base_size = 14) +
    #     labs(title = "BIC vs Lambda",
    #          x = "Lambda", y = "BIC")
    #   
    #   bic_filename <- file.path("plots",
    #                             paste0(n_persons, "n_", iterations, "iter_BIC.png"))
    #   ggsave(bic_filename, plot = p_bic, width = 7, height = 5, dpi = 300)
    # }
    
    #Plot mean out of sample LL
    # if ("mean_LL" %in% names(results_cv$lambda_results)) {
    #   ll_data <- results_cv$lambda_results
    #   p_ll <- ggplot(ll_data, aes(x = lambda, y = mean_LL)) +
    #     geom_line() + geom_point() +
    #     theme_minimal(base_size = 14) +
    #     labs(title = "CV Log-Likelihood vs Lambda",
    #          x = "Lambda", y = "Mean Out-of-Sample Log-Likelihood")
    #   
    #   ll_filename <- file.path("plots",
    #                            paste0(n_persons, "n_", iterations, "iter_LL.png"))
    #   ggsave(ll_filename, plot = p_ll, width = 7, height = 5, dpi = 300)
    # }
    
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
    
    coefs <- as.numeric(coef(final_model))
    names(coefs) <- selected_features
    coefs_full <- setNames(numeric(length(all_features)), all_features)
    coefs_full[selected_features] <- coefs
    #return(coefs_full)
    return(list(coefs_full = coefs_full, best_lambda = results_cv$best_lambda))
    
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
                               paste0(n_persons, "n_", iterations, "iter_lambda_values.csv"))
  write.csv(lambda_df, lambda_filename, row.names = FALSE)
  
  results <- do.call(cbind, results_list)
  results_df <- as.data.frame(t(results))
  # #Evaluation metrics for MCMC
  # # active coefficients only
  # true_act <- true_betas_full[active_features]
  # est_act  <- results_df[, active_features, drop = FALSE]
  # 
  # # sign recovery
  # sign_recovery <- colMeans(sign(est_act) == sign(true_act))
  # 
  # # bias
  # bias <- colMeans(est_act) - true_act
  # 
  # # RMSE
  # rmse <- sqrt(colMeans((t(t(est_act) - true_act))^2))
  # 
  # metrics_active <- data.frame(
  #   Feature = active_features,
  #   True    = true_act,
  #   SignRec = sign_recovery,
  #   Bias    = bias,
  #   RMSE    = rmse
  # )
  # 
  # #zero detection rate
  # # irrelevant coefficients = inactive features with true value = 0
  # est_inact <- results_df[, inactive_features, drop = FALSE]
  # 
  # zero_detect_rate <- colMeans(abs(est_inact) < threshold)
  # 
  # metrics_inactive <- data.frame(
  #   Feature = inactive_features,
  #   ZeroDetect = zero_detect_rate
  # )
  # 
  # lambda_mean <- mean(best_lambdas)
  # lambda_sd   <- sd(best_lambdas)
  # lambda_iqr  <- IQR(best_lambdas)
  # 
  # lambda_metrics <- data.frame(
  #   Mean = lambda_mean,
  #   SD   = lambda_sd,
  #   IQR  = lambda_iqr
  # )
  # 
  # # Save active–feature recovery metrics
  # active_filename <- file.path(
  #   "data",
  #   paste0(n_persons, "n_", iterations, "_active_metrics.csv")
  # )
  # write.csv(metrics_active, active_filename, row.names = FALSE)
  # 
  # # Save inactive–feature sparsity metrics
  # inactive_filename <- file.path(
  #   "data",
  #   paste0(n_persons, "n_", iterations, "_inactive_metrics.csv")
  # )
  # write.csv(metrics_inactive, inactive_filename, row.names = FALSE)
  # 
  # # Save lambda stability metrics
  # lambda_metrics_filename <- file.path(
  #   "data",
  #   paste0(n_persons, "n_", iterations, "_lambda_metrics.csv")
  # )
  # write.csv(lambda_metrics, lambda_metrics_filename, row.names = FALSE)
  
  threshold <- 0.05
  # split active vs inactive for diagnostics
  active_features <- names(true_betas_full[true_betas_full != 0])
  inactive_features <- names(true_betas_full[true_betas_full == 0])
  
  tp_rate <- colMeans(abs(results_df[, active_features]) > threshold, na.rm = TRUE)
  fp_rate <- colMeans(abs(results_df[, inactive_features]) > threshold, na.rm = TRUE)
  
  fn_rate <- 1 - tp_rate
  tn_rate <- 1 - fp_rate
  
  cat("\n=== Detection Summary ===\n")
  
  #cat("True Positives (active recovered):\n"); print(tp_rate)
  #cat("\nFalse Negatives (active missed):\n"); print(fn_rate)
  cat("\nFalse Positives (inactive wrongly picked):\n"); print(head(fp_rate, 10))  # print first few
  cat("\nTrue Negatives (inactive correctly zero):\n"); print(head(tn_rate, 10))
  
  detection_summary <- data.frame(
    Feature = unique(c(names(tp_rate), names(fp_rate))),
    #TP_rate = tp_rate[unique(c(names(tp_rate), names(fp_rate)))],
    #FN_rate = fn_rate[unique(c(names(tp_rate), names(fp_rate)))],
    FP_rate = fp_rate[unique(c(names(tp_rate), names(fp_rate)))],
    TN_rate = tn_rate[unique(c(names(tp_rate), names(fp_rate)))]
  )
  
  summary_filename <- file.path("data",
                                paste0(n_persons, "n_", iterations, "iter_detections.csv"))
  write.csv(detection_summary, summary_filename, row.names = FALSE)
  
  results_active <- results_df[, active_features, drop = FALSE]
  results_inactive <- results_df[, inactive_features, drop = FALSE]
  
  interaction_features <- c("ASC2_male", "ASC3_age", "cost_income",
                            "prot25_child", "prot50_child", "invasive_child")
  results_df_interactions <- results_df[, interaction_features, drop = FALSE]
  #Save the values as csv
  csv_filename <- file.path("data", paste0(n_persons, "n_", iterations, "iter_results.csv"))
  write.csv(results_df_interactions, csv_filename, row.names = FALSE)
  
  true_values <- data.frame(
    Feature = c("ASC2_male", "ASC3_age", "cost_income",
                "prot25_child", "prot50_child", "invasive_child"),
    Forced = c(-0.4, -0.2, 0.3, 0.5, 0.3, -0.4)
  )
  
  plot_data <- tidyr::pivot_longer(results_df_interactions,
                                   cols = everything(),
                                   names_to = "Feature",
                                   values_to = "Coefficient")
  
  p <- ggplot(plot_data, aes(x = Feature, y = Coefficient, fill = Feature)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
    geom_point(data = true_values,
               aes(x = Feature, y = Forced),
               color = "red", shape = 95, size = 6, inherit.aes = FALSE) +
    scale_fill_brewer(palette = "Pastel1") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(title = "Recovery of True Interaction Coefficients",
         subtitle = "Red dash = induced interaction value",
         x = "Interaction Feature",
         y = "Coefficient Estimate")
  
  filename <- file.path("plots", paste0(n_persons, "n_", iterations, "iter.png"))
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  
  print(p)
}


# n_persons = number of persons to consider, chosen randomly based on "id"
# num_covariates: how many features to consider

#start_time <- Sys.time()


# force_feature_interaction(n_persons = 100,
#                           #num_covariates = 24,
#                           #lambda_grid = seq(0.001, 0.02, 0.001),
#                           iterations = 2)

for (i in c(500, 600, 700, 800, 900, 973)) {
  force_feature_interaction(
    n_persons = i,
    iterations = 250
  )
}

# end_time <- Sys.time()
# print(as.numeric(end_time - start_time, units = "mins"))
