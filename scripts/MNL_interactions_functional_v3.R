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
force_feature_interaction <- function(
    n_persons = 400,
    lambda_grid = exp(seq(log(1e-4), log(5e-1), length.out = 10)),
    iterations = 5,
    alpha_grid = c(0.5, 1)
) {

  root <- "data/v3_1"
  out_dir  <- file.path(root,
                        paste0("n_", n_persons, "_iter_", iterations))
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # data
  data <- read.csv("data/doggerbank_full_973_wide.csv")
  set.seed(123)
  selected_ids <- sample(unique(data$id), n_persons)
  data <- data[data$id %in% selected_ids, ]

  # alternatives
  data$ASC21 <- 0; data$ASC22 <- 1; data$ASC23 <- 0
  data$ASC31 <- 0; data$ASC32 <- 0; data$ASC33 <- 1
  
  #scaling
  data$income <- scale(data$income)
  data$age <- scale(data$age)

  alt1 = cbind(0, 0, data$spec101, data$spec251, data$prot251, data$prot501, data$invasive1,
               data$ASC21 * data$male, data$ASC21 * data$age, data$ASC21 * data$income, data$ASC21 * data$child,
               data$ASC31 * data$male, data$ASC31 * data$age, data$ASC31 * data$income, data$ASC31 * data$child,
               data$prot251 * data$male, data$prot251 * data$age, data$prot251 * data$income, data$prot251 * data$child,
               data$prot501 * data$male, data$prot501 * data$age, data$prot501 * data$income, data$prot501 * data$child,
               data$invasive1 * data$male, data$invasive1 * data$age, data$invasive1 * data$income, data$invasive1 * data$child
  )
  alt2 = cbind(1, 0, data$spec102, data$spec252, data$prot252, data$prot502, data$invasive2,
               data$ASC22 * data$male, data$ASC22 * data$age, data$ASC22 * data$income, data$ASC22 * data$child,
               data$ASC32 * data$male, data$ASC32 * data$age, data$ASC32 * data$income, data$ASC32 * data$child,
               data$prot252 * data$male, data$prot252 * data$age, data$prot252 * data$income, data$prot252 * data$child,
               data$prot502 * data$male, data$prot502 * data$age, data$prot502 * data$income, data$prot502 * data$child,
               data$invasive2 * data$male, data$invasive2 * data$age, data$invasive2 * data$income, data$invasive2 * data$child
  )
  alt3 = cbind(0, 1, data$spec103, data$spec253, data$prot253, data$prot503, data$invasive3,
               data$ASC23 * data$male, data$ASC23 * data$age, data$ASC23 * data$income, data$ASC23 * data$child,
               data$ASC33 * data$male, data$ASC33 * data$age, data$ASC33 * data$income, data$ASC33 * data$child,
               data$prot253 * data$male, data$prot253 * data$age, data$prot253 * data$income, data$prot253 * data$child,
               data$prot503 * data$male, data$prot503 * data$age, data$prot503 * data$income, data$prot503 * data$child,
               data$invasive3 * data$male, data$invasive3 * data$age, data$invasive3 * data$income, data$invasive3 * data$child
               )

  true.betas.main <- c(-0.4, -0.4, 0.2, 0.3, 1, 1.3, -0.9) #must match with number of marginals
  true.betas.ia <- rep(0, 20) #placeholder
  true.betas <- c(true.betas.main, true.betas.ia)

  p <- ncol(alt1)
  stopifnot(p == length(true.betas))

  coef_mat <- matrix(NA, nrow = iterations, ncol = p)
  colnames(coef_mat) <- paste0("b", seq_len(p))

  # pval_mat <- matrix(NA, nrow = iterations, ncol = p)
  # colnames(pval_mat) <- colnames(coef_mat)

  K <- 20
  interaction_idx <- (p - K + 1):p
  stopifnot(length(interaction_idx) == K)
  stopifnot(K %% 2 == 0)   # failsafe, must be even
  
  #storing randomized induced beta
  true_beta_mat <- matrix(NA, nrow = iterations, ncol = K)
  
  interaction_names <- c(
    "ASC2_male", "ASC2_age", "ASC2_income", "ASC2_child",
    "ASC3_male", "ASC3_age", "ASC3_income", "ASC3_child",
    "prot25_male", "prot25_age", "prot25_income", "prot25_child",
    "prot50_male", "prot50_age", "prot50_income", "prot50_child",
    "invasive_male", "invasive_age", "invasive_income", "invasive_child"
  )
  
  colnames(true_beta_mat) <- interaction_names
  #to store lambdas
  best_lambda <- numeric(iterations)
  
  #coef_thresh <- 1e-4
  
  #iterate over alpha
  for (a in seq_along(alpha_grid)) {
    
    ## initialize alpha-specific aspects
    alpha <- alpha_grid[a]
    #alpha_sig <- 0.01 #significance not needed anymore
    
    ## REINITIALIZE alpha-specific storage
    coef_mat      <- matrix(NA, iterations, p)
    #pval_mat      <- matrix(NA, iterations, p)
    best_lambda   <- numeric(iterations)
    true_beta_mat <- matrix(NA, iterations, K)
    colnames(true_beta_mat) <- interaction_names

  for (iter in seq_len(iterations)) {
    
    set.seed(123 + iter + 1000 * a)
    
    ## induce oracle truth
    true_betas_iter <- true.betas
    active_idx <- sample(interaction_idx, size = K / 2, replace = FALSE)
    true_betas_iter[interaction_idx] <- 0
    
    #randomized signal induction
    true_betas_iter[active_idx] <-
      round(
        #runif(K / 2, 0.2, 0.6) * 2 * #scale by 2x
        runif(K / 2, 0.1, 1.2) * #broad signals
          sample(c(-1, 1), K / 2, replace = TRUE), 
        digits = 1
      )
    
    true_beta_mat[iter, ] <- true_betas_iter[interaction_idx]

    ## simulate choices
    #Gumbel noise
    g1 <- rgumbel(nrow(alt1))
    g2 <- rgumbel(nrow(alt2))
    g3 <- rgumbel(nrow(alt3))

    #utility functions
    util1 <- alt1 %*% true_betas_iter + g1
    util2 <- alt2 %*% true_betas_iter + g2
    util3 <- alt3 %*% true_betas_iter + g3

    choice <- apply(cbind(util1, util2, util3), 1, which.max)

    #choice indicators
    y1 <- as.integer(choice == 1)
    y2 <- as.integer(choice == 2)
    y3 <- as.integer(choice == 3)

    alt_list    <- list(alt1, alt2, alt3)
    choice_list <- list(y1, y2, y3)

    # tune lambda using alpha
    #lambda tuning, bespoke function at ~/MNL_functions_execution.R
    cv <- tune_lambda_cv_bespoke(
      alt_list     = alt_list,
      choice_list  = choice_list,
      lambda_grid  = lambda_grid,
      person_id    = data$id,
      alpha        = alpha,
      n_folds      = 3,
      th = 0.1
    )
    
    best_lambda[iter] <- cv$best_lambda
    #best_lambda <- rep(0, iterations)
    

    #estimation penalized
    # fit model
    fit <- maxLik(
      function(b) MNL(
        b, alt_list, choice_list,
        lambda = best_lambda[iter],
        alpha = alpha,
        intercept_index = NULL
      ),
      start = rep(0, p),
      method = "BHHH",
      finalHessian = TRUE
    )
    
    # #estimation unpenalized
    # fit <- maxLik(
    #   function(b) MNL_unpenalized(
    #     b, 
    #     alt_list, 
    #     choice_list
    #   ),
    #   start = rep(0, p),
    #   method = "BHHH",
    #   finalHessian = TRUE
    # )

    ## store coef_mat, pval_mat, best_lambda
    coef_mat[iter, ] <- coef(fit)
    
    #fit diagnostics
    if (iter == 1) {
      fit_code  <- rep(NA_integer_, iterations)
      fit_ll    <- rep(NA_real_, iterations)
    }
    
    fit_code[iter] <- tryCatch(as.integer(fit$code), error = function(e) NA_integer_)
    fit_ll[iter]   <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
    
  } #end iterations loop
    
    #save fit diagonstics
    write.csv(
      data.frame(Alpha=alpha, Iteration=seq_len(iterations),
                 BestLambda=best_lambda, FitCode=fit_code, LogLik=fit_ll),
      paste0(out_dir, "/", n_persons, "n_", iterations, "_alpha_", alpha, "_fit_diagnostics.csv"),
      row.names = FALSE
    )
    
    ##coefficients
    coef_interactions <- coef_mat[, interaction_idx, drop = FALSE]
    colnames(coef_interactions) <- interaction_names
    
    write.csv(
      cbind(Alpha = alpha, Iteration = seq_len(iterations), true_beta_mat),
      paste0(out_dir, "/", n_persons, "n_", iterations, "_alpha_", alpha, "_true_interactions.csv"),
      row.names = FALSE
    )
    
    write.csv(
      cbind(Alpha = alpha, Iteration = seq_len(iterations), coef_interactions),
      paste0(out_dir, "/", n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_coefs.csv"),
      row.names = FALSE
    )
    
    ## lambdas
    write.csv(
      data.frame(Alpha = alpha, Iteration = seq_len(iterations), BestLambda = best_lambda),
      paste0(out_dir, "/", n_persons, "n_", iterations, "_alpha_", alpha, "_optimal_lambdas.csv"),
      row.names = FALSE
    )
    
    
    #th_grid <- c(1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5) #threshold grid
    th_grid <- c(0) #turn off thresholding
    
    detected_long_all <- vector("list", length(th_grid))
    rates_all <- vector("list", length(th_grid))
    cm_all <- vector("list", length(th_grid))
    
    truth <- (true_beta_mat != 0) #separate signal from noise
    
    for (i in seq_along(th_grid)) { #no refit
      th <- th_grid[i]
      
      detected <- abs(coef_interactions) >= th
      
      detected_long_all[[i]] <- data.frame(
        Alpha       = alpha,
        Threshold   = th,
        Iteration   = rep(seq_len(iterations), each = K),
        Interaction = rep(interaction_names, times = iterations),
        TrueBeta    = as.vector(t(true_beta_mat)),
        Estimate    = as.vector(t(coef_interactions)),
        Detected    = as.integer(as.vector(t(detected)))
      )
      
      TP <- sum(detected & truth)
      FN <- sum(!detected & truth)
      FP <- sum(detected & !truth)
      TN <- sum(!detected & !truth)
      
      rates_all[[i]] <- data.frame(
        Alpha     = alpha,
        Threshold = th,
        TP = TP, FN = FN, FP = FP, TN = TN,
        TPR = TP / (TP + FN),
        FPR = FP / (FP + TN),
        TNR = TN / (TN + FP),
        FNR = FN / (FN + TP),
        OneMinusFPR = 1 - (FP / (FP + TN))
      )
      
      cm_all[[i]] <- data.frame(
        Alpha = alpha,
        Threshold = th,
        TrueStatus = c("Induced", "Induced", "Not Induced", "Not Induced"),
        Predicted  = c("Detected", "Not Detected", "Detected", "Not Detected"),
        Count      = c(TP, FN, FP, TN)
      )
    } #end threshold loop
    
    detected_long_df <- do.call(rbind, detected_long_all)
    rates_df <- do.call(rbind, rates_all)
    cm_df <- do.call(rbind, cm_all)
    
    write.csv(
      detected_long_df,
      file.path(out_dir,
                paste0(n_persons, "n_", iterations,
                       "_alpha_", alpha, "_interaction_detection_long.csv")),
      row.names = FALSE
    )
    
    write.csv(
      rates_df,
      file.path(out_dir,
                paste0(n_persons, "n_", iterations,
                       "_alpha_", alpha, "_confusion_rates.csv")),
      row.names = FALSE
    )
    
    write.csv(
      cm_df,
      file.path(out_dir,
                paste0(n_persons, "n_", iterations,
                       "_alpha_", alpha, "_confusion_matrix.csv")),
      row.names = FALSE
    )
    
  } #end of alpha loop
} #end function 



## parallel backend
plan(multisession, workers = 59)
options(future.rng.onMisuse = "ignore")

# grid for persons
n_grid <- c(973)

## run experiments in parallel
future_lapply(n_grid, function(n) {
  force_feature_interaction(
    n_persons  = n,
    iterations = 100,
    alpha_grid = c(0.5, 1)
  )
})

#shut down workers
plan(sequential)