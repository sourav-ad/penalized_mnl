## INDUCE AND RECOVER INTERACTIONS WITH COEF-THRESHOLD)

required_packages <- c(
  "maxLik","matrixStats","tidyr","dplyr","glmnet","bgw",
  "Rfast","future.apply","future","evd","ggplot2"
)

install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }
  invisible(lapply(packages, require, character.only = TRUE))
}

install_if_missing(required_packages)

source("functions/utility_functions.R")
source("functions/utility_functions_generalized.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")
source("functions/MNL_functions_execution.R")

force_feature_interaction <- function(
    n_persons = 400,
    lambda_grid = exp(seq(log(1e-4), log(5e-1), length.out = 10)),
    iterations = 100,
    alpha_grid = c(0.5, 1),
    n_folds = 3,
    coef_thresh_grid = c(1e-4)   # add more for threshold sweep
) {

  root <- "data"
  out_dir <- file.path(root, paste0("n_", n_persons, "_iter_", iterations))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  #Data
  data <- read.csv("data/doggerbank_full_973_wide.csv")
  set.seed(123)
  selected_ids <- sample(unique(data$id), n_persons)
  data <- data[data$id %in% selected_ids, ]

  data$ASC21 <- 0; data$ASC22 <- 1; data$ASC23 <- 0
  data$ASC31 <- 0; data$ASC32 <- 0; data$ASC33 <- 1

  # Scaling
  data$income <- scale(data$income)
  data$age    <- scale(data$age)
  data$cost1  <- scale(data$cost1)
  data$cost2  <- scale(data$cost2)
  data$cost3  <- scale(data$cost3)

  #no cost interactions
  col_names <- c(
    "ASC2","ASC3",
    "spec10","spec25","prot25","prot50","invasive","cost",
    "ASC2_male","ASC3_age",
    "prot25_child","prot50_child","invasive_child"
  )

  #alt matrices
  alt1 <- cbind(
    0, 0,
    data$spec101, data$spec251, data$prot251, data$prot501, data$invasive1, data$cost1,
    data$ASC21 * data$male,
    data$ASC31 * data$age,
    data$prot251 * data$child,
    data$prot501 * data$child,
    data$invasive1 * data$child
  )

  alt2 <- cbind(
    1, 0,
    data$spec102, data$spec252, data$prot252, data$prot502, data$invasive2, data$cost2,
    data$ASC22 * data$male,
    data$ASC32 * data$age,
    data$prot252 * data$child,
    data$prot502 * data$child,
    data$invasive2 * data$child
  )

  alt3 <- cbind(
    0, 1,
    data$spec103, data$spec253, data$prot253, data$prot503, data$invasive3, data$cost3,
    data$ASC23 * data$male,
    data$ASC33 * data$age,
    data$prot253 * data$child,
    data$prot503 * data$child,
    data$invasive3 * data$child
  )

  colnames(alt1) <- col_names
  colnames(alt2) <- col_names
  colnames(alt3) <- col_names

  true.betas.main <- c(-0.4, -0.4, 0.2, 0.3, 1, 1.3, -0.9, -0.04)

  interaction_names <- c(
    "ASC2_male",
    "ASC3_age",
    "prot25_child",
    "prot50_child",
    "invasive_child"
  )

  interaction_idx <- match(interaction_names, col_names)
  stopifnot(all(!is.na(interaction_idx)))

  p <- ncol(alt1)
  stopifnot(p == length(col_names))

  #placeholder
  true.betas <- rep(0, p)
  true.betas[match(names(true.betas.main), names(true.betas.main))] <- 0  # no-op safety
  true.betas[1:length(true.betas.main)] <- true.betas.main

  alt_list <- list(alt1, alt2, alt3)
  K <- length(interaction_idx)
  n_active <- ceiling(K / 2)

  #alpha loop
  for (a in seq_along(alpha_grid)) {

    alpha <- alpha_grid[a]

    ## alpha-specific storage objects
    coef_mat      <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(coef_mat) <- col_names

    true_beta_mat <- matrix(NA_real_, nrow = iterations, ncol = K)
    colnames(true_beta_mat) <- interaction_names

    best_lambda <- rep(NA_real_, iterations)
    fit_code    <- rep(NA_integer_, iterations)
    fit_loglik  <- rep(NA_real_, iterations)
    sim_seed    <- rep(NA_integer_, iterations)

    # Iterations
    for (iter in seq_len(iterations)) {

      seed_iter <- 123 + iter + 1000 * a
      sim_seed[iter] <- seed_iter
      set.seed(seed_iter)

      ## induce oracle truth
      true_betas_iter <- true.betas
      true_betas_iter[interaction_idx] <- 0

      active_idx <- sample(interaction_idx, size = n_active, replace = FALSE)

      true_betas_iter[active_idx] <- round(
        runif(n_active, 0.2, 0.6) * 2 * sample(c(-1, 1), n_active, replace = TRUE),
        digits = 1
      )

      true_beta_mat[iter, ] <- true_betas_iter[interaction_idx]

      ## simulate choices (gumbel noise)
      g1 <- rgumbel(nrow(alt1))
      g2 <- rgumbel(nrow(alt2))
      g3 <- rgumbel(nrow(alt3))

      util1 <- alt1 %*% true_betas_iter + g1
      util2 <- alt2 %*% true_betas_iter + g2
      util3 <- alt3 %*% true_betas_iter + g3

      choice <- apply(cbind(util1, util2, util3), 1, which.max)
      y1 <- as.integer(choice == 1)
      y2 <- as.integer(choice == 2)
      y3 <- as.integer(choice == 3)
      choice_list <- list(y1, y2, y3)

      #tune lambda
      cv <- tune_lambda_cv_bespoke(
        alt_list     = alt_list,
        choice_list  = choice_list,
        lambda_grid  = lambda_grid,
        person_id    = data$id,
        alpha        = alpha,
        n_folds      = n_folds
      )

      best_lambda[iter] <- cv$best_lambda

      #penalized fit
      fit <- maxLik(
        function(b) MNL(
          b, alt_list, choice_list,
          lambda = best_lambda[iter],
          alpha  = alpha,
          intercept_index = NULL
        ),
        start = rep(0, p),
        method = "BHHH",
        finalHessian = TRUE
      )

      coef_mat[iter, ] <- coef(fit)
      fit_code[iter]   <- tryCatch(as.integer(fit$code), error = function(e) NA_integer_)
      fit_loglik[iter] <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
    } #end iterations loop

    #save alpha-level outputs

    #true induced interactions
    write.csv(
      data.frame(Alpha = alpha, Iteration = seq_len(iterations), true_beta_mat),
      file.path(out_dir, paste0("alpha_", alpha, "_true_interactions.csv")),
      row.names = FALSE
    )

    ## Estimated interaction coefs
    coef_interactions <- coef_mat[, interaction_idx, drop = FALSE]
    colnames(coef_interactions) <- interaction_names

    write.csv(
      data.frame(Alpha = alpha, Iteration = seq_len(iterations), coef_interactions),
      file.path(out_dir, paste0("alpha_", alpha, "_interaction_coefs.csv")),
      row.names = FALSE
    )

    #full coefficient matris
    write.csv(
      data.frame(Alpha = alpha, Iteration = seq_len(iterations), coef_mat),
      file.path(out_dir, paste0("alpha_", alpha, "_coef_full.csv")),
      row.names = FALSE
    )

    #lambdas
    write.csv(
      data.frame(Alpha = alpha, Iteration = seq_len(iterations), BestLambda = best_lambda),
      file.path(out_dir, paste0("alpha_", alpha, "_optimal_lambdas.csv")),
      row.names = FALSE
    )

    #Fit diagnostics
    write.csv(
      data.frame(
        Alpha = alpha,
        Iteration = seq_len(iterations),
        Seed = sim_seed,
        BestLambda = best_lambda,
        FitCode = fit_code,
        LogLik = fit_loglik
      ),
      file.path(out_dir, paste0("alpha_", alpha, "_fit_diagnostics.csv")),
      row.names = FALSE
    )


    #detection + metrics per threshold (no refitting)
    for (th in coef_thresh_grid) {

      detected <- abs(coef_mat[, interaction_idx, drop = FALSE]) >= th

      detected_long <- data.frame(
        Alpha       = alpha,
        Threshold   = th,
        Iteration   = rep(seq_len(iterations), each = K),
        Interaction = rep(interaction_names, times = iterations),
        TrueBeta    = as.vector(t(true_beta_mat)),
        Estimate    = as.vector(t(coef_interactions)),
        Detected    = as.integer(as.vector(t(detected)))
      )

      write.csv(
        detected_long,
        file.path(out_dir, paste0("alpha_", alpha, "_th_", format(th, scientific = TRUE),
                                  "_interaction_detection_long.csv")),
        row.names = FALSE
      )

      ## Confusion matrix + rates (global across iterations & interactions)
      truth <- (true_beta_mat != 0)

      TP <- sum(detected & truth)
      FN <- sum(!detected & truth)
      FP <- sum(detected & !truth)
      TN <- sum(!detected & !truth)

      write.csv(
        data.frame(
          Alpha = alpha,
          Threshold = th,
          TrueStatus = c("Induced", "Induced", "Not Induced", "Not Induced"),
          Predicted  = c("Detected", "Not Detected", "Detected", "Not Detected"),
          Count      = c(TP, FN, FP, TN)
        ),
        file.path(out_dir, paste0("alpha_", alpha, "_th_", format(th, scientific = TRUE),
                                  "_confusion_matrix.csv")),
        row.names = FALSE
      )

      confusion_rates <- data.frame(
        Alpha = alpha,
        Threshold = th,
        TPR = if ((TP + FN) > 0) TP / (TP + FN) else NA_real_,  # power = 1 - FNR
        FPR = if ((FP + TN) > 0) FP / (FP + TN) else NA_real_,
        TNR = if ((TN + FP) > 0) TN / (TN + FP) else NA_real_,  # 1 - FPR
        FNR = if ((FN + TP) > 0) FN / (FN + TP) else NA_real_
      )

      write.csv(
        confusion_rates,
        file.path(out_dir, paste0("alpha_", alpha, "_th_", format(th, scientific = TRUE),
                                  "_confusion_rates.csv")),
        row.names = FALSE
      )
    } #end coeff threshold loop
  } #end alpha loop

  invisible(TRUE)
}

#Parallel
plan(multisession, workers = 55)
options(future.rng.onMisuse = "ignore")

n_grid <- c(50)

future_lapply(n_grid, function(n) {
  force_feature_interaction(
    n_persons = n,
    iterations = 2,
    alpha_grid = c(0.5, 1),
    lambda_grid = exp(seq(log(1e-4), log(5e-1), length.out = 10)),
    n_folds = 3,
    coef_thresh_grid = c(1e-4)  #extend as needed
  )
})

plan(sequential)
