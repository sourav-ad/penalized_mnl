##   Semi-synthetic analysis for penalized MNL with:
##     (1) respondent-level sample split (selection vs inference)
##     (2) penalized fit on selection split (for screening)
##     (3) screening rule (top-k among interactions; scale-free)
##     (4) unpenalized refit on inference split (for valid Wald p-values)
##     (5) complete bookkeeping: save everything needed to reproduce
##         all metrics/plots without rerunning.
##
## Output:
##   Write CSVs + RDS per (n, alpha):
##     data/<n>n_<iterations>_alpha_<alpha>_ALL_OUTPUTS.rds

#libraries
required_packages <- c(
  "maxLik", "matrixStats", "tidyr", "dplyr", "glmnet", "bgw",
  "Rfast", "future.apply", "future", "evd", "ggplot2"
)

install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]

  if (length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)
  }

  invisible(lapply(packages, require, character.only = TRUE))
}

install_if_missing(required_packages)

#source files
source("functions/utility_functions.R")
source("functions/utility_functions_generalized.R") # has create_alt_matrices2()
source("functions/mnl_function.R")
source("functions/pre_process.R")
source("functions/MNL_functions_execution.R")

#main function: screen + refit + bookkeeping
force_feature_interaction <- function(
    n_persons     = 400,
    lambda_grid   = exp(seq(log(1e-4), log(5e-1), length.out = 10)),
    iterations    = 100,
    alpha_grid    = c(0.5, 1),
    n_folds       = 3,
    alpha_sig     = 0.01,   # detection threshold (on refitted p-values)
    k_select      = NULL,   # if NULL, defaults to K/2 (induced signals)
    out_dir       = "data",
    save_cv_objs  = TRUE    # stores full cv objects in RDS (can be large)
) {

  #output directory
  out_dir <- file.path(out_dir, 
                       paste0("n_", n_persons, "_iter_", iterations))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  #Data 
  data <- read.csv("data/doggerbank_full_973_wide.csv")
  selected_ids <- sample(unique(data$id), n_persons)
  data <- data[data$id %in% selected_ids, ]

  # Alternatives & ASCs 
  data$ASC21 <- 0; data$ASC22 <- 1; data$ASC23 <- 0
  data$ASC31 <- 0; data$ASC32 <- 0; data$ASC33 <- 1

  #scaling 
  data$income <- scale(data$income)
  data$age    <- scale(data$age)
  data$cost2  <- scale(data$cost2)
  data$cost3  <- scale(data$cost3)
  #cost1 is all 0s so interaction detection will be poor

  #alt matrices
  alt1 <- cbind(
    0, 0,
    data$spec101, data$spec251, data$prot251, data$prot501, data$invasive1, data$cost1,
    data$ASC21 * data$male, data$ASC31 * data$age, data$cost1 * data$income,
    data$prot251 * data$child, data$prot501 * data$child, data$invasive1 * data$child
  )

  alt2 <- cbind(
    1, 0,
    data$spec102, data$spec252, data$prot252, data$prot502, data$invasive2, data$cost2,
    data$ASC22 * data$male, data$ASC32 * data$age, data$cost2 * data$income,
    data$prot252 * data$child, data$prot502 * data$child, data$invasive2 * data$child
  )

  alt3 <- cbind(
    0, 1,
    data$spec103, data$spec253, data$prot253, data$prot503, data$invasive3, data$cost3,
    data$ASC23 * data$male, data$ASC33 * data$age, data$cost3 * data$income,
    data$prot253 * data$child, data$prot503 * data$child, data$invasive3 * data$child
  )

  #induced betas
  true.betas.main <- c(-0.4, -0.4, 0.2, 0.3, 1, 1.3, -0.9, -0.04)
  true.betas.ia   <- rep(0, 6) # placeholder
  true.betas      <- c(true.betas.main, true.betas.ia)

  p <- ncol(alt1)
  stopifnot(p == length(true.betas))

  #interaction bookkeeping 
  interaction_idx <- (p - 5):p
  K <- length(interaction_idx)
  stopifnot(K %% 2 == 0)

  interaction_names <- c(
    "ASC2_male",
    "ASC3_age",
    "cost_income",
    "prot25_child",
    "prot50_child",
    "invasive_child"
  )

  #top-k selection, scale-free
  if (is.null(k_select)) k_select <- K / 2

  ## Lists for alt matrices
  alt_list_full_template <- list(alt1, alt2, alt3)

  #loop over alpha grid
  for (a in seq_along(alpha_grid)) {

    alpha <- alpha_grid[a]

    #alpha (0.5, 1) specific storage

    #truth
    true_beta_mat <- matrix(NA_real_, nrow = iterations, ncol = K)
    colnames(true_beta_mat) <- interaction_names

    #tuning
    best_lambda <- numeric(iterations)
    
    ##book keeping

    #penalized fit outputs
    coef_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(coef_mat) <- paste0("b", seq_len(p))

    pval_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(pval_mat) <- paste0("b", seq_len(p))

    #penalized se/z, Wald type
    pen_se_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(pen_se_mat) <- paste0("b", seq_len(p))

    #p vals
    pen_z_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(pen_z_mat) <- paste0("b", seq_len(p))

    #selection tracking over interactions
    selected_mat <- matrix(FALSE, nrow = iterations, ncol = K)
    colnames(selected_mat) <- interaction_names

    #refit coefficients (NA for dropped)
    refit_coef_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(refit_coef_mat) <- paste0("b", seq_len(p))

    refit_se_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(refit_se_mat) <- paste0("b", seq_len(p))

    refit_z_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(refit_z_mat) <- paste0("b", seq_len(p))

    refit_pval_mat <- matrix(NA_real_, nrow = iterations, ncol = p)
    colnames(refit_pval_mat) <- paste0("b", seq_len(p))

    #post refit, p-val criteria on selected interactions
    #dropped implies not detected
    detected_refit_mat <- matrix(FALSE, nrow = iterations, ncol = K)
    colnames(detected_refit_mat) <- interaction_names

    #seeds + splits
    sim_seed_vec <- integer(iterations)
    split_seed_vec <- integer(iterations)

    #exactly which rows selected (iterations x nrow(data))
    #allows reconstructing which rows were used where
    sel_mask_mat <- matrix(FALSE, nrow = iterations, ncol = nrow(data))

    #store selected/kept indices (every data kept)
    selected_int_idx_list <- vector("list", iterations)  
    #Selected Interaction indices in [1..p]
    keep_idx_list <- vector("list", iterations)  
    #full set of column indices used in the unpenalized refit

    #fit diagnostics
    pen_loglik <- rep(NA_real_, iterations)
    refit_loglik <- rep(NA_real_, iterations)
    pen_code  <- integer(iterations)
    refit_code <- integer(iterations)

    #CV objects (complete details)
    #output of tune_lambda_cv_bespoke()
    cv_obj_list <- vector("list", iterations)

    #iteration loop
    for (iter in seq_len(iterations)) {

      #seeds
      sim_seed <- 123 + iter + 1000 * a
      sim_seed_vec[iter] <- sim_seed
      set.seed(sim_seed)

      
      #induce oracle truth 
      true_betas_iter <- true.betas

      #50% signal, 50% noise
      active_idx <- sample(interaction_idx, size = K / 2, replace = FALSE)
      true_betas_iter[interaction_idx] <- 0

      true_betas_iter[active_idx] <-
        round(
          runif(K / 2, 0.2, 0.6) * 2 * sample(c(-1, 1), K / 2, replace = TRUE),
          digits = 1
        )

      true_beta_mat[iter, ] <- true_betas_iter[interaction_idx]

      #gumbel errors
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

      alt_list_full    <- alt_list_full_template
      choice_list_full <- list(y1, y2, y3)

      #to split for selection and unpenalized refit
      #50-50
      split_seed <- 999 + iter + 1000 * a
      split_seed_vec[iter] <- split_seed
      set.seed(split_seed)

      #by person id, no leakage between two
      ids <- unique(data$id)
      sel_ids <- sample(ids, size = floor(length(ids) / 2), replace = FALSE)

      sel_mask <- data$id %in% sel_ids #selection mask
      inf_mask <- !sel_mask #inference mask

      sel_mask_mat[iter, ] <- sel_mask  #all row-level records
      #TRUE and FALSE values for every iteration, per row

      alt_list_sel <- lapply(alt_list_full, function(A) A[sel_mask, , drop = FALSE])
      alt_list_inf <- lapply(alt_list_full, function(A) A[inf_mask, , drop = FALSE])

      #Choices from selection group
      choice_list_sel <- lapply(choice_list_full, function(y) y[sel_mask])
      #Choices from inference group
      choice_list_inf <- lapply(choice_list_full, function(y) y[inf_mask])

      #tune lambda (on selection split, bespoke function)
      cv <- tune_lambda_cv_bespoke(
        alt_list     = alt_list_sel,
        choice_list  = choice_list_sel,
        lambda_grid  = lambda_grid,
        person_id    = data$id[sel_mask],
        alpha        = alpha,
        n_folds      = n_folds
      )

      best_lambda[iter] <- cv$best_lambda
      if (isTRUE(save_cv_objs)) cv_obj_list[[iter]] <- cv

      #penalized MNL fit 
      fit_pen <- maxLik(
        function(b) MNL(
          b, 
          alt_list_sel, 
          choice_list_sel,
          lambda = best_lambda[iter],
          alpha  = alpha,
          intercept_index = NULL
        ),
        start = rep(0, p),
        method = "BHHH", #BHHH with lasso may or may not produce exact 0s
        finalHessian = TRUE
      )

      b_pen <- coef(fit_pen)
      coef_mat[iter, ] <- b_pen

      #penalized pseudo Wald p vals (for heuristics only)
      se_pen <- sqrt(diag(vcov(fit_pen)))
      z_pen  <- b_pen / se_pen
      p_pen  <- 2 * (1 - pnorm(abs(z_pen)))

      pen_se_mat[iter, ] <- se_pen
      pen_z_mat[iter, ]  <- z_pen
      pval_mat[iter, ]   <- p_pen

      #diagnostics, failsafe
      pen_loglik[iter] <- tryCatch(as.numeric(logLik(fit_pen)), error = function(e) NA_real_)
      pen_code[iter]   <- tryCatch(as.integer(fit_pen$code), error = function(e) NA_integer_)

      #screening covariates (top-k among interactions, scale-free)
      ord <- order(abs(b_pen[interaction_idx]), decreasing = TRUE)
      S_int <- interaction_idx[ord[1:k_select]]  # selected interaction indices in [1..p]

      main_idx <- setdiff(seq_len(p), interaction_idx)
      keep_idx <- c(main_idx, S_int)

      selected_int_idx_list[[iter]] <- S_int
      keep_idx_list[[iter]] <- keep_idx

      sel_vec <- rep(FALSE, K)
      sel_vec[match(S_int, interaction_idx)] <- TRUE
      selected_mat[iter, ] <- sel_vec

      #unpenalized refit on inference split (reduced model)
      #filter data from screening, for inference
      #alt matrices size reduces if interactions dropped
      alt_list_inf_red <- lapply(alt_list_inf, function(A) A[, keep_idx, drop = FALSE])

      #no alpha, no lambda
      #unpenalized fit
      #meaningful p vals
      fit_refit <- maxLik(
        function(b) MNL_unpenalized(
          b, #coefficients vec to estimate 
          alt_list_inf_red, 
          choice_list_inf,
          final_eval = FALSE,
          nrep = 6
        ),
        start = b_pen[keep_idx],   #warm start
        #from values already found in penalized fit 
        method = "BHHH",
        finalHessian = TRUE
      )

      #refit p values
      b_refit  <- coef(fit_refit)
      se_refit <- sqrt(diag(vcov(fit_refit)))
      z_refit  <- b_refit / se_refit
      p_refit  <- 2 * (1 - pnorm(abs(z_refit)))

      #map refit outputs to a full length template (NA for dropped)
      #size = covariates considered, always
      refit_coef_full <- rep(NA_real_, p)
      refit_se_full   <- rep(NA_real_, p)
      refit_z_full    <- rep(NA_real_, p)
      refit_p_full    <- rep(NA_real_, p)

      refit_coef_full[keep_idx] <- b_refit
      refit_se_full[keep_idx]   <- se_refit
      refit_z_full[keep_idx]    <- z_refit
      refit_p_full[keep_idx]    <- p_refit

      refit_coef_mat[iter, ] <- refit_coef_full
      refit_se_mat[iter, ]   <- refit_se_full
      refit_z_mat[iter, ]    <- refit_z_full
      refit_pval_mat[iter, ] <- refit_p_full

      #diagnostics, failsafe
      refit_loglik[iter] <- tryCatch(as.numeric(logLik(fit_refit)), error = function(e) NA_real_)
      refit_code[iter]   <- tryCatch(as.integer(fit_refit$code), error = function(e) NA_integer_)

      #refit p-values on selected
      detected_vec <- rep(FALSE, K)

      #positions of selected interactions within the reduced parameter vector
      pos_keep <- match(S_int, keep_idx)
      p_sel <- p_refit[pos_keep]

      #dropped implies not detected (FALSE)
      #selected implies p_val < alpha_sig
      detected_vec[sel_vec] <- (p_sel < alpha_sig)
      detected_refit_mat[iter, ] <- detected_vec
    } # end iterations loop

    #save CSV files

    #penalized interaction coefficients (baseline comparator)
    coef_interactions <- coef_mat[, interaction_idx, drop = FALSE]
    colnames(coef_interactions) <- interaction_names

    write.csv(
      cbind(Alpha = alpha, Iteration = seq_len(iterations), coef_interactions),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_coefs_penalized.csv")),
      row.names = FALSE
    )

    #refit interaction coefficients (primary for post-selection inference; NA if dropped)
    refit_interactions <- refit_coef_mat[, interaction_idx, drop = FALSE]
    colnames(refit_interactions) <- interaction_names

    write.csv(
      cbind(Alpha = alpha, Iteration = seq_len(iterations), refit_interactions),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_coefs_refit.csv")),
      row.names = FALSE
    )

    #Lambdas
    write.csv(
      data.frame(Alpha = alpha, Iteration = seq_len(iterations), BestLambda = best_lambda),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_optimal_lambdas.csv")),
      row.names = FALSE
    )

    #truth (induced betas)
    write.csv(
      cbind(Alpha = alpha, Iteration = seq_len(iterations), true_beta_mat),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_true_beta.csv")),
      row.names = FALSE
    )

    ## Selection df, screening data
    selected_long <- data.frame(
      Alpha       = alpha,
      Iteration   = rep(seq_len(iterations), each = K),
      Interaction = rep(interaction_names, times = iterations),
      TrueBeta    = as.vector(t(true_beta_mat)),
      Selected    = as.integer(as.vector(t(selected_mat)))
    )

    write.csv(
      selected_long,
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_selected_long.csv")),
      row.names = FALSE
    )

    ## Detection summary, using unpenalized refit p-val
    detected_long_refit <- data.frame(
      Alpha       = alpha,
      Iteration   = rep(seq_len(iterations), each = K),
      Interaction = rep(interaction_names, times = iterations),
      TrueBeta    = as.vector(t(true_beta_mat)),
      Selected    = as.integer(as.vector(t(selected_mat))),
      Detected    = as.integer(as.vector(t(detected_refit_mat)))
    )

    write.csv(
      detected_long_refit,
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_detection_refit_long.csv")),
      row.names = FALSE
    )

    #penalized heuristic pval (just in case)
    detected_pen <- (pval_mat[, interaction_idx, drop = FALSE] < alpha_sig)
    detected_long_pen <- data.frame(
      Alpha       = alpha,
      Iteration   = rep(seq_len(iterations), each = K),
      Interaction = rep(interaction_names, times = iterations),
      TrueBeta    = as.vector(t(true_beta_mat)),
      Detected    = as.integer(as.vector(t(detected_pen)))
    )

    write.csv(
      detected_long_pen,
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_interaction_detection_penalized_long.csv")),
      row.names = FALSE
    )

    #confusion matrix (from unpenalized refit p vals)
    truth <- (true_beta_mat != 0)
    detected <- detected_refit_mat

    TP <- sum(detected & truth)
    FN <- sum(!detected & truth)
    FP <- sum(detected & !truth)
    TN <- sum(!detected & !truth)

    write.csv(
      data.frame(
        Alpha = alpha,
        TrueStatus = c("Induced", "Induced", "Not Induced", "Not Induced"),
        Predicted  = c("Detected", "Not Detected", "Detected", "Not Detected"),
        Count      = c(TP, FN, FP, TN)
      ),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_confusion_matrix_refit.csv")),
      row.names = FALSE
    )

    #confusion rates
    confusion_rates <- data.frame(
      Alpha = alpha,
      TPR   = TP /(TP + FN),
      FPR   = FP /(FP + TN),
      TNR   = TN /(TN + FP),
      FNR   = FN /(FN + TP)
    )

    write.csv(
      confusion_rates,
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_confusion_rates_refit.csv")),
      row.names = FALSE
    )

    #every diagnostic
    write.csv(
      data.frame(
        Alpha = alpha,
        Iteration = seq_len(iterations),
        SimSeed = sim_seed_vec,
        SplitSeed = split_seed_vec,
        BestLambda = best_lambda,
        PenLogLik = pen_loglik,
        PenCode = pen_code,
        RefitLogLik = refit_loglik,
        RefitCode = refit_code
      ),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_run_diagnostics.csv")),
      row.names = FALSE
    )

    #more bookkeeping

    #rows, if in selection bunch or inference bunch
    saveRDS(
      sel_mask_mat,
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_sel_mask_mat.rds"))
    )

    #selected/kept indices lists
    #which columns were kept in every iteration
    saveRDS(
      list(selected_int_idx_list = selected_int_idx_list,
           keep_idx_list = keep_idx_list),
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_screening_indices.rds"))
    )

    #Full penalized and refit inference objects
    #comment/uncomment as needed to save CSVs
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), pen_se_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_pen_se.csv")),
              row.names=FALSE)
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), pen_z_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_pen_z.csv")),
              row.names=FALSE)
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), pval_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_pen_p.csv")),
              row.names=FALSE)
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), refit_coef_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_refit_coef_full.csv")),
              row.names=FALSE)
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), refit_se_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_refit_se_full.csv")),
              row.names=FALSE)
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), refit_z_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_refit_z_full.csv")),
              row.names=FALSE)
    write.csv(cbind(Alpha=alpha, Iteration=seq_len(iterations), refit_pval_mat),
              file.path(out_dir, paste0(n_persons,"n_",iterations,"_alpha_",alpha,"_refit_p_full.csv")),
              row.names=FALSE)

    # one big bundle per (n, alpha), all the information
    all_out <- list(
      meta = list(
        n_persons   = n_persons,
        iterations  = iterations,
        alpha       = alpha,
        alpha_sig   = alpha_sig,
        n_folds     = n_folds,
        lambda_grid = lambda_grid,
        k_select    = k_select,
        interaction_idx = interaction_idx,
        interaction_names = interaction_names,
        p = p
      ),
      seeds = list(sim_seed = sim_seed_vec, split_seed = split_seed_vec),
      truth = true_beta_mat,
      tuning = list(best_lambda = best_lambda, cv = if (isTRUE(save_cv_objs)) cv_obj_list else NULL),
      penalized = list(
        coef = coef_mat,
        se   = pen_se_mat,
        z    = pen_z_mat,
        p    = pval_mat,
        logLik = pen_loglik,
        code   = pen_code
      ),
      screening = list(
        selected_mat = selected_mat,
        selected_int_idx_list = selected_int_idx_list,
        keep_idx_list = keep_idx_list,
        sel_mask_mat = sel_mask_mat
      ),
      refit = list(
        coef_full = refit_coef_mat,
        se_full   = refit_se_mat,
        z_full    = refit_z_mat,
        p_full    = refit_pval_mat,
        detected_refit_mat = detected_refit_mat,
        logLik = refit_loglik,
        code   = refit_code
      )
    )

    saveRDS(
      all_out,
      file.path(out_dir, paste0(n_persons, "n_", iterations, "_alpha_", alpha, "_ALL_OUTPUTS.rds"))
    )
  } #end alpha loop

  invisible(TRUE)
}

#parallel backend
plan(multisession, workers = 55)
options(future.rng.onMisuse = "ignore")

#n grid
n_grid <- c(50)

## Run in parallel
future_lapply(n_grid, function(n) {
  force_feature_interaction(
    n_persons   = n,
    iterations  = 2, #iterations x interactions, total elements in confusion matrix
    alpha_grid  = c(0.5, 1),
    lambda_grid = exp(seq(log(1e-4), log(5e-1), length.out = 10)),
    n_folds     = 3,
    alpha_sig   = 0.001,
    k_select    = NULL,      # default K/2
    out_dir     = "data",
    save_cv_objs = TRUE
  )
})

#default back to sequential
plan(sequential)



