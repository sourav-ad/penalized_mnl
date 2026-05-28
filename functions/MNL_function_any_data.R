#Functional MNL with Lasso regularized LL

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
source("functions/utility_functions_generalized.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")

#To change from wide format to long format

prepare_mnl_data <- function(
    data,
    id_col,
    task_col,
    choice_col = NULL,
    choice_indicator_cols = NULL,
    alt_specific_vars,
    demographic_vars = NULL,
    n_alt = 3,
    asc_specs = NULL,
    asc_alts = NULL,
    alt_sep = "",
    drop_missing_over = NULL
) {

  stopifnot(is.data.frame(data))

  required_base <- c(id_col, task_col)
  missing_base <- setdiff(required_base, names(data))
  if (length(missing_base) > 0) {
    stop("Missing required ID/task columns: ", paste(missing_base, collapse = ", "))
  }

  if (is.null(choice_col) && is.null(choice_indicator_cols)) {
    stop("Provide either choice_col or choice_indicator_cols.")
  }

  if (!is.null(choice_col) && !(choice_col %in% names(data))) {
    stop("choice_col not found in data: ", choice_col)
  }

  if (!is.null(choice_indicator_cols)) {
    missing_choice_ind <- setdiff(choice_indicator_cols, names(data))
    if (length(missing_choice_ind) > 0) {
      stop("Missing choice indicator columns: ", paste(missing_choice_ind, collapse = ", "))
    }

    if (length(choice_indicator_cols) != n_alt) {
      stop("choice_indicator_cols must have length equal to n_alt.")
    }
  }

  missing_demo <- setdiff(demographic_vars, names(data))
  if (length(missing_demo) > 0) {
    stop("Missing demographic variables: ", paste(missing_demo, collapse = ", "))
  }

  out <- data

  if (!is.null(drop_missing_over)) {
    out <- out[, colMeans(is.na(out)) <= drop_missing_over, drop = FALSE]
  }

  # Standard internal names
  out$id <- out[[id_col]]
  out$line <- out[[task_col]]

  # Build choice, y1...yJ, choice1...choiceJ
  if (!is.null(choice_col)) {

    chosen <- out[[choice_col]]

    if (!all(chosen %in% seq_len(n_alt))) {
      stop("choice_col must contain alternative numbers 1, ..., n_alt.")
    }

    out$choice <- chosen

    for (j in seq_len(n_alt)) {
      out[[paste0("y", j)]] <- as.integer(chosen == j)
      out[[paste0("choice", j)]] <- as.integer(chosen == j)
    }

  } else {

    for (j in seq_len(n_alt)) {
      out[[paste0("y", j)]] <- as.integer(out[[choice_indicator_cols[j]]])
      out[[paste0("choice", j)]] <- as.integer(out[[choice_indicator_cols[j]]])
    }

    y_mat <- as.matrix(out[paste0("y", seq_len(n_alt))])

    if (any(rowSums(y_mat) != 1)) {
      stop("Each row must have exactly one selected alternative.")
    }

    out$choice <- max.col(y_mat, ties.method = "first")
  }

  # Standardize alternative-specific variables to <base><j>
  for (v in alt_specific_vars) {

    raw_cols <- paste0(v, alt_sep, seq_len(n_alt))
    std_cols <- paste0(v, seq_len(n_alt))

    missing_raw <- setdiff(raw_cols, names(out))
    if (length(missing_raw) > 0) {
      stop(
        "Missing alternative-specific columns for variable '", v, "': ",
        paste(missing_raw, collapse = ", ")
      )
    }

    for (j in seq_len(n_alt)) {
      out[[std_cols[j]]] <- out[[raw_cols[j]]]
    }
  }

  # ASCs: either already present or created by user beforehand.
  # Here we simply standardize/retain the ASC bases supplied.
  if (!is.null(asc_specs)) {
    for (v in asc_specs) {
      raw_cols <- paste0(v, alt_sep, seq_len(n_alt))
      std_cols <- paste0(v, seq_len(n_alt))

      missing_raw <- setdiff(raw_cols, names(out))
      if (length(missing_raw) > 0) {
        stop(
          "Missing ASC columns for '", v, "': ",
          paste(missing_raw, collapse = ", "),
          ". Create them first or remove this ASC from asc_specs."
        )
      }

      for (j in seq_len(n_alt)) {
        out[[std_cols[j]]] <- out[[raw_cols[j]]]
      }
    }
  }

  choice_vars <- c(asc_specs, alt_specific_vars)

  required_cols <- unique(c(
    "id", "line",
    "choice",
    paste0("y", seq_len(n_alt)),
    paste0("choice", seq_len(n_alt)),
    unlist(lapply(choice_vars, function(v) paste0(v, seq_len(n_alt)))),
    demographic_vars
  ))

  missing_required <- setdiff(required_cols, names(out))
  if (length(missing_required) > 0) {
    stop("Missing required standardized columns: ",
         paste(missing_required, collapse = ", "))
  }

  out <- out[, required_cols, drop = FALSE]

  return(list(
    data = out,
    choice_vars = choice_vars,
    demographic_vars = demographic_vars,
    n_alt = n_alt
  ))
}

data_wide_to_long <- function(
    data,
    choice_vars,
    demographic_vars,
    n_alt = 3
) {

  y_cols <- paste0("y", seq_len(n_alt))
  choice_cols <- paste0("choice", seq_len(n_alt))

  required_alt_cols <- unlist(
    lapply(choice_vars, function(v) paste0(v, seq_len(n_alt)))
  )

  required_cols <- unique(c(
    "id", "line",
    "choice",
    y_cols,
    choice_cols,
    required_alt_cols,
    demographic_vars
  ))

  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns for model-ready data: ",
         paste(missing_cols, collapse = ", "))
  }

  df_demo <- data[, required_cols, drop = FALSE]

  long_list <- lapply(seq_len(n_alt), function(j) {

    tmp <- df_demo[, c("id", "line", "choice", demographic_vars), drop = FALSE]

    tmp$choice_option <- j
    tmp$chosen <- df_demo[[paste0("choice", j)]]

    for (v in choice_vars) {
      tmp[[v]] <- df_demo[[paste0(v, j)]]
    }

    tmp
  })

  df_long <- dplyr::bind_rows(long_list) |>
    dplyr::arrange(id, line, choice_option)

  return(list(
    df_demo = df_demo,
    df_long = df_long
  ))
}


#Managing covariate interactions

create_interaction_features <- function(df_long, choice_vars, demographic_vars){
  df_interactions <- df_long[, choice_vars]
  df_interactions_with <- df_long[, demographic_vars]
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


#Optional (can be used for comparision if needed)

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
    #parallel = TRUE
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

#Elastic net parameter tuning using BIC values

lasso_lambda_bic <- function(lambda_grid, alt_matrices, df_long, n = 10,
                             threshold = 0.01, N) {
  stopifnot(!missing(threshold))
  best_lambda <- NULL
  best_BIC <- Inf
  best_res <- NULL

  alt1 <- alt_matrices$alt1
  alt2 <- alt_matrices$alt2
  alt3 <- alt_matrices$alt3

  lambda_results <- data.frame(lambda = lambda_grid, BIC = NA, LL = NA)
  #N <- nrow(df_long)

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    start.values <- rep(0, n)

    res <- maxBFGS(
      function(coeff) MNL(coeff, alt_list, choice_list, lambda, alpha = 0.5, final_eval = FALSE,
                          nrep = nrep, intercept_index = 1),
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

    invisible(MNL(res$estimate, alt_list, choice_list, lambda, alpha = 0.5, final_eval = FALSE,
                  nrep = nrep, intercept_index = 1))
    start.values <- coef(res)

    res <- maxLik(
      function(coeff) MNL(coeff, alt_list, choice_list, lambda, alpha = 0.5, final_eval = FALSE,
                          nrep = nrep, intercept_index = 1),
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

    #unpenalized
    LL_unpenalized <- sum(MNL_unpenalized(res$estimate, alt_list, choice_list, final_eval = FALSE,
                                          nrep = nrep))
    lambda_results$LL[i] <- LL_unpenalized

    active_coeffs <- coef(res)[abs(coef(res)) >= threshold]
    k <- length(active_coeffs)
    lambda_results$k[i] <- k

    BIC_lasso <- -2 * LL_unpenalized + k * log(N)
    lambda_results$BIC[i] <- BIC_lasso

    if (BIC_lasso < best_BIC) {
      best_lambda <- lambda
      best_BIC <- BIC_lasso
      best_res <- res
    }
  }

  # cat("\n====Log-Likelihoods and BIC by Lambda(L1)====\n")
  # print(lambda_results)

  #lambda_results$k[i] <- k

  #unpenalized
  LL_unpenalized_best <- sum(MNL_unpenalized(best_res$estimate, alt_list, choice_list, final_eval = FALSE,
                                             nrep = nrep))

  cat("\n====Lambda tuning(BIC)====\n")
  print(lambda_results)
  cat("\nBest lambda based on BIC:", best_lambda, "\n")
  cat("Corresponding BIC:", best_BIC, "\n")
  cat("Corresponding Log-Likelihood:", LL_unpenalized_best, "\n")

  return(list(
    best_lambda = best_lambda,
    best_BIC = best_BIC,
    best_model = best_res,
    best_LL = LL_unpenalized_best,
    lambda_results = lambda_results
  ))
}


#Elastic net parameter tuning using 5 fold CV on out of sample log likelihood

tune_lambda_cv <- function(df_demo,
                           selected_features,
                           lambda_grid,
                           n_alt = 3,
                           n = 10,
                           n_folds = 5) {
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

      alt_train <- create_alt_matrices2(train_df, selected_features, demographic_vars, n_alt = 3)
      alt_test  <- create_alt_matrices2(test_df, selected_features, demographic_vars, n_alt = 3)

      alt_list_train <- lapply(1:n_alt, function(j) alt_train[[j]])
      alt_list_test  <- lapply(1:n_alt, function(j) alt_test[[j]])

      choice_list_train <- lapply(1:n_alt, function(j) train_df[[paste0("choice", j)]])
      choice_list_test  <- lapply(1:n_alt, function(j) test_df[[paste0("choice", j)]])

      start.values <- rep(0, n)

      res <- maxBFGS(
        function(coeff) MNL(coeff, alt_list_train, choice_list_train, lambda, alpha = 0.5, final_eval = FALSE,
                            nrep = nrep, intercept_index = 1),
        start = start.values,
        print.level = 0,
        iterlim = 200,
        finalHessian = FALSE
      )

      # Evaluate unpenalized LL on test data
      #ll_out_sample <- MNL_cv(res$estimate, alt_list_test, choice_list_test, lambda = 0, alpha = 0.5)
      ll_out_sample <- MNL_unpenalized(res$estimate, alt_list_test, choice_list_test, final_eval = FALSE,
                                       nrep = nrep)
      fold_lls[fold] <- sum(ll_out_sample)
    }

    mean_LL <- mean(fold_lls)
    lambda_results$mean_LL[i] <- mean_LL

    if (mean_LL > best_LL) {
      best_LL <- mean_LL
      best_lambda <- lambda
    }
  }

  cat("\n===== Lambda tuning summary (CV) =====\n")
  print(lambda_results)
  cat("\nBest lambda based on mean out-of-sample LL:", best_lambda, "\n")

  return(list(best_lambda = best_lambda, lambda_results = lambda_results))
}


#Print a final summary table with detailed information about the covariates

# summary_table_mnl <- function(model, selected_features, threshold = 0.01){
#   #stopifnot(!missing(threshold))
#   names(model$estimate) <- selected_features
#   summary_res <- summary(model)
#   coef_df <- as.data.frame(summary_res$estimate)
#   coef_df$Feature <- rownames(coef_df)
#   colnames(coef_df) <- c("Estimate", "Std.Error", "t.value", "p.value", "Feature")
#
#   coef_df$Estimate   <- round(coef_df$Estimate, 4)
#   coef_df$Std.Error  <- round(coef_df$Std.Error, 4)
#   coef_df$t.value    <- round(coef_df$t.value, 3)
#   coef_df$p.value    <- signif(coef_df$p.value, 3)
#
#   coef_df$Shrunk <- ifelse(abs(coef_df$Estimate) < threshold, "Yes", "No")
#   final_table <- coef_df[, c("Feature", "Estimate", "Std.Error", "t.value", "p.value", "Shrunk")]
#   final_table <- final_table[order(-abs(final_table$Estimate)), ]
#   print(final_table, row.names = FALSE)
#   return(final_table)
# }


summary_table_mnl <- function(model,
                              selected_features,
                              threshold = 0.01,
                              method = "MAXLIK"){

  #coefficients depending on method
  if(method == "BGW"){

    beta_hat <- as.numeric(model$estimate)
    se_hat <- model$seBGW
    t_hat  <- model$tstatBGW

    if (is.null(se_hat)) {
      se_hat <- rep(NA_real_, length(beta_hat))
    }

    if (is.null(t_hat)) {
      t_hat <- rep(NA_real_, length(beta_hat))
    }

    p_hat <- 2 * pnorm(-abs(t_hat))

    coef_df <- data.frame(
      Feature  = selected_features,
      Estimate = round(beta_hat, 4),
      Std.Error = round(se_hat, 4),
      t.value   = round(t_hat, 3),
      p.value   = signif(p_hat, 3)
    )

    coef_df$Shrunk <- ifelse(abs(coef_df$Estimate) < threshold, "Yes", "No")

    coef_df <- coef_df[, c("Feature", "Estimate", "Std.Error", "t.value", "p.value", "Shrunk")]
    coef_df <- coef_df[order(-abs(coef_df$Estimate)), ]

    print(coef_df, row.names = FALSE)
    return(coef_df)

  } else {

    #maxLik
    names(model$estimate) <- selected_features

    summary_res <- summary(model)
    coef_df <- as.data.frame(summary_res$estimate)

    coef_df$Feature <- rownames(coef_df)
    colnames(coef_df) <- c("Estimate", "Std.Error", "t.value", "p.value", "Feature")

    coef_df$Estimate   <- round(coef_df$Estimate, 4)
    coef_df$Std.Error  <- round(coef_df$Std.Error, 4)
    coef_df$t.value    <- round(coef_df$t.value, 3)
    coef_df$p.value    <- signif(coef_df$p.value, 3)

    coef_df$Shrunk <- ifelse(abs(coef_df$Estimate) < threshold, "Yes", "No")

    final_table <- coef_df[, c("Feature", "Estimate", "Std.Error", "t.value", "p.value", "Shrunk")]
    final_table <- final_table[order(-abs(final_table$Estimate)), ]

    print(final_table, row.names = FALSE)
    return(final_table)
  }
}



##Parallelization function

lasso_lambda_bic_parallel <- function(lambda_grid, alt_list, choice_list, n = 10,
                                      threshold = 0.01, N, alpha = 0.5) {
  stopifnot(!missing(threshold))
  library(future.apply)

  # Make sure the plan is set before running this
  # e.g., plan(multisession)

  results_list <- future_lapply(lambda_grid, function(lambda) {
    start.values <- rep(0, n)
    cat("Running for lambda =", lambda, "on PID", Sys.getpid(), "\n")
    # First optimization with maxBFGS
    res <- maxBFGS(
      function(coeff) MNL(coeff, alt_list, choice_list, lambda, alpha = alpha,
                          final_eval = FALSE, nrep = nrep, intercept_index = 1),
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

    # Refined optimization with maxLik
    start.values <- coef(res)

    res <- maxLik(
      function(coeff) MNL(coeff, alt_list, choice_list, lambda, alpha = alpha,
                          final_eval = FALSE, nrep = nrep, intercept_index = 1),
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

    # Compute unpenalized LL for BIC
    LL_unpenalized <- sum(MNL_unpenalized(res$estimate, alt_list, choice_list,
                                          final_eval = FALSE, nrep = nrep))

    #active_coeffs <- coef(res)[abs(coef(res)) >= threshold]
    active_coeffs <- coef(res)[abs(coef(res)) >= 1e-6]
    k <- length(active_coeffs)
    BIC_lasso <- -2 * LL_unpenalized + k * log(N)

    return(list(lambda = lambda, BIC = BIC_lasso, LL = LL_unpenalized, k = k, res = res))
  })

  # Bind results into dataframe
  lambda_results <- do.call(rbind, lapply(results_list, function(r) {
    data.frame(lambda = r$lambda, BIC = r$BIC, LL = r$LL, k = r$k)
  }))

  # Find the best result
  best_idx <- which.min(lambda_results$BIC)
  best_lambda <- lambda_results$lambda[best_idx]
  best_BIC <- lambda_results$BIC[best_idx]
  best_LL <- lambda_results$LL[best_idx]
  best_res <- results_list[[best_idx]]$res

  cat("\n====Lambda tuning (BIC)====\n")
  print(lambda_results)
  cat("\nBest lambda based on BIC:", best_lambda, "\n")
  cat("Corresponding BIC:", best_BIC, "\n")
  cat("Corresponding Log-Likelihood:", best_LL, "\n")

  return(list(
    best_lambda = best_lambda,
    best_BIC = best_BIC,
    best_model = best_res,
    best_LL = best_LL,
    lambda_results = lambda_results
  ))
}

#Tune using CV

# tune_lambda_cv_parallel <- function(df_demo,
#                                     selected_features,
#                                     lambda_grid = NULL,
#                                     demographic_vars,
#                                     n_alt = 3,
#                                     n = 10,
#                                     n_folds = 5,
#                                     n_lambda = 30,
#                                     lambda_min = 1e-4,
#                                     lambda_max = 1,
#                                     patience = 3) {
#   library(future.apply)
#
#   if (is.null(lambda_grid)) {
#     lambda_grid <- exp(seq(log(lambda_min), log(lambda_max), length.out = n_lambda))
#   }
#
#   # Create folds (respondent-wise split)
#   set.seed(123)
#   id_list <- unique(df_demo$id)
#   folds <- cut(seq_along(id_list), breaks = n_folds, labels = FALSE)
#   id_folds <- split(id_list, folds)
#
#   # Parallelized outer loop
#   lambda_results_list <- future_lapply(lambda_grid, function(lambda) {
#
#     cat("Running lambda =", lambda, "on PID", Sys.getpid(), "\n")
#
#     fold_lls <- numeric(n_folds)
#
#     for (fold in 1:n_folds) {
#       test_ids <- id_folds[[fold]]
#       train_ids <- setdiff(id_list, test_ids)
#
#       train_df <- df_demo[df_demo$id %in% train_ids, ]
#       test_df <- df_demo[df_demo$id %in% test_ids, ]
#
#       alt_train <- create_alt_matrices2(train_df, selected_features, demographic_vars, n_alt)
#       alt_test  <- create_alt_matrices2(test_df, selected_features, demographic_vars, n_alt)
#
#       alt_list_train <- lapply(1:n_alt, function(j) alt_train[[j]])
#       alt_list_test  <- lapply(1:n_alt, function(j) alt_test[[j]])
#
#       choice_list_train <- lapply(1:n_alt, function(j) train_df[[paste0("choice", j)]])
#       choice_list_test  <- lapply(1:n_alt, function(j) test_df[[paste0("choice", j)]])
#
#       start.values <- rep(0, n)
#
#       res <- maxBFGS(
#         function(coeff) MNL(coeff, alt_list_train, choice_list_train, lambda, alpha = 0.5, final_eval = FALSE,
#                             nrep = nrep, intercept_index = 1),
#         start = start.values,
#         print.level = 0,
#         iterlim = 200,
#         finalHessian = FALSE
#       )
#
#       ll_out_sample <- MNL_unpenalized(res$estimate, alt_list_test, choice_list_test, final_eval = FALSE,
#                                        nrep = nrep)
#       fold_lls[fold] <- sum(ll_out_sample)
#     }
#
#     mean_LL <- mean(fold_lls)
#     return(data.frame(lambda = lambda, mean_LL = mean_LL))
#   })
#
#   # Combine results
#   lambda_results <- do.call(rbind, lambda_results_list)
#
#   best_idx <- which.max(lambda_results$mean_LL)
#   best_lambda <- lambda_results$lambda[best_idx]
#   best_LL <- lambda_results$mean_LL[best_idx]
#
#   cat("\n===== Lambda tuning summary (CV, Parallel) =====\n")
#   print(lambda_results)
#   cat("\nBest lambda based on mean out-of-sample LL:", best_lambda, "\n")
#
#   return(list(best_lambda = best_lambda, lambda_results = lambda_results))
# }


tune_lambda_cv_parallel <- function(df_demo,
                                    selected_features,
                                    lambda_grid = NULL,
                                    demographic_vars,
                                    n_alt = 3,
                                    n = 10,
                                    n_folds = 5,
                                    n_lambda = 30,
                                    lambda_min = 1e-4,
                                    lambda_max = 1,
                                    patience = 3,
                                    alpha = 1,
                                    optimizer = "BFGS",
                                    early_stop = TRUE,
                                    nrep = NULL) {
  library(future.apply)

  # --- Log-spaced λ grid if not provided ---
  if (is.null(lambda_grid)) {
    lambda_grid <- exp(seq(log(lambda_min), log(lambda_max), length.out = n_lambda))
  }

  # Create folds
  set.seed(123)
  id_list <- unique(df_demo$id)
  folds <- cut(seq_along(id_list), breaks = n_folds, labels = FALSE)
  id_folds <- split(id_list, folds)

  results <- data.frame(lambda = numeric(0), mean_LL = numeric(0))
  best_LL <- -Inf
  no_improve_count <- 0

  for (lambda in lambda_grid) {
    fold_lls <- numeric(n_folds)

    for (fold in 1:n_folds) {
      test_ids <- id_folds[[fold]]
      train_ids <- setdiff(id_list, test_ids)

      train_df <- df_demo[df_demo$id %in% train_ids, ]
      test_df  <- df_demo[df_demo$id %in% test_ids, ]

      alt_train <- create_alt_matrices2(train_df, selected_features, demographic_vars, n_alt)
      alt_test  <- create_alt_matrices2(test_df,  selected_features, demographic_vars, n_alt)

      alt_list_train <- lapply(1:n_alt, function(j) alt_train[[j]])
      alt_list_test  <- lapply(1:n_alt, function(j) alt_test[[j]])

      choice_list_train <- lapply(1:n_alt, function(j) train_df[[paste0("choice", j)]])
      choice_list_test  <- lapply(1:n_alt, function(j) test_df[[paste0("choice", j)]])

      start.values <- rep(0, n)

      # res <- maxBFGS(
      # function(coeff) MNL(coeff, alt_list_train, choice_list_train,
      #                       lambda,
      #                       #alpha = 0.5,
      #                       alpha,
      #                       final_eval = FALSE, nrep = nrep, intercept_index = 1),
      #   start = start.values,
      #   print.level = 0,
      #   iterlim = 200,
      #   finalHessian = FALSE
      # )

      if (optimizer == "BFGS") {

        res <- maxBFGS(
          function(coeff) MNL(coeff,
                              alt_list_train,
                              choice_list_train,
                              lambda,
                              alpha,
                              final_eval = FALSE,
                              nrep = nrep,
                              intercept_index = 1,
                              out = "logprobs"),
          start = start.values,
          print.level = 0,
          iterlim = 200,
          finalHessian = FALSE
        )

        beta_hat <- res$estimate

      } else if (optimizer == "BGW") {

        calcR <- function(coeff) {

          coeff <- as.numeric(coeff)

          MNL(coeff,
              alt_list_train,
              choice_list_train,
              lambda,
              alpha,
              final_eval = FALSE,
              nrep = nrep,
              intercept_index = 1,
              out = "choiceprobs")
        }

        res <- bgw::bgw_mle(
          calcR = calcR,
          betaStart = start.values,
          bgw_settings = list(
            printLevel = 0
          )
        )

        #print(res$estimate)
        beta_hat <- as.numeric(res$estimate)


      } else {

        stop("Unsupported optimizer")

      }


      ll_out_sample <- MNL_unpenalized(beta_hat,
                                       alt_list_test,
                                       choice_list_test,
                                       final_eval = FALSE, nrep = nrep)


      fold_lls[fold] <- sum(ll_out_sample)
      #fold_lls[fold] <- mean(ll_out_sample)
    }

    mean_LL <- mean(fold_lls)
    results <- rbind(results, data.frame(lambda = lambda, mean_LL = mean_LL))

    # Early stopping, if LL doesn't improve it won't go further after 3 steps
    if(early_stop){

      #best_LL <- Inf
      if (mean_LL > best_LL) {
        best_LL <- mean_LL
        no_improve_count <- 0
      } else {
        no_improve_count <- no_improve_count + 1
      }

      if (no_improve_count >= patience) {
        message("Early stopping at lambda = ", lambda, " after ", patience, " declines.")
        break
      }}
  }

  best_idx <- which.max(results$mean_LL)
  best_lambda <- results$lambda[best_idx]

  cat("\n===== Lambda tuning summary (CV + Log Grid + Early stopping) =====\n")
  print(results)
  cat("\nBest lambda based on mean out-of-sample LL:", best_lambda, "\n")

  return(list(best_lambda = best_lambda, lambda_results = results))
}




####bespoke lambda tuning for semi synthetic data ONLY

tune_lambda_cv_bespoke <- function(
    alt_list,
    choice_list,
    lambda_grid,
    person_id,
    alpha = 1,
    n_folds = 5,
    start = NULL,
    method = "BHHH",
    th = NULL,
    optimizer = "BFGS"
) {

  ## NOT for the df_demo pipeline

  #number of alternatives and observations
  n_alt <- length(alt_list)
  n_obs <- nrow(alt_list[[1]])
  p     <- ncol(alt_list[[1]])

  #sanity checks
  stopifnot(length(choice_list) == n_alt)
  stopifnot(all(sapply(alt_list, nrow) == n_obs))
  stopifnot(length(lambda_grid) > 1)

  #default starting values
  if (is.null(start)) {
    start <- rep(0, p)
  }

  #CV folds with unique respondents in each fold, no leakage

  persons <- unique(person_id)
  fold_id_person <- sample(
    rep(seq_len(n_folds), length.out = length(persons))
  )
  names(fold_id_person) <- persons

  fold_id <- fold_id_person[as.character(person_id)]

  # build folds as row indices
  folds <- split(seq_len(n_obs), fold_id)

  #CV scores
  cv_ll <- numeric(length(lambda_grid))

  #loop over lambda values
  for (l in seq_along(lambda_grid)) {

    lambda <- lambda_grid[l]
    ll_fold <- numeric(n_folds)

    #loop over k folds

    for (k in seq_len(n_folds)) {

      idx_test  <- folds[[k]]
      idx_train <- setdiff(seq_len(n_obs), idx_test)


      alt_train <- lapply(alt_list, function(A)
        A[idx_train, , drop = FALSE]
      )
      alt_test <- lapply(alt_list, function(A)
        A[idx_test, , drop = FALSE]
      )

      #split choices
      choice_train <- lapply(choice_list, function(y)
        y[idx_train]
      )
      choice_test <- lapply(choice_list, function(y)
        y[idx_test]
      )

      if(optimizer == "BFGS" || optimizer == "BHHH"){

        #fit penalized model
        fit <- maxLik(
          function(b) MNL(
            b,
            alt_train,
            choice_train,
            lambda = lambda,
            alpha = alpha,
            intercept_index = NULL
          ),
          start = start,
          method = method,
          finalHessian = FALSE
        )
        bets_hat <- coef(fit)

      }  else if (optimizer == "BGW"){

        calcR <- function(coeff) {

          MNL(
            coeff,
            alt_train,
            choice_train,
            lambda = lambda,
            alpha = alpha,
            intercept_index = NULL,
            out = "choiceprobs"
          )
        }

        fit <- bgw::bgw_mle(
          calcR = calcR,
          betaStart = start,
          bgw_settings = list(printLevel = 0)
        )


        beta_hat <- as.numeric(fit$estimate)

      }

      #unpenalized log-likelihood on testing fold
      ll_test <- MNL(
        beta_hat,
        alt_test,
        choice_test,
        lambda = 0,          #no penalty on test data
        alpha = alpha,
        final_eval = TRUE
      )

      # #unpenalized log-likelihood on testing fold
      # ll_test <- MNL(
      #   coef(fit),
      #   #b_hat, #pass thresholded coefficients
      #   alt_test,
      #   choice_test,
      #   lambda = 0,          #no penalty on test data
      #   alpha = alpha,
      #   final_eval = TRUE
      # )

      #store sum log-likelihood
      ll_fold[k] <- sum(ll_test)
    }

    # average CV log-likelihood across folds
    cv_ll[l] <- mean(ll_fold)
  }

  #return optimal lambda

  list(
    best_lambda = lambda_grid[which.max(cv_ll)],
    cv_ll       = cv_ll,
    lambda_grid = lambda_grid
  )
}


#Threshold tuning to find the optimal threshold
tune_threshold_cv <- function(df_demo,
                              selected_features,
                              demographic_vars,
                              lambda,
                              alpha = 0.5,
                              threshold_grid = c(1e-4, 1e-3, 1e-2, 5e-2, 1e-1),
                              n_alt = 3,
                              n = length(selected_features),
                              n_folds = 5,
                              optimizer = "BFGS",
                              nrep = nrep,
                              rule = "best_ll",
                              ll_tolerance = 0) {

  set.seed(123)

  id_list <- unique(df_demo$id)
  folds <- cut(seq_along(id_list), breaks = n_folds, labels = FALSE)
  id_folds <- split(id_list, folds)

  threshold_results <- data.frame(
    threshold = threshold_grid,
    mean_LL = NA_real_,
    mean_n_retained = NA_real_
  )

  all_fold_lls <- matrix(
    NA_real_,
    nrow = length(threshold_grid),
    ncol = n_folds
  )

  all_fold_retained <- matrix(
    NA_real_,
    nrow = length(threshold_grid),
    ncol = n_folds
  )

  for (fold in seq_len(n_folds)) {

    test_ids <- id_folds[[fold]]
    train_ids <- setdiff(id_list, test_ids)

    train_df <- df_demo[df_demo$id %in% train_ids, ]
    test_df  <- df_demo[df_demo$id %in% test_ids, ]

    alt_train <- create_alt_matrices2(
      train_df,
      selected_features,
      demographic_vars,
      n_alt
    )

    alt_test <- create_alt_matrices2(
      test_df,
      selected_features,
      demographic_vars,
      n_alt
    )

    alt_list_train <- lapply(seq_len(n_alt), function(j) alt_train[[j]])
    alt_list_test  <- lapply(seq_len(n_alt), function(j) alt_test[[j]])

    choice_list_train <- lapply(seq_len(n_alt), function(j) train_df[[paste0("choice", j)]])
    choice_list_test  <- lapply(seq_len(n_alt), function(j) test_df[[paste0("choice", j)]])

    start.values <- rep(0, n)

    if (optimizer %in% c("BFGS", "BHHH")) {

      fit <- maxLik::maxLik(
        function(coeff)
          MNL(
            coeff,
            alt_list_train,
            choice_list_train,
            lambda = lambda,
            alpha = alpha,
            final_eval = FALSE,
            nrep = nrep,
            intercept_index = 1,
            out = "logprobs"
          ),
        start = start.values,
        method = optimizer,
        iterlim = 200,
        print.level = 0,
        finalHessian = FALSE
      )

      beta_hat <- as.numeric(fit$estimate)

    } else if (optimizer == "BGW") {

      calcR <- function(coeff) {
        MNL(
          as.numeric(coeff),
          alt_list_train,
          choice_list_train,
          lambda = lambda,
          alpha = alpha,
          final_eval = FALSE,
          nrep = nrep,
          intercept_index = 1,
          out = "choiceprobs"
        )
      }

      fit <- bgw::bgw_mle(
        calcR = calcR,
        betaStart = start.values,
        bgw_settings = list(printLevel = 0)
      )

      beta_hat <- as.numeric(fit$estimate)

    } else {
      stop("Unsupported optimizer")
    }

    for (i in seq_along(threshold_grid)) {

      th <- threshold_grid[i]
      beta_th <- beta_hat

      beta_th[abs(beta_th) < th] <- 0

      ll_test <- MNL_unpenalized(
        beta_th,
        alt_list_test,
        choice_list_test,
        final_eval = FALSE,
        nrep = nrep
      )

      all_fold_lls[i, fold] <- sum(ll_test)
      all_fold_retained[i, fold] <- sum(beta_th != 0)
    }
  }

  threshold_results$mean_LL <- rowMeans(all_fold_lls)
  threshold_results$mean_n_retained <- rowMeans(all_fold_retained)

  best_LL <- max(threshold_results$mean_LL)

  if (rule == "best_ll") {

    best_idx <- which.max(threshold_results$mean_LL)

  } else if (rule == "sparsest_within_tolerance") {

    eligible <- threshold_results[
      threshold_results$mean_LL >= best_LL - ll_tolerance,
      ,
      drop = FALSE
    ]

    best_idx <- as.integer(rownames(
      eligible[which.min(eligible$mean_n_retained), , drop = FALSE]
    ))

  } else {
    stop("Unsupported threshold rule. Use 'best_ll' or 'sparsest_within_tolerance'.")
  }

  best_threshold <- threshold_results$threshold[best_idx]

  cat("\n===== Threshold tuning summary =====\n")
  print(threshold_results)
  cat("\nSelected threshold:", best_threshold, "\n")

  return(list(
    best_threshold = best_threshold,
    threshold_results = threshold_results
  ))
}
