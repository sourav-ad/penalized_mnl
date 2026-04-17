
run_penalized_mnl <- function(
    data,
    choice_vars,
    demographic_vars,
    lambda_grid = exp(seq(log(6e-4), log(3e-2), length.out = 10)),
    n_alt = 3,
    n_folds = 5,
    alpha = 0.5,
    threshold = 0.01,
    optimizer = "BFGS",
    method = "MAXLIK"
){
#browser()
  required_packages <- c(
    "maxLik","matrixStats","tidyr","dplyr",
    "glmnet","bgw","Rfast","future.apply",
    "future","ggplot2"
  )

  lapply(required_packages, require, character.only = TRUE)

  source("functions/utility_functions.R")
  source("functions/utility_functions_generalized.R")
  source("functions/mnl_function.R")
  source("functions/pre_process.R")
  source("functions/MNL_functions_execution.R")

  #preprocessing

  output <- data_wide_to_long(data, n_alt)

  df_demo <- output$df_demo
  df_long <- output$df_long

  #interactions

  final_df <- create_interaction_features(
    df_long,
    choice_vars,
    demographic_vars
  )

  selected_features <- colnames(final_df)

  num_covariates <- ncol(final_df)

  #alternative matrices

  alt_matrices <- create_alt_matrices2(
    df_demo,
    selected_features = selected_features,
    demographic_vars = demographic_vars,
    n_alt = n_alt
  )

  alt_list <- lapply(1:n_alt, function(j) alt_matrices[[j]])
  choice_list <- lapply(1:n_alt, function(j) df_demo[[paste0("choice", j)]])

  #cv

  plan(multisession)
  options(future.rng.onMisuse = "ignore")

  results_cv <- tune_lambda_cv_parallel(
    df_demo,
    selected_features,
    lambda_grid,
    demographic_vars,
    n_alt = n_alt,
    n = num_covariates,
    n_folds = n_folds,
    optimizer = optimizer,
    early_stop = TRUE
  )

  best_lambda <- results_cv$best_lambda

  #final model

  start.values <- rep(0, length(selected_features))

  if (optimizer %in% c("BHHH", "BFGS")) {
    
    final_model <- maxLik(
      function(coeff)
        MNL(
          coeff,
          alt_list,
          choice_list,
          lambda = best_lambda,
          alpha = alpha,
          final_eval = FALSE,
          nrep = 6,
          intercept_index = 1,
          out = "logprobs"
        ),
      start = start.values,
      method = optimizer,
      iterlim = 200,
      print.level = 0,
      finalHessian = TRUE
    )
    
    coefficients_table <- summary_table_mnl(
      final_model,
      selected_features,
      threshold = threshold,
      method = "maxLik"
    )
    
  } else if (optimizer == "BGW") {
    
    calcR <- function(coeff){
      MNL(
        coeff,
        alt_list,
        choice_list,
        lambda = best_lambda,
        alpha = alpha,
        final_eval = FALSE,
        nrep = 6,
        intercept_index = 1,
        out = "choiceprobs"
      )
    }
    
    
    final_model <- bgw::bgw_mle(
      calcR = calcR,
      betaStart = start.values,
      bgw_settings = list(printLevel = 0)
    )
    
    # coefficients_table <- summary_table_mnl(
    #   final_model,
    #   selected_features,
    #   threshold = threshold,
    #   method = method
    # )
    
  } else {
    stop("Unsupported optimizer")
  }

  coefficients_table <- summary_table_mnl(
    final_model,
    selected_features,
    threshold = threshold,
    method = method
  )

  plan(sequential)

  #outputs

  list(
    model = final_model,
    coefficients = coefficients_table,
    best_lambda = best_lambda,
    cv_results = results_cv$lambda_results
  )

}
