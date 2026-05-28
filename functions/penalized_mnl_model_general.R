run_penalized_mnl <- function(
    data,
    choice_vars = NULL,
    demographic_vars = NULL,
    #lambda_grid = exp(seq(log(6e-4), log(3e-2), length.out = 10)),
    lambda_grid = NULL,
    n_lambda = 20,
    lambda_min = 1e-8,
    lambda_max = 1e-1,
    n_alt = 3,
    n_folds = 5,
    alpha = 0.5,
    threshold = 0.01,
    threshold_grid = NULL,
    tune_threshold = FALSE,
    threshold_rule = "best_ll",
    threshold_ll_tolerance = 0,
    optimizer = "BFGS",
    method = "MAXLIK",
    data_schema = NULL,
    n_workers = NULL,
    nrep = NULL
){
#browser()
  required_packages <- c(
    "maxLik","matrixStats","tidyr","dplyr",
    "glmnet","bgw","Rfast","future.apply",
    "future","ggplot2", "parallelly"
  )

  install_if_missing <- function(packages) {
    missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
    
    if(length(missing_packages) > 0) {
      install.packages(missing_packages, dependencies = TRUE)  # Install missing packages
    }
    
    # Load all packages
    lapply(packages, require, character.only = TRUE)
  }
  
  install_if_missing(required_packages)

  source("functions/utility_functions.R")
  source("functions/utility_functions_generalized.R")
  source("functions/mnl_function.R")
  source("functions/pre_process.R")
  source("functions/MNL_function_any_data.R") #Handles threshold detection

  #preprocessing

  # Optional schema-driven adapter for non-Dogger-Bank datasets
  if (!is.null(data_schema)) {
    
    prepared <- do.call(
      prepare_mnl_data,
      c(list(data = data), data_schema)
    )
    
    data <- prepared$data
    choice_vars <- prepared$choice_vars
    demographic_vars <- prepared$demographic_vars
    n_alt <- prepared$n_alt
    
    cat("\nPrepared choice vars:\n")
    print(choice_vars)
    
    cat("\nPrepared demographic vars:\n")
    print(demographic_vars)
    
    cat("\nPrepared data columns:\n")
    print(names(data))
  }
  
  if (is.null(choice_vars) || is.null(demographic_vars)) {
    stop("Provide choice_vars and demographic_vars, or provide data_schema.")
  }
  
  # If the user does not provide a lambda grid, create a broad log-spaced grid.
  # n_lambda controls how fine the search is.
  if (is.null(lambda_grid)) {
    lambda_grid <- exp(seq(
      log(lambda_min),
      log(lambda_max),
      length.out = n_lambda
    ))
  }
  
  cat("\nLambda grid used for CV:\n")
  print(lambda_grid)
  
  #nrep detection
  if (is.null(nrep)) {
    
    if (!("id" %in% names(data))) {
      stop(
        "Cannot auto-detect nrep because column 'id' is not present. ",
        "Provide nrep manually or make sure prepare_mnl_data() creates an 'id' column."
      )
    }
    
    reps_per_id <- table(data$id)
    
    if (length(unique(reps_per_id)) != 1) {
      stop(
        "Cannot auto-detect nrep because respondents have unequal numbers of choice tasks. ",
        "Please provide nrep manually or balance the data first."
      )
    }
    
    nrep <- as.integer(unique(reps_per_id))
  }
  
  cat("\nUsing nrep =", nrep, "\n")
  
  # Convert standardized wide data to long data for interaction construction
  output <- data_wide_to_long(
    data = data,
    choice_vars = choice_vars,
    demographic_vars = demographic_vars,
    n_alt = n_alt
  )
  
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
  
  cat("\nSelected features entering model:\n")
  print(selected_features)
  
  cat("\nNumber of features entering model:", num_covariates, "\n")

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

  if (is.null(n_workers)) {
    n_workers <- max(1, min(8, parallelly::availableCores() - 1))
  }
  
  plan(multisession, workers = n_workers)
  
  options(future.rng.onMisuse = "ignore")

  results_cv <- tune_lambda_cv_parallel(
    df_demo,
    selected_features,
    lambda_grid,
    demographic_vars,
    n_alt = n_alt,
    n = num_covariates,
    n_folds = n_folds,
    alpha = alpha,
    optimizer = optimizer,
    early_stop = TRUE,
    nrep = nrep
  )

  best_lambda <- results_cv$best_lambda
  
  if (best_lambda == min(lambda_grid)) {
    warning(
      "Best lambda is at the lower boundary of the lambda grid. ",
      "Consider decreasing lambda_min."
    )
  }
  
  if (best_lambda == max(lambda_grid)) {
    warning(
      "Best lambda is at the upper boundary of the lambda grid. ",
      "Consider increasing lambda_max."
    )
  }
  
  #threshold tuning
  threshold_results <- NULL
  
  if (tune_threshold) {
    
    if (is.null(threshold_grid)) {
      threshold_grid <- c(1e-4, 1e-3, 1e-2, 5e-2, 1e-1)
    }
    
    threshold_cv <- tune_threshold_cv(
      df_demo = df_demo,
      selected_features = selected_features,
      demographic_vars = demographic_vars,
      lambda = best_lambda,
      alpha = alpha,
      threshold_grid = threshold_grid,
      n_alt = n_alt,
      n = num_covariates,
      n_folds = n_folds,
      optimizer = optimizer,
      nrep = nrep,
      rule = threshold_rule,
      ll_tolerance = threshold_ll_tolerance
    )
    
    threshold <- threshold_cv$best_threshold
    threshold_results <- threshold_cv$threshold_results
  }

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
          nrep = nrep,
          intercept_index = 1,
          out = "logprobs"
        ),
      start = start.values,
      method = optimizer,
      iterlim = 200,
      print.level = 0,
      finalHessian = TRUE
    )
    
    # coefficients_table <- summary_table_mnl(
    #   final_model,
    #   selected_features,
    #   threshold = threshold,
    #   method = "maxLik"
    # )
    
  } else if (optimizer == "BGW") {
    
    calcR <- function(coeff){
      MNL(
        coeff,
        alt_list,
        choice_list,
        lambda = best_lambda,
        alpha = alpha,
        final_eval = FALSE,
        nrep = nrep,
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
    best_threshold = threshold,
    cv_results = results_cv$lambda_results,
    threshold_results = threshold_results
  )

}
