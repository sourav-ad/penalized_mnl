###MNL function to obtain Log likelihood

#Generalized code for any number of utility function is provided below

#Used in colprods()
nrep <- 6

intercept_index <- 1

MNL <- function(coeff, alt1, alt2, alt3, lambda, alpha = 0.5, final_eval = FALSE) {

  # util1 = (alt1 %*% coeff[1:n])
  # util2 = (alt2 %*% coeff[1:n])
  # util3 = (alt3 %*% coeff[1:n])
  
  util1 = (alt1 %*% coeff[1:ncol(alt1)])
  util2 = (alt2 %*% coeff[1:ncol(alt2)])
  util3 = (alt3 %*% coeff[1:ncol(alt3)])

  innerArg = ((exp(util1) * df_demo$choice1) + (exp(util2) * df_demo$choice2) + (exp(util3) * df_demo$choice3)) / ( exp(util1) + exp(util2) + exp(util3))
  #returns a vector of probabilities
  choice_probs = colprods(matrix(innerArg, nrow=nrep))
  #choice_probs = ()

  # log-likelihood
  LL = log(choice_probs)

  #We need a vector of log likelihoods
  if (length(LL) == 1) {
    stop("Error: LL is returning a single scalar instead of a vector.")
  }

  #L1 norm
  # penalty <- if(!is.null(intercept_index)){
  #   lambda * sum(abs(coeff[-intercept_index]))
  # } else {
  #   lambda * sum(abs(coeff))
  # }

  #L1 + L2
  penalty <- if (!is.null(intercept_index)) {
    coeff_excl_intercept <- coeff[-intercept_index]
    lambda * (alpha * sum(abs(coeff_excl_intercept)) + ((1 - alpha) / 2) * sum(coeff_excl_intercept^2))
  } else {
    lambda * (alpha * sum(abs(coeff)) + ((1 - alpha) / 2) * sum(coeff^2))
  }
  
  LL_elastic_net <- LL - penalty
  return(LL_elastic_net)
}



nrep <- 6
intercept_index<- 1

#Generalized function
MNL_cv <- function(coeff, alt_list, choice_list, lambda, alpha = 0.5, final_eval = FALSE) {
  n_alt <- length(alt_list)  # Number of alternatives
  
  # Compute utilities for each alternative
  utils <- lapply(alt_list, function(alt) alt %*% coeff[1:ncol(alt)])
  
  # Exponentiate utilities
  exp_utils <- lapply(utils, exp)
  
  # Compute numerator: sum(exp(util) * choice)
  numerator <- Reduce(`+`, mapply(`*`, exp_utils, choice_list, SIMPLIFY = FALSE))
  
  # Compute denominator: sum(exp(util))
  denominator <- Reduce(`+`, exp_utils)
  
  # Compute choice probabilities
  choice_probs <- colprods(matrix(numerator / denominator, nrow = nrep))
  
  # Log-likelihood
  LL <- log(choice_probs)
  
  if (length(LL) == 1) {
    stop("Error: LL is returning a single scalar instead of a vector.")
  }
  
  #L1 norm
  # penalty <- if(!is.null(intercept_index)){
  #   lambda * sum(abs(coeff[-intercept_index]))
  # } else {
  #   lambda * sum(abs(coeff))
  # }
  
  #L1 + L2, elastic net
  penalty <- if (!is.null(intercept_index)) {
    coeff_excl_intercept <- coeff[-intercept_index]
    lambda * (alpha * sum(abs(coeff_excl_intercept)) + ((1 - alpha) / 2) * sum(coeff_excl_intercept^2))
  } else {
    lambda * (alpha * sum(abs(coeff)) + ((1 - alpha) / 2) * sum(coeff^2))
  }
  
  LL_elastic_net <- LL - penalty
  return(LL_elastic_net)
}
