nrep <- 6
intercept_index<- 1

#Generalized function
MNL <- function(coeff, alt_list, choice_list, lambda, final_eval = FALSE) {
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
  penalty <- if(!is.null(intercept_index)){
    lambda * sum(abs(coeff[-intercept_index]))
  } else {
    lambda * sum(abs(coeff))
  }
  
  LL_lasso <- LL - penalty
  #return the penalized log likelihood
  return(LL_lasso)
}
