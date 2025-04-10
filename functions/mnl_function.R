###MNL function to obtain Log likelihood

#Generalized code for any number of utility function is provided below

#Used in colprods()
nrep <- 6

# MNL <- function(coeff, alt_list, choice_list, nrep) {
#   n_alt <- length(alt_list)  # Number of alternatives
# 
#   # Compute utilities for each alternative
#   utils <- lapply(alt_list, function(alt) alt %*% coeff[1:ncol(alt)])
# 
#   # Exponentiate utilities
#   exp_utils <- lapply(utils, exp)
# 
#   # Compute numerator: sum(exp(util) * choice)
#   numerator <- Reduce(`+`, mapply(`*`, exp_utils, choice_list, SIMPLIFY = FALSE))
# 
#   # Compute denominator: sum(exp(util))
#   denominator <- Reduce(`+`, exp_utils)
# 
#   # Compute choice probabilities
#   choice_probs <- colprods(matrix(numerator / denominator, nrow = nrep))
# 
#   # Log-likelihood
#   LL <- log(choice_probs)
# 
#   if (length(LL) == 1) {
#     stop("Error: LL is returning a single scalar instead of a vector.")
#   }
# 
#   return(LL)
# }


MNL <- function(coeff, alt1, alt2, alt3) {
  
  util1 = (alt1 %*% coeff[1:n])
  util2 = (alt2 %*% coeff[1:n])
  util3 = (alt3 %*% coeff[1:n])
  
  innerArg = ((exp(util1) * df_demo$choice1) + (exp(util2) * df_demo$choice2) + (exp(util3) * df_demo$choice3)) / ( exp(util1) + exp(util2) + exp(util3))
  #returns a vector of probabilities
  choice_probs = colprods(matrix(innerArg,nrow=nrep))
  #choice_probs = ()
  
  # log-likelihood  
  LL = log(choice_probs)
  
  #We need a vector of log likelihoods
  if (length(LL) == 1) {
    stop("Error: LL is returning a single scalar instead of a vector.")
  }
  
  return(LL)
}




#Generalized function
# MNL <- function(coeff, alt_list, choice_list) {
#   n_alt <- length(alt_list)  # Number of alternatives
#   
#   # Compute utilities for each alternative
#   utils <- lapply(alt_list, function(alt) alt %*% coeff[1:ncol(alt)])
#   
#   # Exponentiate utilities
#   exp_utils <- lapply(utils, exp)
#   
#   # Compute numerator: sum(exp(util) * choice)
#   numerator <- Reduce(`+`, mapply(`*`, exp_utils, choice_list, SIMPLIFY = FALSE))
#   
#   # Compute denominator: sum(exp(util))
#   denominator <- Reduce(`+`, exp_utils)
#   
#   # Compute choice probabilities
#   choice_probs <- colprods(matrix(numerator / denominator, nrow = nrep))
#   
#   # Log-likelihood
#   LL <- log(choice_probs)
#   
#   if (length(LL) == 1) {
#     stop("Error: LL is returning a single scalar instead of a vector.")
#   }
#   
#   return(LL)
# }


# #compact test:
# if (0){
#   save(start.values, alt1, alt2, alt3, df_demo, file="data/tmp_elasticnet.Rda")
#   load("data/tmp_elasticnet.Rda")
#   LL = MNL(start.values, alt1, alt2, alt3, nrep=6)
# }