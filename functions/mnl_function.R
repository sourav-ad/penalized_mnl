###MNL function to obtain Log likelihood

#Generalized code for any number of utility function is provided below

#Used in colprods()
nrep <- 6

intercept_index <- 1

MNL <- function(coeff, alt1, alt2, alt3, lambda, final_eval = FALSE) {

  util1 = (alt1 %*% coeff[1:n])
  util2 = (alt2 %*% coeff[1:n])
  util3 = (alt3 %*% coeff[1:n])

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
  penalty <- if(!is.null(intercept_index)){
    lambda * sum(abs(coeff[-intercept_index]))
  } else {
    lambda * sum(abs(coeff))
  }

  LL_lasso <- LL - penalty
  #return the penalized log likelihood
  return(LL_lasso)
}


# #compact test:
# if (0){
#   save(start.values, alt1, alt2, alt3, df_demo, file="data/tmp_elasticnet.Rda")
#   load("data/tmp_elasticnet.Rda")
#   LL = MNL(start.values, alt1, alt2, alt3, nrep=6)
# }