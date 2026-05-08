###############################################################################
# Penalized Multinomial Logit (MNL) example

# Estimate a MNL model that can automatically screen
# many possible interaction effects between choice attributes and respondent
# characteristics.

# Practical interpretation:
# - Large coefficients indicate stronger effects.
# - Very small coefficients may be shrunk toward zero and can be treated as
#   weak or unimportant effects.

# This file contains the penalized MNL workflow
# Also installs/loads missing packages
source("functions/penalized_mnl_model.R")

### CAUTION: SOME PACKAGES MAY TAKE A LONG TIME TO INSTALL ###
#Libraries

required_packages <- c("maxLik","matrixStats","tidyr","dplyr",
                       "glmnet","bgw","Rfast","future.apply",
                       "future","ggplot2", "parallelly")

install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  if(length(missing_packages) > 0) {
    install.packages(missing_packages, dependencies = TRUE)  # Install missing packages
  }
  
  # Load all packages
  lapply(packages, require, character.only = TRUE)
}

install_if_missing(required_packages)

# Dogger Bank choice experiment data formatted like Apollo
data <- read.csv("data/doggerbank_full_regularization.csv")

# "main effects" of the choice model.
choice_vars <- c(
  "ASC2","ASC3",
  "cost","spec10","spec25",
  "prot25","prot50","invasive"
) #main effects

# Respondent characteristics or demographic variables
# add/change as needed 
demographic_vars <- c(
  "male", "edu", "job1")

####DEFAULT VALUES

# lambda_grid <- exp(seq(log(6e-4), log(3e-2), length.out = 10))
# n_alt <- 3
# n_folds <- 5
# alpha <- 0.5
# threshold <- 0.01 #screening threshold

# Estimate the penalized MNL model

# optimizer:
#
# "BGW" is faster.
# "BFGS" is the general-purpose optimizer.

# method:
# Controls how the final results are handled.
# Use "BGW" if estimating with BGW.
# Use "MAXLIK" if estimating with BFGS through maxLik.

# model <- run_penalized_mnl(
#   data,
#   choice_vars,
#   demographic_vars,
#   lambda_grid,
#   n_alt,
#   n_folds,
#   alpha,
#   threshold,
#   optimizer = "BGW", #Choose BGW or BFGS
#   method = "BGW" #Choose BGW or MAXLIK (if optimizer = "BFGS")
# )

model <- run_penalized_mnl(
  data = data,
  choice_vars = choice_vars,
  demographic_vars = demographic_vars,
  optimizer = "BGW",
  method = "BGW",
  n_workers = 30 #change as needed
)

# estimated coefficients
model$coefficients

####PLOT RESULTS (uncomment if needed)

#to get the LL plot, plot the output from model$cv_results

# plot_cv <- ggplot(model$cv_results, aes(x = lambda, y = mean_LL)) +
#   geom_line(linewidth = 1.5, colour = "grey") +
#   geom_point(size = 3) +
#   scale_x_log10() +
#   theme_bw(base_size = 16) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"),
#     text = element_text(family = "Helvetica"),
#     axis.text = element_text(size = 20),
#     axis.title = element_text(size = 30)
#   ) +
#   labs(
#     x = "Lambda",
#     y = "Mean Out-of-Sample LL"
#   )
# 
# plot_cv

#Save plot

# fname <- paste0(
#   "plots/ll_vs_lambda",
#   "_alpha_", alpha,
#   "_demvars_", length(demographic_vars),
#   ".png"
# )
# ggsave(fname, plot_cv, width = 10, height = 7, dpi = 400)
