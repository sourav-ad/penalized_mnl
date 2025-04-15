###########MULTINOMIAL LOGIT WITH ELASTIC NET################

###To run the code fully utilizing readline(), use source("mnl_extended.R", echo = TRUE)

# Uses the complete Dogger Bank data in wide format
# Turns a wide format survey data into long form
# Performs an elastic net regression based on the response variable
# Displays top candidates by coefficient if needed
# Uses the top n candidates in an MNL model 
# MNL model is implemented from scratch
# Displays statistical significance
# USER CAN INTERACT TO CHANGE THE COVARIATES
# Slower, since more covariates would lead to longer runtime for elastic net regression

# -------------------------------------------------------------------------


#Libraries

required_packages <- c("maxLik", "matrixStats", "tidyr", "dplyr", "glmnet", "bgw", "Rfast")

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
#Set appropriate working directory
source("functions/utility_functions.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")

#Data (change path as needed)
data <- read.csv("data/doggerbank_full_973_wide.csv")
#person id
data$person <- rep(1:(dim(data)[1]/6), each=6)

#Remove columns 'choice1' to 'choice 6' and create new columns instead
#We modify the data to have an idea about each respondent's choice among the 3 options provided 
data <- data[, !names(data) %in% c('choice1', 'choice2', 'choice3', 'choice4', 'choice5', 'choice6')]

n_alt <- 3

#Automated for any number of choices. Response variable must be named y1, y2,...

# choice variable
data$choice <- 0

# data$choice[data$y1==1] <- 1
# data$choice[data$y2==1] <- 2
# data$choice[data$y3==1] <- 3

#As many as user wants

for (i in 1:n_alt) {
  data$choice[data[[paste0("y", i)]] == 1] <- i
}

#Choice options for compatibility with given code
#Since the utility function is defined differently

# data$choice1 <- ifelse(data$choice == 1, 1, 0)
# data$choice2 <- ifelse(data$choice == 2, 1, 0)
# data$choice3 <- ifelse(data$choice == 3, 1, 0)

# choice1, choice2,...

for (i in 1:n_alt) {
  data[[paste0("choice", i)]] <- ifelse(data$choice == i, 1, 0)
}


#pre-processing
#remove columns with more than 50% empty rows
#data <- data[, colMeans(is.na(data)) <= 0.5]

#rename columns q22_7, q22_9 for easier recognition of interaction

data <- data %>%
  rename('q227' = 'q22_7', 'q229' = 'q22_9')

#View(data)

#Automate this

#Take a subset of the data
df_demo <- data[, c('id', 'line',
                    'invasive1', 'cost1', 'spec101', 'spec251', 'prot251', 'prot501',
                    'invasive2','cost2', 'spec102', 'spec252', 'prot252', 'prot502',
                    'invasive3', 'cost3', 'spec103', 'spec253', 'prot253', 'prot503', 'q227', 'q229',
                    'edu', 'male', 'job', 'age', 'choice', 'y1', 'y2', 'y3',
                    'choice1', 'choice2', 'choice3', 'job1', 'job2', 
                    'job3', 'job4', 'job5', 'job6', 'job7', 'job8', 'q1', 'q2', 'q6', 'q7', 'q10')]

#If adding columns that are constant across responses, change 'constant_vars' in utility_functions.R


# Convert wide data format to long
df_long <- df_demo %>%
  pivot_longer(
    cols = starts_with('cost')|starts_with('spec10')|starts_with('spec25')|starts_with('prot25')|starts_with('prot50')|starts_with('invasive')|starts_with('protest'),
    names_to = c(".value", "choice_option"),
    names_pattern = "(.*)([1-3])"
  ) %>%
  mutate(
    choice_option = as.integer(choice_option),  # Convert choice_option to integer
    chosen = as.integer(choice_option == choice)  # Create binary indicator for chosen option
  ) %>%
  arrange(id, line, choice_option)  # Arrange the data for clarity


df_long <- df_long %>%
  mutate(
    chosen = ifelse(choice_option == choice, 1, 0)  # 1 if choice_option matches choice, else 0
  )

#We now use Elastic net to uncover which covariates and their interactions are more influential
#Helps us choose candidates for utility functions in the choice model

#Interested in demographic feature interaction


#User interaction to determine interactions


#Without user interaction

#Columns that vary with response
df_interactions <- df_long[, c('cost', 'spec10', 'spec25', 'prot25', 'prot50',
                               'invasive')]

#Columns that do not vary with response
df_interactions_with <- df_long[, c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
                                    'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
                                    'job5', 'job6', 'job7', 'job8')]

#Create an empty dataframe with same rows as df_interactions
interaction_df <- data.frame(matrix(nrow = nrow(df_interactions), ncol = 0))

# Create interactions **ONLY between df_interactions and df_interactions_with**
for (col1 in colnames(df_interactions)) {
  for (col2 in colnames(df_interactions_with)) {
    interaction_name <- paste0(col1, "_", col2)
    
    # Create interaction term as the product
    interaction_df[[interaction_name]] <- df_interactions[[col1]] * df_interactions_with[[col2]]
  }
}

#View(interaction_df)

#Combine long df with the interaction df
#Will have columns of original features + interaction between features
final_df <- cbind(df_interactions, interaction_df)

#Scale and remove NA if any, before Elastic Net
final_df_scaled <- final_df %>%
  mutate(across(where(is.numeric), scale))

final_df_scaled <- final_df_scaled[, colSums(is.na(final_df_scaled)) == 0]

#Elastic Net on columns of original features + interaction between features
#Uncomment if a different df_interactions is chosen (a different subset of columns)
X <- as.matrix(final_df_scaled)
y <- as.numeric(df_long$chosen) #Column containing information if a response was chosen or not

set.seed(123)

#Run the elastic net if a different df_demo is used
#Save the best_lambda accordingly
elastic_net_model <- cv.glmnet(
  x = X,
  y = y,
  alpha = 0.5,            # Elastic Net (0.5 = balance between Lasso and Ridge)
  family = "binomial",    # For binary classification
  nfolds = 5,             # 5-fold cross-validation
  maxit = 5000,          # High number of iterations for convergence
  type.measure = "class"  # Classification accuracy
)

best_lambda <- elastic_net_model$lambda.min
cat("Best Lambda:", best_lambda, "\n")

#use the best lambda as obtained to fit the final model

final_model <- glmnet(
  x = X,
  y = y,
  alpha = 0.5,
  family = "binomial",
  lambda = best_lambda
)

#Save the coefficients
#Will contain coefficients of covariates and their interactions in descending order
#Top n candidates can be chosenf or the MNL
coefficients <- coef(final_model)
coefficients_df <- as.data.frame(as.matrix(coefficients))
colnames(coefficients_df)[1] <- "coefficient"
coefficients_df$feature <- rownames(coefficients_df)

sorted_coefficients <- coefficients_df %>%
  filter(feature != "(Intercept)") %>%
  arrange(desc(abs(coefficient)))

########CHANGE AS NEEDED###############
#top n coefficients from Elastic Net
n <- 20
########CHANGE AS NEEDED###############
selected_features <- sorted_coefficients$feature[1:n]
#print(selected_features)

#Use the function create_alt_matrices() on the df_demo subset of data
#Use the wide form

alt_matrices <- create_alt_matrices(df_demo, selected_features, n_alt)

#create_alt_matrices() to be generalized to any number of utility functions first

#utility function matrices (applicable when 3 choices are presented)
alt1 <- alt_matrices$alt1
alt2 <- alt_matrices$alt2
alt3 <- alt_matrices$alt3

#alt_list <- mget(paste0("alt", 1:n_alt))
alt_list <- list(alt1, alt2, alt3)


#A required parameter in the code 
nset <- nrow(df_demo)
L1_lambda_grid <- c(seq(0.005, 0.009, 0.001))

best_L1_lambda <- NULL
best_BIC <- Inf
best_res <- NULL

# Store all BICs
lambda_results <- data.frame(lambda = L1_lambda_grid, BIC = NA)

#initial values of betas
start.values <- rep(0, n)

#MNL() function to be generalized to any number as well

choice_list <- lapply(1:n_alt, function(i) df_demo[[paste0("choice", i)]])

for (i in seq_along(L1_lambda_grid)) {
  
  lambda <- L1_lambda_grid[i]
  #initial values of betas
  start.values <- rep(0, n)
  # estimate model (without inverting Hessian)
  res = maxBFGS(
        #MNL,
        function(coeff) MNL(coeff, alt_list, choice_list, lambda, final_eval = FALSE),
        grad=NULL,
        hess=NULL,
        start=start.values,
        fixed=NULL,
        print.level=0,
        iterlim=200,
        constraints=NULL,
        tol=1e-25, reltol=1e-25,
        finalHessian=FALSE,
        parscale=rep(1, length=length(start))
      )
  
  invisible(MNL(res$estimate, alt_list, choice_list, lambda, final_eval = TRUE))
  
  # estimate model (with inverting Hessian)
  #new start values are the values obtained above as starting values
  start.values = coef(res)
  
  res = maxLik(
        #MNL,
        function(coeff) MNL(coeff, alt_list, choice_list, lambda, final_eval = FALSE),  
        grad=NULL, 
        hess=NULL, 
        start=start.values, 
        fixed=NULL, 
        print.level=1, 
        method="BHHH", 
        iterlim=2,
        #iterlim=0, 
        constraints=NULL, 
        tol=1e-04, reltol=1e-04,
        finalHessian=TRUE
      )
  
  # Compute penalized BIC (based on penalized LL)
  N <- nrow(df_long)
  LL_unpenalized <- MNL(res$estimate, alt_list, choice_list, lambda = 0, final_eval = FALSE) 
  LL_unpenalized<- sum(LL_unpenalized)
  threshold <- 1e-3
  active_coeffs <- coef(res)[abs(coef(res)) >= threshold]
  k <- length(active_coeffs)
  print(k)
  BIC_lasso <- -2 * LL_unpenalized + k * log(N)
  lambda_results$BIC[i] <- BIC_lasso
  
  # Check for best
  if (BIC_lasso < best_BIC) {
    best_L1_lambda <- lambda
    best_BIC <- BIC_lasso
    best_res <- res
  }
  
}

plot(lambda_results$lambda, lambda_results$BIC, type = "b", 
     xlab = "Lambda (L1 Penalty)", ylab = "BIC", 
     main = "Model Selection using BIC")

cat("\n===== Lambda tuning summary =====\n")
print(lambda_results)
cat("\nBest lambda based on BIC:", best_L1_lambda, "\n")
cat("Best BIC:", best_BIC, "\n")

names(best_res$estimate) <- selected_features
MNL.res <- best_res
print(summary(MNL.res))

final_coeff <- best_res$estimate
threshold <- 1e-3
zero_indices <- which(abs(final_coeff) < threshold)
zero_coeffs <- names(final_coeff)[zero_indices]

cat("==Coefficients shrunk after Lasso regularization==\n")

if(length(zero_coeffs) == 0){
  cat("None of the coeffcients were shrunk \n")
} else {
  cat(paste(zero_coeffs, collapse = ", "), "\n")
}