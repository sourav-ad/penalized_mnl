###########MULTINOMIAL LOGIT WITH ELASTIC NET################

# Turns a wide format survey data into long form
# Performs an elastic net regression based on the response variable
# Displays top candidates by coefficient if needed
# Uses the top n candidates in an MNL model 
# Displays statistical significance

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
source("functions/utility_functions.R")
source("functions/mnl_function.R")
source("functions/pre_process.R")

#Data (change path as needed)
data <- read.csv("data/Doggerbank_data_973_wide.csv")
#person id
data$person <- rep(1:(dim(data)[1]/6), each=6)

# choice variable
data$choice <- 0
data$choice[data$y1==1] <- 1
data$choice[data$y2==1] <- 2
data$choice[data$y3==1] <- 3

#Choice options for compatibility with given code 
#Since the utility function is defined differently
data$choice1 <- ifelse(data$choice == 1, 1, 0)
data$choice2 <- ifelse(data$choice == 2, 1, 0)
data$choice3 <- ifelse(data$choice == 3, 1, 0)

#pre-processing

#remove columns with more than 50% empty rows
data <- data[, colMeans(is.na(data)) <= 0.5]

#fill country columns using the mode value
country_cols <- c('england', 'scotland', 'wales', 'nireland')

fill_na_with_mode <- function(column) {
  #mode_val <- names(which.max(table(column, useNA = "no")))  # Calculate mode (optional)
  column[is.na(column)] <- 2                        # Replace NA with its own category
  return(column)
}

data[country_cols] <- lapply(data[country_cols], fill_na_with_mode)

#rename columns q22_7, q22_9 for easier recognition of interaction

data <- data %>%
  rename('q227' = 'q22_7', 'q229' = 'q22_9')

#Take a subset of the data 
df_demo <- data[, c('id', 'line',
                  'invasive1', 'cost1', 'spec101', 'spec251', 'prot251', 'prot501',
                  'invasive2','cost2', 'spec102', 'spec252', 'prot252', 'prot502',
                  'invasive3', 'cost3', 'spec103', 'spec253', 'prot253', 'prot503', 'q227', 'q229',
                  'edu', 'male', 'job', 'age', 'choice','protest1', 'protest2', 'protest3', 'y1', 'y2', 'y3',
                  'choice1', 'choice2', 'choice3')]

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

#Columns that vary with response
df_interactions <- df_long[, c('cost', 'spec10', 'spec25', 'prot25', 'prot50',
                               'invasive')]

#Columns that do not vary with response
df_interactions_with <- df_long[, c('male', 'edu', 'job', 'age', 'q227', 'q229')]

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
# elastic_net_model <- cv.glmnet(
#   x = X,
#   y = y,
#   alpha = 0.5,            # Elastic Net (0.5 = balance between Lasso and Ridge)
#   family = "binomial",    # For binary classification
#   nfolds = 5,             # 5-fold cross-validation
#   maxit = 5000,          # High number of iterations for convergence
#   type.measure = "class"  # Classification accuracy
# )
# 
# best_lambda <- elastic_net_model$lambda.min
# cat("Best Lambda:", best_lambda, "\n")

#use the best lambda as obtained to fit the final model
best_lambda <- 0.0005874177

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

#To view the top candidates
#print(sorted_coefficients$feature[1:10])

#Total number of covariates (marginal + interaction)
#length(sorted_coefficients$feature)


########CHANGE AS NEEDED###############
#top n coefficients from Elastic Net
n <- 20
########CHANGE AS NEEDED###############
selected_features <- sorted_coefficients$feature[1:n]
#print(selected_features)

#Use the function create_alt_matrices() on the df_demo subset of data
#Use the wide form
alt_matrices <- create_alt_matrices(df_demo, selected_features)

#utility function matrices
alt1 <- alt_matrices$alt1
alt2 <- alt_matrices$alt2
alt3 <- alt_matrices$alt3


#A required parameter in the code 
nset <- nrow(df_demo)

#initial values of betas
start.values <- rep(0, n)

# # estimate model (without inverting Hessian)
res = maxBFGS(
              #MNL,
              function(coeff) MNL(coeff, alt1, alt2, alt3, final_eval = FALSE),
              grad=NULL,
              hess=NULL,
              start=start.values,
              fixed=NULL,
              print.level=1,
              iterlim=200,
              constraints=NULL,
              tol=1e-25, reltol=1e-25,
              finalHessian=FALSE,
              parscale=rep(1, length=length(start))
              )

invisible(MNL(res$estimate, alt1, alt2, alt3, final_eval = TRUE))

# estimate model (with inverting Hessian)
#new start values are the values obtained above as starting values
start.values = coef(res)

res = maxLik(
            #MNL,
            function(coeff) MNL(coeff, alt1, alt2, alt3),  
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

#Give names to the betas, top n candidates
names(res$estimate) = selected_features

cat(paste(",",res$estimate))

MNL.res = res

summary(MNL.res)

#Bayesian Information Criteria
N <- nrow(df_long)
BIC_value <- -2 * logLik(res) + length(coef(res)) * log(N)
cat("BIC value for elastic-net + MNL: ", BIC_value, "\n")

final_coeff <- res$estimate
threshold <- 1e-4
zero_indices <- which(abs(final_coeff) < threshold)
zero_coeffs <- names(final_coeff)[zero_indices]

cat("===== Coefficients shrunk after Lasso regularization =====\n")

if(length(zero_coeffs) == 0){
  cat("None of the coeffcients were shrunk \n")
} else {
  cat(paste(zero_coeffs, collapse = ", "), "\n")
}
