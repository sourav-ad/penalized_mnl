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

source("functions/penalized_mnl_model_general.R")
#data <- read.csv("data/plastics-survey-uk-us.csv")
data <- read.csv("data/GoCoase3_peanlised.csv")
names(data) <- gsub("_", "", names(data)) #remove underscores in column names
#data <- read.csv("data/doggerbank_full_regularization.csv")

#Example GoCoase
data_schema <- list(
  #identify each respondent/individual
  id_col = "id",
  #identify the choice task/question faced by the respondent
  task_col = "idII",
  #column that records which alternative was chosen
  choice_col = "choice",

  #Do not include ASC variables
  #choice specific
  #needs to be mentioned without any underscores 
  alt_specific_vars = c(
    "sandnu",
    "deichnu",
    "bepfnu"
    # "ufernu","rucknu","steilnu"
    )
  ,

  #respondent specific
  demographic_vars = c(
    "genderf",
    "balticseaSH"
  ),

  #number of alternatives in each choice task
  n_alt = 3,

  #Specify the ASC columns
  asc_specs = c("ASC1", "ASCsq"),

  #separator between the base variable name and the alternative number
  alt_sep = ""
)


 
# Example plastics
# data_schema <- list(
#   #identify each respondent/individual
#   id_col = "person",
#   #identify the choice task/question faced by the respondent
#   task_col = "set",
#   #column that records which alternative was chosen
#   choice_col = "choice",
# 
#   #Do not include ASC variables
#   #choice specific
#   alt_specific_vars = c(
#     "beach90alt",
#     "foreign20alt",
#     "intl15alt",
#     "coastal15alt"
#   )
#   ,
# 
#   #respondent specific
#   demographic_vars = c(
#     "household",
#     "sectoroilgas",
#     "environorg",
#     "usavisit"
#   ),
# 
#   #number of alternatives in each choice task
#   n_alt = 3,
# 
#   #Specify the ASC columns
#   #asc_specs = c("ASC1", "ASCsq"),
# 
#   #separator between the base variable name and the alternative number
#   alt_sep = ""
# )



#Doggerbank example
# data_schema <- list(
#   id_col = "id",
#   task_col = "line",
#   
#   # Dogger Bank already has y1, y2, y3
#   choice_indicator_cols = c("y1", "y2", "y3"),
#   
#   alt_specific_vars = c(
#     "cost",
#     "spec10",
#     "spec25"
#   ),
#   
#   demographic_vars = c(
#     "male",
#     "edu",
#     "job1"
#   ),
#   
#   n_alt = 3,
#   
#   # Dogger Bank already has:
#   # ASC21 ASC22 ASC23
#   # ASC31 ASC32 ASC33
#   asc_specs = c("ASC2", "ASC3"),
#   
#   alt_sep = ""
# )

#number of times an experiment is repeated
nrep <- as.integer(unique(table(data[[data_schema$id_col]]))) #extracted from data

model <- run_penalized_mnl(
  data = data,
  data_schema = data_schema,
  lambda_grid = exp(seq(log(6e-4), log(3e-2), length.out = 10)), #should be chosen appropriately
  n_folds = 5,
  alpha = 0.5, #can be adjusted
  threshold = 0.01, #for determining "shrunk" criteria
  optimizer = "BFGS",
  method = "MAXLIK",
  n_workers = 25 #change as per system capabilities
  #n_workers = parallelly::availableCores() #uses all the cores available in the system
)

model$coefficients