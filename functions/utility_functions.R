# Function to dynamically create alternative-specific variables with interactions
#Generalized code for any number of utility function is provided below

#To handle e.g like cost_spec10
#it will find cost, spec10 and multiply cost*spec10 for alt1, alt2, alt3
create_alt_matrices <- function(df, selected_features) {
  alt1_list <- list()
  alt2_list <- list()
  alt3_list <- list()

  #Variables that do not need a suffix i.e are constant for a respondent
  #Change as per requirement
  constant_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229') #reduced Dogger bank data

  # constant_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
  #                    'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
  #                    'job5', 'job6', 'job7', 'job8') #full Dogger bank data

  for (feature in selected_features) {
    if (grepl("_", feature)) {  # Detect interaction terms by looking for "_"
      #This will need renaming columns like q22_7 to q227
      components <- unlist(strsplit(feature, "_"))  # Split into component terms

      #Accounting for if to use suffix or not
      alt1_comp1 <- if (components[1] %in% constant_vars) df[[components[1]]] else df[[paste0(components[1], "1")]]
      alt1_comp2 <- if (components[2] %in% constant_vars) df[[components[2]]] else df[[paste0(components[2], "1")]]
      alt2_comp1 <- if (components[1] %in% constant_vars) df[[components[1]]] else df[[paste0(components[1], "2")]]
      alt2_comp2 <- if (components[2] %in% constant_vars) df[[components[2]]] else df[[paste0(components[2], "2")]]
      alt3_comp1 <- if (components[1] %in% constant_vars) df[[components[1]]] else df[[paste0(components[1], "3")]]
      alt3_comp2 <- if (components[2] %in% constant_vars) df[[components[2]]] else df[[paste0(components[2], "3")]]

      # Construct interaction terms for each alternative
      alt1_list[[feature]] <- alt1_comp1 * alt1_comp2
      alt2_list[[feature]] <- alt2_comp1 * alt2_comp2
      alt3_list[[feature]] <- alt3_comp1 * alt3_comp2
    } else {
      # Directly assign non-interaction terms
      alt1_list[[feature]] <- if (feature %in% constant_vars) df[[feature]] else df[[paste0(feature, "1")]]
      alt2_list[[feature]] <- if (feature %in% constant_vars) df[[feature]] else df[[paste0(feature, "2")]]
      alt3_list[[feature]] <- if (feature %in% constant_vars) df[[feature]] else df[[paste0(feature, "3")]]
    }
  }
  # Combine into matrices
  alt1 <- do.call(cbind, alt1_list)
  alt2 <- do.call(cbind, alt2_list)
  alt3 <- do.call(cbind, alt3_list)

  return(list(alt1 = alt1, alt2 = alt2, alt3 = alt3))
}




# # Generalized function to create alternative-specific utility matrices with interactions
# create_alt_matrices <- function(df, selected_features, n_alt) {
#   # Store output 
#   alt_matrices <- vector("list", n_alt)
#   names(alt_matrices) <- paste0("alt", 1:n_alt)
#   
#   # Variables that are respondent-specific (don't vary by alternative)
#   constant_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2', 
#                      'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4', 
#                      'job5', 'job6', 'job7', 'job8')  # Adjust as needed
#   
#   # Create an empty list of lists
#   #one sublist for each alternative
#   for (i in 1:n_alt) {
#     alt_matrices[[i]] <- list()
#   }
#   
#   for (feature in selected_features) {
#     
#     if (grepl("_", feature)) {
#       components <- unlist(strsplit(feature, "_"))
#       
#       for (i in 1:n_alt) {
#         comp1 <- if (components[1] %in% constant_vars) df[[components[1]]] else df[[paste0(components[1], i)]]
#         comp2 <- if (components[2] %in% constant_vars) df[[components[2]]] else df[[paste0(components[2], i)]]
#         alt_matrices[[i]][[feature]] <- comp1 * comp2
#       }
#       
#     } else {
#       for (i in 1:n_alt) {
#         alt_matrices[[i]][[feature]] <- if (feature %in% constant_vars) df[[feature]] else df[[paste0(feature, i)]]
#       }
#     }
#   }
#   
#   # Collapse each alternative's list of columns into a matrix
#   for (i in 1:n_alt) {
#     alt_matrices[[i]] <- do.call(cbind, alt_matrices[[i]])
#   }
#   
#   return(alt_matrices)
# }



#Helper function to take in integer
num_variables <- function(prompt) {
  cat(prompt)
  as.integer(readline())
}

name_variables <- function(prompt) {
  cat(prompt)
  vars <- strsplit(readline(), ",")[[1]]
  vars <- trimws(vars)  # remove extra spaces
  return(vars)
}
