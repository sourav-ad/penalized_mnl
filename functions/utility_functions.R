# Function to dynamically create alternative-specific variables with interactions


#To handle e.g like cost_spec10
#it will find cost, spec10 and multiply cost*spec10 for alt1, alt2, alt3
create_alt_matrices <- function(df, selected_features, demographic_vars) {
  alt1_list <- list()
  alt2_list <- list()
  alt3_list <- list()

  #Variables that do not need a suffix i.e are constant for a respondent
  #Change as per requirement
  #demographic_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229') #reduced Dogger bank data

  # demographic_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
  #                     'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
  #                     'job5', 'job6', 'job7', 'job8') #full Dogger bank data

  for (feature in selected_features) {
    if (grepl("_", feature)) {  # Detect interaction terms by looking for "_"
      #This will need renaming columns like q22_7 to q227
      components <- unlist(strsplit(feature, "_"))  # Split into component terms

      #Accounting for if to use suffix or not
      alt1_comp1 <- if (components[1] %in% demographic_vars) df[[components[1]]] else df[[paste0(components[1], "1")]]
      alt1_comp2 <- if (components[2] %in% demographic_vars) df[[components[2]]] else df[[paste0(components[2], "1")]]
      alt2_comp1 <- if (components[1] %in% demographic_vars) df[[components[1]]] else df[[paste0(components[1], "2")]]
      alt2_comp2 <- if (components[2] %in% demographic_vars) df[[components[2]]] else df[[paste0(components[2], "2")]]
      alt3_comp1 <- if (components[1] %in% demographic_vars) df[[components[1]]] else df[[paste0(components[1], "3")]]
      alt3_comp2 <- if (components[2] %in% demographic_vars) df[[components[2]]] else df[[paste0(components[2], "3")]]

      # Construct interaction terms for each alternative
      alt1_list[[feature]] <- alt1_comp1 * alt1_comp2
      alt2_list[[feature]] <- alt2_comp1 * alt2_comp2
      alt3_list[[feature]] <- alt3_comp1 * alt3_comp2
    } else {
      # Directly assign non-interaction terms
      alt1_list[[feature]] <- if (feature %in% demographic_vars) df[[feature]] else df[[paste0(feature, "1")]]
      alt2_list[[feature]] <- if (feature %in% demographic_vars) df[[feature]] else df[[paste0(feature, "2")]]
      alt3_list[[feature]] <- if (feature %in% demographic_vars) df[[feature]] else df[[paste0(feature, "3")]]
    }
  }
  # Combine into matrices
  alt1 <- do.call(cbind, alt1_list)
  alt2 <- do.call(cbind, alt2_list)
  alt3 <- do.call(cbind, alt3_list)
  
  # alt1 <- as.matrix(do.call(cbind, alt1_list)) * 1.0
  # alt2 <- as.matrix(do.call(cbind, alt2_list)) * 1.0
  # alt3 <- as.matrix(do.call(cbind, alt3_list)) * 1.0

  return(list(alt1 = alt1, alt2 = alt2, alt3 = alt3))
}