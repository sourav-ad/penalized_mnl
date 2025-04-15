# Generalized function to create alternative-specific utility matrices with interactions


create_alt_matrices <- function(df, selected_features, n_alt) {
  # Store output
  alt_matrices <- vector("list", n_alt)
  names(alt_matrices) <- paste0("alt", 1:n_alt)
  
  # Variables that are respondent-specific (don't vary by alternative)
  constant_vars <- c('male', 'edu', 'job', 'age', 'q227', 'q229', 'q1', 'q2',
                     'q6', 'q7', 'q10', 'job1', 'job2', 'job3', 'job4',
                     'job5', 'job6', 'job7', 'job8')  # Adjust as needed
  
  # Create an empty list of lists
  #one sublist for each alternative
  for (i in 1:n_alt) {
    alt_matrices[[i]] <- list()
  }
  
  for (feature in selected_features) {
    
    if (grepl("_", feature)) {
      components <- unlist(strsplit(feature, "_"))
      
      for (i in 1:n_alt) {
        comp1 <- if (components[1] %in% constant_vars) df[[components[1]]] else df[[paste0(components[1], i)]]
        comp2 <- if (components[2] %in% constant_vars) df[[components[2]]] else df[[paste0(components[2], i)]]
        alt_matrices[[i]][[feature]] <- comp1 * comp2
      }
      
    } else {
      for (i in 1:n_alt) {
        alt_matrices[[i]][[feature]] <- if (feature %in% constant_vars) df[[feature]] else df[[paste0(feature, i)]]
      }
    }
  }
  
  # Collapse each alternative's list of columns into a matrix
  for (i in 1:n_alt) {
    alt_matrices[[i]] <- do.call(cbind, alt_matrices[[i]])
  }
  
  return(alt_matrices)
}