# Generalized function to create alternative-specific utility matrices with interactions


create_alt_matrices2 <- function(df, selected_features, demographic_vars, n_alt = 3) {
  alt_matrices <- vector("list", n_alt)
  names(alt_matrices) <- paste0("alt", 1:n_alt)

  for (i in 1:n_alt) {
    alt_matrices[[i]] <- list()
  }

  for (feature in selected_features) {

    if (grepl("_", feature)) {
      components <- unlist(strsplit(feature, "_"))

      for (i in 1:n_alt) {
        comp1 <- if (components[1] %in% demographic_vars) df[[components[1]]] else df[[paste0(components[1], i)]]
        comp2 <- if (components[2] %in% demographic_vars) df[[components[2]]] else df[[paste0(components[2], i)]]
        alt_matrices[[i]][[feature]] <- comp1 * comp2
      }

    } else {
      for (i in 1:n_alt) {
        alt_matrices[[i]][[feature]] <- if (feature %in% demographic_vars) df[[feature]] else df[[paste0(feature, i)]]
      }
    }
  }

  for (i in 1:n_alt) {
    alt_matrices[[i]] <- do.call(cbind, alt_matrices[[i]])
  }

  return(alt_matrices)
}
