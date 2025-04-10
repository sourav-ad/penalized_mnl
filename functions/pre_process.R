##preprocessing

#Suggestion: Missing country entries to be their own category as in not_country 0, country 1 and unknown 2

preprocess_data <- function(data, column){
  #remove columns with more than 50% empty rows
  data <- data[, colMeans(is.na(data)) <= 0.5]
  
  #fill country columns using the mode value
  country_cols <- c('country', 'england', 'scotland', 'wales', 'nireland')
  
  fill_na_with_mode <- function(column) {
    mode_val <- names(which.max(table(column, useNA = "no")))  # Calculate mode
    column[is.na(column)] <- mode_val                         # Replace NA with mode
    return(column)
  }
}

#Fill up NA entries with a new category
# 2 if entries 0/1, n+1 if entries 0 to n or 1 to n
#Will eventually go to pre_process.R

fill_na <- function(column){
  #Get the unique values 
  unique_values <- sort(unique(na.omit(column)))
  
  #Checking numerical entries 
  if (is.numeric(unique_values)) {
    min_val <- min(unique_values)
    max_val <- max(unique_values)
    
    if (length(unique_values) == 2 && all(unique_values %in% c(0, 1))) {
      # Binary column (0/1) → Replace NA with 2
      column[is.na(column)] <- 2  
    } else if (all(unique_values %in% seq(min_val, max_val))) {
      # Ordinal column (0:n or 1:n) → Replace NA with n+1
      column[is.na(column)] <- max_val + 1  
    } 
    
  }
  
  return(column)
}