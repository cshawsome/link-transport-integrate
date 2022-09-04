hotdeck <- function(dataset_to_impute, hotdeck_dataset, imputation_mat){
  for(var in rownames(imputation_mat)){
    contingency_cell_vars <- 
      names(imputation_mat[var, which(imputation_mat[var, ] == 1)])
    
    dataset_to_impute %<>% 
      unite(col = "contingency_cell", all_of(contingency_cell_vars), sep = " | ", 
            remove = FALSE)
    
    hotdeck_dataset %<>% 
      unite(col = "contingency_cell", all_of(contingency_cell_vars), sep = " | ", 
            remove = FALSE)
    
    for(cell in unique(dataset_to_impute$contingency_cell)){
      #identify indices for imputation 
      indices <- which(is.na(dataset_to_impute[, var]))[
        which(is.na(dataset_to_impute[, var])) %in% 
          which(dataset_to_impute$contingency_cell == cell)]
      
      if(length(indices) == 0){
        next
      } else{
        #identify hotdeck set for indices
        hot_deck_set <- hotdeck_dataset %>% filter(contingency_cell == cell) %>%
          drop_na(all_of(var))
        
        #track size of hotdeck pool
        dataset_to_impute[indices, paste0(var, "_pool")] <- nrow(hot_deck_set)
        
        #impute values
        dataset_to_impute[indices, c(var, paste0(var, "_cat"))] <- 
          sample_n(hot_deck_set[, c(var, paste0(var, "_cat"))], 
                   size = length(indices), replace = TRUE)
      }
    }
  }
  
  return(dataset_to_impute)
}

#test function
dataset_to_impute <- HCAP
hotdeck_dataset <- HCAP
imputation_mat <- hotdeck_vars_mat
