generate_synthetic_neuropsych <- function(dataset, neuropsych_dataset, 
                                          sample_size, path_to_results){
  #---- bootstrap dataset ----
  superpop <- sample_n(dataset, size = sample_size, replace = TRUE) %>% 
    mutate("HHIDPN" = seq(1, sample_size))
  
  #---- hot decking ----
  contingency_cell_vars <- c("female","black", "hispanic", "age_cat", 
                             "edyrs_cat") 
  
  hotdeck_vars <- c("mmse_norm", "delrc", "animal_naming", "wrc_yes", "wrc_no", 
                    "imm_story", "del_story", "imm_cp", "del_cp", "trailsA")
  
  for(var in hotdeck_vars){
    superpop %<>% unite(col = "contingency_cell", all_of(contingency_cell_vars), 
                        sep = " | ", remove = FALSE)
    
    neuropsych_dataset %<>% 
      unite(col = "contingency_cell", all_of(contingency_cell_vars), sep = " | ", 
            remove = FALSE)
    
    for(cell in unique(superpop$contingency_cell)){
      #identify indices for imputation in superpop
      indices <- which(superpop$contingency_cell == cell)
      
      #identify hotdeck set for indices
      hot_deck_set <- neuropsych_dataset %>% filter(contingency_cell == cell)
      
      #track size of hotdeck pool
      superpop[indices, paste0(var, "_pool")] <- nrow(hot_deck_set)
      
      #impute values
      superpop[indices, c(var, paste0(var, "_cat"))] <- 
        sample_n(hot_deck_set[, c(var, paste0(var, "_cat"))], 
                 size = length(indices), replace = TRUE)
    }
    
    #use imputed var in next contingency cell iteration
    contingency_cell_vars <- c(contingency_cell_vars, paste0(var, "_cat"))
  }
  
  #---- save results ----
  write_csv(superpop, paste0(path_to_results, "superpop_", sample_size, ".csv"))
  
}

#test function
dataset <- HRS_analytic
neuropsych_dataset <- HCAP_analytic
sample_size = 10000
path_to_results <- paste0(path_to_box, "data/superpopulations/")
