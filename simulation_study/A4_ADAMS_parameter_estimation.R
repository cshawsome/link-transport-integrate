#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **ADAMS imputed data ----
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned")) 

#---- **variable selection results ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv")) %>% 
  dplyr::select("Variable") %>% unlist()

#---- **stack data and add cell IDs ----
#Merge Black, Hispanic, and Stroke (ever/never)
# Ex: 001 is a white participant with a history of stroke
ADAMS_imputed_stacked <- do.call(rbind, ADAMS_imputed_clean) %>% 
  unite("cell_ID", c("Black", "Hispanic", "Astroke"), sep = "", remove = FALSE)

#Sanity check
nrow(ADAMS_imputed_stacked) == 25*nrow(ADAMS_imputed_clean[[1]])
table(ADAMS_imputed_stacked$cell_ID, useNA = "ifany")

#---- parameter estimation ----
continuous_vars <- selected_vars[str_detect(selected_vars, "_Z")] %>% 
  str_remove_all(., "_Z")

normal_parameter_list <- list() 

#---- **Normal distribution ----
for(cell in cell_IDs){
  filtered_data <- 
    lapply(ADAMS_imputed_clean, function(x) x %>% filter(cell_ID == cell))
  
  #counts are going to differ across datasets because stroke is imputed
  min_count <- lapply(filtered_data, nrow) %>% unlist() %>% min()
  filtered_data <- 
    lapply(filtered_data, function(x) 
      sample_n(x, size = min_count, replace = FALSE) %>% 
        dplyr::select(all_of(continuous_vars)))
  
  #---- **mean matrix ----
  M = Reduce("+", filtered_data)/length(filtered_data)
  
  #---- **row covariance ----
  #independent because independent observations
  U = diag(1, nrow = min_count, ncol = min_count)
  
  #---- **column covariance ----
  #center matrices and take (matrix)^T(matrix)
  centered_data <- lapply(filtered_data, function(x) x - M)
  prod_data <- lapply(centered_data, function(x) 
    t(as.matrix(x)) %*% as.matrix(x))
  V = Reduce("+", prod_data)/(min_count*length(prod_data))
  
  #---- **store in list ----
  normal_parameter_list[[cell]] <- list("M" = M, "U" = U, "V" = V)
}

#---- **save parameters ----
saveRDS(normal_parameter_list, paste0(path_to_box, "analyses/simulation_study/", 
                                      "continuous_distribution_parameters/", 
                                      "normal_parameters"))
