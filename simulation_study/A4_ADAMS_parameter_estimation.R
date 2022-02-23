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
for(group in c("Unimpaired", "Other", "MCI")){
  filtered_data <- 
    ADAMS_imputed_stacked[which(ADAMS_imputed_stacked[, group] == 1), ]
  
  #---- ****predictors and outcomes ----
  X <- as.matrix(filtered_data[, c("Black", "Hispanic", "Astroke")] %>% 
                   mutate("Intercept" = 1)) 
  Y <- as.matrix(filtered_data[, continuous_vars])
  
  #---- ****beta_hat ----
  beta_hat <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  #---- ****row covariance ----
  U = diag(1, nrow = min_count, ncol = min_count)
  
  #---- ****column covariance ----
  #center matrices and take (matrix)^T(matrix)
  centered_data <- lapply(filtered_data, function(x) x - M)
  prod_data <- lapply(centered_data, function(x) 
    t(as.matrix(x)) %*% as.matrix(x))
  V = Reduce("+", prod_data)/(min_count*length(prod_data))
  
  #---- ****store in list ----
  normal_parameter_list[[group]] <- 
    list("beta_center" = beta_hat, "row_cov" = XtX_inv, 
         "sigma_center" = sigma_hat, "sigma_dof" = dof)
}
#---- **save parameters ----
saveRDS(normal_parameter_list, paste0(path_to_box, "analyses/simulation_study/", 
                                      "continuous_distribution_parameters/", 
                                      "normal_parameters"))
