#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))

#---- **ADAMS imputed data ----
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned")) 

#---- **variable selection results ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv")) %>% 
  dplyr::select("data_label") %>% unlist()

#---- **stack data and add cell IDs ----
#Merge Black, Hispanic, and Stroke (ever/never)
# Ex: 001 is a white participant with a history of stroke
ADAMS_imputed_stacked <- do.call(rbind, ADAMS_imputed_clean) %>% 
  unite("cell_ID", c("Black", "Hispanic", "Astroke"), sep = "", 
        remove = FALSE) %>% 
  rename_at(vars(variable_labels$ADAMS[-1]), ~ variable_labels$data_label[-1])
  
# #Sanity check
# nrow(ADAMS_imputed_stacked) == 25*nrow(ADAMS_imputed_clean[[1]])
# table(ADAMS_imputed_stacked$cell_ID, useNA = "ifany")

#---- parameter estimation ----
continuous_vars <- selected_vars[str_detect(selected_vars, "_Z")] %>% 
  str_remove_all(., "_Z")

normal_parameter_list <- list() 

#---- **Normal distribution ----
for(group in c("Unimpaired", "Other", "MCI")){
  filtered_data <- 
    ADAMS_imputed_stacked[which(ADAMS_imputed_stacked[, group] == 1), ]
  
  #---- ****predictors and outcomes ----
  X <- as.matrix(filtered_data[, c("black", "hispanic", "stroke")] %>% 
                   mutate("Intercept" = 1)) 
  Y <- as.matrix(filtered_data[, continuous_vars])
  
  #---- ****row covariance ----
  XtX_inv <- solve(t(X)%*%X)
  
  #---- ****beta_hat ----
  beta_hat <- XtX_inv%*%t(X)%*%Y
  
  #---- ****sigma_hat ----
  #center matrices and take (matrix)^T(matrix)
  centered_data <- Y - X%*%beta_hat
  sigma_hat <- t(centered_data)%*%centered_data
  
  #---- ****dof ----
  #dof = n + p + 1 - q
  dof = nrow(Y) + ncol(beta_hat) + 1 - nrow(beta_hat)
  
  #---- ****store in list ----
  normal_parameter_list[[group]] <- 
    list("beta_center" = beta_hat, "row_cov" = XtX_inv, 
         "sigma_center" = sigma_hat, "sigma_dof" = dof)
}
#---- **save parameters ----
saveRDS(normal_parameter_list, paste0(path_to_box, "analyses/simulation_study/", 
                                      "continuous_distribution_parameters/", 
                                      "normal_parameters"))
