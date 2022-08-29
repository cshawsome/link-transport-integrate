#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **ADAMS imputed data ----
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/chunk_1/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(ADAMS_imputed_clean[[1]]))

#---- **variable selection results ----
selected_vars <- 
  read_csv(paste0(path_to_box, 
                  "data/variable_selection/model_coefficients.csv")) %>% 
  dplyr::select("data_label") %>% unlist()

#---- **stack data and add cell IDs ----
#Merge Black, Hispanic, and Stroke (ever/never)
# Ex: 001 is a white participant with a history of stroke
ADAMS_imputed_stacked <- do.call(rbind, ADAMS_imputed_clean) %>% 
  unite("cell_ID", c("Black", "Hispanic", "Astroke"), sep = "", 
        remove = FALSE) %>% 
  rename_at(vars(variable_labels$ADAMS), ~ variable_labels$data_label)
  
# #Sanity check
# nrow(ADAMS_imputed_stacked) == 25*nrow(ADAMS_imputed_clean[[1]])
# table(ADAMS_imputed_stacked$cell_ID, useNA = "ifany")

#---- **summary stats ----
#---- ****means/medians continous vars ----
test_subset <- ADAMS_imputed_clean[[1]]
continuous_vars <- 
  colnames(test_subset)[str_detect(colnames(test_subset), "_Z")] 

subset <- test_subset %>% filter(Dementia == 1)

cbind(summarize_at(subset, continuous_vars, .funs = "mean") %>% t() %>% 
  set_colnames("mean"), 
  summarize_at(subset, continuous_vars, .funs = "median") %>% t() %>% 
    set_colnames("median"))

#---- ****race x dem counts ----
test_subset <- ADAMS_imputed_clean[[1]]

table(test_subset$White, test_subset$Black, test_subset$Hispanic, 
      test_subset$Dementia) %>% 
  as.data.frame %>% filter(!Freq == 0) %>% 
  set_colnames(c("White", "Black", "Hispanic", "Dementia", "Freq"))

#---- scaled parameter estimation ----
continuous_vars <- selected_vars[str_detect(selected_vars, "_Z")] 

normal_parameter_list <- list() 

#---- **Normal distribution ----
for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
  filtered_data <- 
    ADAMS_imputed_stacked[which(ADAMS_imputed_stacked[, group] == 1), ]
  
  #---- ****predictors and outcomes ----
  X <- as.matrix(filtered_data[, c("black", "hispanic", "stroke")] %>% 
                   mutate("Intercept" = 1) %>% 
                   dplyr::select("Intercept", everything())) 
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
  #dof = n/25 + p + 1 - q
  # dividing by 25 because sample size is artificially inflated due to stacking
  dof = nrow(Y)/25 + ncol(beta_hat) + 1 - nrow(beta_hat)
  
  #---- ****store in list ----
  normal_parameter_list[[group]] <- 
    list("beta_center" = beta_hat, "row_cov" = XtX_inv, 
         "sigma_center" = sigma_hat, "sigma_dof" = dof)
}
#---- **save parameters ----
saveRDS(normal_parameter_list, paste0(path_to_box, "analyses/simulation_study/", 
                                      "continuous_distribution_parameters/", 
                                      "normal_parameters"))

#---- raw data parameter estimation ----
continuous_vars <- selected_vars[str_detect(selected_vars, "_Z")] %>% 
  str_remove("_Z")

normal_parameter_list <- list() 

#---- **Normal distribution ----
for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
  filtered_data <- 
    ADAMS_imputed_stacked[which(ADAMS_imputed_stacked[, group] == 1), ]
  
  #---- ****predictors and outcomes ----
  X <- as.matrix(filtered_data[, c("black", "hispanic", "stroke")] %>% 
                   mutate("Intercept" = 1) %>% 
                   dplyr::select("Intercept", everything())) 
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
  #dof = n/25 + p + 1 - q
  # dividing by 25 because sample size is artificially inflated due to stacking
  dof = nrow(Y)/25 + ncol(beta_hat) + 1 - nrow(beta_hat)
  
  #---- ****store in list ----
  normal_parameter_list[[group]] <- 
    list("beta_center" = beta_hat, "row_cov" = XtX_inv, 
         "sigma_center" = sigma_hat, "sigma_dof" = dof)
}
#---- **save parameters ----
saveRDS(normal_parameter_list, paste0(path_to_box, "analyses/simulation_study/", 
                                      "continuous_distribution_parameters/", 
                                      "normal_parameters_raw_data"))
