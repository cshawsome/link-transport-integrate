#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **prior imputed clean ----
prior_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/prior_data/MI/MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric)) 

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(prior_imputed_clean[[1]]))

#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/model_coefficients.csv")) 

#---- relabel columns ----
prior_imputed_clean <- 
  lapply(prior_imputed_clean, 
         function(x) rename_at(x, vars(variable_labels$ADAMS), ~ 
                                 variable_labels$data_label)) 

#---- predictors ----
#unimpaired model predictors
unimpaired_preds <- selected_vars %>% 
  filter(data_label != "Intercept" & Unimpaired != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#other model predictors
other_preds <- selected_vars %>% 
  filter(data_label != "Intercept" & Other != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#mci model predictors
mci_preds <- selected_vars %>% 
  filter(data_label != "Intercept" & MCI != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#---- models ----
model_function <- function(data, unimpaired_pred, other_preds, mci_preds){
  unimpaired_model <- 
    glm(formula(paste("Unimpaired ~ ", 
                      paste(unimpaired_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = data)
  
  other_model <- 
    glm(formula(paste("Other ~ ", paste(other_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = data %>% 
          filter(Unimpaired == 0))
  
  mci_model <- 
    glm(formula(paste("MCI ~ ", paste(mci_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = data %>% 
          filter(Unimpaired == 0 & Other == 0)) 
  
  return(list("unimpaired_betas" = coefficients(unimpaired_model),
              "other_betas" = coefficients(other_model), 
              "mci_betas" = coefficients(mci_model), 
              "unimpaired_cov" = vcov(unimpaired_model), 
              "other_cov" = vcov(other_model), 
              "mci_cov" = vcov(mci_model)))
}

estimates <- 
  lapply(prior_imputed_clean, model_function, unimpaired_preds, other_preds, 
         mci_preds) 

# #---- check distributions ----
# for(group in c("unimpaired", "other", "mci")){
#   data <- lapply(estimates, "[[", paste0(group, "_betas")) %>%
#     do.call(rbind, .) %>% t() %>% as.data.frame()
# 
#   for(var in rownames(data)){
#     show(hist(as.numeric(data[var, ]), main = paste0(group, " ", var)))
#   }
# }

#---- format output ----
for(est in c("betas", "cov")){
  for(group in c("unimpaired", "other", "mci")){
    data <- lapply(estimates, "[[", paste0(group, "_", est)) 
    
    if(est == "betas"){
      data %<>% 
        do.call(rbind, .) %>% t() %>% as.data.frame()
      
      data %<>% set_colnames(seq(1, ncol(data))) %>% 
        mutate("preds" = c("(Intercept)", get(paste0(group, "_preds"))))
      
      write_csv(data, paste0(path_to_box, "data/prior_data/latent_class_", group, 
                             "_", est, ".csv"))
    } else{
      saveRDS(data, paste0(path_to_box, "data/prior_data/latent_class_", group, 
                           "_", est))
    }
  }
}
