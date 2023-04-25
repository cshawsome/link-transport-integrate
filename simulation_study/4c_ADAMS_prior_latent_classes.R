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

# JZ: check dimension
length(prior_imputed_clean)
dim(prior_imputed_clean[[1]])

#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/model_coefficients.csv")) 

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

#dementia model predictors
dementia_preds <- selected_vars %>% 
  filter(data_label != "Intercept" & Dementia != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#---- models ----
# JZ -- for my own understanding: 
# fit the glm model on each imputed copy of the dataset (10000 copies in total)
# for each Dx (4 Dx's) based on the variables determined in the 
# variable selection process, and then collect coef and vcov to use as prior

model_function <- 
  function(data, unimpaired_preds, other_preds, mci_preds, dementia_preds){
    for(group in c("Unimpaired", "Other", "MCI", "Dementia")){
      assign(paste0(tolower(group), "_model"), 
             glm(formula(paste(group, " ~ ", 
                               paste(get(paste0(tolower(group), "_preds")), 
                                     collapse = " + "), collapse = "")), 
                 family = "binomial", data = data))
    }
    
    return(list("unimpaired_betas" = coefficients(unimpaired_model),
                "other_betas" = coefficients(other_model), 
                "mci_betas" = coefficients(mci_model), 
                "dementia_betas" = coefficients(dementia_model),
                "unimpaired_cov" = vcov(unimpaired_model), 
                "other_cov" = vcov(other_model), 
                "mci_cov" = vcov(mci_model), 
                "dementia_cov" = vcov(dementia_model)))
  }

estimates <- 
  lapply(prior_imputed_clean, model_function, unimpaired_preds, other_preds, 
         mci_preds, dementia_preds) 

# #---- check distributions ----
# for(group in c("unimpaired", "other", "mci", "dementia")){
#   data <- lapply(estimates, "[[", paste0(group, "_betas")) %>%
#     do.call(rbind, .) %>% t() %>% as.data.frame()
# 
#   for(var in rownames(data)){
#     show(hist(as.numeric(data[var, ]), main = paste0(group, " ", var)))
#   }
# }

#---- format output ----
for(est in c("betas", "cov")){
  for(group in c("unimpaired", "other", "mci", "dementia")){
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
