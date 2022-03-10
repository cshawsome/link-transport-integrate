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
  readRDS(paste0(path_to_box, 
                 "data/ADAMS/prior_data/MI/MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(prior_imputed_clean[[1]]))

#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv")) 

#---- relabel columns ----
prior_imputed_clean <- 
  lapply(prior_imputed_clean, 
         function(x) rename_at(x, vars(variable_labels$ADAMS), ~ 
                                 variable_labels$data_label)) 

#---- predictors ----
#unimpaired model predictors
unimpaired_preds <- selected_vars %>% 
  filter(data_label != "(Intercept)" & Unimpaired != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#other model predictors
other_preds <- selected_vars %>% 
  filter(data_label != "(Intercept)" & Other != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#mci model predictors
mci_preds <- selected_vars %>% 
  filter(data_label != "(Intercept)" & MCI != 0) %>% 
  dplyr::select(data_label) %>% unlist()

#---- model ----
bootstrap_models <- function(prop){
  subsample <- sample_frac(ADAMS_subset, size = prop, replace = TRUE)
  
  unimpaired_model <- 
    glm(formula(paste("AUnimpaired ~ ", 
                      paste(unimpaired_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample)
  
  other_model <- 
    glm(formula(paste("AOther ~ ", paste(other_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(AUnimpaired == 0))
  
  mci_model <- 
    glm(formula(paste("AMCI ~ ", paste(mci_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(AUnimpaired == 0 & AOther == 0)) 
  
  return(list("unimpaired_betas" = coefficients(unimpaired_model),
              "other_betas" = coefficients(other_model), 
              "mci_betas" = coefficients(mci_model), 
              "unimpaired_cov" = as.vector(vcov(unimpaired_model)), 
              "other_cov" = as.vector(vcov(other_model)), 
              "mci_cov" = as.vector(vcov(mci_model))))
}

bootstrap_runs <- replicate(10000, bootstrap_models(prop = 1), simplify = FALSE)

#---- check distributions ----
for(group in c("unimpaired", "other", "mci")){
  data <- lapply(bootstrap_runs, "[[", paste0(group, "_betas")) %>% 
    do.call(rbind, .) %>% t() %>% as.data.frame() 
  
  for(var in rownames(data)){
    show(hist(as.numeric(data[var, ]), main = paste0(group, " ", var)))
  }
}

#---- format output ----
for(est in c("betas", "cov")){
  for(group in c("unimpaired", "other", "mci")){
    data <- lapply(bootstrap_runs, "[[", paste0(group, "_", est)) %>% 
      do.call(rbind, .) %>% t() %>% as.data.frame() 
    
    if(est == "betas"){
      data %<>% mutate("preds" = c("(Intercept)", get(paste0(group, "_preds"))))
    }
    
    data %>% 
      write_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                       "ADAMS_test/latent_class_", group, "_", est, ".csv"))
  }
}
