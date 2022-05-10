#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "broom", "sjPlot", "gridExtra", "magrittr", "glmnet")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- variable list ----
#variables and groups in order of priority for model inclusion 
# (based on conceptual model and prioritizing variables with least missingness
# in unimputed datasets)
#sociodemographics: "AAGE_Z", "Black", "Hispanic", "Female",  "EDYRS_Z", 
# "Not working", "Retired", "Married/partnered"
# 
# neuropsych_gen_cog: "ANMSETOT_norm_Z", "ANBNTTOT_Z", "ANIMMCR_Z", "ANDELCOR_Z",
# "ANSER7T_Z", "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z",
# "ANWM2TOT_Z", "ANBWC20_Z", "ANBWC86_Z", "ANCPTOT_Z", "ANRCPTOT_Z",
# "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", "ANSCISOR",  "ANPRES", "ANSMEM2_Better",
# "ANSMEM2_Worse"
# 
# functional: "Aadla_Z", "Aiadla_Z"
# 
# health: "Abmi_derived_Z", "Astroke", "Adiabe", "Ahearte", "Ahibpe", "Asmoken",
#   "Amoderate_drinking", "Aheavy_drinking"

#Because we Z-scored continuous variables, the interpretations for beta are for 
# a 1 SD change in the continuous predictors

var_list <- c("AAGE_Z", "Black", "Hispanic", "Female",  "EDYRS_Z", "Not working", 
              "Retired", "Married/partnered", "ANMSETOT_norm_Z", "ANIMMCR_Z", 
              "ANDELCOR_Z", "ANSER7T_Z", "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", 
              "ANWM1TOT_Z","ANWM2TOT_Z", "ANBWC20_Z", "ANBWC86_Z", "ANCPTOT_Z", 
              "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", "ANSCISOR",  
              "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse", "Aadla_Z", "Aiadla_Z", 
              "Abmi_derived_Z", "Astroke", "Adiabe", "Ahearte", "Ahibpe", 
              "Asmoken", "Amoderate_drinking", "Aheavy_drinking")

#---- stack data and add weights ----
ADAMS_imputed_stacked <- do.call(rbind, ADAMS_imputed_clean) %>% 
  mutate("weights" = 1/length(ADAMS_imputed_clean))

# #Sanity check
# nrow(ADAMS_imputed_stacked) == 25*nrow(ADAMS_imputed_clean[[1]])
# 1/length(ADAMS_imputed_clean)
# table(ADAMS_imputed_stacked$weights, useNA = "ifany")

#---- variable selection ----
#---- **best lambda function ----
#choose the best lambda for lasso regression based on 1000 runs of 10-fold CV
best_lambda <- function(x, y, weights){
  cv_model <- cv.glmnet(x, y, alpha = 1, weights = weights, family = "binomial")
  min_lambda <- cv_model$lambda.min
  return(min_lambda)
}

#---- **lasso regression function ----
lasso_reg <- function(data, var_list){
  #---- preallocate for results ----
  model_list <- list()
  lambda_vec <- vector(length = 3) %>% 
    set_names(c("Unimpaired", "Other", "MCI"))
  
  for(class in c("Unimpaired", "Other", "MCI")){
    if(class == "Unimpaired"){
      model_data <- data
    } else if(class == "Other"){
      model_data <- data %>% filter(Unimpaired == 0)
    } else{
      model_data <- data %>% filter(Unimpaired == 0 & Other == 0)
    }
    
    #---- define outcome and predictor variables ----
    y <- model_data[, class]
    x <- data.matrix(model_data[, var_list])
    weights <- model_data[, "weights"]
    
    #---- choose best lambda ----
    lambda <- median(replicate(1000, best_lambda(x, y, weights)))
    
    #---- store results ----
    model_list[[class]] <- 
      glmnet(x, y, alpha = 1, lambda = lambda, weights = weights, 
             family = "binomial")
    
    lambda_vec[class] <- lambda
  }
  return(list("models" = model_list, "lambdas" = lambda_vec))
}

#About 3.5 hours
start <- Sys.time()
variable_selection <- lasso_reg(ADAMS_imputed_stacked, var_list)
end <- Sys.time() - start

#---- **list predictors ----
selected_vars <- 
  matrix(nrow = (length(var_list) + 1), ncol = (1 + 3)) %>% 
  as.data.frame() %>%
  set_colnames(c("Variable", "Unimpaired", "Other", "MCI")) %>% 
  mutate("Variable" = c("(Intercept)", var_list))

for(model in c("Unimpaired", "Other", "MCI")){
  selected_vars[, model] <- 
    as.vector(coef(variable_selection$models[[model]]))
}  

#---- **truncate decimals ----
#truncate at 2 decimal places
selected_vars[, 2:4] <- trunc(selected_vars[, 2:4]*10^2)/10^2

#remove variables that are never selected
selected_vars %<>% mutate("times_selected" = rowSums(selected_vars[, -1])) %>% 
  filter(times_selected != 0)

#----**relabel vars ----
selected_vars %<>% left_join(., variable_labels, by = c("Variable" = "ADAMS")) 

#---- **save results ----
saveRDS(variable_selection, paste0(path_to_box, "analyses/simulation_study/", 
                                   "variable_selection/ADAMS_lasso_models"))

write_csv(selected_vars %>% 
            dplyr::select(c("data_label", "Unimpaired", "Other", "MCI")), 
          paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                 "model_coefficients.csv"))
