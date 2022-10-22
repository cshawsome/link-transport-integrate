#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "broom", "sjPlot", "gridExtra", "magrittr", "glmnet")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/chunk_1/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))

# #---- add interactions ----
# add_interactions <- function(dataset, cog_vars){
#   for(var in cog_vars){
#     #race interaction
#     dataset %<>% mutate(!!sym(paste0("black*", var)) := Black*!!sym(var))
#     dataset %<>% mutate(!!sym(paste0("hispanic*", var)) := Hispanic*!!sym(var))
#     
#     #edyrs interaction
#     dataset %<>% mutate(!!sym(paste0("edyrs_Z*", var)) := EDYRS_Z*!!sym(var))
#   }
#   
#   return(dataset)
# }
# 
# cog_vars <-
#   c("ANMSETOT_norm_Z", "ANIMMCR_Z", "ANDELCOR_Z", "ANSER7T_Z", "ANAFTOT_Z", 
#     "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", "ANWM2TOT_Z", "ANBWC20", 
#     "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", "ANSCISOR",  
#     "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse")
# 
# ADAMS_imputed_clean %<>% lapply(., function(x) add_interactions(x, cog_vars))

#---- variable list ----
#variables are grouped in order of priority for model inclusion 
# (based on conceptual model and prioritizing variables with least missingness
# in unimputed datasets)
#sociodemographics: "AAGE_Z", "Black", "Hispanic", "Female",  "EDYRS_Z", 
# "Not working", "Retired", "Married/partnered"
# 
# neuropsych_gen_cog: "ANMSETOT_norm_Z", "ANIMMCR_Z", "ANDELCOR_Z", "ANSER7T_Z", 
# "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", "ANWM2TOT_Z", 
# "ANBWC20", "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", 
# "ANSCISOR",  "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse"
# 
# functional: "Aadla_Z", "Aiadla_Z"
# 
# health: "Abmi_derived_Z", "Astroke", "Adiabe", "Ahearte", "Ahibpe", "Asmoken",
#   "Amoderate_drinking", "Aheavy_drinking"
# 
# race x cognitive variables interactions: black x 
#  ("ANMSETOT_norm_Z", "ANIMMCR_Z", "ANDELCOR_Z", "ANSER7T_Z", 
#   "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", "ANWM2TOT_Z", 
#   "ANBWC20", "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", 
#   "ANSCISOR",  "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse"), 
#   hispanic x 
#  ("ANMSETOT_norm_Z", "ANIMMCR_Z", "ANDELCOR_Z", "ANSER7T_Z", 
#   "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", "ANWM2TOT_Z", 
#   "ANBWC20", "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", 
#   "ANSCISOR",  "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse")
# 
# education x cognitive variables interactions: EDYRS_Z x 
#  ("ANMSETOT_norm_Z", "ANIMMCR_Z", "ANDELCOR_Z", "ANSER7T_Z", 
#   "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", "ANWM2TOT_Z", 
#   "ANBWC20", "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", 
#   "ANSCISOR",  "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse")

#Because we Z-scored continuous variables, the interpretations for beta are for 
# a 1 SD change in the continuous predictors

var_list <- 
  c("AAGE_Z", "Black", "Hispanic", "Female",  "EDYRS_Z", "Not working", "Retired", 
    "Married/partnered", "ANMSETOT_norm_Z", "ANIMMCR_Z", "ANDELCOR_Z", 
    "ANSER7T_Z", "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z",
    "ANWM2TOT_Z", "ANBWC20", "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", 
    "ANCACTUS", "ANSCISOR", "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse", 
    "Aadla_Z", "Aiadla_Z", "Abmi_derived_Z", "Astroke", "Adiabe", "Ahearte", 
    "Ahibpe", "Asmoken", "Amoderate_drinking", "Aheavy_drinking") 
    
    # paste0("black*", cog_vars), paste0("hispanic*", cog_vars), 
    # paste0("edyrs_Z*", cog_vars))

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
  lambda_vec <- vector(length = 4) %>%
    set_names(c("Unimpaired", "Other", "MCI", "Dementia"))
  
  for(class in c("Unimpaired", "Other", "MCI", "Dementia")){
    #Try unconditional models
    model_data <- data
    
    # #Impaired models conditional on being impaired
    # if(class == "Unimpaired"){
    #   model_data <- data
    # } else {
    #   model_data <- data %>% filter(Unimpaired == 0)
    # } 
    
    # #Conditional models
    # if(class == "Unimpaired"){
    #   model_data <- data
    # } else if(class == "Other"){
    #   model_data <- data %>% filter(Unimpaired == 0)
    # } else{
    #   model_data <- data %>% filter(Unimpaired == 0 & Other == 0)
    # }
    
    #---- define outcome and predictor variables ----
    y <- model_data[, class]
    x <- data.matrix(model_data[, var_list])
    weights <- model_data[, "weights"]
    
    #---- choose best lambda ----
    lambda <- median(replicate(1, best_lambda(x, y, weights)))
    
    #---- store results ----
    model_list[[class]] <- 
      glmnet(x, y, alpha = 1, weights = weights, lambda = lambda, 
             family = "binomial")
    
    lambda_vec[class] <- lambda
  }
  return(list("models" = model_list, "lambdas" = lambda_vec))
}

#About 1.5 mins
start <- Sys.time()
variable_selection <- lasso_reg(ADAMS_imputed_stacked, var_list)
end <- Sys.time() - start

#---- **list predictors ----
selected_vars <- 
  matrix(nrow = (length(var_list) + 1), ncol = (1 + 4)) %>% 
  as.data.frame() %>%
  set_colnames(c("Variable", "Unimpaired", "Other", "MCI", "Dementia")) %>% 
  mutate("Variable" = c("(Intercept)", var_list))

for(model in c("Unimpaired", "Other", "MCI", "Dementia")){
  selected_vars[, model] <- 
    as.vector(coef(variable_selection$models[[model]]))
}  

#---- **truncate decimals ----
#truncate at 2 decimal places
selected_vars[, 2:5] <- trunc(selected_vars[, 2:5]*10^2)/10^2

#remove variables that are never selected
selected_vars %<>% mutate("times_selected" = rowSums(selected_vars[, -1])) %>% 
  filter(times_selected != 0)

#----**relabel vars ----
selected_vars %<>% left_join(., variable_labels, by = c("Variable" = "ADAMS")) 

#---- **save results ----
saveRDS(variable_selection, paste0(path_to_box, "data/variable_selection/", 
                                   "ADAMS_lasso_models"))

write_csv(selected_vars %>% 
            dplyr::select(c("data_label", "Unimpaired", "Other", "MCI", "Dementia")), 
          paste0(path_to_box, "data/variable_selection/model_coefficients.csv"))

# #---- test ----
# lasso_models <- readRDS(paste0(path_to_box, "data/variable_selection/", 
#                                "ADAMS_lasso_models"))
