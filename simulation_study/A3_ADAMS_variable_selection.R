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

#---- variable list ----
#variables and groups in order of priority for model inclusion 
# (based on conceptual model and prioritizing variables with least missingness
# in unimputed datasets)
#sociodemographics: "AAGE_Z", "Black", "Hispanic", "Female",  "EDYRS", 
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

var_list <- c("AAGE_Z", "Black", "Hispanic", "Female",  "EDYRS", "Not working", 
              "Retired", "Married/partnered", "ANMSETOT_norm_Z", "ANBNTTOT_Z", 
              "ANIMMCR_Z", "ANDELCOR_Z", "ANSER7T_Z", "ANAFTOT_Z", "ANRECYES_Z", 
              "ANRECNO_Z", "ANWM1TOT_Z","ANWM2TOT_Z", "ANBWC20_Z", "ANBWC86_Z", 
              "ANCPTOT_Z", "ANRCPTOT_Z", "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", 
              "ANSCISOR",  "ANPRES", "ANSMEM2_Better", "ANSMEM2_Worse", 
              "Aadla_Z", "Aiadla_Z", "Abmi_derived_Z", "Astroke", "Adiabe", 
              "Ahearte", "Ahibpe", "Asmoken", "Amoderate_drinking", 
              "Aheavy_drinking")

#---- stack data and add weights ----
ADAMS_imputed_stacked <- do.call(rbind, ADAMS_imputed_clean) %>% 
  mutate("weights" = 1/nrow(ADAMS_imputed_stacked))

# #Sanity check
# nrow(ADAMS_imputed_stacked) == 25*nrow(ADAMS_imputed_clean[[1]])
# 1/nrow(ADAMS_imputed_stacked)
# table(ADAMS_imputed_stacked$weights, useNA = "ifany")

#---- variable selection ----
#---- **best lambda function ----
#choose the best lambda for lasso regression based on 1000 runs of 10-fold CV
best_lambda <- function(x, y, weights){
  cv_model <- cv.glmnet(x, y, alpha = 1, weights = weights)
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
    weights <- data.matrix(model_data[, "weights"])
    
    #---- choose best lambda ----
    lambda <- median(replicate(1000, best_lambda(x, y, weights)))
    
    #---- store results ----
    model_list[[class]] <- 
      glmnet(x, y, alpha = 1, lambda = lambda, weights = weights)
    
    lambda_vec[class] <- lambda
  }
  return(list("models" = model_list, "lambdas" = lambda_vec))
}

start <- Sys.time()
variable_selection <- lasso_reg(ADAMS_imputed_stacked, var_list)
end <- Sys.time() - start

#---- **list predictors ----
selected_vars <- 
  matrix(nrow = (length(var_list) + 1), ncol = (1 + 3*25)) %>% 
  as.data.frame() %>%
  set_colnames(c("Variable", 
                 str_replace(apply(expand_grid(c("Unimpaired", "Other", "MCI"), 
                                               seq(1, 25)) %>% 
                                     set_colnames(c("model", "number")) %>% 
                                     arrange(number), 1, paste0, 
                                   collapse = "_"), " ", ""))) %>% 
  mutate("Variable" = c("(Intercept)", var_list))


for(m in 1:25){
  for(model in c("Unimpaired", "Other", "MCI")){
    selected_vars[, paste0(model, "_", m)] <- 
      as.vector(coef(variable_selection[[m]]$models[[model]]))
  }  
}

#remove variables that are never selected
selected_vars %<>% mutate("times_selected" = rowSums(selected_vars[, -1])) %>% 
  filter(times_selected != 0)

#---- **save results ----
saveRDS(variable_selection, paste0(path_to_box, "analyses/variable_selection/", 
                                   "ADAMS_lasso_models"))

write_csv(selected_vars, paste0(path_to_box, "analyses/variable_selection/", 
                                "model_coefficients.csv"))

#---- OLD ----
# #---- models ----
# 
# 
# #---- **Unimpaired vs. Impaired ----
# unimpaired_v_impaired <- glm(Unimpaired ~ AAGE_Z + Black + Hispanic + 
#                                `Not working` + Retired + ANMSETOT_norm_Z + 
#                                ANIMMCR_Z + ANAFTOT_Z + ANRECYES_Z + ANRECNO_Z + 
#                                SELFCOG_Z + Aadla_Z + Astroke + Adiabe,
#                              family = "binomial", data = avg_ADAMS_imputed)
# 
# #check sample size (n = 826): Null deviance dof + 1
# summary(unimpaired_v_impaired)
# 
# #H_0 is that model fits
# p_val_unimpaired_v_impaired <- 
#   1 - pchisq(unimpaired_v_impaired$deviance, unimpaired_v_impaired$df.residual)
# 
# unimpaired_v_impaired_summary <- 
#   tidy(unimpaired_v_impaired, exponentiate = TRUE, conf.int = TRUE, 
#        conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
# unimpaired_v_impaired_summary
# 
# unimpaired_preds <- unimpaired_v_impaired_summary$term
# 
# #---- **Other vs. MCI or Dementia ----
# #conditional on being classified as impaired
# other_v_MCI_dem <- glm(Other ~ AAGE_Z + ANMSETOT_norm_Z + ANDELCOR_Z + 
#                          ANRECNO_Z, 
#                        family = "binomial", 
#                        data = avg_ADAMS_imputed %>% filter(Unimpaired == 0))
# 
# #check sample size (n = 519): Null deviance dof + 1
# summary(other_v_MCI_dem)
# 
# #H_0 is that model fits
# p_val_other_v_MCI_dem <- 
#   1 - pchisq(other_v_MCI_dem$deviance, other_v_MCI_dem$df.residual)
# 
# other_v_MCI_dem_summary <- 
#   tidy(other_v_MCI_dem, exponentiate = TRUE, conf.int = TRUE, 
#        conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
# other_v_MCI_dem_summary
# 
# other_preds <- other_v_MCI_dem_summary$term
# 
# #---- **MCI vs. Dementia ----
# #conditional as being classified as being impaired but not having other impairment
# MCI_v_dem <- glm(MCI ~ Black + Hispanic + ANMSETOT_norm_Z + ANAFTOT_Z + 
#                    SELFCOG_Z + Astroke + Asmoken, 
#                  family = "binomial", 
#                  data = avg_ADAMS_imputed %>% 
#                    filter(Unimpaired == 0 & Other == 0))
# 
# #check sample size (n = 371): Null deviance dof + 1
# summary(MCI_v_dem)
# 
# #H_0 is that model fits
# p_val_MCI_v_dem <- 
#   1 - pchisq(MCI_v_dem$deviance, MCI_v_dem$df.residual)
# 
# MCI_v_dem_summary <- 
#   tidy(MCI_v_dem, exponentiate = TRUE, conf.int = TRUE, 
#        conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
# MCI_v_dem_summary
# 
# MCI_preds <- MCI_v_dem_summary$term
# 
# #---- **variable selection table ----
# tab_model(unimpaired_v_impaired, other_v_MCI_dem, MCI_v_dem, digits = 3, 
#           title = "", show.loglik = TRUE,
#           show.dev = TRUE,
#           file = paste0(path_to_box, "analyses/variable_selection/",
#                         "dem_class_multi_part_models.html"))
# 
# #---- distribution of predictors ----
# predictor_dists <- function(data, outcome, predictors){
#   if(outcome == "Unimpaired"){
#     model <- tidy(glm(as.formula(
#       paste(outcome, paste0(predictors[-1], collapse = " + "), sep = " ~ ")), 
#       family = "binomial", data = data)) 
#   } else if(outcome == "Other"){
#     model <- tidy(glm(as.formula(
#       paste(outcome, paste0(predictors[-1], collapse = " + "), sep = " ~ ")), 
#       family = "binomial", data = data %>% filter(Unimpaired == 0)))
#   } else{
#     model <- tidy(glm(as.formula(
#       paste(outcome, paste0(predictors[-1], collapse = " + "), sep = " ~ ")), 
#       family = "binomial", data = data %>% filter(Unimpaired == 0 & Other == 0)))
#   }
#   
#   return(model[, c("term", "estimate", "p.value")])
# }
# 
# unimpaired_pred_dists <- 
#   lapply(ADAMS_imputed_clean, predictor_dists, "Unimpaired", 
#          unimpaired_preds) %>% do.call(rbind, .) %>% 
#   set_colnames(c("Predictor", "beta", "p-value")) %>% 
#   mutate("OR" = exp(beta))
# 
# other_pred_dists <- 
#   lapply(ADAMS_imputed_clean, predictor_dists, "Other", 
#          other_preds) %>% do.call(rbind, .) %>% 
#   set_colnames(c("Predictor", "beta", "p-value")) %>% 
#   mutate("OR" = exp(beta))
# 
# MCI_pred_dists <- 
#   lapply(ADAMS_imputed_clean, predictor_dists, "MCI", 
#          MCI_preds) %>% do.call(rbind, .) %>% 
#   set_colnames(c("Predictor", "beta", "p-value")) %>% 
#   mutate("OR" = exp(beta))
# 
# #---- **plot dists ----
# for(model in c("unimpaired", "other", "MCI")){
#   data <- get(paste0(model, "_pred_dists"))
#   
#   pdf(paste0(path_to_box, 
#              "analyses/variable_selection/", model, "_dists.pdf"), 
#       paper = "letter", height = 10.5, width = 8, onefile = TRUE)
#   
#   print(ggplot(data = data, aes(x = OR)) + 
#           geom_density() + theme_minimal() + 
#           ggtitle(paste0(model, " Model ORs")) + 
#           facet_wrap(vars(Predictor), scales = "free"))
#   
#   print(ggplot(data = data, aes(x = `p-value`)) + 
#           geom_density() + theme_minimal() + 
#           ggtitle(paste0(model, " Model p-values")) + 
#           facet_wrap(vars(Predictor), scales = "free"))
#   
#   dev.off()
# }
# 
# #---- contingency cell counts ----
# #full table
# table(avg_ADAMS_imputed$Black, avg_ADAMS_imputed$Hispanic, 
#       avg_ADAMS_imputed$`Not working`, avg_ADAMS_imputed$Retired,
#       avg_ADAMS_imputed$Astroke, avg_ADAMS_imputed$Adiabe, 
#       avg_ADAMS_imputed$Asmoken) %>% 
#   as.data.frame() %>% 
#   set_colnames(c("Black", "Hispanic", "Not Working", "Retired", "Stroke", 
#                  "Diabetes", "Smoker", "Count")) %>% 
#   filter(!(Black == 1 & Hispanic == 1)) %>% 
#   filter(!(`Not Working` == 1 & Retired == 1)) %>%
#   write_csv(paste0(path_to_box, 
#                    "analyses/variable_selection/full_contingency.csv"))
# 
# #simplified table
# table(avg_ADAMS_imputed$Black, avg_ADAMS_imputed$Hispanic, 
#       avg_ADAMS_imputed$Astroke) %>% 
#   as.data.frame() %>% 
#   set_colnames(c("Black", "Hispanic", "Stroke", "Count")) %>% 
#   filter(!(Black == 1 & Hispanic == 1)) %>% 
#   write_csv(paste0(path_to_box, 
#                    "analyses/variable_selection/simplified_contingency.csv"))
# 
# #---- save estimates ----
# #---- **LB preds ----
# #lower bound = least dementia
# LB_preds <- 
#   list("unimpaired" = unimpaired_pred_dists %>% group_by(Predictor) %>% 
#          summarise_at(.vars = "beta", max), 
#        "other" = other_pred_dists %>% group_by(Predictor) %>% 
#          summarise_at(.vars = "beta", max), 
#        "MCI" = MCI_pred_dists %>% group_by(Predictor) %>% 
#          summarise_at(.vars = "beta", max))
# 
# #---- **UB preds ----
# #upper bound = most dementia
# UB_preds <- 
#   list("unimpaired" = unimpaired_pred_dists %>% group_by(Predictor) %>% 
#          summarise_at(.vars = "beta", min), 
#        "other" = other_pred_dists %>% group_by(Predictor) %>% 
#          summarise_at(.vars = "beta", min), 
#        "MCI" = MCI_pred_dists %>% group_by(Predictor) %>% 
#          summarise_at(.vars = "beta", min))
# 
# saveRDS(LB_preds, paste0(path_to_box, "analyses/variable_selection/LB_preds"))
# saveRDS(UB_preds, paste0(path_to_box, "analyses/variable_selection/UB_preds"))
