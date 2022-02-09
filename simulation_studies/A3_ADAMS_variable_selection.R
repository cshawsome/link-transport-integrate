#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "broom", "sjPlot", "gridExtra", "magrittr")

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
#neuropsych_gen_cog: "ANMSETOT_norm_Z", "ANBNTTOT_Z", "ANIMMCR_Z", "ANDELCOR_Z", 
# "ANSER7T_Z", "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", 
# "ANWM2TOT_Z", "ANBWC20_Z", "ANBWC86_Z", "ANCPTOT_Z", "ANRCPTOT_Z", 
# "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", "ANSCISOR",  "ANPRES", "ANSMEM2_Better", 
# "ANSMEM2_Worse", "avg_proxy_cog_Better", "avg_proxy_cog_Worse"
# 
#functional: "Aadla_Z", "Aiadla_Z"
# 
#health: "Abmi_derived_Z", "Astroke", "Adiabe", "Ahearte", "Ahibpe", "Asmoken", 
#   "Amoderate_drinking", "Aheavy_drinking"

#---- average MI datasets ----
avg_ADAMS_imputed <- 
  Reduce("+", ADAMS_imputed_clean)/length(ADAMS_imputed_clean)

indicator_vars <- c("Working", "White", "ANSMEM2_Same", "avg_proxy_cog_Same", 
                    "Ano_drinking", "Female", "Married/partnered", "Not working", 
                    "Retired", "Black", "Hispanic", "ANPRES", "ANSCISOR",
                    "ANSMEM2_Better", "ANSMEM2_Worse", "avg_proxy_cog_Better", 
                    "avg_proxy_cog_Worse", "Adiabe", "Ahearte", "Ahibpe", 
                    "Asmoken", "Astroke", "Amoderate_drinking", "Aheavy_drinking", 
                    "Unimpaired", "MCI", "Dementia", "Other", "ANCACTUS")
avg_ADAMS_imputed[, indicator_vars] <- round(avg_ADAMS_imputed[, indicator_vars])

#---- models ----
#Because we Z-scored continuous variables, the interpretations for beta are for 
# a 1 SD change in the continuous predictors

#---- **Unimpaired vs. Impaired ----
unimpaired_v_impaired <- glm(Unimpaired ~ AAGE_Z + Black + Hispanic + 
                               `Not working` + Retired + ANMSETOT_norm_Z + 
                               ANIMMCR_Z + ANAFTOT_Z + ANRECYES_Z + ANRECNO_Z + 
                               SELFCOG_Z + Aadla_Z + Astroke + Adiabe,
                             family = "binomial", data = avg_ADAMS_imputed)

#check sample size (n = 826): Null deviance dof + 1
summary(unimpaired_v_impaired)

#H_0 is that model fits
p_val_unimpaired_v_impaired <- 
  1 - pchisq(unimpaired_v_impaired$deviance, unimpaired_v_impaired$df.residual)

unimpaired_v_impaired_summary <- 
  tidy(unimpaired_v_impaired, exponentiate = TRUE, conf.int = TRUE, 
       conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
unimpaired_v_impaired_summary

unimpaired_preds <- unimpaired_v_impaired_summary$term

#---- **Other vs. MCI or Dementia ----
#conditional on being classified as impaired
other_v_MCI_dem <- glm(Other ~ AAGE_Z + ANMSETOT_norm_Z + ANDELCOR_Z + 
                         ANRECNO_Z, 
                       family = "binomial", 
                       data = avg_ADAMS_imputed %>% filter(Unimpaired == 0))

#check sample size (n = 519): Null deviance dof + 1
summary(other_v_MCI_dem)

#H_0 is that model fits
p_val_other_v_MCI_dem <- 
  1 - pchisq(other_v_MCI_dem$deviance, other_v_MCI_dem$df.residual)

other_v_MCI_dem_summary <- 
  tidy(other_v_MCI_dem, exponentiate = TRUE, conf.int = TRUE, 
       conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
other_v_MCI_dem_summary

other_preds <- other_v_MCI_dem_summary$term

#---- **MCI vs. Dementia ----
#conditional as being classified as being impaired but not having other impairment
MCI_v_dem <- glm(MCI ~ Black + Hispanic + ANMSETOT_norm_Z + ANAFTOT_Z + 
                   SELFCOG_Z + avg_proxy_cog_Better + avg_proxy_cog_Worse + 
                   Astroke + Asmoken, 
                 family = "binomial", 
                 data = avg_ADAMS_imputed %>% 
                   filter(Unimpaired == 0 & Other == 0))

#check sample size (n = 371): Null deviance dof + 1
summary(MCI_v_dem)

#H_0 is that model fits
p_val_MCI_v_dem <- 
  1 - pchisq(MCI_v_dem$deviance, MCI_v_dem$df.residual)

MCI_v_dem_summary <- 
  tidy(MCI_v_dem, exponentiate = TRUE, conf.int = TRUE, 
       conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
MCI_v_dem_summary

MCI_preds <- MCI_v_dem_summary$term

#---- **variable selection table ----
tab_model(unimpaired_v_impaired, other_v_MCI_dem, MCI_v_dem, digits = 3, 
          title = "", show.loglik = TRUE,
          show.dev = TRUE,
          file = paste0(path_to_box, "analyses/variable_selection/",
                        "dem_class_multi_part_models.html"))

#---- distribution of predictors ----
predictor_dists <- function(data, outcome, predictors){
  if(outcome == "Unimpaired"){
    model <- tidy(glm(as.formula(
      paste(outcome, paste0(predictors[-1], collapse = " + "), sep = " ~ ")), 
      family = "binomial", data = data)) 
  } else if(outcome == "Other"){
    model <- tidy(glm(as.formula(
      paste(outcome, paste0(predictors[-1], collapse = " + "), sep = " ~ ")), 
      family = "binomial", data = data %>% filter(Unimpaired == 0)))
  } else{
    model <- tidy(glm(as.formula(
      paste(outcome, paste0(predictors[-1], collapse = " + "), sep = " ~ ")), 
      family = "binomial", data = data %>% filter(Unimpaired == 0 & Other == 0)))
  }
  
  return(model[, c("term", "estimate", "p.value")])
}

unimpaired_pred_dists <- 
  lapply(ADAMS_imputed_clean, predictor_dists, "Unimpaired", 
         unimpaired_preds) %>% do.call(rbind, .) %>% 
  set_colnames(c("Predictor", "beta", "p-value")) %>% 
  mutate("OR" = exp(beta))

other_pred_dists <- 
  lapply(ADAMS_imputed_clean, predictor_dists, "Other", 
         other_preds) %>% do.call(rbind, .) %>% 
  set_colnames(c("Predictor", "beta", "p-value")) %>% 
  mutate("OR" = exp(beta))

MCI_pred_dists <- 
  lapply(ADAMS_imputed_clean, predictor_dists, "MCI", 
         MCI_preds) %>% do.call(rbind, .) %>% 
  set_colnames(c("Predictor", "beta", "p-value")) %>% 
  mutate("OR" = exp(beta))

#---- **plot dists ----
for(model in c("unimpaired", "other", "MCI")){
  data <- get(paste0(model, "_pred_dists"))
  
  pdf(paste0(path_to_box, 
             "analyses/variable_selection/", model, "_dists.pdf"), 
      paper = "letter", height = 10.5, width = 8, onefile = TRUE)
  
  print(ggplot(data = data, aes(x = OR)) + 
          geom_density() + theme_minimal() + 
          ggtitle(paste0(model, " Model ORs")) + 
          facet_wrap(vars(Predictor), scales = "free"))
  
  print(ggplot(data = data, aes(x = `p-value`)) + 
          geom_density() + theme_minimal() + 
          ggtitle(paste0(model, " Model p-values")) + 
          facet_wrap(vars(Predictor), scales = "free"))
  
  dev.off()
}

#---- save estimates ----
#---- **LB preds ----
#lower bound = least dementia
LB_preds <- 
  list("unimpaired" = unimpaired_pred_dists %>% group_by(Predictor) %>% 
         summarise_at(.vars = "beta", max), 
       "other" = other_pred_dists %>% group_by(Predictor) %>% 
         summarise_at(.vars = "beta", max), 
       "MCI" = MCI_pred_dists %>% group_by(Predictor) %>% 
         summarise_at(.vars = "beta", max))

#---- **UB preds ----
#upper bound = most dementia
LB_preds <- 
  list("unimpaired" = unimpaired_pred_dists %>% group_by(Predictor) %>% 
         summarise_at(.vars = "beta", min), 
       "other" = other_pred_dists %>% group_by(Predictor) %>% 
         summarise_at(.vars = "beta", min), 
       "MCI" = MCI_pred_dists %>% group_by(Predictor) %>% 
         summarise_at(.vars = "beta", min))


