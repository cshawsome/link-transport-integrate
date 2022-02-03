#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "broom")

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
#neuropsych_gen_cog: "ANMSETOT_norm", "ANBNTTOT_Z", "ANIMMCR_Z", "ANDELCOR_Z", 
# "ANSER7T_Z", "ANAFTOT_Z", "ANRECYES_Z", "ANRECNO_Z", "ANWM1TOT_Z", 
# "ANWM2TOT_Z", "ANBWC20_Z", "ANBWC86_Z", "ANCPTOT_Z", "ANRCPTOT_Z", 
# "ANTMASEC_Z", "SELFCOG_Z", "ANCACTUS", "ANSCISOR",  "ANPRES", "ANSMEM2_Better", 
# "ANSMEM2_Worse", "avg_proxy_cog_Better", "avg_proxy_cog_Worse"
# 
# functional: "Aadla", "Aiadla"
# 
# health: "Abmi_derived_Z", "Astroke", "Adiabe", "Ahearte", "Ahibpe", "Asmoken", 
#   "Amoderate_drinking", "Aheavy_drinking"

#---- average MI datasets ----
avg_ADAMS_imputed <- 
  Reduce("+", ADAMS_imputed_clean)/length(ADAMS_imputed_clean)
avg_ADAMS_imputed <- round(avg_ADAMS_imputed)

#---- models ----
#Because we Z-scored continuous variables, the interpretations for beta are for 
# a 1 SD change in the continuous predictors

#---- **Unimpaired vs. Impaired ----
unimpaired_v_impaired <- glm(Unimpaired ~ AAGE_Z + Black + Hispanic +
                               `Not working` + Retired + ANMSETOT_norm + 
                               ANBNTTOT_Z + ANIMMCR_Z + ANAFTOT_Z + ANRECYES_Z + 
                               ANBWC86_Z,
                             family = "binomial", data = avg_ADAMS_imputed)

#H_0 is that model fits
p_val_unimpaired_v_impaired <- 
  1 - pchisq(unimpaired_v_impaired$deviance, unimpaired_v_impaired$df.residual)

unimpaired_v_impaired_summary <- 
  tidy(unimpaired_v_impaired, exponentiate = TRUE, conf.int = TRUE, 
     conf.level = 0.95) %>% mutate_if(is.numeric, round, 3) %>% as.data.frame()
unimpaired_v_impaired_summary

unimpaired_v_impaired_preds <- unimpaired_v_impaired_summary$term












