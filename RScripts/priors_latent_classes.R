#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "openxlsx")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv")) %>% 
  mutate_if(is.character, as.factor)

#Set factor levels
sapply(ADAMS_subset, class)
ADAMS_subset %<>% 
  mutate("GENDER_label" = fct_relevel(GENDER_label, "Male"), 
         "ETHNIC_label" = fct_relevel(ETHNIC_label, "White"), 
         "EDYRScat_label" = fct_relevel(EDYRScat_label, "No School Completed", 
                                        "1st-8th Grade", "Some High School", 
                                        "High School Diploma", "Some College"), 
         "AAMARRD_label" = fct_relevel(AAMARRD_label, "Married/partnered"), 
         "AACURRWK_label" = fct_relevel(AACURRWK_label, "Retired"), 
         "DRINKcat_label" = fct_relevel(DRINKcat_label, "No Drinking", 
                                        "Moderate Drinking"))

#Sanity check
levels(ADAMS_subset$GENDER_label)
levels(ADAMS_subset$ETHNIC_label)
levels(ADAMS_subset$EDYRScat_label)
levels(ADAMS_subset$AAMARRD_label)
levels(ADAMS_subset$AACURRWK_label)
levels(ADAMS_subset$DRINKcat_label)

#---- dem class indicators ----
ADAMS_subset %<>% 
  #Normal vs. Impaired
  mutate("ANormal" = ifelse(Adem_dx_cat == "Normal", 1, 0),
         #Other vs. Dementia/MCI
         "AOther" = ifelse(Adem_dx_cat == "Other", 1, 0), 
         #MCI vs. Dementia
         "AMCI" = ifelse(Adem_dx_cat == "MCI", 1, 0))

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$ANormal, useNA = "ifany")
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AOther, useNA = "ifany")
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AMCI, useNA = "ifany")

#---- var list ----
vars <- 
  #Demographics
  c("AAGE", "GENDER_label", "ETHNIC_label", "EDYRScat_label", "AAMARRD_label", 
    "AACURRWK_label", 
    #Neuropsych
    "ANMSETOT", "ANBWC20", "ANBWC86", "ANSER7T", "ANSCISOR", "ANCACTUS", 
    "ANPRES", "ANVCPRES", "ANAFTOT", "ANBNTTOT", "ANCPTOT", "ANDCPTOT", 
    "ANIMMCR", "ANDELCOR", "ANRECYES", "ANRECNO", "ANWM1TOT", "ANWM2TOT", 
    "ANTMASEC", "ANTMBSEC", "ANSDMTOT", 
    #Health and health behaviors
    "Astroke", "Ahibpe", "Adiabe", "Ahearte", "Acancre", "Abmi", "Aiadla", 
    "Aadla", "Agrossa", "Afinea", "Acesd", "Asmoken", "DRINKcat_label",
    #Cognition
    "Apstmem", "proxy_cog", 
    #Classes
    "ANormal", "AOther", "AMCI")

#---- var select: Normal vs. Impaired ----
normal_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  dplyr::select(-c("AOther", "AMCI")) 

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 214 people)
sort(colSums(is.na(normal_model_data)))

#---- **model ----
normal_model <- glm(ANormal ~ AAGE + ETHNIC_label + AACURRWK_label + ANMSETOT + 
                      ANSER7T + ANIMMCR + ANRECYES + ANWM1TOT + proxy_cog,
                    family = "binomial", data = normal_model_data)

summary(normal_model)

#H_0 is that model fits
p_val_normal_model <- 
  1 - pchisq(normal_model$deviance, normal_model$df.residual)

tidy_normal_model <- tidy(normal_model, exponentiate = TRUE, conf.int = TRUE, 
                          conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_normal_model

#---- var select: Other vs. Dementia/MCI ----
other_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(ANormal == 0) %>% dplyr::select(-c("ANormal", "AMCI")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 137 people)
sort(colSums(is.na(other_model_data)))

#---- **model ----
other_model <- glm(AOther ~ AAGE + GENDER_label + EDYRScat_label + ANMSETOT + 
                     ANIMMCR + ANDELCOR + proxy_cog,
                   family = "binomial", data = other_model_data)

summary(other_model)

#H_0 is that model fits
p_val_other_model <- 
  1 - pchisq(other_model$deviance, other_model$df.residual)

tidy_other_model <- tidy(other_model, exponentiate = TRUE, conf.int = TRUE, 
                         conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_other_model

#---- var select: MCI vs. Dementia ----
MCI_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(ANormal == 0 & AOther == 0) %>% dplyr::select(-c("ANormal", "AOther")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 100 people)
sort(colSums(is.na(MCI_model_data)))

#---- **model ----
MCI_model <- glm(AMCI ~ AAGE + ETHNIC_label + EDYRScat_label + Aiadla + 
                   ANMSETOT + Astroke + Abmi + ANIMMCR + proxy_cog,
                 family = "binomial", data = MCI_model_data)

summary(MCI_model)

#H_0 is that model fits
p_val_MCI_model <- 
  1 - pchisq(MCI_model$deviance, MCI_model$df.residual)

tidy_MCI_model <- tidy(MCI_model, exponentiate = TRUE, conf.int = TRUE, 
                       conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_MCI_model

#---- save output ----
table_list <- list("Normal vs. Impaired" = tidy_normal_model, 
                   "Other vs. MCI or Dementia" = tidy_other_model, 
                   "MCI vs. Dementia" = tidy_MCI_model)
write.xlsx(table_list, file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                                     "/tables/priors/", 
                                     "dem_class_nested_regressions.xlsx"))
