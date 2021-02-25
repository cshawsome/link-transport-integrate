#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "openxlsx", "sjPlot")

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
#try not to lose more than 25% of the sample (about 214 people; n = 642)
#sensitivity: do not lose more than 50% of the sample (about n = 428)
#---- **var missingness ----
sort(colSums(is.na(normal_model_data)))

#---- **model 25 ----
normal_model_25 <- glm(ANormal ~ AAGE + ETHNIC_label +  ANMSETOT + ANSER7T +
                         ANIMMCR + ANRECYES + ANWM1TOT + proxy_cog,
                       family = "binomial", data = normal_model_data)

summary(normal_model_25)

#H_0 is that model fits
p_val_normal_model_25 <- 
  1 - pchisq(normal_model_25$deviance, normal_model_25$df.residual)

tidy_normal_model_25 <- tidy(normal_model_25, exponentiate = TRUE, 
                             conf.int = TRUE, conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_normal_model_25

#---- **model 50 ----
normal_model_50 <- glm(ANormal ~ AAGE + ETHNIC_label +  ANMSETOT + ANSER7T +
                         ANIMMCR + ANRECYES + ANWM1TOT + proxy_cog,
                       family = "binomial", data = normal_model_data)

summary(normal_model_50)

#H_0 is that model fits
p_val_normal_model_50 <- 
  1 - pchisq(normal_model_25$deviance, normal_model_25$df.residual)

tidy_normal_model_50 <- tidy(normal_model_50, exponentiate = TRUE, 
                             conf.int = TRUE, conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_normal_model_50

#---- var select: Other vs. Dementia/MCI ----
other_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(ANormal == 0) %>% dplyr::select(-c("ANormal", "AMCI")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 137 people; n = 412)
#sensitivity: try not to lose more than 50% of the sample (about n = 275)
#---- **var missingness ----
sort(colSums(is.na(other_model_data)))

#---- **model 25 ----
other_model_25 <- glm(AOther ~ AAGE + ANMSETOT + ANIMMCR + ANDELCOR,
                      family = "binomial", data = other_model_data)

summary(other_model_25)

#H_0 is that model fits
p_val_other_model_25 <- 
  1 - pchisq(other_model_25$deviance, other_model_25$df.residual)

tidy_other_model_25 <- tidy(other_model_25, exponentiate = TRUE, 
                            conf.int = TRUE, conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_other_model_25

#---- **model 50 ----
other_model_50 <- glm(AOther ~ AAGE + ANMSETOT + ANIMMCR + ANDELCOR + 
                        proxy_cog,
                      family = "binomial", data = other_model_data)

summary(other_model_50)

#H_0 is that model fits
p_val_other_model_50 <- 
  1 - pchisq(other_model_50$deviance, other_model_50$df.residual)

tidy_other_model_50 <- tidy(other_model_50, exponentiate = TRUE, 
                            conf.int = TRUE, conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_other_model_50

#---- var select: MCI vs. Dementia ----
MCI_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(ANormal == 0 & AOther == 0) %>% dplyr::select(-c("ANormal", "AOther")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 100 people)
#sensitivity: try not to lose more than 50% of the sample (about 200 people)
#---- **var missingness ----
sort(colSums(is.na(MCI_model_data)))

#---- **model 25 ----
MCI_model_25 <- glm(AMCI ~ Aiadla + ANMSETOT + Astroke + Abmi + ANIMMCR, 
                    family = "binomial", data = MCI_model_data)

summary(MCI_model_25)

#H_0 is that model fits
p_val_MCI_model_25 <- 
  1 - pchisq(MCI_model_25$deviance, MCI_model_25$df.residual)

tidy_MCI_model_25 <- tidy(MCI_model_25, exponentiate = TRUE, conf.int = TRUE, 
                          conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_MCI_model_25

#---- **model 50 ----
MCI_model_50 <- glm(AMCI ~ Aiadla + ANMSETOT + Astroke + Abmi + ANIMMCR + 
                      ANWM2TOT, 
                    family = "binomial", data = MCI_model_data)

summary(MCI_model_50)

#H_0 is that model fits
p_val_MCI_model_50 <- 
  1 - pchisq(MCI_model_50$deviance, MCI_model_50$df.residual)

tidy_MCI_model_50 <- tidy(MCI_model_50, exponentiate = TRUE, conf.int = TRUE, 
                          conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_MCI_model_50

#---- sjPlot ----
tab_model(normal_model_25, other_model_25, MCI_model_25, digits = 3, 
          title = "Up to 25% sample dropped", show.loglik = TRUE,
          show.dev = TRUE,
          file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                        "/tables/priors/", 
                        "dem_class_nested_regressions_25.html")) 


tab_model(normal_model_50, other_model_50, MCI_model_50, digits = 3, 
          title = "Up to 50% sample dropped", show.loglik = TRUE,
          show.dev = TRUE,
          file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                        "/tables/priors/", 
                        "dem_class_nested_regressions_50.html")) 

#---- save output .xlsx ----
# table_list <- list("Normal vs. Impaired" = tidy_normal_model, 
#                    "Other vs. MCI or Dementia" = tidy_other_model, 
#                    "MCI vs. Dementia" = tidy_MCI_model)
# write.xlsx(table_list, file = paste0("/Users/CrystalShaw/Box/Dissertation/",
#                                      "/tables/priors/", 
#                                      "dem_class_nested_regressions.xlsx"))
