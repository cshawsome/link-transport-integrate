#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "sjPlot")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/ADAMS/cleaned/ADAMS_subset_mixed.csv"), 
           col_types = cols("AYEAR" = col_character(), 
                            "Astroke" = col_character(), 
                            "Ahibpe" = col_character(), 
                            "Adiabe" = col_character(), 
                            "Ahearte" = col_character(), 
                            "Acancre" = col_character(), 
                            "Asmoken" = col_character(), 
                            "White" = col_character(), 
                            "Black" = col_character(), 
                            "Hispanic" = col_character())) %>% 
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

# #Sanity check
# levels(ADAMS_subset$GENDER_label)
# levels(ADAMS_subset$ETHNIC_label)
# levels(ADAMS_subset$EDYRScat_label)
# levels(ADAMS_subset$AAMARRD_label)
# levels(ADAMS_subset$AACURRWK_label)
# levels(ADAMS_subset$DRINKcat_label)

#---- Z-score ----
ADAMS_subset %<>% mutate_if(is.numeric, scale)

#---- dem class indicators ----
ADAMS_subset %<>% 
  #Unimpaired vs. Impaired
  mutate("AUnimpaired" = ifelse(Adem_dx_cat == "Unimpaired", 1, 0),
         #Other vs. Dementia/MCI
         "AOther" = ifelse(Adem_dx_cat == "Other", 1, 0), 
         #MCI vs. Dementia
         "AMCI" = ifelse(Adem_dx_cat == "MCI", 1, 0))

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AUnimpaired, useNA = "ifany")
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AOther, useNA = "ifany")
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AMCI, useNA = "ifany")

#---- var list ----
vars <- 
  #Demographics
  c("AAGE", "GENDER_label", "ETHNIC_label", "EDYRScat_label", "AAMARRD_label", 
    "AACURRWK_label", 
    #Neuropsych
    "ANMSETOT", "ANMSETOT_norm", "ANBWC20", "ANBWC86", "ANSER7T", "ANSCISOR", 
    "ANCACTUS", "ANPRES", "ANVCPRES", "ANAFTOT", "ANBNTTOT", "ANCPTOT", 
    "ANDCPTOT", "ANIMMCR", "ANDELCOR", "ANRECYES", "ANRECNO", "ANWM1TOT", 
    "ANWM2TOT", "ANTMASEC", "ANTMBSEC", "ANSDMTOT", 
    #Health and health behaviors
    "Astroke", "Ahibpe", "Adiabe", "Ahearte", "Acancre", "Abmi", "Aiadla", 
    "Aadla", "Agrossa", "Afinea", "Acesd", "Asmoken", "DRINKcat_label",
    #Cognition
    "Apstmem", "proxy_cog", 
    #Classes
    "AUnimpaired", "AOther", "AMCI")

#---- NORMED MMSE ----
#---- var select: Unimpaired vs. Impaired ----
unimpaired_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  dplyr::select(-c("AOther", "AMCI")) 

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 214 people; n = 642)
#sensitivity: do not lose more than 50% of the sample (about n = 428)
#---- **var missingness ----
sort(colSums(is.na(unimpaired_model_data) %>% 
               set_colnames(colnames(unimpaired_model_data))))

#---- **model 25 ----
unimpaired_model_25 <- glm(AUnimpaired ~ AAGE + ETHNIC_label + ANMSETOT_norm + 
                             ANSER7T + ANIMMCR + ANRECYES + ANWM1TOT + proxy_cog, 
                           family = "binomial", data = unimpaired_model_data)

summary(unimpaired_model_25)

#H_0 is that model fits
p_val_unimpaired_model_25 <- 
  1 - pchisq(unimpaired_model_25$deviance, unimpaired_model_25$df.residual)

tidy_unimpaired_model_25 <- tidy(unimpaired_model_25, exponentiate = TRUE, 
                                 conf.int = TRUE, conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_unimpaired_model_25

# #---- ****Sanity checking ORs for race (race only) ----
# unimpaired_model_25_check1 <- glm(AUnimpaired ~ AAGE + ETHNIC_label,
#                               family = "binomial", data = unimpaired_model_data)
# 
# summary(unimpaired_model_25_check1)
# 
# #H_0 is that model fits
# p_val_unimpaired_model_25_check1 <- 
#   1 - pchisq(unimpaired_model_25_check1$deviance, 
#              unimpaired_model_25_check1$df.residual)
# 
# tidy_unimpaired_model_25_check1 <- tidy(unimpaired_model_25_check1, exponentiate = TRUE, 
#                                     conf.int = TRUE, conf.level = 0.95) %>% 
#   mutate_if(is.numeric, round, 4) %>% as.data.frame()
# #show results
# tidy_unimpaired_model_25_check1
# 
# #---- ****Sanity checking ORs for race (race + MMSE) ----
# unimpaired_model_25_check2 <- glm(AUnimpaired ~ AAGE + ETHNIC_label + ANMSETOT,
#                               family = "binomial", data = unimpaired_model_data)
# 
# summary(unimpaired_model_25_check2)
# 
# #H_0 is that model fits
# p_val_unimpaired_model_25_check2 <- 
#   1 - pchisq(unimpaired_model_25_check2$deviance, 
#              unimpaired_model_25_check2$df.residual)
# 
# tidy_unimpaired_model_25_check2 <- tidy(unimpaired_model_25_check2, exponentiate = TRUE, 
#                                     conf.int = TRUE, conf.level = 0.95) %>% 
#   mutate_if(is.numeric, round, 4) %>% as.data.frame()
# #show results
# tidy_unimpaired_model_25_check2

# #---- sjPlot ----
# tab_model(unimpaired_model_25_check1, unimpaired_model_25_check2, digits = 3, 
#           title = "Checking Race/Ethnicity --> Impairment Classification", 
#           show.loglik = TRUE, show.dev = TRUE,
#           file = paste0("/Users/CrystalShaw/Box/Dissertation/",
#                         "/tables/priors/", 
#                         "dem_class_nested_regressions_25_check.html")) 

#---- var select: Other vs. Dementia/MCI ----
other_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(AUnimpaired == 0) %>% dplyr::select(-c("AUnimpaired", "AMCI")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 137 people; n = 412)
#---- **var missingness ----
sort(colSums(is.na(other_model_data)) %>% set_names(colnames(other_model_data)))

#---- **model 25 ----
other_model_25 <- glm(AOther ~ AAGE + ANMSETOT_norm + ANDELCOR + ANIMMCR, 
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

#---- var select: MCI vs. Dementia ----
MCI_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(AUnimpaired == 0 & AOther == 0) %>% dplyr::select(-c("AUnimpaired", "AOther")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 100 people)
#---- **var missingness ----
sort(colSums(is.na(MCI_model_data)) %>% set_names(colnames(MCI_model_data)))

#---- **model 25 ----
MCI_model_25 <- glm(AMCI ~ Aiadla + ANMSETOT_norm + Astroke + Abmi + ANIMMCR, 
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

#---- sjPlot ----
tab_model(unimpaired_model_25, other_model_25, MCI_model_25, digits = 3, 
          title = "Up to 25% sample dropped", show.loglik = TRUE,
          show.dev = TRUE,
          file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                        "/tables/priors/", 
                        "dem_class_nested_regressions_25.html"))

#---- select important predictors ----
ADAMS_subset %<>% mutate("(Intercept)" = 1)
vars <- c("HHIDPN", c(names(coefficients(unimpaired_model_25)), 
                      names(coefficients(other_model_25)), 
                      names(coefficients(MCI_model_25))), 
          "White", "AUnimpaired", "AOther", "AMCI", "Adem_dx_cat", "ETHNIC_label")
vars[which(vars == "ETHNIC_labelBlack")] <- "Black"
vars[which(vars == "ETHNIC_labelHispanic")] <- "Hispanic"
vars[which(vars == "Astroke1")] <- "Astroke"

ADAMS_subset %<>% dplyr::select(all_of(vars)) %>% na.omit()

#---- split sample ----
set.seed(20210623)
train <- sample_frac(ADAMS_subset, size = 0.70, replace = FALSE)
test <- ADAMS_subset %>% filter(!HHIDPN %in% train$HHIDPN)

write_csv(train, file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                               "/data/ADAMS/cleaned/ADAMS_train.csv"))
write_csv(test, file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                              "/data/ADAMS/cleaned/ADAMS_test.csv"))

#---- OLD ----
#---- RAW MMSE ----
#---- var select: Unimpaired vs. Impaired ----
unimpaired_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  dplyr::select(-c("AOther", "AMCI")) 

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 214 people; n = 642)
#sensitivity: do not lose more than 50% of the sample (about n = 428)
#---- **var missingness ----
sort(colSums(is.na(unimpaired_model_data) %>% 
               set_colnames(colnames(unimpaired_model_data))))

#---- **model 25 ----
unimpaired_model_25 <- glm(AUnimpaired ~ AAGE + ETHNIC_label + ANMSETOT + 
                             ANSER7T + ANIMMCR + ANRECYES + ANWM1TOT + proxy_cog, 
                       family = "binomial", data = unimpaired_model_data)

summary(unimpaired_model_25)

#H_0 is that model fits
p_val_unimpaired_model_25 <- 
  1 - pchisq(unimpaired_model_25$deviance, unimpaired_model_25$df.residual)

tidy_unimpaired_model_25 <- tidy(unimpaired_model_25, exponentiate = TRUE, 
                             conf.int = TRUE, conf.level = 0.95) %>% 
  mutate_if(is.numeric, round, 4) %>% as.data.frame()
#show results
tidy_unimpaired_model_25

# #---- ****Sanity checking ORs for race (race only) ----
# unimpaired_model_25_check1 <- glm(AUnimpaired ~ AAGE + ETHNIC_label,
#                               family = "binomial", data = unimpaired_model_data)
# 
# summary(unimpaired_model_25_check1)
# 
# #H_0 is that model fits
# p_val_unimpaired_model_25_check1 <- 
#   1 - pchisq(unimpaired_model_25_check1$deviance, 
#              unimpaired_model_25_check1$df.residual)
# 
# tidy_unimpaired_model_25_check1 <- tidy(unimpaired_model_25_check1, exponentiate = TRUE, 
#                                     conf.int = TRUE, conf.level = 0.95) %>% 
#   mutate_if(is.numeric, round, 4) %>% as.data.frame()
# #show results
# tidy_unimpaired_model_25_check1
# 
# #---- ****Sanity checking ORs for race (race + MMSE) ----
# unimpaired_model_25_check2 <- glm(AUnimpaired ~ AAGE + ETHNIC_label + ANMSETOT,
#                               family = "binomial", data = unimpaired_model_data)
# 
# summary(unimpaired_model_25_check2)
# 
# #H_0 is that model fits
# p_val_unimpaired_model_25_check2 <- 
#   1 - pchisq(unimpaired_model_25_check2$deviance, 
#              unimpaired_model_25_check2$df.residual)
# 
# tidy_unimpaired_model_25_check2 <- tidy(unimpaired_model_25_check2, exponentiate = TRUE, 
#                                     conf.int = TRUE, conf.level = 0.95) %>% 
#   mutate_if(is.numeric, round, 4) %>% as.data.frame()
# #show results
# tidy_unimpaired_model_25_check2

# #---- sjPlot ----
# tab_model(unimpaired_model_25_check1, unimpaired_model_25_check2, digits = 3, 
#           title = "Checking Race/Ethnicity --> Impairment Classification", 
#           show.loglik = TRUE, show.dev = TRUE,
#           file = paste0("/Users/CrystalShaw/Box/Dissertation/",
#                         "/tables/priors/", 
#                         "dem_class_nested_regressions_25_check.html")) 

#---- var select: Other vs. Dementia/MCI ----
other_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(AUnimpaired == 0) %>% dplyr::select(-c("AUnimpaired", "AMCI")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 137 people; n = 412)
#---- **var missingness ----
sort(colSums(is.na(other_model_data)) %>% set_names(colnames(other_model_data)))

#---- **model 25 ----
other_model_25 <- glm(AOther ~ AAGE + ANMSETOT + ANDELCOR + ANIMMCR, 
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

#---- var select: MCI vs. Dementia ----
MCI_model_data <- ADAMS_subset %>% dplyr::select(all_of(vars)) %>% 
  filter(AUnimpaired == 0 & AOther == 0) %>% 
  dplyr::select(-c("AUnimpaired", "AOther")) 

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat)

#check missingness-- add variables to model by level of missingness
#try not to lose more than 25% of the sample (about 100 people)
#---- **var missingness ----
sort(colSums(is.na(MCI_model_data)) %>% set_names(colnames(MCI_model_data)))

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

#---- sjPlot ----
tab_model(unimpaired_model_25, other_model_25, MCI_model_25, digits = 3, 
          title = "Up to 25% sample dropped", show.loglik = TRUE,
          show.dev = TRUE,
          file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                        "/tables/priors/", 
                        "dem_class_nested_regressions_25.html"))

#---- select important predictors ----
ADAMS_subset %<>% mutate("(Intercept)" = 1)
vars <- c("HHIDPN", c(names(coefficients(unimpaired_model_25)), 
                      names(coefficients(other_model_25)), 
                      names(coefficients(MCI_model_25))), 
          "White", "AUnimpaired", "AOther", "AMCI", "Adem_dx_cat", 
          "ETHNIC_label")
vars[which(vars == "ETHNIC_labelBlack")] <- "Black"
vars[which(vars == "ETHNIC_labelHispanic")] <- "Hispanic"
vars[which(vars == "Astroke1")] <- "Astroke"

ADAMS_subset %<>% dplyr::select(all_of(vars)) %>% na.omit()




