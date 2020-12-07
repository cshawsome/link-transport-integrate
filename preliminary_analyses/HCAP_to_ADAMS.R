#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- color palette ----
green = "#66a182"

#---- import data ----
#---- **neuropsych ----
#ADAMS
neuropsych_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
neuropsych_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                 "ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych_A <- read_da_dct(neuropsych_data_path_A, 
                                  neuropsych_dict_path_A, HHIDPN = "TRUE")

#HCAP
HCAP_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16da/HC16HP_R.da")
HCAP_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16sta/HC16HP_R.dct")

HCAP <- read_da_dct(HCAP_data_path, HCAP_dict_path, HHIDPN = "TRUE")
cog_test_labels <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                   "cog_test_meaningful_labels.csv"))

#---- **sociodemographic ----
#ADAMS
ADAMS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                  "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                  "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") 

#HCAP
HRS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                "tracker/trk2018_3/TRK2018TR_R.da")
HRS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                "tracker/trk2018_3/TRK2018TR_R.dct")

HRS_tracker <- read_da_dct(HCAP_tracker_data_path, HCAP_tracker_dict_path, 
                           HHIDPN = "TRUE") 

#---- join data ----
HCAP <- left_join(HCAP, HRS_tracker, by = "HHIDPN")
ADAMS <- left_join(ADAMS_neuropsych_A, ADAMS_tracker, by = "HHIDPN")

#---- select variables ----
HCAP_vars <- c("HHIDPN", "GENDER", "HISPANIC", "RACE", "SCHLYRS",
               #for age calculation
               "BIRTHYR", "H1RIWYEAR",
               #MMSE
               "H1RMSESCORE")

ADAMS_vars <- c("HHIDPN", 
                #MMSE
                "ANMSETOT")

HCAP_subset <- HCAP %>% dplyr::select(all_of(HCAP_vars))
ADAMS_subset <- ADAMS %>% dplyr::select(all_of(ADAMS_vars))

# #Variable check
# colSums(is.na(HCAP_subset))

#Drop people without an interview year
HCAP_subset %<>% filter(!is.na(H1RIWYEAR))
#Drop people with missing MMSE data in ADAMS
ADAMS_subset %<>% filter(ANMSETOT <= 30)

#---- format data ----
#---- **age ----
HCAP_subset %<>% mutate("HCAP_age" = H1RIWYEAR - BIRTHYR)

# #Sanity check
# ggplot(data = HCAP_subset) + 
#   geom_histogram(aes(HCAP_age), color = green, fill = green) + 
#   theme_minimal() + xlab("Age") + theme(text = element_text(size = 14))
# ggsave(filename = "HCAP_age.jpeg", plot = last_plot(), device = "jpeg", 
#        path = paste0("/Users/CrystalShaw/Box/Dissertation/", 
#                      "preliminary_analyses/", 
#                      "HCAP_synthetic/figures/"), width = 7, height = 5, 
#        units = "in")
# 
# summary(HCAP_subset$HCAP_age)

#Drop original variables
HCAP_subset %<>% dplyr::select(-c("H1RIWYEAR", "BIRTHYR"))

# #Sanity check
# colnames(HCAP_subset)

#---- **race/ethnicity ----
# #Variable check
# #0 = "Not Obtained"; 1-3 = diff categories Hispanic; 5 = "Non-Hispanic"
# table(HCAP_subset$HISPANIC) 
# #0 = "Not Obtained"; 1 = "White"; 2 = "Black"; 7 = "Other"
# table(HCAP_subset$RACE)
# table(HCAP_subset$RACE, HCAP_subset$HISPANIC)

#Drop missing race
HCAP_subset %<>% filter(RACE != 0)

HCAP_subset %<>% 
  mutate("race_ethnic_cat" = 
           case_when(RACE == 1 & HISPANIC == 5 ~ "Non-hispanic White", 
                     HISPANIC %in% c(1, 2, 3) ~ "Hispanic", 
                     RACE == 2 & (HISPANIC == 5 | HISPANIC == 0) ~ 
                       "Non-hispanic Black", 
                     RACE == 7 & HISPANIC == 5 ~ "Other"))

# #Sanity check
# table(HCAP_subset$race_ethnic_cat, useNA = "ifany")
# table(HCAP_subset$race_ethnic_cat, useNA = "ifany")/nrow(HCAP_subset)

#Drop original variables
HCAP_subset %<>% dplyr::select(-c("RACE", "HISPANIC"))

#Sanity check
colnames(HCAP_subset)

#---- **sex/gender ----
# #Variable check
# #1 = "Male"; 2 = "Female"
# table(HCAP_subset$GENDER)

HCAP_subset %<>% mutate("female" = ifelse(GENDER == 2, 1, 0))

# #Sanity check
# table(HCAP_subset$female, useNA = "ifany")
# table(HCAP_subset$female, useNA = "ifany")/nrow(HCAP_subset)

#Drop original variable
HCAP_subset %<>% dplyr::select(-c("GENDER"))

# #Sanity check
# colnames(HCAP_subset)

#---- education distributions ----
# #Variable Check
# #0-17 = years of education; 99 = "Not Ascertained"
# table(HCAP_subset$SCHLYRS)

#Remove missing education
HCAP_subset %<>% filter(SCHLYRS != 99)

ggplot(data = HCAP_subset) + 
  geom_histogram(aes(SCHLYRS), color = green, fill = green) +
  theme_minimal() + xlab("Years of Education") + 
  theme(text = element_text(size = 14))

ggsave(filename = "HCAP_edyrs.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

# #Variable check
# summary(HCAP_subset$SCHLYRS)

#---- neuropsych distributions ----
#HCAP
ggplot(data = HCAP_subset) + 
  geom_histogram(aes(H1RMSESCORE), color = green, fill = green) +
  theme_minimal() + xlab("MMSE Score") + theme(text = element_text(size = 14))

ggsave(filename = "HCAP_MMSE_score.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

#ADAMS
ggplot(data = ADAMS_subset) + 
  geom_histogram(aes(ANMSETOT), color = green, fill = green) +
  theme_minimal() + xlab("MMSE Score") + theme(text = element_text(size = 14))

ggsave(filename = "ADAMS_MMSE_score.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

# #Variable check
# summary(HCAP_subset$H1RMSESCORE)
summary(ADAMS_subset$ANMSETOT)

#---- modeling HCAP data ----
#Formatting race/ethnicity factors
HCAP_subset %<>% mutate_at("race_ethnic_cat", as.factor) %>% 
  mutate_at("race_ethnic_cat", function(x) fct_relevel(x, "Non-hispanic White"))

#---- **main effects models ----
MMSE_main_effects <- lm(H1RMSESCORE ~ HCAP_age + female + 
                          race_ethnic_cat + 
                          SCHLYRS, data = HCAP_subset)

summary(MMSE_main_effects)
plot(MMSE_main_effects$residuals)
hist(MMSE_main_effects$fitted.values)
summary(MMSE_main_effects$fitted.values)

#---- **2-way interactions ----
MMSE_2way_intx <- lm(H1RMSESCORE ~ HCAP_age*female + HCAP_age*race_ethnic_cat + 
                       female*race_ethnic_cat + female*SCHLYRS + 
                       race_ethnic_cat*SCHLYRS, 
                     data = HCAP_subset)

summary(MMSE_2way_intx)
tidy(MMSE_2way_intx, conf.int = TRUE)
plot(MMSE_2way_intx$residuals)
hist(MMSE_2way_intx$fitted.values)
summary(MMSE_2way_intx$fitted.values)

#---- **3-way interactions ----
MMSE_3way_intx <- lm(H1RMSESCORE ~ HCAP_age*female*race_ethnic_cat + 
                       SCHLYRS*female*race_ethnic_cat, 
                     data = HCAP_subset)

summary(MMSE_3way_intx)
plot(MMSE_3way_intx$residuals)
hist(MMSE_3way_intx$fitted.values)
summary(MMSE_3way_intx$fitted.values)

#---- mixture distributions ----
#Try to recover the HCAP shape distribution without taking covariates into 
# consideration... yet
samp_size = 10000

high_score_prop = 0.60; high_score_mean = 27; high_score_SD = 2
med_score_prop = 0.35; med_score_mean = 22; med_score_SD = 3.5
low_score_prop = 0.05; low_score_mean = 15; low_score_SD = 5

high_scores <- rnorm(n = floor(high_score_prop*samp_size), 
                     mean = high_score_mean, sd = high_score_SD)
high_scores[high_scores > 30] <- 30
high_scores[high_scores < 0] <- 0

med_scores <- rnorm(n = floor(med_score_prop*samp_size), 
                     mean = med_score_mean, sd = med_score_SD)
med_scores[med_scores > 30] <- 30
med_scores[med_scores < 0] <- 0

low_scores <- rnorm(n = floor(low_score_prop*samp_size), 
                    mean = low_score_mean, sd = low_score_SD)
low_scores[low_scores > 30] <- 30
low_scores[low_scores < 0] <- 0

HCAP_plot_data <- data.frame("values" = HCAP_subset$H1RMSESCORE, 
                             "data_type" = "HCAP")
synthetic_HCAP_plot_data <- 
  data.frame("values" = c(high_scores, med_scores, low_scores), 
             "data_type" = "Synthetic HCAP") 

plot_data <- rbind(HCAP_plot_data, synthetic_HCAP_plot_data)

#Comparison plot
ggplot(data = plot_data, aes(x = values, color = data_type, fill = data_type)) + 
  geom_histogram(alpha = 0.5, position = "identity") + theme_minimal() + 
  xlab("MMSE Total Scores") + theme(text = element_text(size = 14))

ggsave(filename = "synthetic_HCAP_MMSE_score.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

  
