#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

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

#---- select variables ----
HCAP_vars <- c("HHIDPN", "GENDER", "HISPANIC", "RACE", "SCHLYRS",
               #for age calculation
               "BIRTHYR", "H1RIWYEAR",
               #MMSE
               "H1RMSESCORE")

HCAP_subset <- HCAP %>% dplyr::select(all_of(HCAP_vars))

# #Variable check
# colSums(is.na(HCAP_subset))

#Drop people without an interview year
HCAP_subset %<>% filter(!is.na(H1RIWYEAR))

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
ggplot(data = HCAP_subset) + 
  geom_histogram(aes(H1RMSESCORE), color = green, fill = green) +
  theme_minimal() + xlab("MMSE Score") + theme(text = element_text(size = 14))

ggsave(filename = "HCAP_MMSE_score.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

# #Variable check
# summary(HCAP_subset$H1RMSESCORE)

#---- Modeling HCAP data ----
  
