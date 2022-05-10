#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "haven", "labelled", "forcats", 
       "NormPsy")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- import data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **RAND ----
rand_waves <- 13 #Corresponding to HCAP

rand_variables = 
  #HHIDPN, 
  c("hhidpn",
    #Health and health behaviors: BMI, ADLs, IADLs, stroke (ever/never), 
    # hypertension (ever/never), diabetes (ever/never), CVD (ever/never), 
    # smoker (yes/no), days drinking, number drinks/day,
    paste0("r", rand_waves, "bmi"), paste0("r", rand_waves, "adla"), 
    paste0("r", rand_waves, "iadla"), paste0("r", rand_waves, "stroke"), 
    paste0("r", rand_waves, "hibpe"), paste0("r", rand_waves, "diabe"),
    paste0("r", rand_waves, "hearte"), paste0("r", rand_waves, "smoken"),
    paste0("r", rand_waves, "drinkd"), paste0("r", rand_waves, "drinkn"),
    #Cognitive variables: serial 7s, total MMSE score, BWC (20),
    # subjective cognitive change
    paste0("r", rand_waves, "ser7"), paste0("r", rand_waves, "cogtot"), 
    paste0("r", rand_waves, "bwc20"), paste0("r", rand_waves, "pstmem"))

RAND <- read_dta(paste0(path_to_box, "data/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **HRS tracker ----
HRS_tracker_data_path <- 
  paste0(path_to_box, "data/HRS_tracker/trk2018_3/TRK2018TR_R.da")

HRS_tracker_dict_path <- 
  paste0(path_to_box, "data/HRS_tracker/trk2018_3/TRK2018TR_R.dct")

HRS_tracker <- read_da_dct(HRS_tracker_data_path, HRS_tracker_dict_path, 
                           HHIDPN = "TRUE") %>% 
  #select variables: ID, HCAP selection, married/partnered status, sex/gender,
  # age, race, ethnicity, years of education, 2016 interview type
  dplyr::select("HHIDPN", "HCAP_SELECT", "PCOUPLE", "GENDER", "PAGE", "RACE", 
                "HISPANIC", "SCHLYRS", "PIWTYPE") %>%
  #Participated in 2016 HRS
  filter(PIWTYPE == 1) %>%  
  #get rid of leading 0
  mutate_at("HHIDPN", as.numeric) %>% 
  mutate_at("HHIDPN", as.character)

#---- **HRS Core ----
HRS_core_data_path <- 
  paste0(path_to_box, "data/HRS/Core_files/h16core/h16da/H16J_R.da")
HRS_core_dict_path <- 
  paste0(path_to_box, "data/HRS/Core_files/h16core/h16sta/H16J_R.dct")

HRS_core <- read_da_dct(HRS_core_data_path, HRS_core_dict_path, 
                        HHIDPN = "TRUE") %>% 
  #select variables: ID, employment status 
  dplyr::select("HHIDPN", "PJ005M1") %>% 
  #get rid of leading 0
  mutate_at("HHIDPN", as.numeric) %>% 
  mutate_at("HHIDPN", as.character)

#---- **HCAP ----
HCAP_data_path <- paste0(path_to_box, "data/HCAP/HC16/HC16da/HC16HP_R.da")

HCAP_dict_path <- paste0(path_to_box, "data/HCAP/HC16/HC16sta/HC16HP_R.dct")

HCAP_2016 <- read_da_dct(HCAP_data_path, HCAP_dict_path, HHIDPN = "TRUE") %>% 
  #select variables: ID, total MMSE score, object naming (scissor, cactus), 
  # president naming, animal naming, word recall (immediate), 
  # word recall (delayed), word list recall (yes, no), 
  # brave man recall (immediate), story recall (immediate), 
  # constructional praxis (immediate), constructional praxis (delayed), trails A
  dplyr::select("HHIDPN", "H1RMSESCORE", "H1RTICSSCISOR", "H1RTICSCACTUS", 
                "H1RTICSPRES", "H1RAFSCORE", "H1RWLIMM1SCORE", "H1RWLIMM2SCORE", 
                "H1RWLIMM3SCORE", "H1RWLDELSCORE", "H1RWLRECYSCORE", 
                "H1RWLRECNSCORE", "H1RBMIMMSCORE", "H1RLMIMMSCORE", 
                "H1RCPIMMSCORE", "H1RCPDELSCORE", "H1RTMASCORE") %>% 
  #get rid of leading 0
  mutate_at("HHIDPN", as.numeric) %>% 
  mutate_at("HHIDPN", as.character)

#---- join data ----
HCAP <- left_join(HCAP_2016, HRS_tracker, by = "HHIDPN") %>% 
  left_join(., HRS_core, by = "HHIDPN") %>%
  left_join(., RAND, by = "HHIDPN")

#---- clean data ----
#---- **HCAP select ----
# #these should all be 1, but looked in HCAP documentation and the sample size 
# # for HCAP is correct even though there is one person with a 2 (not selected)
# table(HCAP$HCAP_SELECT, useNA = "ifany")

#---- **2016 couple status ----
# table(HCAP$PCOUPLE, useNA = "ifany")

#change 5 = no to 0 = no
HCAP %<>% mutate_at(.vars = "PCOUPLE", function(x) ifelse(x == 5, 0, 1))

# #Sanity check
# table(HCAP$PCOUPLE, useNA = "ifany")

#---- **sex/gender ----
#table(HCAP$GENDER, useNA = "ifany")
HCAP %<>% 
  mutate(GENDER_label = as.factor(ifelse(GENDER == 1, "Male", "Female"))) %>% 
  mutate("Female" = ifelse(GENDER_label == "Female", 1, 0))

# #Sanity check
# table(HCAP$GENDER, HCAP$GENDER_label)
# table(HCAP$Female, useNA = "ifany")

#---- **age ----
# table(HCAP$PAGE, useNA = "ifany")

#restrict to 70+ (N = 2610, dropped n = 13,534)
HCAP %<>% filter(PAGE >= 886)

#---- **race ----
#table(HCAP$RACE, useNA = "ifany")
HCAP %<>% 
  mutate(RACE_label = as.factor(case_when(RACE == 1 ~ "White", 
                                          RACE == 2 ~ "Black", 
                                          RACE == 7 ~ "Other"))) %>% 
  mutate("RACE_White" = ifelse(RACE_label == "White", 1, 0), 
         "RACE_Black" = ifelse(RACE_label == "Black", 1, 0), 
         "RACE_Other" = ifelse(RACE_label == "Other", 1, 0))

# #Sanity check
# table(HCAP$RACE_label, HCAP$RACE_White, useNA = "ifany")
# table(HCAP$RACE_label, HCAP$RACE_Black, useNA = "ifany")
# table(HCAP$RACE_label, HCAP$RACE_Other, useNA = "ifany")

#---- **hispanic ----
#table(HCAP$HISPANIC, useNA = "ifany")

#change 5 = no to 0 = no, combine all hispanic types and set 0 = not obtained to NA
HCAP %<>% mutate("HISPANIC_indicator" = case_when(HISPANIC %in% c(1, 2, 3) ~ 1, 
                                                  HISPANIC %in% c(0, 5) ~ 0))

# #Sanity check
# table(HCAP$HISPANIC, HCAP$HISPANIC_indicator, useNA = "ifany")

#---- **race/ethnicity ----
HCAP %<>% 
  mutate(ETHNIC_label = 
           case_when(HISPANIC_indicator == 1 ~ "Hispanic",
                     HISPANIC_indicator == 0 & RACE_label == "White" ~ "White",
                     HISPANIC_indicator == 0 & RACE_label == "Black" ~ "Black", 
                     HISPANIC_indicator == 0 & RACE_label == "Other" ~ "Other")) %>% 
  mutate("White" = ifelse(ETHNIC_label == "White", 1, 0), 
         "Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0), 
         "Other" = ifelse(ETHNIC_label == "Other", 1, 0))

# #Sanity check
# table(HCAP$ETHNIC_label, HCAP$White, useNA = "ifany")
# table(HCAP$ETHNIC_label, HCAP$Black, useNA = "ifany")
# table(HCAP$ETHNIC_label, HCAP$Hispanic, useNA = "ifany")
# table(HCAP$ETHNIC_label, HCAP$Other, useNA = "ifany")

#restrict to White, Black, Hispanic (to align with ADAMS data) 
# (N = 2551, dropped n = 59)

HCAP %<>% filter(Other == 0)
