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
  #get rid of leading 0
  mutate_at("HHIDPN", as.numeric) %>% 
  mutate_at("HHIDPN", as.character)

#---- join data ----
HCAP <- left_join(HCAP_2016, HRS_tracker, by = "HHIDPN") %>% 
  left_join(., HRS_core, by = "HHIDPN") %>%
  left_join(., RAND, by = "HHIDPN")

