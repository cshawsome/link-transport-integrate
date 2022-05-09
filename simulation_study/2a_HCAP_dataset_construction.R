#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "haven", "labelled", "forcats", 
       "NormPsy")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- import data ----
#---- **RAND ----
rand_waves <- 13 #Corresponding to HCAP

rand_variables = 
  #HHIDPN, gender, race, ethnicity
  c("hhidpn", "ragender", "raracem", "rahispan", 
    #Sampling variables: marital status, coupled household status,
    paste0("r", rand_waves, "mstat"), paste0("h", hrs_waves, "cpl"), 
    #Sociodemographics: age at beginning interview
    paste0("r", hrs_waves, "agey_b"),
    #Health and health behaviors: BMI, IADLs, stroke (ever/never)
    paste0("r", hrs_waves, "bmi"), paste0("r", hrs_waves, "iadla"),
    paste0("r", hrs_waves, "stroke"),
    #Cognitive variables: serial 7s, immediate word recall, delayed word recall, 
    # total MMSE score
    paste0("r", hrs_waves, "ser7"),
    paste0("r", hrs_waves, "imrc"), paste0("r", hrs_waves, "dlrc"), 
    paste0("r", hrs_waves, "cogtot"))

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **HRS tracker ----
HRS_tracker_data_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS_tracker/trk2018_3/", 
         "TRK2018TR_R.da")
HRS_tracker_dict_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS_tracker/trk2018_3/", 
         "TRK2018TR_R.dct")

HRS_tracker <- read_da_dct(HRS_tracker_data_path, HRS_tracker_dict_path, 
                           HHIDPN = "TRUE") %>% 
  #Participated in 2016 HRS
  filter(PIWTYPE == 1)  

#---- **HCAP ----
HCAP_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                         "HCAP/HC16/HC16da/HC16HP_R.da")

HCAP_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                         "HCAP/HC16/HC16sta/HC16HP_R.dct")

HCAP_2016 <- read_da_dct(HCAP_data_path, HCAP_dict_path, HHIDPN = "TRUE") 

