#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "haven")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- ADAMS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- ****ADAMS tracker ----
ADAMS_tracker_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path,
                             HHIDPN = "TRUE") %>% 
  #select variables: ID, Wave A participation, HRS cognitive test at sampling, 
  #                  interview year, marital status, age, race/ethnicity, 
  #                  years of education, employment status, HRS proxy cognition
  dplyr::select("HHIDPN", "AASSESS", "GENDER", "SELFCOG", "AYEAR", "AAMARRD", 
                "AAGE", "ETHNIC", "EDYRS", "AACURRWK", "PROXCOG") %>%
  #N = 1170
  #filter to those who completed Wave A assessment (N = 856; dropped n = 314)
  filter(AASSESS == 1)

#---- ****neuropsych measures ----
ADAMS_neuropsych_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
ADAMS_neuropsych_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych <- 
  read_da_dct(ADAMS_neuropsych_data_path, ADAMS_neuropsych_dict_path,
              HHIDPN = "TRUE") %>% 
  #select variables: ID, MMSE, backwards count (20), backwards count (86), 
  #                  serial 7s, item naming (scissors), item naming (cactus), 
  #                  President naming, VP naming, animal naming, 
  #                  Boston naming test, immediate word recall, 
  #                  delayed word recall, word list recognition (yes), 
  #                  word list recognition (no), immediate story recall, 
  #                  delayed story recall, immediate constructional praxis, 
  #                  delayed constructional praxis, symbol/digit substitution, 
  #                  trails A, trails B, subjective cognitive change
  dplyr::select("HHIDPN", "ANMSETOT", "ANBWC201", "ANBWC202", "ANBWC861", 
                "ANBWC862", "ANSER7T", "ANSCISOR", "ANCACTUS", "ANPRES", 
                "ANVCPRES", "ANAFTOT", "ANBNTTOT", "ANIMMCR1", "ANIMMCR2", 
                "ANIMMCR3", "ANDELCOR", "ANRECYES", "ANRECNO", "ANWM1TOT", 
                "ANWM2TOT", "ANCPTOT", "ANRCPTOT", "ANSDMTOT", "ANTMASEC", 
                "ANTMBSEC", "ANSMEM2")

#---- ****dementia dx ----
ADAMS_demdx_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1ada/ADAMS1AD_R.da")
ADAMS_demdx_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1asta/ADAMS1AD_R.dct")

ADAMS_demdx <- read_da_dct(ADAMS_demdx_data_path, ADAMS_demdx_dict_path,
                           HHIDPN = "TRUE") %>% 
  dplyr::select("HHIDPN", "ADCCDX1")

#---- ****proxy measures ----
ADAMS_proxy_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1ada/ADAMS1AG_R.da")
ADAMS_proxy_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1asta/ADAMS1AG_R.dct")

ADAMS_proxy <- read_da_dct(ADAMS_proxy_data_path, ADAMS_proxy_dict_path,
                           HHIDPN = "TRUE") %>% 
  #select variables: IQCODE items, proxy type
  dplyr::select("HHIDPN", paste0("AGQ", c(seq(14, 29), 101)))

#---- ****RAND variables ----
rand_waves <- seq(5, 6, by = 1) #Corresponding to ADAMS
rand_variables <- 
  c("hhidpn", 
    #Health and health behaviors (ever/never stroke, ever/never
    # hypertension, ever/never diabetes, ever/never cvd, BMI, IADLs, ADLs, 
    # depressive symptoms, smokes now, number days drinking per week, 
    # number drinks/day) 
    paste0("r", rand_waves, "stroke"), paste0("r", rand_waves, "hibpe"), 
    paste0("r", rand_waves, "diabe"), paste0("r", rand_waves, "hearte"),
    paste0("r", rand_waves, "bmi"), paste0("r", rand_waves, "iadla"), 
    paste0("r", rand_waves, "adla"), paste0("r", rand_waves, "cesd"), 
    paste0("r", rand_waves, "smoken"), paste0("r", rand_waves, "drinkd"),
    paste0("r", rand_waves, "drinkn"))

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **join data ----

#---- **data cleaning ----

#---- OLD ----

#---- **RAND ----
rand_waves <- seq(5, 9, by = 1) #Corresponding to ADAMS +- 1 wave (for imputation)
rand_variables <- 
  c("hhidpn", 
    #Health and health behaviors (ever/never stroke, ever/never
    # hypertension, ever/never diabetes, ever/never cvd, ever/never cancer, 
    # BMI, IADLs, ADLs, gross motor, fine motor, depressive symptoms, 
    # self-reported change in memory, smokes now, number days drinking per week, 
    # number drinks/day) 
    paste0("r", rand_waves, "stroke"), paste0("r", rand_waves, "hibpe"), 
    paste0("r", rand_waves, "diabe"), paste0("r", rand_waves, "hearte"),
    paste0("r", rand_waves, "cancre"), paste0("r", rand_waves, "bmi"),
    paste0("r", rand_waves, "iadla"), paste0("r", rand_waves, "adla"),
    paste0("r", rand_waves, "grossa"), paste0("r", rand_waves, "finea"),
    paste0("r", rand_waves, "cesd"), paste0("r", rand_waves, "pstmem"), 
    paste0("r", rand_waves, "smoken"), paste0("r", rand_waves, "drinkd"),
    paste0("r", rand_waves, "drinkn"))

RAND <- read_dta(paste0(path_to_box, "Dissertation/data/HRS/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **ADAMS ----
ADAMS <- rbind(read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_test.csv")), 
               read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_train.csv")))






