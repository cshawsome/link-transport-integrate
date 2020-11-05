#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "stringr", "ggcorrplot", "psych")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
neuropsych_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
neuropsych_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                 "ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych_A <- read_da_dct(neuropsych_data_path_A, 
                                  neuropsych_dict_path_A, HHIDPN = "TRUE")

cog_test_labels <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                   "cog_test_meaningful_labels.csv"))

#---- variables of interest ----
ADAMS_vars <- c("HHIDPN", paste0("MSE", seq(1, 22, by = 1)), "BWC86", "SER7T", 
                "SCISOR", "CACTUS", "ANPRES", "ANAFTOT", "CPTOT", "IMMCR", 
                "DELCOR", "RECYES", "RECNO", "WM1TOT", "WM2TOT", "MASEC", 
                "MBSEC", "ANSDMTOT")

not_these <- c("ANMSE11T", "ANRCPTOT")

ADAMSA_assessment <- ADAMS_neuropsych_A %>% 
  dplyr::select(contains(ADAMS_vars)) %>% 
  dplyr::select(-c(contains(not_these))) 
