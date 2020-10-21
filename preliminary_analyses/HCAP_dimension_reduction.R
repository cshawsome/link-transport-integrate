#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16da/HC16HP_R.da")
dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16sta/HC16HP_R.dct")

HCAP <- read_da_dct(data_path, dict_path, HHIDPN = "TRUE")

#---- variables of interest ----
HCAP_vars <- c("HHIDPN", "MSE", "TICS", "WLIMM1", "WLIMM2", "WLIMM3", "1066", 
               "SCORE")
not_these <- c("TEST", "INTRO", "H1RMSESCORE", "H1RTICSSCORE", "H1R11066SCORE", 
               "SCORESE", "CESDSCORE", "STSCORE")

HCAP_assessment <- HCAP %>% dplyr::select(contains(HCAP_vars)) %>% 
  dplyr::select(-c(contains(not_these), "H1RMSE11T")) #H1RMSE11T1 needs to stay

#---- sum scores and averages of repeated trials ----


#---- make correlation matrix ----


