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

#---- coding correct/incorrect answers ----
#Need to recode 2 = correct (different level) to 1
recode_correct <- c("ANMSE17", "ANMSE19", "ANMSE21")

#Need to recode 6 = start over to 0
recode_incorrect <- 
  colnames(ADAMSA_assessment)[as.logical(
    str_detect(colnames(ADAMSA_assessment), "BWC"))] 

for(var in recode_incorrect){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] == 6), var] <- 0
}

for(var in recode_correct){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] == 2), var] <- 1
}

# #Sanity check
# for(var in colnames(ADAMSA_assessment)){
#   print(var)
#   print(table(ADAMSA_assessment[, var], useNA = "ifany"))
# }

