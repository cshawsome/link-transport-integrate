#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "stringr", "broom", "openxlsx")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data: neuropsych ----
neuropsych_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
neuropsych_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                 "ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych_A <- read_da_dct(neuropsych_data_path_A, 
                                  neuropsych_dict_path_A, HHIDPN = "TRUE")

#---- import data: sociodemographic ----
ADAMS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                  "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                  "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") 

#---- import data: demdx ----
demdx_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1a/adams1ada/ADAMS1AD_R.da")
demdx_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                            "ADAMS/adams1a/adams1asta/ADAMS1AD_R.dct")

ADAMS_demdx_A <- read_da_dct(demdx_data_path_A, demdx_dict_path_A, 
                             HHIDPN = "TRUE") 

#---- format data: neuropsych ----

#---- format data: sociodemographics ----
ADAMS_tracker %<>% 
  dplyr::select("HHIDPN", "AAGE", "GENDER", "ETHNIC", "EDYRS") %>% 
  mutate("female" = ifelse(GENDER == 1, 0, 1), 
         "ethnic_cat" = case_when(ETHNIC == 1 ~ "Non-hispanic White", 
                                  ETHNIC == 2 ~ "Non-hispanic Black", 
                                  ETHNIC == 3 ~ "Hispanic"))

#---- format data: dem dx ----
ADAMS_demdx_A %<>% 
  dplyr::select("HHIDPN", "ADFDX1") %>% 
  mutate("dem_dx" = 
           case_when(ADFDX1 %in% c(1, 2) ~ "Probable/Possible AD", 
                     ADFDX1 %in% c(3, 4) ~ 
                       "Probable/Possible Vascular Dementia", 
                     ADFDX1 %in% 
                       c(5, 8, 14, 23, 24, 25, 26, 27, 21, 28, 29, 30) ~ 
                       "Other",
                     ADFDX1 %in% c(18) ~ "Probable Dementia",
                     ADFDX1 %in% c(10, 13, 15) ~ "Dementia", 
                     ADFDX1 %in% c(20, 22) ~ "MCI", 
                     ADFDX1 == 31 ~ "Normal"), 
         "impaired" = 
           ifelse(dem_dx %in% c("Normal", "Other"), 0, 1))

#---- generate sythetic data ----