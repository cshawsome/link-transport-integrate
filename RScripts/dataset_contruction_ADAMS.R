#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
#---- **neuropsych ----
for(wave in c("a", "b", "c", "d")){
  neuropsych_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                                 "ADAMS1", str_to_upper(wave), "N_R.da")
  
  neuropsych_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                                 "ADAMS1", str_to_upper(wave), "N_R.dct")
  
  if(wave == "a"){
    ADAMS_neuropsych <- read_da_dct(neuropsych_data_path, neuropsych_dict_path, 
                                    HHIDPN = "TRUE") 
  } else{
    ADAMS_neuropsych %<>% 
      left_join(., read_da_dct(neuropsych_data_path, neuropsych_dict_path, 
                               HHIDPN = "TRUE"), 
                by = "HHIDPN")
  }
}

#---- **sociodemographic ----
ADAMS_tracker_data_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
         "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/",
         "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") 
#---- **dem dx ----
for(wave in c("a", "b", "c", "d")){
  demdx_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                            "ADAMS1", str_to_upper(wave), "D_R.da")
  demdx_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                            "ADAMS1", str_to_upper(wave), "D_R.dct")
  if(wave == "a"){
    ADAMS_demdx <- read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
      dplyr::select("HHIDPN", paste0(str_to_upper(wave), "DFDX1")) %>% 
      set_colnames(c("HHIDPN", "dem_dx")) %>% 
      mutate(new_col = 
               case_when(dem_dx %in% c(1, 2) ~ "Probable/Possible AD", 
                         dem_dx %in% c(3, 4) ~ 
                           "Probable/Possible Vascular Dementia", 
                         dem_dx %in% 
                           c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 21, 28, 29, 
                             30, 33) ~ "Other",
                         dem_dx %in% c(18, 32) ~ "Probable Dementia",
                         dem_dx %in% c(10, 13, 15, 16, 17, 19) ~ "Dementia", 
                         dem_dx %in% c(20, 22) ~ "MCI", 
                         dem_dx == 31 ~ "Normal")) %>% 
      dplyr::select(-c("dem_dx")) %>%
      set_colnames(c("HHIDPN", 
                     paste0(str_to_upper(wave), "dem_dx_cat")))
  } else{
    ADAMS_demdx %<>% 
      left_join(., read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
                  dplyr::select("HHIDPN", 
                                paste0(str_to_upper(wave), "DFDX1")) %>% 
                  set_colnames(c("HHIDPN", "dem_dx")) %>% 
                  mutate(new_col = 
                           case_when(dem_dx %in% c(1, 2) ~ 
                                       "Probable/Possible AD", 
                                     dem_dx %in% c(3, 4) ~ 
                                       "Probable/Possible Vascular Dementia", 
                                     dem_dx %in% 
                                       c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 
                                         21, 28, 29, 30, 33) ~ "Other",
                                     dem_dx %in% c(18, 32) ~ 
                                       "Probable Dementia",
                                     dem_dx %in% c(10, 13, 15, 16, 17, 19) ~ 
                                       "Dementia", 
                                     dem_dx %in% c(20, 22) ~ "MCI", 
                                     dem_dx == 31 ~ "Normal")) %>% 
                  dplyr::select(-c("dem_dx")) %>% 
                  set_colnames(c("HHIDPN", 
                                 paste0(str_to_upper(wave), "dem_dx_cat"))), 
                by = "HHIDPN")
  }
}

#---- join data ----
ADAMS <- left_join(ADAMS_tracker, ADAMS_neuropsych, by = "HHIDPN") %>% 
  left_join(., ADAMS_demdx, by = "HHIDPN")

#---- select vars ----
vars <- c("HHIDPN", "NMSETOT", )
