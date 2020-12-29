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

#filter to those who completed Wave A assessment
ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") %>% filter(AASSESS == 1) 
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
      dplyr::select("HHIDPN", paste0(str_to_upper(wave), "DFDX1")) 
  } else{
    ADAMS_demdx %<>% 
      left_join(., read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
                  dplyr::select("HHIDPN", 
                                paste0(str_to_upper(wave), "DFDX1")))
  }
}

for(wave in c("A", "B", "C", "D")){
  dem_dx_var <- paste0(wave, "DFDX1")
  ADAMS_demdx[, paste0(wave, "dem_dx_cat")] <- 
    case_when(ADAMS_demdx[, dem_dx_var] %in% c(1, 2) ~ "Probable/Possible AD", 
              ADAMS_demdx[, dem_dx_var] %in% c(3, 4) ~ 
                "Probable/Possible Vascular Dementia", 
              ADAMS_demdx[, dem_dx_var] %in% 
                c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 21, 28, 29, 
                  30, 33) ~ "Other",
              ADAMS_demdx[, dem_dx_var] %in% c(18, 32) ~ "Probable Dementia",
              ADAMS_demdx[, dem_dx_var] %in% c(10, 13, 15, 16, 17, 19) ~ 
                "Dementia", 
              ADAMS_demdx[, dem_dx_var] %in% c(20, 22) ~ "MCI", 
              ADAMS_demdx[, dem_dx_var] == 31 ~ "Normal")
}

# #Sanity check
# for(wave in c("A", "B", "C", "D")){
#   dem_dx_var <- paste0(wave, "DFDX1")
#   print(paste0("Wave ", wave))
#   print(table(ADAMS_demdx[, dem_dx_var], 
#               ADAMS_demdx[, paste0(wave, "dem_dx_cat")], 
#               useNA = "ifany"))
# }

#Remove original variables
ADAMS_demdx %<>% dplyr::select(-c(paste0(c("A", "B", "C", "D"), "DFDX1")))

#---- join data ----
ADAMS <- left_join(ADAMS_tracker, ADAMS_neuropsych, by = "HHIDPN") %>% 
  left_join(., ADAMS_demdx, by = "HHIDPN")

#---- select vars ----
vars <- c("HHIDPN",
          #Sociodemographics
          "AGE", "GENDER", "ETHNIC", "EDYRS",
          #Cognition
          "NMSETOT", 
          #Dem dx
          "dem_dx_cat")

ADAMS_subset <- ADAMS %>% dplyr::select(contains(all_of(vars))) %>% 
  #Drop age bracket, age at selection, and language variables
  dplyr::select(-contains(c("BKT", "AGESEL", "LANGUAGE")))

ADAMS_waves <- c("A", "B", "C", "D")

#---- clean: AGE ----
for(wave in ADAMS_waves){
  col <- paste0(wave, "AGE")
  new_col <- paste0(col, "_cat")
  
  ADAMS_subset[, new_col] <- ADAMS_subset[, col] - 
    (min(ADAMS_subset[, col], na.rm = TRUE) - 1)
  
  ADAMS_subset[which(ADAMS_subset[, new_col] > 900), new_col] <- NA
}

#Sanity check
for(wave in ADAMS_waves){
  col <- paste0(wave, "AGE")
  new_col <- paste0(col, "_cat")
  
  print(paste0("Wave", wave))
  print(table(ADAMS_subset[, new_col], ADAMS_subset[, col], useNA = "ifany"))
}

#---- clean: MMSE ----
#Variable check
for(wave in c("A", "B", "C", "D")){
  MMSE_var <- paste0(wave, "NMSETOT")
  print(paste0("Wave ", wave))
  print(table(ADAMS_subset[, MMSE_var], useNA = "ifany"))
}

ADAMS_subset %<>% 
  mutate_at(.vars = paste0(c("A", "B", "C", "D"), "NMSETOT"), 
            function(x) ifelse(x > 30, NA, x))

#---- save dataset ----
write_csv(ADAMS_subset, path = paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                      "data/cleaned/ADAMS_subset.csv"))


