#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "haven", "labelled")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
#---- **RAND ----
rand_waves <- seq(5, 9, by = 1) #Corresponding to ADAMS
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

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **ADAMS tracker ----
ADAMS_tracker_data_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
         "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/",
         "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

#filter to those who completed Wave A assessment
ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") %>% filter(AASSESS == 1) 
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
  left_join(., ADAMS_demdx, by = "HHIDPN") %>% 
  left_join(., RAND_subset, by = "HHIDPN")

#---- select vars ----
vars <- c("HHIDPN",
          #ADAMS interview year
          "AYEAR",
          #Sociodemographics (age, sex/gender, ethnicity, years of education, 
          # marital status, current working status)
          "AGE", "GENDER", "ETHNIC", "EDYRS", "AMARRD", "ACURRWK",
          #Health and health behaviors-- relevant to Wave A for now 
          # (ever stroke)
          paste0("r", seq(5, 7), "stroke"),
          #Cognition (total score MMSE)
          "NMSETOT", 
          #Dem dx
          "dem_dx_cat")

ADAMS_subset <- ADAMS %>% dplyr::select(contains(all_of(vars))) %>% 
  #Drop age bracket, age at selection, and language variables
  dplyr::select(-contains(c("BKT", "AGESEL", "LANGUAGE")))

ADAMS_waves <- c("A", "B", "C", "D")

#---- variable check ----
table(ADAMS_subset$GENDER, useNA = "ifany")
table(ADAMS_subset$ETHNIC, useNA = "ifany")
table(ADAMS_subset$EDYRS, useNA = "ifany")
table(ADAMS_subset$AAGE, useNA = "ifany")

#---- clean: neuropsych ----
#---- **MMSE ----
#Variable check
for(wave in ADAMS_waves){
  MMSE_var <- paste0(wave, "NMSETOT")
  print(paste0("Wave ", wave))
  print(table(ADAMS_subset[, MMSE_var], useNA = "ifany"))
}

ADAMS_subset %<>% 
  mutate_at(.vars = paste0(c("A", "B", "C", "D"), "NMSETOT"), 
            function(x) ifelse(x > 30, NA, x))

#---- transform: sociodemographics ----
#We want to use normal approximations to these variables 

#---- **age ----
#Variable check
hist(ADAMS_subset$AAGE)

ADAMS_subset %<>% mutate("log_AAGE" = log(AAGE))

#Sanity check
hist(ADAMS_subset$log_AAGE)

#---- **education ----
#Variable check
hist(ADAMS_subset$EDYRS)

#---- transform: neuropsych ----
#We want to use normal approximations to these variables 

#---- **MMSE ----
#Variable check-- this one is going to be skewed no matter what and it's one of 
# the variables that is going to help define out dementia classes
hist(ADAMS_subset$ANMSETOT)

#---- save dataset ----
write_csv(ADAMS_subset, path = paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                      "data/cleaned/ADAMS_subset_mixed.csv"))


