#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "haven", "labelled", "forcats")

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
#---- **ADAMS neuropsych ----
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

#---- **ADAMS informant questionnaire ----
for(wave in c("a", "b", "c", "d")){
  proxy_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                            "ADAMS1", str_to_upper(wave), "G_R.da")
  proxy_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                            "ADAMS1", str_to_upper(wave), "G_R.dct")
  if(wave == "a"){
    ADAMS_proxy <- read_da_dct(proxy_data_path, proxy_dict_path, 
                               HHIDPN = "TRUE") %>% 
      dplyr::select("HHIDPN", paste0(str_to_upper(wave), "GQ", seq(14, 29))) 
  } else{
    ADAMS_proxy %<>% 
      left_join(., read_da_dct(proxy_data_path, proxy_dict_path, 
                               HHIDPN = "TRUE") %>% 
                  dplyr::select("HHIDPN", 
                                paste0(str_to_upper(wave), "GQ", seq(14, 29))))
  }
}

#---- **ADAMS dem dx ----
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
  left_join(., ADAMS_proxy, by = "HHIDPN") %>% 
  left_join(., ADAMS_demdx, by = "HHIDPN") %>% 
  left_join(., RAND, by = "HHIDPN")

#---- select vars ----
vars <- c("HHIDPN",
          #ADAMS interview year
          "AYEAR",
          #Sociodemographics (age, sex/gender, ethnicity, years of education, 
          # marital status, current working status)
          "AGE", "GENDER", "ETHNIC", "EDYRS", "AMARRD", "ACURRWK",
          #Health and health behaviors-- relevant to Wave A for now 
          # (ever stroke, ever hypertension, ever diabetes, ever CVD, 
          # ever cancer, BMI, IADLs, ADLs, gross motor, fine motor, depressive
          # symptoms, subjective cognitive complaints, smoking status, 
          # number days drinking, number drinks per day) 
          paste0("r", seq(5, 7), "stroke"), paste0("r", seq(5, 7), "hibp"), 
          paste0("r", seq(5, 7), "diab"), paste0("r", seq(5, 7), "heart"),
          paste0("r", seq(5, 7), "cancr"), paste0("r", seq(5, 7), "bmi"), 
          paste0("r", seq(5, 7), "iadla"), paste0("r", seq(5, 7), "adla"), 
          paste0("r", seq(5, 7), "grossa"), paste0("r", seq(5, 7), "finea"),
          paste0("r", seq(5, 7), "cesd"), paste0("r", seq(5, 7), "pstmem"), 
          paste0("r", seq(5, 7), "smoken"), paste0("r", seq(5, 7), "drinkd"),
          paste0("r", seq(5, 7), "drinkn"),
          #WAVE A only
          #Cognition (total score MMSE, backwards count 20 (>1 trials), 
          #backwards count 86 (>1 trials), total score serial 7s, 
          #object naming scissors, object naming cactus, President, VP, 
          #total score animal naming, total score Boston naming, 
          #total score immediate constructional praxis, total score delayed 
          #constructional praxis, immediate word recall (>1 trials), delayed
          #word recall, word recall (yes), word recall (no), 
          #total score Wechsler logical memory I, 
          #total score Wechsler logical memory II, Trails A, Trails B, 
          #total score symbol digit modalities)
          "ANMSETOT", "ANBWC20", "ANBWC86", "ANSER7T", "ANSCISOR", "ANCACTUS", 
          "ANPRES", "ANVCPRES", "ANAFTOT", "ANBNTTOT", "ANCPTOT", "ANDCPTOT", 
          "ANIMMCR", "ANDELCOR", "ANRECYES", "ANRECNO", "ANWM1TOT", "ANWM2TOT", 
          "ANTMASEC", "ANTMBSEC", "ANSDMTOT",
          #proxy cognition: IQCODE
          "AGQ",
          #Dem dx
          "dem_dx_cat")

ADAMS_subset <- ADAMS %>% dplyr::select(contains(all_of(vars))) %>% 
  #Drop age bracket, age at selection, and language variables
  dplyr::select(-contains(c("BKT", "AGESEL", "LANGUAGE")))

ADAMS_waves <- c("A", "B", "C", "D")

#---- clean: sociodemographics ----
#---- **variable check ----
#table(ADAMS_subset$AAGE, useNA = "ifany")

#---- **sex/gender ----
#table(ADAMS_subset$GENDER, useNA = "ifany")
ADAMS_subset %<>% 
  mutate(GENDER_label = as.factor(ifelse(GENDER == 1, "Male", "Female"))) %>% 
  mutate(GENDER_label = fct_relevel(GENDER_label, "Female"))

# #Sanity check
# levels(ADAMS_subset$GENDER_label)

#---- **race/ethnicity ----
#table(ADAMS_subset$ETHNIC, useNA = "ifany")
ADAMS_subset %<>% 
  mutate(ETHNIC_label = as.factor(case_when(ETHNIC == 1 ~ "White", 
                                            ETHNIC == 2 ~ "Black", 
                                            ETHNIC == 3 ~ "Hispanic"))) %>% 
  mutate(ETHNIC_label = fct_relevel(ETHNIC_label, "White"))

# #Sanity check
# levels(ADAMS_subset$ETHNIC_label)

#---- **education ----
#table(ADAMS_subset$EDYRS, useNA = "ifany")
ADAMS_subset %<>% 
  mutate(EDYRScat = case_when(EDYRS == 0 ~ 1, 
                              EDYRS %in% seq(1, 8) ~ 2, 
                              EDYRS %in% seq(9, 11) ~ 3, 
                              EDYRS == 12 ~ 4, 
                              EDYRS %in% seq(13, 15) ~ 5, 
                              EDYRS >= 16 ~ 6)) %>% 
  mutate(EDYRScat_label = 
           as.factor(case_when(EDYRScat == 1 ~ "No School Completed", 
                               EDYRScat == 2 ~ "1st-8th Grade", 
                               EDYRScat == 3 ~ "Some High School", 
                               EDYRScat == 4 ~ "High School Diploma",
                               EDYRScat == 5 ~ "Some College", 
                               EDYRScat == 6 ~ "Bachelor's Degree or Higher"))
  ) %>% 
  mutate(EDYRScat_label = 
           fct_relevel(EDYRScat_label, "No School Completed", "1st-8th Grade", 
                       "Some High School", "High School Diploma", 
                       "Some College"))

# #Sanity check
# table(ADAMS_subset$EDYRS, ADAMS_subset$EDYRScat, useNA = "ifany")
# table(ADAMS_subset$EDYRScat, ADAMS_subset$EDYRScat_label, useNA = "ifany")
# levels(ADAMS_subset$EDYRScat_label)

#---- **marital status ----
#table(ADAMS_subset$AAMARRD, useNA = "ifany")
for(wave in ADAMS_waves){
  ADAMS_subset %<>% 
    mutate(!!paste0(wave, "AMARRD_label") := 
             case_when(!!sym(paste0(wave, "AMARRD")) %in% c(1, 3, 4) ~ 
                         "Not married/partnered",
                       !!sym(paste0(wave, "AMARRD"))  == 2 ~ 
                         "Married/partnered", 
                       !!sym(paste0(wave, "AMARRD")) == 5 ~ "Widowed", 
                       !!sym(paste0(wave, "AMARRD")) == 8 ~ "Unknown")) %>% 
    mutate(!!paste0(wave, "AMARRD_label") := 
             fct_relevel(!!sym(paste0(wave, "AMARRD_label")), 
                         "Married/partnered"))
}

# #Sanity check
# for(wave in ADAMS_waves){
#   print(paste0("Wave ", wave))
#   show(table(ADAMS_subset[, paste0(wave, "AMARRD")],
#         ADAMS_subset[, paste0(wave, "AMARRD_label")], useNA = "ifany"))
#   show(levels(ADAMS_subset[, paste0(wave, "AMARRD_label")]))
# }

#---- **retirement status ----
#table(ADAMS_subset$AACURRWK, useNA = "ifany")

#looks like there's no data on working status past wave B
for(wave in c("A", "B")){
  ADAMS_subset %<>% 
    mutate(!!paste0(wave, "ACURRWK_label") := 
             case_when(!!sym(paste0(wave, "ACURRWK")) == 1 ~ "Working",
                       !!sym(paste0(wave, "ACURRWK"))  == 2 ~ "Retired", 
                       !!sym(paste0(wave, "ACURRWK")) == 3 ~ "Semi-retired", 
                       !!sym(paste0(wave, "ACURRWK")) == 4 ~ "Disabled", 
                       !!sym(paste0(wave, "ACURRWK")) == 5 ~ "Unemployed")) %>% 
    mutate(!!paste0(wave, "ACURRWK_label") := 
             fct_relevel(!!sym(paste0(wave, "ACURRWK_label")), "Retired"))
}

# #Sanity check
# for(wave in c("A", "B")){
#   print(paste0("Wave ", wave))
#   show(table(ADAMS_subset[, paste0(wave, "ACURRWK")],
#         ADAMS_subset[, paste0(wave, "ACURRWK_label")], useNA = "ifany"))
#   show(levels(ADAMS_subset[, paste0(wave, "ACURRWK_label")]))
# }

#---- clean: health and health behaviors ----
#---- **wave-updated vars ----
#For repeated measures, want to take the wave most representative of ADAMS wave
wave_updated_vars <- c("stroke", "hibpe", "diabe", "hearte", "cancre", "bmi", 
                       "iadla", "adla", "grossa", "finea", "cesd", "pstmem", 
                       "smoken", "drinkd", "drinkn")

for(var in wave_updated_vars){
  ADAMS_subset %<>% 
    mutate(!!paste0("A", var) := 
             case_when(AYEAR == 2001 ~ !!sym(paste0("r5", var)), 
                       AYEAR %in% c(2002, 2003) ~ !!sym(paste0("r6", var)), 
                       AYEAR == 2004 ~ !!sym(paste0("r7", var))))
}

# #Sanity check
# View(ADAMS_subset[, c("AYEAR", paste0("r", seq(5, 7), "cesd"), "Acesd")] %>%
#        filter(!is.na(Acesd)))
# colnames(ADAMS_subset)
# 
#Variable check
for(var in wave_updated_vars){
  print(var)
  show(table(ADAMS_subset[, paste0("A", var)], useNA = "ifany"))
}

#---- **drinking ----
ADAMS_subset %<>% 
  mutate("drinks_per_week" = Adrinkd*Adrinkn) %>%
  mutate("DRINKcat" = 
           case_when(drinks_per_week == 0 ~ 1,
                     GENDER_label == "Male" &
                       (drinks_per_week >= 1 & drinks_per_week < 14) ~ 2,
                     GENDER_label == "Female" &
                       (drinks_per_week >= 1 & drinks_per_week < 7) ~ 2,
                     GENDER_label == "Male" &
                       (drinks_per_week >= 14 | Adrinkn >= 4) ~ 3,
                     GENDER_label == "Female" &
                       (drinks_per_week >= 7 | Adrinkn >= 3) ~ 3)) %>% 
  mutate("DRINKcat_label" = 
           as.factor(case_when(DRINKcat == 1 ~ "No Drinking", 
                               DRINKcat == 2 ~ "Moderate Drinking", 
                               DRINKcat == 3 ~ "Heavy Drinking"))) %>% 
  mutate(DRINKcat_label = fct_relevel(DRINKcat_label, "No Drinking", 
                                      "Moderate Drinking"))

# #Sanity check
# View(ADAMS_subset[, c("GENDER_label", "drinks_per_week", "Adrinkn", 
#                       "DRINKcat_label")])
# levels(ADAMS_subset$DRINKcat_label)

#---- clean: neuropsych ----
#---- **MMSE ----
# #Variable check
# table(ADAMS_subset$ANMSETOT, useNA = "ifany")

ADAMS_subset %<>% 
  mutate_at(.vars = "ANMSETOT", function(x) ifelse(x > 30, NA, x))

# #Sanity check
# table(ADAMS_subset$ANMSETOT, useNA = "ifany")

#---- **BWC 20 and 86 ----
# #Variable check
# table(ADAMS_subset$ANBWC201, useNA = "ifany")
# table(ADAMS_subset$ANBWC202, useNA = "ifany")
# table(ADAMS_subset$ANBWC861, useNA = "ifany")
# table(ADAMS_subset$ANBWC862, useNA = "ifany")

#Recode
ADAMS_subset %<>% 
  mutate_at(.vars = c("ANBWC201", "ANBWC202", "ANBWC861", "ANBWC862"), 
            #Missing/refused  
            function(x) ifelse(x > 6, NA, x)) %>% 
  mutate_at(.vars = c("ANBWC201", "ANBWC202", "ANBWC861", "ANBWC862"), 
            #restart
            function(x) ifelse(x == 6, 0, x)) 

# #Sanity check
# table(ADAMS_subset$ANBWC201, useNA = "ifany")
# table(ADAMS_subset$ANBWC202, useNA = "ifany")
# table(ADAMS_subset$ANBWC861, useNA = "ifany")
# table(ADAMS_subset$ANBWC862, useNA = "ifany")

#Take the higher score
ADAMS_subset %<>% mutate("ANBWC20" = pmax(ANBWC201, ANBWC202, na.rm = TRUE), 
                         "ANBWC86" = pmax(ANBWC861, ANBWC862, na.rm = TRUE))

# #Sanity check
# View(ADAMS_subset[, c("ANBWC201", "ANBWC202", "ANBWC20")])
# View(ADAMS_subset[, c("ANBWC861", "ANBWC862", "ANBWC86")])

#---- **serial 7s ----
# #Variable check
# table(ADAMS_subset$ANSER7T, useNA = "ifany")

#Recode 
ADAMS_subset %<>% 
  mutate_at(.vars = c("ANSER7T"), 
            #Missing/refused  
            function(x) ifelse(x > 5, NA, x))

# #Sanity check
# table(ADAMS_subset$ANSER7T, useNA = "ifany")

#---- **object naming: cactus, scissors ----
# #Variable check
# table(ADAMS_subset$ANCACTUS, useNA = "ifany")
# table(ADAMS_subset$ANSCISOR, useNA = "ifany")

#Recode
ADAMS_subset %<>% 
  mutate_at(.vars = c("ANCACTUS", "ANSCISOR"), 
            #Missing/refused  
            function(x) ifelse(x > 1, NA, x)) 

# #Sanity check
# table(ADAMS_subset$ANCACTUS, useNA = "ifany")
# table(ADAMS_subset$ANSCISOR, useNA = "ifany")
                
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


