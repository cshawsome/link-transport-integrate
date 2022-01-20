#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "haven", "magrittr", "NormPsy")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- ADAMS ----
#---- **read data ----
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

RAND <- read_dta(paste0(path_to_box, "data/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character) %>% rename("HHIDPN" = "hhidpn")

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **join data ----
ADAMS <- left_join(ADAMS_tracker, ADAMS_neuropsych, by = "HHIDPN") %>% 
  left_join(., ADAMS_proxy, by = "HHIDPN") %>% 
  left_join(., ADAMS_demdx, by = "HHIDPN") %>% 
  left_join(., RAND, by = "HHIDPN")

#---- **clean data ----
#---- ****sex/gender ----
#table(ADAMS$GENDER, useNA = "ifany")
ADAMS %<>% 
  mutate(GENDER_label = as.factor(ifelse(GENDER == 1, "Male", "Female"))) %>% 
  mutate("Female" = ifelse(GENDER_label == "Female", 1, 0))

# #Sanity check
# table(ADAMS$GENDER, ADAMS$GENDER_label)
# table(ADAMS$Female, useNA = "ifany")

#---- ****race/ethnicity ----
#table(ADAMS$ETHNIC, useNA = "ifany")
ADAMS %<>% 
  mutate(ETHNIC_label = as.factor(case_when(ETHNIC == 1 ~ "White", 
                                            ETHNIC == 2 ~ "Black", 
                                            ETHNIC == 3 ~ "Hispanic"))) %>% 
  mutate("White" = ifelse(ETHNIC_label == "White", 1, 0), 
         "Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0))

# #Sanity check
# table(ADAMS$ETHNIC_label, ADAMS$White, useNA = "ifany")
# table(ADAMS$ETHNIC_label, ADAMS$Black, useNA = "ifany")
# table(ADAMS$ETHNIC_label, ADAMS$Hispanic, useNA = "ifany")

#---- ****interview year ----
#table(ADAMS$AYEAR, useNA = "ifany")

#---- ****HRS cognition ----
# table(ADAMS$SELFCOG, useNA = "ifany")
# table(ADAMS$SELFCOG, useNA = "ifany")/nrow(ADAMS)

#---- ****marital status ----
#table(ADAMS$AAMARRD, useNA = "ifany")
ADAMS %<>% 
  mutate("AAMARRD_label" = case_when(AAMARRD == 1 ~ "Single", 
                                     AAMARRD == 2 ~ "Married or common law", 
                                     AAMARRD == 3 ~ "Divorced", 
                                     AAMARRD == 4 ~ "Separated", 
                                     AAMARRD == 5 ~ "Widow", 
                                     AAMARRD == 8 ~ "Don't Know")) %>%
  mutate("AAMARRD_collapsed_label" = 
           case_when(AAMARRD %in% c(1, 3, 4, 5) ~ "Not married/partnered", 
                     AAMARRD == 2 ~ "Married/partnered")) %>% 
  mutate("Married/partnered" = 
           ifelse(AAMARRD_collapsed_label == "Married/partnered", 1, 0))

# #Sanity check
# table(ADAMS$AAMARRD_label, ADAMS$AAMARRD_collapsed_label, "useNA" = "ifany")
# table(ADAMS$`Married/partnered`, useNA = "ifany")
# table(ADAMS$`Married/partnered`, useNA = "ifany")/nrow(ADAMS)

#---- ****age ----
#table(ADAMS$AAGE, useNA = "ifany")

#---- ****education ----
#table(ADAMS$EDYRS, useNA = "ifany")

#---- ****employment status ----
#table(ADAMS$AACURRWK, useNA = "ifany")
ADAMS %<>% 
  mutate("AACURRWK_label" = case_when(AACURRWK == 1 ~ "Working", 
                                      AACURRWK == 2 ~ "Retired", 
                                      AACURRWK == 3 ~ "Semi-retired", 
                                      AACURRWK == 4 ~ "Disabled", 
                                      AACURRWK == 5 ~ "Unemployed")) %>%
  mutate("AACURRWK_collapsed_label" = 
           case_when(AACURRWK %in% c(1, 3) ~ "Working", 
                     AACURRWK == 2 ~ "Retired", 
                     AACURRWK %in% c(4, 5) ~ "Not working")) %>% 
  mutate("Working" = ifelse(AACURRWK_collapsed_label == "Working", 1, 0), 
         "Retired" = ifelse(AACURRWK_collapsed_label == "Retired", 1, 0), 
         "Not working" = ifelse(AACURRWK_collapsed_label == "Not working", 1, 0))

# #Sanity check
# table(ADAMS$AACURRWK_label, ADAMS$AACURRWK_collapsed_label, "useNA" = "ifany")
# table(ADAMS$AACURRWK_collapsed_label, ADAMS$Working, useNA = "ifany")
# table(ADAMS$AACURRWK_collapsed_label, ADAMS$Retired, useNA = "ifany")
# table(ADAMS$AACURRWK_collapsed_label, ADAMS$`Not working`, useNA = "ifany")
# table(ADAMS$AACURRWK_collapsed_label, useNA = "ifany")/nrow(ADAMS)

#---- ****MMSE ----
#table(ADAMS$ANMSETOT, useNA = "ifany")
ADAMS %<>% 
  mutate_at(.vars = "ANMSETOT", function(x) ifelse(x > 30, NA, x)) %>% 
  #normalized MMSE
  mutate("ANMSETOT_norm" = normMMSE(ANMSETOT))

# #Sanity check
# hist(ADAMS$ANMSETOT)
# hist(ADAMS$ANMSETOT_norm)
# table(ADAMS$ANMSETOT_norm, useNA = "ifany")
# table(ADAMS$ANMSETOT_norm, useNA = "ifany")/nrow(ADAMS)

#---- ****BWC 20 and 86 ----
# table(ADAMS$ANBWC201, useNA = "ifany")
# table(ADAMS$ANBWC202, useNA = "ifany")
# table(ADAMS$ANBWC861, useNA = "ifany")
# table(ADAMS$ANBWC862, useNA = "ifany")

#Recode
ADAMS %<>% 
  mutate_at(.vars = c("ANBWC201", "ANBWC202", "ANBWC861", "ANBWC862"), 
            #Missing/refused  
            function(x) ifelse(x > 6, NA, x)) %>% 
  mutate_at(.vars = c("ANBWC201", "ANBWC202", "ANBWC861", "ANBWC862"), 
            #restart
            function(x) ifelse(x == 6, 0, x)) 

#Take the higher score
ADAMS %<>% mutate("ANBWC20" = pmax(ANBWC201, ANBWC202, na.rm = TRUE), 
                  "ANBWC86" = pmax(ANBWC861, ANBWC862, na.rm = TRUE))

# #Sanity check
# View(ADAMS[, c("ANBWC201", "ANBWC202", "ANBWC20")])
# View(ADAMS[, c("ANBWC861", "ANBWC862", "ANBWC86")])
# table(ADAMS$ANBWC20, useNA = "ifany")
# table(ADAMS$ANBWC86, useNA = "ifany")
# table(ADAMS$ANBWC20, useNA = "ifany")/nrow(ADAMS)
# table(ADAMS$ANBWC86, useNA = "ifany")/nrow(ADAMS)

#---- ****serial 7s ----
#table(ADAMS$ANSER7T, useNA = "ifany")

#Recode 
ADAMS %<>% mutate_at(.vars = c("ANSER7T"), 
                     #Missing/refused  
                     function(x) ifelse(x > 5, NA, x))

# #Sanity check
# table(ADAMS$ANSER7T, useNA = "ifany")
# table(ADAMS$ANSER7T, useNA = "ifany")/nrow(ADAMS)

#---- HCAP ----







