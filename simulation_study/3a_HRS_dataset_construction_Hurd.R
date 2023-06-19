#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "haven", "labelled", "magrittr", "tidyr")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- read data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS tracker ----
HRS_tracker_data_path <- 
  paste0(path_to_box, "data/HRS_tracker/trk2018_3/TRK2018TR_R.da")
HRS_tracker_dict_path <- 
  paste0(path_to_box, "data/HRS_tracker/trk2018_3/TRK2018TR_R.dct")

HRS_tracker <- read_da_dct(HRS_tracker_data_path, HRS_tracker_dict_path,
                           HHIDPN = "TRUE") %>% 
  #select variables: ID, 2016 Wave participation, HCAP selection, 
  # 2016 married/partnered status, sex/gender, age, race, ethnicity, 
  # years of education
  dplyr::select("HHIDPN", "PIWTYPE", "HCAP_SELECT", "PCOUPLE", "GENDER", "PAGE", 
                "RACE", "HISPANIC", "SCHLYRS") %>%
  #N = 43398
  #filter to those who completed 2016 Wave interview (N = 20911; dropped n = 22487)
  filter(PIWTYPE == 1) %>% 
  #get rid of leading 0
  mutate_at("HHIDPN", as.numeric) %>% 
  mutate_at("HHIDPN", as.character)

#---- **HRS Core ----
HRS_core_data_path <- 
  paste0(path_to_box, "data/HRS/Core_files/h16core/h16da/H16J_R.da")
HRS_core_dict_path <- 
  paste0(path_to_box, "data/HRS/Core_files/h16core/h16sta/H16J_R.dct")

HRS_core <- read_da_dct(HRS_core_data_path, HRS_core_dict_path, 
                        HHIDPN = "TRUE") %>% 
  #select variables: ID, employment status 
  dplyr::select("HHIDPN", "PJ005M1") %>% 
  #get rid of leading 0
  mutate_at("HHIDPN", as.numeric) %>% 
  mutate_at("HHIDPN", as.character)

#---- **RAND variables ----
hurd_waves <- 12
rand_waves <- 13

rand_variables <- 
  c("hhidpn",
    #Cognition (immediate word recall, delayed word recall, serial 7s, 
    # backwards count (20), object naming (scissors and cactus), 
    # president naming, subjective cognitive decline, date recall (day, month, 
    # year, day of the week))
    #Health and health behaviors (ever/never stroke, ever/never diabetes, 
    # ever/never CVD, ever/never hypertension, smokes now, 
    # drinking days per week, number of drinks per day, adl, iadl, iadlz, bmi)
    #Proxy indicator-- this is always 1 when cognitive test items are missing
    paste0("r", c(hurd_waves, rand_waves), "imrc"), 
    paste0("r", c(hurd_waves, rand_waves), "dlrc"), 
    paste0("r", c(hurd_waves, rand_waves), "ser7"), 
    paste0("r", c(hurd_waves, rand_waves), "bwc20"), 
    paste0("r", rand_waves, "cogtot"), 
    paste0("r", c(hurd_waves, rand_waves), "cact"), 
    paste0("r", c(hurd_waves, rand_waves), "scis"), 
    paste0("r", c(hurd_waves, rand_waves), "pres"), 
    paste0("r", rand_waves, "pstmem"), 
    paste0("r", c(hurd_waves, rand_waves), "dy"),
    paste0("r", c(hurd_waves, rand_waves), "mo"),
    paste0("r", c(hurd_waves, rand_waves), "yr"),
    paste0("r", c(hurd_waves, rand_waves), "dw"),
    paste0("r", rand_waves, "stroke"), 
    paste0("r", rand_waves, "diabe"), paste0("r", rand_waves, "hearte"), 
    paste0("r", rand_waves, "hibpe"), paste0("r", rand_waves, "smoken"), 
    paste0("r", rand_waves, "drinkd"), paste0("r", rand_waves, "drinkn"), 
    paste0("r", c(hurd_waves, rand_waves), "adla"), 
    paste0("r", rand_waves, "iadla"), 
    paste0("r", c(hurd_waves, rand_waves), "iadlza"),
    paste0("r", rand_waves, "bmi"), paste0("r", rand_waves, "proxy"))

RAND <- read_dta(paste0(path_to_box, "data/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character) %>% rename("HHIDPN" = "hhidpn")

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- join data ----
HRS <- left_join(HRS_tracker, HRS_core, by = "HHIDPN") %>% 
  left_join(., RAND, by = "HHIDPN")

#---- clean data ----
#---- **HCAP selection ----
# table(HRS$HCAP_SELECT, useNA = "ifany")

#change 2 = no to 0 = no
HRS %<>% mutate_at(.vars = "HCAP_SELECT", function(x) ifelse(x == 2, 0, 1))

#---- **age ----
# table(HRS$PAGE, useNA = "ifany")

#restrict to 70+ (N = 7337, dropped n = 13,574)
HRS %<>% filter(PAGE >= 70)

#Hurd model: create age categories
HRS %<>% 
  mutate("PAGE_cat" = 
           cut(PAGE, breaks = c(min(HRS$PAGE), 75, 80, 85, 90, max(HRS$PAGE)), 
               include.lowest = TRUE, right = FALSE))

# #Sanity check
# table(HRS$PAGE_cat, useNA = "ifany")

#---- **2016 couple status ----
# table(HRS$PCOUPLE, useNA = "ifany")

#change 5 = no to 0 = no
HRS %<>% mutate_at(.vars = "PCOUPLE", function(x) ifelse(x == 5, 0, 1))

# #Sanity check
# table(HRS$PCOUPLE, useNA = "ifany")

#---- **sex/gender ----
#table(HRS$GENDER, useNA = "ifany")
HRS %<>% 
  mutate(GENDER_label = as.factor(ifelse(GENDER == 1, "Male", "Female"))) %>% 
  mutate("Female" = ifelse(GENDER_label == "Female", 1, 0))

# #Sanity check
# table(HRS$GENDER, HRS$GENDER_label)
# table(HRS$Female, useNA = "ifany")

#---- **race ----
#table(HRS$RACE, useNA = "ifany")
HRS %<>% 
  mutate(RACE_label = as.factor(case_when(RACE == 1 ~ "White", 
                                          RACE == 2 ~ "Black", 
                                          RACE == 7 ~ "Other"))) %>% 
  mutate("RACE_White" = ifelse(RACE_label == "White", 1, 0), 
         "RACE_Black" = ifelse(RACE_label == "Black", 1, 0), 
         "RACE_Other" = ifelse(RACE_label == "Other", 1, 0))

# #Sanity check
# table(HRS$RACE_label, HRS$RACE_White, useNA = "ifany")
# table(HRS$RACE_label, HRS$RACE_Black, useNA = "ifany")
# table(HRS$RACE_label, HRS$RACE_Other, useNA = "ifany")

#---- **hispanic ----
#table(HRS$HISPANIC, useNA = "ifany")

#change 5 = no to 0 = no, combine all hispanic types and set 0 = not obtained to NA
HRS %<>% mutate("HISPANIC_indicator" = case_when(HISPANIC %in% c(1, 2, 3) ~ 1, 
                                                 HISPANIC == 5 ~ 0))

# #Sanity check
# table(HRS$HISPANIC, HRS$HISPANIC_indicator, useNA = "ifany")

#---- **race/ethnicity ----
HRS %<>% 
  mutate(ETHNIC_label = 
           case_when(HISPANIC_indicator == 1 ~ "Hispanic",
                     HISPANIC_indicator == 0 & RACE_label == "White" ~ "White",
                     HISPANIC_indicator == 0 & RACE_label == "Black" ~ "Black", 
                     HISPANIC_indicator == 0 & RACE_label == "Other" ~ "Other")) %>% 
  mutate("White" = ifelse(ETHNIC_label == "White", 1, 0), 
         "Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0), 
         "Other" = ifelse(ETHNIC_label == "Other", 1, 0))

# #Sanity check
# table(HRS$ETHNIC_label, HRS$White, useNA = "ifany")
# table(HRS$ETHNIC_label, HRS$Black, useNA = "ifany")
# table(HRS$ETHNIC_label, HRS$Hispanic, useNA = "ifany")
# table(HRS$ETHNIC_label, HRS$Other, useNA = "ifany")

#drop people missing race/ethnicity (N = 7335, dropped n = 2)
HRS %<>% drop_na(ETHNIC_label)

#restrict to White, Black, Hispanic (to align with ADAMS data) 
# (N = 7166, dropped n = 169)
HRS %<>% filter(Other == 0)

#---- **years of education ----
HRS %<>% mutate("SCHLYRS" = ifelse(SCHLYRS > 17, NA, SCHLYRS))

# #Sanity check
# table(HRS$SCHLYRS, useNA = "ifany")

#drop people missing education data (N = 7162, dropped n = 4)
HRS %<>% drop_na(SCHLYRS)

#education categories
HRS %<>% 
  mutate("edu_cat" = 
           cut(SCHLYRS, breaks = c(min(HRS$SCHLYRS), 11, 12, max(HRS$SCHLYRS)), 
               include.lowest = TRUE, right = TRUE, 
               labels = c("LTHS", "HSGED", "GTHS")))

# #Sanity Check
# table(HRS$edu_cat)
# table(HRS$SCHLYRS)
# sum(HRS$SCHLYRS > 12)
# sum(HRS$SCHLYRS < 12)

#---- **employment status ----
#table(HRS$PJ005M1, useNA = "ifany")
HRS %<>% 
  mutate("PJ005M1_label" = 
           case_when(PJ005M1 == 1 ~ "Working now", 
                     PJ005M1 == 2 ~ "Unemployed and looking for work", 
                     PJ005M1 == 3 ~ "Temporarily laid off", 
                     PJ005M1 == 4 ~ "Disabled", 
                     PJ005M1 == 5 ~ "Retired", 
                     PJ005M1 == 6 ~ "Homemaker", 
                     PJ005M1 == 7 ~ "Other", 
                     PJ005M1 == 8 ~ "On leave")) %>%
  mutate("PJ005M1_collapsed_label" = 
           case_when(PJ005M1 == 1 ~ "Working", 
                     PJ005M1 == 5 ~ "Retired", 
                     PJ005M1 %in% c(2, 3, 4, 6, 7, 8) ~ "Not working")) %>% 
  mutate("Working" = ifelse(PJ005M1_collapsed_label == "Working", 1, 0), 
         "Retired" = ifelse(PJ005M1_collapsed_label == "Retired", 1, 0), 
         "Not working" = ifelse(PJ005M1_collapsed_label == "Not working", 1, 0))

# #Sanity check
# table(HRS$PJ005M1, HRS$PJ005M1_label, "useNA" = "ifany")
# table(HRS$PJ005M1_label, HRS$PJ005M1_collapsed_label, "useNA" = "ifany")
# table(HRS$PJ005M1_collapsed_label, HRS$Working, useNA = "ifany")
# table(HRS$PJ005M1_collapsed_label, HRS$Retired, useNA = "ifany")
# table(HRS$PJ005M1_collapsed_label, HRS$`Not working`, useNA = "ifany")

#drop people missing employment status (N = 7138, dropped n = 24)
HRS %<>% drop_na(Working)

#---- **subjective cognitive change ----
HRS %<>% mutate("subj_cog_better" = ifelse(r13pstmem == 1, 1, 0), 
                "subj_cog_same" = ifelse(r13pstmem == 2, 1, 0), 
                "subj_cog_worse" = ifelse(r13pstmem == 3, 1, 0))

#---- **immediate word recall ----
# table(HRS$r13imrc, useNA = "ifany")

#---- **delayed word recall ----
# table(HRS$r13dlrc, useNA = "ifany")

#---- **serial 7s ----
# table(HRS$r13ser7, useNA = "ifany")

#---- **backwards count 20 ----
# table(HRS$r13bwc20, useNA = "ifany")

#Hurd model: raw score
HRS %<>% mutate("r12bwc20_raw" = r12bwc20)
HRS %<>% mutate("r13bwc20_raw" = r13bwc20)

#count corrects on second try as correct
HRS %<>% mutate_at(.vars = "r13bwc20", function(x) ifelse(x >= 1, 1, 0))

# #Sanity check
# table(HRS$r13bwc20_raw)

#---- **total cognition ----
# table(HRS$r13cogtot, useNA = "ifany")

#---- **ADLs ----
# table(HRS$r13adla, useNA = "ifany")

#---- **IADLs ----
# table(HRS$r13iadla, useNA = "ifany")

#---- **bmi ----
# table(HRS$r13bmi, useNA = "ifany")

#---- derived variables ----
#---- **drinking behavior ----
HRS %<>% 
  mutate("drinks_per_week" = r13drinkd*r13drinkn) %>%
  mutate("drink_cat" = 
           case_when(drinks_per_week == 0 ~ 1,
                     Female == 0 & 
                       (drinks_per_week >= 1 & drinks_per_week < 14) ~ 2,
                     Female == 1 &
                       (drinks_per_week >= 1 & drinks_per_week < 7) ~ 2,
                     Female == 0 &
                       (drinks_per_week >= 14 | r13drinkn >= 4) ~ 3,
                     Female == 1 &
                       (drinks_per_week >= 7 | r13drinkn >= 3) ~ 3)) %>% 
  mutate("drink_cat_label" = 
           case_when(drink_cat == 1 ~ "No Drinking", 
                     drink_cat == 2 ~ "Moderate Drinking", 
                     drink_cat == 3 ~ "Heavy Drinking")) %>% 
  mutate("no_drinking" = ifelse(drink_cat_label == "No Drinking", 1, 0), 
         "moderate_drinking" = 
           ifelse(drink_cat_label == "Moderate Drinking", 1, 0), 
         "heavy_drinking" = 
           ifelse(drink_cat_label == "Heavy Drinking", 1, 0))

#remove people missing alcohol consumption (N = 7093, dropped n = 45)
HRS %<>% drop_na(no_drinking)

#---- summarize missingness ----
colMeans(is.na(HRS))[which(colMeans(is.na(HRS)) > 0)]

#---- **filter participants missing health variables ----
health_vars <- c(paste0("r13", c("bmi", "smoken", "hibpe", "diabe", "hearte", 
                                 "stroke")))

#remove people missing any health variables (N = 6886, dropped n = 207)
HRS %<>% drop_na(all_of(health_vars))

# #Sanity check
# colMeans(is.na(HRS))[which(colMeans(is.na(HRS)) > 0)]

#---- **summarize missingness on any cognitive assessment ----
#There are 752 participants missing at least one cognitive assessment in HRS
subset <- names(colMeans(is.na(HRS))[which(colMeans(is.na(HRS)) > 0)])
remove_vars <- c("r13drinkd", "r13pstmem", "memimp16", "RACE_label", 
                 "RACE_White", "RACE_Black", "RACE_Other", "drinks_per_week")
cog_vars <- subset[-which(subset %in% remove_vars)]

nrow(HRS) - nrow(HRS %>% drop_na(all_of(cog_vars)))

#---- **CC cognitive assessments ----
HRS %<>% drop_na(all_of(cog_vars))

# #Sanity check
# colMeans(is.na(HRS))[which(colMeans(is.na(HRS)) > 0)]

#---- Hurd analysis vars ----
#---- **date summary ----
for(wave in c(12, 13)){
  HRS %<>% 
    mutate(!!paste0("r", wave, "date_recall") := 
             rowSums(across(!!paste0("r", wave, c("mo", "dy", "yr", "dw")))))
}

# #Sanity Check
# head(HRS$r12dw + HRS$r12dy + HRS$r12yr + HRS$r12mo)
# head(HRS$r12date_recall)
# 
# head(HRS$r13dw + HRS$r13dy + HRS$r13yr + HRS$r13mo)
# head(HRS$r13date_recall)

#---- **change scores ----
change_score_vars <- c("adla", "iadlza", "date_recall", "bwc20_raw", "ser7", 
                       "scis", "cact", "pres", "imrc", "dlrc")

for (var in change_score_vars){
  HRS %<>% mutate(!!paste0(var, "_change") := 
                    !!sym(paste0("r13", var)) - !!sym(paste0("r12", var)))
}

# #Sanity check
# table(HRS$r13adla - HRS$r12adla)
# table(HRS$adla_change)
# 
# table(HRS$r13bwc20_raw - HRS$r12bwc20_raw)
# table(HRS$bwc20_raw_change)

#---- algorithms ----
#---- **LKW ----
#LWK sum score
HRS %<>% 
  mutate("LKW_sum_score" =  
           rowSums(across(paste0("r13", c("imrc", "dlrc", "ser7", "bwc20_raw")))))

# #Sanity check
# table(HRS$r13imrc + HRS$r13dlrc + HRS$r13ser7 + HRS$r13bwc20_raw)
# table(HRS$LKW_sum_score)

#LWK classification
HRS %<>% mutate("dem_LKW" = ifelse(LKW_sum_score <= 6, 1, 0))

# #Sanity check
# sum(HRS$LKW_sum_score <= 6)
# table(HRS$dem_LKW)

#---- rename columns ----
variable_labels <- 
  variable_labels[which(variable_labels$HRS %in% colnames(HRS)), ]

HRS %<>% 
  rename_at(vars(variable_labels$HRS), ~ variable_labels$data_label)

#---- save dataset ----
HRS %>% 
  dplyr::select(-one_of("GENDER", "RACE", "HISPANIC", "PIWTYPE", "PJ005M1", 
                        "r13drinkd", "r13drinkn", "r13pstmem", "GENDER_label", 
                        "RACE_label", "RACE_White", "RACE_Black", "RACE_Other", 
                        "HISPANIC_indicator", "ETHNIC_label", "PJ005M1_label", 
                        "PJ005M1_collapsed_label", "drinks_per_week", 
                        "drink_cat", "drink_cat_label")) %>% 
  write_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv"))
