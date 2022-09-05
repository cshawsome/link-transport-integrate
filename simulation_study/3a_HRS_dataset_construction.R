#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "haven", "labelled", "magrittr")

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
  #select variables: ID, 2016 Wave participation, 2016 married/partnered status, 
  # sex/gender, age, race, ethnicity, years of education
  dplyr::select("HHIDPN", "PIWTYPE", "PCOUPLE", "GENDER", "PAGE", "RACE", 
                "HISPANIC", "SCHLYRS") %>%
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
rand_waves <- 13
rand_variables <- 
  c("hhidpn",
    #Cognition (immediate word recall, delayed word recall, serial 7s, 
    # backwards count (20), object naming (scissors and cactus), 
    # president naming, subjective cognitive decline)
    #Health and health behaviors (ever/never stroke, ever/never diabetes, 
    # ever/never CVD, ever/never hypertension, smokes now, 
    # drinking days per week, number of drinks per day, adl, iadl, bmi)
    #Proxy indicator-- this is always 1 when cognitive test items are missing
    paste0("r", rand_waves, "imrc"), paste0("r", rand_waves, "dlrc"), 
    paste0("r", rand_waves, "ser7"), paste0("r", rand_waves, "bwc20"), 
    paste0("r", rand_waves, "cogtot"), paste0("r", rand_waves, "cact"), 
    paste0("r", rand_waves, "scis"), paste0("r", rand_waves, "pres"), 
    paste0("r", rand_waves, "pstmem"), paste0("r", rand_waves, "stroke"), 
    paste0("r", rand_waves, "diabe"), paste0("r", rand_waves, "hearte"), 
    paste0("r", rand_waves, "hibpe"), paste0("r", rand_waves, "smoken"), 
    paste0("r", rand_waves, "drinkd"), paste0("r", rand_waves, "drinkn"), 
    paste0("r", rand_waves, "adla"), paste0("r", rand_waves, "iadla"), 
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
#---- **age ----
# table(HRS$PAGE, useNA = "ifany")

#restrict to 70+ (N = 7377, dropped n = 13,534)
HRS %<>% filter(PAGE >= 70)

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

#restrict to White, Black, Hispanic (to align with ADAMS data) 
# (N = 7166, dropped n = 171)

HRS %<>% filter(Other == 0)

#---- **years of education ----
HRS %<>% mutate("SCHLYRS" = ifelse(SCHLYRS > 17, NA, SCHLYRS))

# #Sanity check
# table(HRS$SCHLYRS, useNA = "ifany")

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

#count corrects on second try as correct
HRS %<>% mutate_at(.vars = "r13bwc20", function(x) ifelse(x >= 1, 1, 0))

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

#---- summarize missingness ----
colMeans(is.na(HRS))

#---- CC HRS ----
#keep those missing co
HRS_CC <- na.omit(HRS)

#---- rename columns ----
variable_labels <- 
  variable_labels[which(variable_labels$HRS %in% colnames(HRS_CC)), ]

HRS_CC %<>% 
  rename_at(vars(variable_labels$HRS), ~ variable_labels$data_label)

#---- derived variable bins ----
HRS_CC %<>% 
  mutate(
    #---- **age ----
    "age_cat" = cut(age, breaks = c(min(HRS_CC$age), 75, 80, 85, 90, 
                                    max(HRS_CC$age)), 
                    include.lowest = TRUE, right = FALSE), 
    #---- **education ----
    "edyrs_cat" = cut(edyrs, 
                      breaks = c(0, 0.9, 11, 12, 15, 16, 
                                 max(HRS_CC$edyrs, na.rm = TRUE)),
                      labels = c("none", "less than HS", "HS", "some college", 
                                 "college", "graduate studies"),
                      include.lowest = TRUE, right = TRUE), 
    #---- **immediate word recall ----
    "immrc_cat" = cut_number(immrc, n = 5), 
    #---- **serial 7s ----
    "ser7_cat" = cut_number(ser7, n = 2), 
    #---- **hrs cognition ----
    "hrs_cog_cat" = cut_number(hrs_cog, n = 5), 
    #---- **adl ----
    "adl_cat" = ifelse(adl == 0, "none", "any"), 
    #---- **iadl ----
    "iadl_cat" = ifelse(iadl == 0, "none", "any"), 
    #---- **bmi ----
    "bmi_cat" = cut_number(bmi, n = 5))

# #Sanity check
# #sum should be 6314
# check_vars <- c("age_cat", "edyrs_cat", "immrc_cat", "ser7_cat", "hrs_cog_cat", 
#                 "adl_cat", "iadl_cat", "bmi_cat")
# 
# for(var in check_vars){
#   print(sum(table(HRS_CC[, var])))
# }

#---- save datasets ----
#clean data
HRS %>% write_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv"))

#pared-down analytic data
remove <- c("PIWTYPE", "RACE", "RACE_label", "RACE_White", "RACE_Black", 
            "RACE_Other", "HISPANIC", "HISPANIC_indicator", "ETHNIC_label", 
            "Other", "GENDER", "GENDER_label", "PJ005M1", "PJ005M1_label", 
            "PJ005M1_collapsed_label", "drinks_per_week", "drink_cat", 
            "drink_cat_label", "r13drinkd", "r13drinkn", "r13pstmem")

HRS_CC %>% dplyr::select(-all_of(remove)) %>%
  write_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))
