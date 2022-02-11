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
  #select variables: ID, 2016 Wave participation, 2016 HCAP selection, 
  # 2016 married/partnered status, sex/gender, age, race, ethnicity, 
  dplyr::select("HHIDPN", "PIWTYPE", "HCAP_SELECT", "PCOUPLE", "GENDER", "PAGE", 
                "RACE", "HISPANIC", "SCHLYRS") %>%
  #N = 43398
  #filter to those who completed 2016 Wave interview (N = 20911; dropped n = 22487)
  filter(PIWTYPE == 1)

#---- **HRS Core ----
HRS_core_data_path <- 
  paste0(path_to_box, "data/HRS/Core_files/h16core/h16da/H16J_R.da")
HRS_core_dict_path <- 
  paste0(path_to_box, "data/HRS/Core_files/h16core/h16sta/H16J_R.dct")

HRS_core <- read_da_dct(HRS_core_data_path, HRS_core_dict_path, 
                        HHIDPN = "TRUE") %>% 
  #select variables: ID, employment status 
  dplyr::select("HHIDPN", "PJ005M1")

#---- **RAND variables ----
rand_waves <- 13 #Corresponding to 2016 HRS
rand_variables <- 
  c("hhidpn",
    #Cognition (object naming (scissors and cactus), president naming, 
    # subjective cognitive decline)
    #Health and health behaviors (ever/never stroke, ever/never diabetes, 
    # ever/never CVD, ever/never hypertension, smokes now, 
    # drinking days per week, number of drinks per day)
    paste0("r", rand_waves, "cact"), paste0("r", rand_waves, "scis"), 
    paste0("r", rand_waves, "pres"), paste0("r", rand_waves, "pstmem"),
    paste0("r", rand_waves, "stroke"), paste0("r", rand_waves, "diabe"), 
    paste0("r", rand_waves, "hearte"), paste0("r", rand_waves, "hibpe"), 
    paste0("r", rand_waves, "smoken"), paste0("r", rand_waves, "drinkd"), 
    paste0("r", rand_waves, "drinkn"))

RAND <- read_dta(paste0(path_to_box, "data/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character) %>% rename("HHIDPN" = "hhidpn")

#Remove labeled data format
val_labels(RAND) <- NULL

#---- join data ----
HRS <- left_join(HRS_tracker, HRS_core, by = "HHIDPN") %>% 
  left_join(., RAND, by = "HHIDPN")

#---- clean data ----
#---- **2016 HCAP selection ----
# table(HRS$HCAP_SELECT, useNA = "ifany")
# table(HRS$HCAP_SELECT, useNA = "ifany")/nrow(HRS)

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

#---- **age ----
# table(HRS$PAGE, useNA = "ifany")

#restrict to 70+ (N = 7377, dropped n = 13,534)
HRS %<>% filter(PAGE >= 70)

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

#---- **employment status ----
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

#---- **health and health behaviors ----
# table(ADAMS$AYEAR, useNA = "ifany")
#For repeated measures, want to take the wave most representative of ADAMS wave
wave_updated_vars <- c("stroke", "hibpe", "diabe", "hearte", "bmi", 
                       "iadla", "adla", "smoken", "drinkd", "drinkn")

for(var in wave_updated_vars){
  ADAMS %<>% 
    mutate(!!paste0("A", var) := 
             case_when(AYEAR %in% c(2001, 2002) ~ !!sym(paste0("r5", var)), 
                       AYEAR %in% c(2003, 2004) ~ !!sym(paste0("r6", var))))
}

# #Sanity check
# View(ADAMS[, c("AYEAR", paste0("r", seq(5, 6), "stroke"), "Astroke")] %>% 
#        filter(!is.na(Astroke)))
# colnames(ADAMS)

#---- **filter: missing all neurospych + general cognitive measures ----
neuro_cog_measures <- c("SELFCOG", "ANMSETOT_norm", "ANBWC20", "ANBWC86", 
                        "ANSER7T", "ANSCISOR", "ANCACTUS", "ANPRES", "ANAFTOT", 
                        "ANBNTTOT", "ANIMMCR", "ANDELCOR", "ANRECYES", "ANRECNO", 
                        "ANWM1TOT", "ANWM2TOT", "ANCPTOT", "ANRCPTOT", 
                        "ANTMASEC")

ADAMS %<>% 
  mutate("num_cog_measures" = 
           rowSums(!is.na(ADAMS %>% 
                            dplyr::select(all_of(neuro_cog_measures))))) %>% 
  #N = 826; dropped n = 30
  filter(num_cog_measures > 0)

# #Sanity check
# table(ADAMS$num_cog_measures, useNA = "ifany")

#---- **summarize missingness ----
colMeans(is.na(ADAMS))

#---- imputation-specific variables ----
#---- **ADAMS proxy type ----
# table(ADAMS$AGQ101, useNA = "ifany")
ADAMS %<>% 
  mutate("proxy_type_label" = case_when(AGQ101 == 1 ~ "Spouse", 
                                        AGQ101 == 2 ~ "Child", 
                                        AGQ101 == 3 ~ "Grandchild", 
                                        AGQ101 == 4 ~ "Professional", 
                                        AGQ101 == 5 ~ "Child-in-law", 
                                        AGQ101 == 6 ~ "Sibling",
                                        AGQ101 == 7 ~ "Niece/Nephew",
                                        AGQ101 == 8 ~ "Sibling of Spouse", 
                                        AGQ101 == 9 ~ "Parent/Parent-in-law",
                                        AGQ101 == 10 ~ "Other Relative",
                                        AGQ101 == 13 ~ "Other", 
                                        AGQ101 == 15 ~ "Blank")) %>%
  mutate("proxy_type_collapsed_label" = 
           case_when(AGQ101 == 1 ~ "Spouse", 
                     AGQ101 == 2 ~ "Child", 
                     AGQ101 %in% c(3, 5, 6, 7, 8, 9, 10) ~ "Other Relative", 
                     AGQ101 %in% c(4, 13) ~ "Other")) %>% 
  mutate("proxy_Spouse" = 
           ifelse(proxy_type_collapsed_label == "Spouse", 1, 0), 
         "proxy_Child" = 
           ifelse(proxy_type_collapsed_label == "Child", 1, 0), 
         "proxy_Other_Relative" = 
           ifelse(proxy_type_collapsed_label == "Other Relative", 1, 0),
         "proxy_Other" = 
           ifelse(proxy_type_collapsed_label == "Other", 1, 0))

# #Sanity check
# table(ADAMS$proxy_type_label, ADAMS$proxy_type_collapsed_label, useNA = "ifany")
# table(ADAMS$proxy_type_collapsed_label, ADAMS$proxy_Spouse, useNA = "ifany")
# table(ADAMS$proxy_type_collapsed_label, ADAMS$proxy_Child, useNA = "ifany")
# table(ADAMS$proxy_type_collapsed_label, ADAMS$proxy_Other, useNA = "ifany")
# table(ADAMS$proxy_type_collapsed_label, ADAMS$proxy_Other_Relative, 
#       useNA = "ifany")

#---- **HRS married/partnered status ----
# table(ADAMS$r5mpart, useNA = "ifany")

#---- **HRS BWC 20 and 86 ----
# table(ADAMS$r5bwc20, useNA = "ifany")
# table(ADAMS$r6bwc86, useNA = "ifany")
# table(ADAMS$ANBWC20, useNA = "ifany")

bwc_vars <- c(paste0("r", cog_test_waves, "bwc20"), 
              paste0("r", seq(5, 6), "bwc86"))
ADAMS %<>% 
  mutate_at(.vars = all_of(bwc_vars), 
            #Correct on 2nd try = correct
            function(x) ifelse(x == 2, 1, x))

# #Sanity check
# table(ADAMS$r5bwc20, useNA = "ifany")
# table(ADAMS$r6bwc86, useNA = "ifany")

#---- **HRS object naming: cactus, scissors ----
# table(ADAMS$ANCACTUS, useNA = "ifany")
# table(ADAMS$ANSCISOR, useNA = "ifany")
# table(ADAMS$r5scis, useNA = "ifany")
# table(ADAMS$r5cact, useNA = "ifany")

#---- **HRS total cognition ----
# table(ADAMS$SELFCOG, useNA = "ifany")
# table(ADAMS$r5cogtot, useNA = "ifany")

#---- **HRS 10-word recall (immediate and delayed) ----
# table(ADAMS$ANIMMCR1, useNA = "ifany")
# table(ADAMS$ANDELCOR, useNA = "ifany")
# table(ADAMS$r5imrc, useNA = "ifany")
# table(ADAMS$r5dlrc, useNA = "ifany")

#---- **HRS President naming ----
# table(ADAMS$ANPRES, useNA = "ifany")
# table(ADAMS$r5pres, useNA = "ifany")

#---- **HRS serial 7s ----
# table(ADAMS$ANSER7T, useNA = "ifany")
# table(ADAMS$r5ser7, useNA = "ifany")

#---- **HRS subjective cognitive decline ----
# table(ADAMS$r5pstmem, useNA = "ifany")
for(wave in cog_test_waves){
  ADAMS %<>%  
    mutate(!!paste0("r", wave, "pstmem_Better") := 
             ifelse(!!sym(paste0("r", wave, "pstmem")) == 1, 1, 0), 
           !!paste0("r", wave, "pstmem_Same") := 
             ifelse(!!sym(paste0("r", wave, "pstmem")) == 2, 1, 0), 
           !!paste0("r", wave, "pstmem_Worse") := 
             ifelse(!!sym(paste0("r", wave, "pstmem")) == 3, 1, 0))
}

# #Sanity check
# table(ADAMS$r5pstmem, ADAMS$r5pstmem_Better, useNA = "ifany")
# table(ADAMS$r5pstmem, ADAMS$r5pstmem_Same, useNA = "ifany")
# table(ADAMS$r5pstmem, ADAMS$r5pstmem_Worse, useNA = "ifany")
# table(ADAMS$r6pstmem, ADAMS$r6pstmem_Better, useNA = "ifany")

#---- sanity check health vars ----
#---- **bmi ----
bmi_vars <- paste0("r", rand_waves, "bmi")
# lapply(ADAMS[, bmi_vars], FUN = hist)
# lapply(ADAMS[, bmi_vars], FUN = table)
# lapply(ADAMS[, bmi_vars], FUN = which.max)
# 
# #obs 297 is an outlier, so check height and weight data
height_vars <- paste0("r", rand_waves, "height")
# weight_vars <- paste0("r", rand_waves, "weight")
# ADAMS[297, height_vars] # 3.08 ft
# mean(unlist(ADAMS[297, weight_vars])) # 139.5 lbs

#set this person to missing for all their bmi data and all height data
ADAMS[297, bmi_vars] <- NA
ADAMS[297, height_vars] <- NA

#---- **height ----
# height_vars <- paste0("r", rand_waves, "height")
# lapply(ADAMS[, height_vars], FUN = hist)
# lapply(ADAMS[, height_vars], FUN = table)
# lapply(ADAMS[, height_vars], FUN = which.min)
# lapply(ADAMS[, height_vars], FUN = which.max)
# 
# #obs 574 may be a low outlier, so check height data
# ADAMS[574, height_vars] # 4.4 ft, seems fine
# 
# #obs 701 may be a high outlier, so check height data
# ADAMS[701, height_vars] # 6.3 ft, seems fine

#---- **weight ----
weight_vars <- paste0("r", rand_waves, "weight")
# lapply(ADAMS[, weight_vars], FUN = hist)
# lapply(ADAMS[, weight_vars], FUN = table)
# lapply(ADAMS[, weight_vars], FUN = function(x) min(x, na.rm = TRUE))
# lapply(ADAMS[, weight_vars], FUN = function(x) max(x, na.rm = TRUE))
# lapply(ADAMS[, weight_vars], FUN = which.min)
# 
# #obs 409 may be a low outlier in wave 7
# ADAMS[409, weight_vars] # 57 lbs, seems like an error

#set this person to missing for all their bmi data and all weight data
ADAMS[409, c(bmi_vars, weight_vars)] <- NA

#---- **drinking days ----
# drinkd_vars <- paste0("r", rand_waves, "drinkd")
# lapply(ADAMS[, drinkd_vars], FUN = table)

#---- **drinks per day ----
# drinkn_vars <- paste0("r", rand_waves, "drinkn")
# lapply(ADAMS[, drinkn_vars], FUN = table)
# 
# #who is this person with 16 drinks per day?? 
# # I'm going to leave them in, cuz who's to say **shrug**
# lapply(ADAMS[, drinkn_vars], FUN = which.max)
# ADAMS[651, c(drinkd_vars, drinkn_vars)]

#---- save datasets ----
#clean data
ADAMS %>% write_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_clean.csv"))

#pared-down analytic data
remove <- c("AASSESS", "AACURRWK", "AACURRWK_collapsed_label", "AACURRWK_label", 
            "AAMARRD", "AAMARRD_collapsed_label", "AAMARRD_label", "Adem_cat", 
            "ADFDX1", "ANSMEM2", "ANSMEM2_collapsed_label", "ANSMEM2_label", 
            paste0("AGQ", c(seq(14, 29), 101)), "avg_proxy_cog", "ETHNIC",
            "avg_proxy_cog_label", "avg_proxy_cog_collapsed_label", "GENDER",
            paste0("ANBWC", c("201", "202", "861", "862")), "GENDER_label",
            paste0("ANIMMCR", seq(1, 3)), "ETHNIC_label", "num_cog_measures", 
            "proxy_type_label", "proxy_type_collapsed_label", 
            paste0("r", cog_test_waves, "pstmem"), 
            paste0("r", rand_waves, "bmi"))

ADAMS %>% dplyr::select(-all_of(remove)) %>% 
  write_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))
