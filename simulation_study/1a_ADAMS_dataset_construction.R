#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "haven", "labelled", "magrittr", "NormPsy")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- read data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **ADAMS tracker ----
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
                "AAGE", "ETHNIC", "EDYRS", "AACURRWK") %>%
  #N = 1170
  #filter to those who completed Wave A assessment (N = 856; dropped n = 314)
  filter(AASSESS == 1)

#---- **neuropsych measures ----
ADAMS_neuropsych_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
ADAMS_neuropsych_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych <- 
  read_da_dct(ADAMS_neuropsych_data_path, ADAMS_neuropsych_dict_path,
              HHIDPN = "TRUE") %>% 
  #select variables: ID, MMSE, backwards count (20), serial 7s, 
  # item naming (scissors), item naming (cactus), President naming, 
  # animal naming, immediate word recall, delayed word recall, 
  # word list recognition (yes), word list recognition (no), 
  # immediate story recall, delayed story recall, 
  # immediate constructional praxis, delayed constructional praxis, trails A,
  # subjective cognitive change
  dplyr::select("HHIDPN", "ANMSETOT", "ANBWC201", "ANBWC202", "ANSER7T", 
                "ANSCISOR", "ANCACTUS", "ANPRES", "ANAFTOT", "ANIMMCR1", 
                "ANIMMCR2", "ANIMMCR3", "ANDELCOR", "ANRECYES", "ANRECNO", 
                "ANWM1TOT", "ANWM2TOT", "ANCPTOT", "ANRCPTOT", "ANTMASEC", 
                "ANSMEM2")

#---- **dementia dx ----
ADAMS_demdx_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1ada/ADAMS1AD_R.da")
ADAMS_demdx_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1asta/ADAMS1AD_R.dct")

ADAMS_demdx <- read_da_dct(ADAMS_demdx_data_path, ADAMS_demdx_dict_path,
                           HHIDPN = "TRUE") %>% 
  dplyr::select("HHIDPN", "ADFDX1")

#---- **proxy measures ----
ADAMS_proxy_data_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1ada/ADAMS1AG_R.da")
ADAMS_proxy_dict_path <- 
  paste0(path_to_box, "data/ADAMS/adams1a/adams1asta/ADAMS1AG_R.dct")

ADAMS_proxy <- read_da_dct(ADAMS_proxy_data_path, ADAMS_proxy_dict_path,
                           HHIDPN = "TRUE") %>% 
  #select variables: IQCODE items, proxy type
  dplyr::select("HHIDPN", paste0("AGQ", c(seq(14, 29), 101)))

#---- **RAND variables ----
rand_waves <- seq(4, 7, by = 1) #Corresponding to ADAMS + imputation
cog_test_waves <- seq(5, 7, by = 1) #Corresponding to ADAMS + imputation
rand_variables <- 
  c("hhidpn",
    #Sociodemographics (marital status)
    #Health and health behaviors (ever/never stroke, ever/never
    #   hypertension, ever/never diabetes, ever/never cvd, BMI, height, weight, 
    #   IADLs, ADLs, smokes now, number days drinking per week, number drinks/day) 
    #Cognitive tests (backwards count 20, item naming (scissors and 
    #   cactus), immediate and delayed word recall, serial 7s, President naming, 
    #   subjective cognitive decline, total cognition score)
    paste0("r", cog_test_waves, "mpart"),
    paste0("r", rand_waves, "stroke"), paste0("r", rand_waves, "hibpe"), 
    paste0("r", rand_waves, "diabe"), paste0("r", rand_waves, "hearte"),
    paste0("r", rand_waves, "bmi"), paste0("r", rand_waves, "height"), 
    paste0("r", rand_waves, "weight"), paste0("r", rand_waves, "iadla"), 
    paste0("r", rand_waves, "adla"), paste0("r", rand_waves, "smoken"), 
    paste0("r", rand_waves, "drinkd"), paste0("r", rand_waves, "drinkn"), 
    paste0("r", cog_test_waves, "bwc20"), paste0("r", cog_test_waves, "scis"), 
    paste0("r", cog_test_waves, "cact"), paste0("r", cog_test_waves, "imrc"), 
    paste0("r", cog_test_waves, "dlrc"), paste0("r", cog_test_waves, "ser7"), 
    paste0("r", cog_test_waves, "pres"), paste0("r", cog_test_waves, "pstmem"), 
    paste0("r", cog_test_waves, "cogtot"))

RAND <- read_dta(paste0(path_to_box, "data/HRS/RAND_longitudinal/STATA/", 
                        "randhrs1992_2018v1.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character) %>% rename("HHIDPN" = "hhidpn")

#Remove labeled data format
val_labels(RAND) <- NULL

#---- join data ----
ADAMS <- left_join(ADAMS_tracker, ADAMS_neuropsych, by = "HHIDPN") %>% 
  left_join(., ADAMS_proxy, by = "HHIDPN") %>% 
  left_join(., ADAMS_demdx, by = "HHIDPN") %>% 
  left_join(., RAND, by = "HHIDPN")

#---- clean data ----
#---- **sex/gender ----
#table(ADAMS$GENDER, useNA = "ifany")
ADAMS %<>% 
  mutate(GENDER_label = as.factor(ifelse(GENDER == 1, "Male", "Female"))) %>% 
  mutate("Female" = ifelse(GENDER_label == "Female", 1, 0))

# #Sanity check
# table(ADAMS$GENDER, ADAMS$GENDER_label)
# table(ADAMS$Female, useNA = "ifany")

#---- **race/ethnicity ----
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

#---- **interview year ----
#table(ADAMS$AYEAR, useNA = "ifany")

#---- **HRS cognition ----
# table(ADAMS$SELFCOG, useNA = "ifany")

#---- **marital status ----
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

#---- **age ----
#table(ADAMS$AAGE, useNA = "ifany")

#---- **education ----
#table(ADAMS$EDYRS, useNA = "ifany")

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

#---- **MMSE ----
#table(ADAMS$ANMSETOT, useNA = "ifany")
ADAMS %<>% 
  mutate_at(.vars = "ANMSETOT", function(x) ifelse(x > 30, NA, x)) %>% 
  #normalized MMSE
  mutate("ANMSETOT_norm" = normMMSE(ANMSETOT))

# #Sanity check
# hist(ADAMS$ANMSETOT)
# hist(ADAMS$ANMSETOT_norm)
# table(ADAMS$ANMSETOT_norm, useNA = "ifany")

#---- **BWC 20 ----
# table(ADAMS$ANBWC201, useNA = "ifany")
# table(ADAMS$ANBWC202, useNA = "ifany")
# table(ADAMS$ANBWC861, useNA = "ifany")
# table(ADAMS$ANBWC862, useNA = "ifany")
ADAMS %<>% 
  mutate_at(.vars = c("ANBWC201", "ANBWC202"), 
            #Missing/refused  
            function(x) ifelse(x > 6, NA, x)) %>% 
  mutate_at(.vars = c("ANBWC201", "ANBWC202"), 
            #restart
            function(x) ifelse(x == 6, 0, x)) 

#Take the higher score
ADAMS %<>% mutate("ANBWC20" = pmax(ANBWC201, ANBWC202, na.rm = TRUE))

# #Sanity check
# View(ADAMS[, c("ANBWC201", "ANBWC202", "ANBWC20")])
# table(ADAMS$ANBWC20, useNA = "ifany")

#---- **serial 7s ----
#table(ADAMS$ANSER7T, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANSER7T"), 
                     #Missing/refused  
                     function(x) ifelse(x > 5, NA, x))

# #Sanity check
# table(ADAMS$ANSER7T, useNA = "ifany")

#---- **object naming: cactus, scissors ----
# table(ADAMS$ANCACTUS, useNA = "ifany")
# table(ADAMS$ANSCISOR, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANCACTUS", "ANSCISOR"), 
                     #Missing/refused  
                     function(x) ifelse(x > 1, NA, x)) 

# #Sanity check
# table(ADAMS$ANCACTUS, useNA = "ifany")
# table(ADAMS$ANSCISOR, useNA = "ifany")

#---- **President naming ----
# table(ADAMS$ANPRES, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANPRES"), 
                     #Missing/refused  
                     function(x) ifelse(x > 1, NA, x)) 

# #Sanity check
# table(ADAMS$ANPRES, useNA = "ifany")

#---- **animal naming ----
#table(ADAMS$ANAFTOT, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANAFTOT"), 
                     #Missing/refused  
                     function(x) ifelse(x > 33, NA, x)) 

# #Sanity check
# table(ADAMS$ANAFTOT, useNA = "ifany")

#---- **10-word recall (immediate and delayed) ----
# table(ADAMS$ANIMMCR1, useNA = "ifany")
# table(ADAMS$ANDELCOR, useNA = "ifany")
ADAMS %<>% 
  mutate_at(.vars = c("ANIMMCR1", "ANIMMCR2", "ANIMMCR3", "ANDELCOR"), 
            #Missing/refused  
            function(x) ifelse(x > 10, NA, x)) %>% 
  #Best of 3 immediate recall trials
  mutate("ANIMMCR" = pmax(ANIMMCR1, ANIMMCR2, ANIMMCR3, na.rm = TRUE))

# #Sanity check
# View(ADAMS[, c("ANIMMCR1", "ANIMMCR2", "ANIMMCR3", "ANIMMCR")])
# table(ADAMS$ANIMMCR, useNA = "ifany")
# table(ADAMS$ANDELCOR, useNA = "ifany")

#---- **word list recognition (yes/no) ----
# table(ADAMS$ANRECNO, useNA = "ifany")
# table(ADAMS$ANRECYES, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANRECNO", "ANRECYES"), 
                     #Missing/refused  
                     function(x) ifelse(x > 10, NA, x))

# #Sanity check
# table(ADAMS$ANRECNO, useNA = "ifany")
# table(ADAMS$ANRECYES, useNA = "ifany")

#---- **story recall (immediate and delayed) ----
# table(ADAMS$ANWM1TOT, useNA = "ifany")
# table(ADAMS$ANWM2TOT, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANWM1TOT", "ANWM2TOT"), 
                     #Missing/refused  
                     function(x) ifelse(x > 37, NA, x))

# #Sanity check
# table(ADAMS$ANWM1TOT, useNA = "ifany")
# table(ADAMS$ANWM2TOT, useNA = "ifany")

#---- **constructional praxis (immediate and delayed) ----
# table(ADAMS$ANCPTOT, useNA = "ifany")
# table(ADAMS$ANRCPTOT, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANCPTOT", "ANRCPTOT"), 
                     #Missing/refused  
                     function(x) ifelse(x > 11, NA, x))

# #Sanity check
# table(ADAMS$ANCPTOT, useNA = "ifany")
# table(ADAMS$ANRCPTOT, useNA = "ifany")

#---- **trails A ----
# table(ADAMS$ANTMASEC, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANTMASEC"),
                     #Missing/refused
                     function(x) ifelse(x > 900, NA, x))

# #Sanity check
# table(ADAMS$ANTMASEC, useNA = "ifany")

#---- **subjective cognitive change ----
# table(ADAMS$ANSMEM2, useNA = "ifany")
ADAMS %<>% mutate_at(.vars = c("ANSMEM2"),
                     #Missing/refused
                     function(x) ifelse(x > 5, NA, x)) %>% 
  mutate("ANSMEM2_label" = case_when(ANSMEM2 == 1 ~ "Much Better", 
                                     ANSMEM2 == 2 ~ "Better", 
                                     ANSMEM2 == 3 ~ "Same", 
                                     ANSMEM2 == 4 ~ "Worse", 
                                     ANSMEM2 == 5 ~ "Much Worse")) %>% 
  mutate("ANSMEM2_collapsed_label" = 
           case_when(ANSMEM2_label %in% c("Much Better", "Better") ~ "Better", 
                     ANSMEM2_label == "Same" ~ "Same", 
                     ANSMEM2_label %in% c("Worse", "Much Worse") ~ "Worse")) %>% 
  mutate("ANSMEM2_Better" = ifelse(ANSMEM2_collapsed_label == "Better", 1, 0), 
         "ANSMEM2_Same" = ifelse(ANSMEM2_collapsed_label == "Same", 1, 0), 
         "ANSMEM2_Worse" = ifelse(ANSMEM2_collapsed_label == "Worse", 1, 0))

# #Sanity check
# table(ADAMS$ANSMEM2, useNA = "ifany")
# table(ADAMS$ANSMEM2, ADAMS$ANSMEM2_label, useNA = "ifany")
# table(ADAMS$ANSMEM2_label, ADAMS$ANSMEM2_collapsed_label, useNA = "ifany")
# table(ADAMS$ANSMEM2_Better, ADAMS$ANSMEM2_collapsed_label, useNA = "ifany")
# table(ADAMS$ANSMEM2_Same, ADAMS$ANSMEM2_collapsed_label, useNA = "ifany")
# table(ADAMS$ANSMEM2_Worse, ADAMS$ANSMEM2_collapsed_label, useNA = "ifany")

#---- **proxy cognition ----
# table(ADAMS$AGQ14, useNA = "ifany")
# table(ADAMS$AGQ29, useNA = "ifany")
ADAMS %<>% mutate("avg_proxy_cog" = ADAMS %>% 
                    dplyr::select(paste0("AGQ", seq(14, 29))) %>% 
                    rowMeans(., na.rm = TRUE)) %>% 
  mutate("avg_proxy_cog" = ifelse(is.nan(avg_proxy_cog), NA, avg_proxy_cog)) %>%
  #floor to get whole numbers
  mutate_at(.vars = c("avg_proxy_cog"), floor) %>% 
  mutate("avg_proxy_cog_label" = 
           case_when(avg_proxy_cog == 1 ~ "Much Better", 
                     avg_proxy_cog == 2 ~ "Better", 
                     avg_proxy_cog == 3 ~ "Same", 
                     avg_proxy_cog == 4 ~ "Worse", 
                     avg_proxy_cog == 5 ~ "Much Worse")) %>% 
  mutate("avg_proxy_cog_collapsed_label" = 
           case_when(avg_proxy_cog_label %in% 
                       c("Much Better", "Better") ~ "Better", 
                     avg_proxy_cog_label == "Same" ~ "Same", 
                     avg_proxy_cog_label %in% 
                       c("Worse", "Much Worse") ~ "Worse")) %>% 
  mutate("avg_proxy_cog_Better" = 
           ifelse(avg_proxy_cog_collapsed_label == "Better", 1, 0), 
         "avg_proxy_cog_Same" = 
           ifelse(avg_proxy_cog_collapsed_label == "Same", 1, 0), 
         "avg_proxy_cog_Worse" = 
           ifelse(avg_proxy_cog_collapsed_label == "Worse", 1, 0))

# #Sanity check
# View(ADAMS %>% dplyr::select(paste0("AGQ", seq(14, 29)), "avg_proxy_cog"))
# table(ADAMS$avg_proxy_cog, useNA = "ifany")
# table(ADAMS$avg_proxy_cog, ADAMS$avg_proxy_cog_label, useNA = "ifany")
# table(ADAMS$avg_proxy_cog_label, ADAMS$avg_proxy_cog_collapsed_label,
#       useNA = "ifany")
# table(ADAMS$avg_proxy_cog_Better, ADAMS$avg_proxy_cog_collapsed_label,
#       useNA = "ifany")
# table(ADAMS$avg_proxy_cog_Same, ADAMS$avg_proxy_cog_collapsed_label,
#       useNA = "ifany")
# table(ADAMS$avg_proxy_cog_Worse, ADAMS$avg_proxy_cog_collapsed_label,
#       useNA = "ifany")

# #Distribution check: it is not the case that all those missing proxies are high
# # functioning
# ADAMS %>% filter(is.na(avg_proxy_cog)) %>% dplyr::select(ANMSETOT) %>% 
#   table(., useNA = "ifany")

#---- **dementia dx ----
# table(ADAMS$ADFDX1, useNA = "ifany")
ADAMS %<>% 
  #collapse categories
  mutate("Adem_dx_label_collapsed" = 
           case_when(ADFDX1 %in% c(1, 2) ~ "Probable/Possible AD", 
                     ADFDX1 %in% c(3, 4) ~ 
                       "Probable/Possible Vascular Dementia", 
                     ADFDX1 %in% c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 21, 
                                   28, 29, 30, 33) ~ "Other",
                     ADFDX1 %in% c(18, 32) ~ "Probable Dementia",
                     ADFDX1 %in% c(10, 13, 15, 16, 17, 19) ~ "Dementia", 
                     ADFDX1 %in% c(20, 22) ~ "MCI", 
                     ADFDX1 == 31 ~ "Normal")) %>% 
  #further collapsing categories
  mutate("Adem_cat" = 
           case_when(Adem_dx_label_collapsed %in% 
                       c("Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia", 
                         "Probable Dementia", "Dementia") ~ "Dementia",
                     Adem_dx_label_collapsed == "Other" ~ "Other", 
                     Adem_dx_label_collapsed == "MCI" ~ "MCI", 
                     Adem_dx_label_collapsed == "Normal" ~ "Unimpaired")) %>% 
  #dummy vars
  mutate("Dementia" = ifelse(Adem_cat == "Dementia", 1, 0), 
         "Other" = ifelse(Adem_cat == "Other", 1, 0), 
         "MCI" = ifelse(Adem_cat == "MCI", 1, 0), 
         "Unimpaired" = ifelse(Adem_cat == "Unimpaired", 1, 0))

# #Sanity check
# table(ADAMS$Adem_dx_label_collapsed, ADAMS$Adem_cat, useNA = "ifany")
# table(ADAMS$Adem_cat, ADAMS$Dementia, useNA = "ifany")
# table(ADAMS$Adem_cat, ADAMS$Other, useNA = "ifany")
# table(ADAMS$Adem_cat, ADAMS$MCI, useNA = "ifany")
# table(ADAMS$Adem_cat, ADAMS$Unimpaired, useNA = "ifany")

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
neuro_cog_measures <- c("SELFCOG", "ANMSETOT_norm", "ANBWC20", "ANSER7T", 
                        "ANSCISOR", "ANCACTUS", "ANPRES", "ANAFTOT", "ANIMMCR", 
                        "ANDELCOR", "ANRECYES", "ANRECNO", "ANWM1TOT", 
                        "ANWM2TOT", "ANCPTOT", "ANRCPTOT", "ANTMASEC")

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

bwc_vars <- c(paste0("r", cog_test_waves, "bwc20"))
ADAMS %<>% 
  mutate_at(.vars = all_of(bwc_vars), 
            #Correct on 2nd try = correct
            function(x) ifelse(x == 2, 1, x))

# #Sanity check
# table(ADAMS$r5bwc20, useNA = "ifany")

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
            paste0("ANBWC", c("201", "202")), "GENDER_label",
            paste0("ANIMMCR", seq(1, 3)), "ETHNIC_label", "num_cog_measures", 
            "proxy_type_label", "proxy_type_collapsed_label", 
            paste0("r", cog_test_waves, "pstmem"), 
            paste0("r", rand_waves, "bmi"))

ADAMS %>% dplyr::select(-all_of(remove)) %>% 
  write_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))
