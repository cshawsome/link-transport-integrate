#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy", "labelled", 
       "gtsummary", "writexl", "plyr")

options(scipen = 999)

#---- source scripts ----
source(paste0(here::here("functions", "read_da_dct.R")))

#---- Table 4.1: HRS, HCAP, ADAMS, Superpop characteristics ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- ****variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- ****ADAMS ----
ADAMS <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv")) 

ADAMS_labels <- 
  variable_labels[which(variable_labels$ADAMS %in% colnames(ADAMS)), ]

ADAMS %<>% rename_at(vars(ADAMS_labels$ADAMS), ~ ADAMS_labels$data_label)

#---- ****HCAP ----
HCAP <- read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_clean.csv"))

HCAP_labels <- 
  variable_labels[which(variable_labels$HCAP %in% colnames(HCAP)), ]

HCAP %<>% rename_at(vars(HCAP_labels$HCAP), ~ HCAP_labels$data_label)

#---- ****HRS ----
HRS <- read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv"))

HRS_labels <- 
  variable_labels[which(variable_labels$HRS %in% colnames(HRS)), ]

HRS %<>% rename_at(vars(HRS_labels$HRS), ~ HRS_labels$data_label)

#---- **superpopulation ----
superpop <- 
  read_csv(paste0(path_to_box, "data/superpopulations/superpop_1000000.csv"))

#---- **variable selection ----
selected_vars <- 
  read_csv(paste0(path_to_box, 
                  "data/variable_selection/model_coefficients.csv")) %>% 
  filter(data_label != "Intercept") %>% dplyr::select("data_label") %>% 
  unlist() %>% unname() %>% str_remove(., "_Z")

ADAMS %<>% dplyr::select(c(all_of(selected_vars), 
                           "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  mutate("dataset" = "ADAMS")

HCAP %<>% dplyr::select(all_of(selected_vars)) %>% mutate("dataset" = "HCAP") 

HRS %<>% 
  dplyr::select(selected_vars[all_of(selected_vars %in% colnames(HRS))]) %>% 
  mutate("dataset" = "HRS") 

superpop %<>% dplyr::select(c(all_of(selected_vars), 
                              "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  mutate("dataset" = "Superpopulation")

#---- **join data ----
table_data <- rbind.fill(ADAMS, HCAP) %>% rbind.fill(., HRS) %>% 
  rbind.fill(superpop)

#---- **create label vars ----
table_data %<>% 
  mutate("race_label" = 
           case_when(black == 1 ~ "Black", 
                     hispanic == 1 ~ "Hispanic", 
                     TRUE ~ "White"), 
         "employment_label" = 
           case_when(not_working == 1 ~ "Not working", 
                     is.na(not_working) ~ "Missing",
                     retired == 1 ~ "Retired", 
                     TRUE ~ "Working"), 
         "subjective_cognition_label" = 
           case_when(subj_cog_better == 1 ~ "Better than 2 years ago",
                     is.na(subj_cog_better) ~ "Missing",
                     subj_cog_worse == 1 ~ "Worse than 2 years ago", 
                     TRUE ~ "Same as 2 years ago"), 
         "alcohol_consumption_label" = 
           case_when(moderate_drinking == 1 ~ "Moderate drinking", 
                     is.na(moderate_drinking) ~ "Missing",
                     heavy_drinking == 1 ~ "Heavy drinking", 
                     TRUE ~ "No drinking"), 
         "impairment_label" = 
           case_when(Unimpaired == 1 ~ "Unimpaired", 
                     MCI == 1 ~ "MCI", 
                     Dementia == 1 ~ "Dementia", 
                     Other == 1 ~ "Other"), 
         "married_partnered_label" = 
           case_when(married_partnered == 1 ~ "Married/Partnered", 
                     is.na(married_partnered) ~ "Missing",
                     TRUE ~ "Not Married/Partnered"), 
         "female_label" = 
           case_when(female == 1 ~ "Female", 
                     is.na(female) ~ "Missing", 
                     TRUE ~ "Male"), 
         "bwc20_label" = 
           case_when(bwc20 == 1 ~ "Backwards count (20): correct",
                     is.na(bwc20) ~ "Missing",
                     TRUE ~ "Backwards count (20): incorrect"), 
         "cactus_label" = 
           case_when(cactus == 1 ~ "Item naming (cactus): correct", 
                     is.na(cactus) ~ "Missing",
                     cactus == 0 ~ "Item naming (cactus): incorrect"), 
         "scissor_label" = 
           case_when(scissor == 1 ~ "Item naming (scissor): correct", 
                     is.na(scissor) ~ "Missing",
                     scissor == 0 ~ "Item naming (scissor): incorrect"), 
         "pres_label" = 
           case_when(pres == 1 ~ "President naming: correct", 
                     is.na(pres) ~ "Missing",
                     pres == 0 ~ "President naming: incorrect"), 
         "stroke_label" = 
           case_when(stroke == 1 ~ "History of stroke",
                     is.na(stroke) ~ "Missing",
                     stroke == 0 ~ "No History of stroke"), 
         "diabe_label" = 
           case_when(diabe == 1 ~ "History of diabetes",
                     is.na(diabe) ~ "Missing",
                     diabe == 0 ~ "No History of diabetes"), 
         "hearte_label" = 
           case_when(hearte == 1 ~ "History of heart disease",
                     is.na(hearte) ~ "Missing",
                     hearte == 0 ~ "No History of heart disease"), 
         "hibpe_label" = 
           case_when(hibpe == 1 ~ "History of hypertension",
                     is.na(hibpe) ~ "Missing",
                     hibpe == 0 ~ "No History of hypertension"), 
         "smoken_label" = 
           case_when(smoken == 1 ~ "Current smoker",
                     is.na(smoken) ~ "Missing",
                     smoken == 0 ~ "Not Current smoker"))

#---- **set factor levels ----
table_data %<>% 
  mutate_at("race_label", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate_at("employment_label", function(x) 
    factor(x, levels = c("Working", "Not working", "Retired", "Missing"))) %>%
  mutate_at("subjective_cognition_label", function(x) 
    factor(x, levels = c("Same as 2 years ago", "Better than 2 years ago", 
                         "Worse than 2 years ago", "Missing"))) %>%
  mutate_at("alcohol_consumption_label", function(x) 
    factor(x, levels = c("No drinking", "Moderate drinking", 
                         "Heavy drinking", "Missing"))) %>%
  mutate_at("impairment_label", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("married_partnered_label", function(x)
    factor(x, levels = c("Married/Partnered", "Not Married/Partnered", 
                         "Missing"))) %>% 
  mutate_at("female_label", function(x)
    factor(x, levels = c("Female", "Male", "Missing"))) %>% 
  mutate_at("bwc20_label", function(x)
    factor(x, levels = c("Backwards count (20): correct", 
                         "Backwards count (20): incorrect", "Missing"))) %>% 
  mutate_at("cactus_label", function(x)
    factor(x, levels = c("Item naming (cactus): correct", 
                         "Item naming (cactus): incorrect", "Missing"))) %>% 
  mutate_at("scissor_label", function(x)
    factor(x, levels = c("Item naming (scissor): correct", 
                         "Item naming (scissor): incorrect", "Missing"))) %>% 
  mutate_at("pres_label", function(x)
    factor(x, levels = c("President naming: correct",
                         "President naming: incorrect", "Missing"))) %>% 
  mutate_at("stroke_label", function(x)
    factor(x, levels = c("History of stroke", "No History of stroke", 
                         "Missing"))) %>% 
  mutate_at("diabe_label", function(x)
    factor(x, levels = c("History of diabetes", "No History of diabetes", 
                         "Missing"))) %>% 
  mutate_at("hearte_label", function(x)
    factor(x, levels = c("History of heart disease", 
                         "No History of heart disease", "Missing"))) %>% 
  mutate_at("hibpe_label", function(x)
    factor(x, levels = c("History of hypertension", 
                         "No History of hypertension", "Missing"))) %>% 
  mutate_at("smoken_label", function(x)
    factor(x, levels = c("Current smoker", "Not Current smoker", "Missing"))) 

#---- **label variables ----
table_data %<>%
  labelled::set_variable_labels(age = "Age",
                                female_label = "Female",
                                edyrs = "Years of Education",
                                married_partnered_label = "Married/Partnered",
                                mmse_norm = "Total MMSE (normalized)",
                                immrc = "Immediate word recall",
                                delrc = "Delayed word recall",
                                ser7 = "Serial 7s",
                                animal_naming = "Animal naming",
                                wrc_yes = "Word recall (yes)",
                                wrc_no = "Word recall (no)",
                                imm_story = "Immediate story recall",
                                del_story = "Delayed story recall",
                                bwc20_label = "Backwards count (20): correct",
                                imm_cp = "Immediate constructional praxis",
                                del_cp = "Delayed constructional praxis",
                                trailsA = "Trails A",
                                hrs_cog = "HRS total cognition",
                                cactus_label = "Item naming (cactus): correct",
                                scissor_label = "Item naming (scissor): correct",
                                pres_label = "President naming: correct",
                                adl = "ADLs",
                                iadl = "IADLs",
                                bmi = "BMI",
                                stroke_label = "History of stroke",
                                diabe_label = "History of diabetes",
                                hearte_label = "History of heart disease",
                                hibpe_label = "History of hypertension",
                                smoken_label = "Current smoker",
                                race_label = "Race/Ethnicity",
                                employment_label = "Employment status",
                                subjective_cognition_label =
                                  "Subjective cognitive status",
                                alcohol_consumption_label = "Alcohol consumption",
                                impairment_label = "Impairment group")

#---- **label categories ----
table_data %<>% 
  labelled::set_value_labels(dataset = 
                               c("ADAMS\nBaseline (2002)\n" = "ADAMS", 
                                 "HCAP 70+\nBaseline (2016)\n" = "HCAP", 
                                 "HRS 70+\n(2016)\n" = "HRS", 
                                 "Superpopulation" = "Superpopulation\n"))

#---- **make table ----
table <- table_data %>%
  #set labelled variables as factors 
  purrr::modify_if(labelled::is.labelled, labelled::to_factor) %>%
  gtsummary::tbl_summary(missing = "no", 
                         missing_text = "Missing",
                         by = dataset,
                         statistic = list(all_continuous() ~ 
                                            c("{mean} ({sd})", 
                                              "{N_miss} ({p_miss}%)")),
                         #If we don't specify these variables types, all levels 
                         #  will be summarized
                         type = list(all_continuous() ~ "continuous2",
                                     iadl ~ "continuous2",
                                     adl ~ "continuous2",
                                     ser7 ~ "continuous2"),
                         
                         #Specifying the number of decimal places for 
                         #  categorical vars
                         digits = list(all_continuous() ~ 1, 
                                       female_label ~ c(0, 1),
                                       married_partnered_label ~ c(0, 1), 
                                       bwc20_label ~ c(0, 1), 
                                       cactus_label ~ c(0, 1), 
                                       scissor_label ~ c(0, 1), 
                                       pres_label ~ c(0, 1), 
                                       stroke_label ~ c(0, 1), 
                                       diabe_label ~ c(0, 1), 
                                       hearte_label ~ c(0, 1), 
                                       hibpe_label ~ c(0, 1), 
                                       smoken_label ~ c(0, 1), 
                                       race_label ~ c(0, 1), 
                                       employment_label ~ c(0, 1), 
                                       subjective_cognition_label ~ c(0, 1), 
                                       alcohol_consumption_label ~ c(0, 1), 
                                       impairment_label ~ c(0, 1)), 
                         #order variables
                         include = 
                           c(age, female_label, race_label, edyrs, 
                             employment_label, married_partnered_label, bmi, 
                             stroke_label, diabe_label, hearte_label, 
                             hibpe_label, smoken_label, 
                             alcohol_consumption_label, adl, iadl, immrc, delrc, 
                             ser7, cactus_label, scissor_label, pres_label, 
                             bwc20_label, hrs_cog, subjective_cognition_label, 
                             mmse_norm, animal_naming, wrc_yes, wrc_no, 
                             imm_story, del_story, imm_cp, del_cp, trailsA, 
                             impairment_label)) %>%
  
  #Renaming the header
  gtsummary::modify_header(label = "Variable") %>%
  
  #Moving labels from the bottom to next to each of the variables
  gtsummary::add_stat_label(location = "row") %>%
  
  #Modify missing label
  modify_table_body(dplyr::mutate,
                    label = ifelse(label == "N missing (% missing)",
                                   "Missing, n (%)", label)) %>%
  
  #setting as a tibble so that it can be output to excel
  gtsummary::as_tibble()

#---- **save output ----
#Going to need to clean up this output in excel
writexl::write_xlsx(gtsummary::as_tibble(table), path = paste0(
  path_to_box, "tables/chapter_4/table4.1_sample_characteristics.xlsx"))

#---- Appendix Table E.1: Imputed HRS, HCAP, ADAMS, Superpop characteristics ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- ****variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- ****ADAMS ----
ADAMS_imputed_list <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/chunk_1/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

ADAMS <- Reduce(`+`, ADAMS_imputed_list) / length(ADAMS_imputed_list)

ADAMS_labels <- 
  variable_labels[which(variable_labels$ADAMS %in% colnames(ADAMS)), ]

ADAMS %<>% rename_at(vars(ADAMS_labels$ADAMS), ~ ADAMS_labels$data_label)

#---- ****HCAP ----
HCAP <- read_csv(paste0(path_to_box, 
                        "data/HCAP/cleaned/HCAP_analytic_for_sim.csv")) 

HCAP_labels <- 
  variable_labels[which(variable_labels$HCAP %in% colnames(HCAP)), ]

HCAP %<>% rename_at(vars(HCAP_labels$HCAP), ~ HCAP_labels$data_label)

#---- ****HRS ----
HRS <- read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv")) 

HRS_labels <- 
  variable_labels[which(variable_labels$HRS %in% colnames(HRS)), ]

HRS %<>% rename_at(vars(HRS_labels$HRS), ~ HRS_labels$data_label)

#---- **superpopulation ----
superpop <- 
  read_csv(paste0(path_to_box, "data/superpopulations/superpop_1000000.csv")) 

#---- **variable selection ----
selected_vars <- 
  read_csv(paste0(path_to_box, 
                  "data/variable_selection/model_coefficients.csv")) %>% 
  filter(data_label != "Intercept") %>% dplyr::select("data_label") %>% 
  unlist() %>% unname() %>% str_remove(., "_Z")

ADAMS %<>% 
  dplyr::select(c(all_of(selected_vars), 
                  "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  mutate("dataset" = "ADAMS")

HCAP %<>% dplyr::select(all_of(selected_vars)) %>% mutate("dataset" = "HCAP") 

HRS %<>% 
  dplyr::select(selected_vars[all_of(selected_vars %in% colnames(HRS))]) %>% 
  mutate("dataset" = "HRS") 

superpop %<>% dplyr::select(c(all_of(selected_vars), 
                              "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  mutate("dataset" = "Superpopulation")

#---- **format ADAMS indicator vars ----
indicator_vars <- 
  c("retired", "not_working", "married_partnered", "stroke", "diabe", "hearte", 
    "hibpe", "moderate_drinking", "heavy_drinking", "cactus", "scissor", "pres", 
    "bwc20", "subj_cog_better", "subj_cog_worse")

ADAMS %<>% 
  mutate_at(all_of(indicator_vars), function(x) ifelse(x > 0.5, 1, 0))


#---- **join data ----
table_data <- rbind.fill(ADAMS, HCAP) %>% rbind.fill(., HRS) %>% 
  rbind.fill(superpop)

#---- **create label vars ----
table_data %<>% 
  mutate("race_label" = 
           case_when(black == 1 ~ "Black", 
                     hispanic == 1 ~ "Hispanic", 
                     TRUE ~ "White"), 
         "employment_label" = 
           case_when(not_working == 1 ~ "Not working", 
                     is.na(not_working) ~ "Missing",
                     retired == 1 ~ "Retired", 
                     TRUE ~ "Working"), 
         "subjective_cognition_label" = 
           case_when(subj_cog_better == 1 ~ "Better than 2 years ago",
                     is.na(subj_cog_better) ~ "Missing",
                     subj_cog_worse == 1 ~ "Worse than 2 years ago", 
                     TRUE ~ "Same as 2 years ago"), 
         "alcohol_consumption_label" = 
           case_when(moderate_drinking == 1 ~ "Moderate drinking", 
                     is.na(moderate_drinking) ~ "Missing",
                     heavy_drinking == 1 ~ "Heavy drinking", 
                     TRUE ~ "No drinking"), 
         "impairment_label" = 
           case_when(Unimpaired == 1 ~ "Unimpaired", 
                     MCI == 1 ~ "MCI", 
                     Dementia == 1 ~ "Dementia", 
                     Other == 1 ~ "Other"), 
         "married_partnered_label" = 
           case_when(married_partnered == 1 ~ "Married/Partnered", 
                     is.na(married_partnered) ~ "Missing",
                     TRUE ~ "Not Married/Partnered"), 
         "female_label" = 
           case_when(female == 1 ~ "Female", 
                     is.na(female) ~ "Missing", 
                     TRUE ~ "Male"), 
         "bwc20_label" = 
           case_when(bwc20 == 1 ~ "Backwards count (20): correct",
                     is.na(bwc20) ~ "Missing",
                     TRUE ~ "Backwards count (20): incorrect"), 
         "cactus_label" = 
           case_when(cactus == 1 ~ "Item naming (cactus): correct", 
                     is.na(cactus) ~ "Missing",
                     cactus == 0 ~ "Item naming (cactus): incorrect"), 
         "scissor_label" = 
           case_when(scissor == 1 ~ "Item naming (scissor): correct", 
                     is.na(scissor) ~ "Missing",
                     scissor == 0 ~ "Item naming (scissor): incorrect"), 
         "pres_label" = 
           case_when(pres == 1 ~ "President naming: correct", 
                     is.na(pres) ~ "Missing",
                     pres == 0 ~ "President naming: incorrect"), 
         "stroke_label" = 
           case_when(stroke == 1 ~ "History of stroke",
                     is.na(stroke) ~ "Missing",
                     stroke == 0 ~ "No History of stroke"), 
         "diabe_label" = 
           case_when(diabe == 1 ~ "History of diabetes",
                     is.na(diabe) ~ "Missing",
                     diabe == 0 ~ "No History of diabetes"), 
         "hearte_label" = 
           case_when(hearte == 1 ~ "History of heart disease",
                     is.na(hearte) ~ "Missing",
                     hearte == 0 ~ "No History of heart disease"), 
         "hibpe_label" = 
           case_when(hibpe == 1 ~ "History of hypertension",
                     is.na(hibpe) ~ "Missing",
                     hibpe == 0 ~ "No History of hypertension"), 
         "smoken_label" = 
           case_when(smoken == 1 ~ "Current smoker",
                     is.na(smoken) ~ "Missing",
                     smoken == 0 ~ "Not Current smoker"))

#---- **set factor levels ----
table_data %<>% 
  mutate_at("race_label", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate_at("employment_label", function(x) 
    factor(x, levels = c("Working", "Not working", "Retired"))) %>%
  mutate_at("subjective_cognition_label", function(x) 
    factor(x, levels = c("Same as 2 years ago", "Better than 2 years ago", 
                         "Worse than 2 years ago"))) %>%
  mutate_at("alcohol_consumption_label", function(x) 
    factor(x, levels = c("No drinking", "Moderate drinking", 
                         "Heavy drinking"))) %>%
  mutate_at("impairment_label", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("married_partnered_label", function(x)
    factor(x, levels = c("Married/Partnered", "Not Married/Partnered"))) %>% 
  mutate_at("female_label", function(x)
    factor(x, levels = c("Female", "Male"))) %>% 
  mutate_at("bwc20_label", function(x)
    factor(x, levels = c("Backwards count (20): correct", 
                         "Backwards count (20): incorrect"))) %>% 
  mutate_at("cactus_label", function(x)
    factor(x, levels = c("Item naming (cactus): correct", 
                         "Item naming (cactus): incorrect"))) %>% 
  mutate_at("scissor_label", function(x)
    factor(x, levels = c("Item naming (scissor): correct", 
                         "Item naming (scissor): incorrect"))) %>% 
  mutate_at("pres_label", function(x)
    factor(x, levels = c("President naming: correct",
                         "President naming: incorrect"))) %>% 
  mutate_at("stroke_label", function(x)
    factor(x, levels = c("History of stroke", "No History of stroke"))) %>% 
  mutate_at("diabe_label", function(x)
    factor(x, levels = c("History of diabetes", "No History of diabetes"))) %>% 
  mutate_at("hearte_label", function(x)
    factor(x, levels = c("History of heart disease", 
                         "No History of heart disease"))) %>% 
  mutate_at("hibpe_label", function(x)
    factor(x, levels = c("History of hypertension", 
                         "No History of hypertension"))) %>% 
  mutate_at("smoken_label", function(x)
    factor(x, levels = c("Current smoker", "Not Current smoker"))) 

#---- **label variables ----
table_data %<>%
  labelled::set_variable_labels(age = "Age",
                                female_label = "Female",
                                edyrs = "Years of Education",
                                married_partnered_label = "Married/Partnered",
                                mmse_norm = "Total MMSE (normalized)",
                                immrc = "Immediate word recall",
                                delrc = "Delayed word recall",
                                ser7 = "Serial 7s",
                                animal_naming = "Animal naming",
                                wrc_yes = "Word recall (yes)",
                                wrc_no = "Word recall (no)",
                                imm_story = "Immediate story recall",
                                del_story = "Delayed story recall",
                                bwc20_label = "Backwards count (20): correct",
                                imm_cp = "Immediate constructional praxis",
                                del_cp = "Delayed constructional praxis",
                                trailsA = "Trails A",
                                hrs_cog = "HRS total cognition",
                                cactus_label = "Item naming (cactus): correct",
                                scissor_label = "Item naming (scissor): correct",
                                pres_label = "President naming: correct",
                                adl = "ADLs",
                                iadl = "IADLs",
                                bmi = "BMI",
                                stroke_label = "History of stroke",
                                diabe_label = "History of diabetes",
                                hearte_label = "History of heart disease",
                                hibpe_label = "History of hypertension",
                                smoken_label = "Current smoker",
                                race_label = "Race/Ethnicity",
                                employment_label = "Employment status",
                                subjective_cognition_label =
                                  "Subjective cognitive status",
                                alcohol_consumption_label = "Alcohol consumption",
                                impairment_label = "Impairment group")

#---- **label categories ----
table_data %<>% 
  labelled::set_value_labels(dataset = 
                               c("ADAMS\nBaseline (2002)\n" = "ADAMS", 
                                 "HCAP 70+\nBaseline (2016)\n" = "HCAP", 
                                 "HRS 70+\n(2016)\n" = "HRS", 
                                 "Superpopulation" = "Superpopulation\n"))

#---- **make table ----
table <- table_data %>%
  #set labelled variables as factors 
  purrr::modify_if(labelled::is.labelled, labelled::to_factor) %>%
  gtsummary::tbl_summary(missing = "no", 
                         by = dataset,
                         statistic = list(all_continuous() ~ 
                                            c("{mean} ({sd})")),
                         #If we don't specify these variables types, all levels 
                         #  will be summarized
                         type = list(all_continuous() ~ "continuous",
                                     iadl ~ "continuous",
                                     adl ~ "continuous",
                                     ser7 ~ "continuous", 
                                     female_label ~ "dichotomous",
                                     married_partnered_label ~ "dichotomous",
                                     bwc20_label ~ "dichotomous",
                                     cactus_label ~ "dichotomous",
                                     scissor_label ~ "dichotomous",
                                     pres_label ~ "dichotomous",
                                     stroke_label ~ "dichotomous",
                                     diabe_label ~ "dichotomous",
                                     hearte_label ~ "dichotomous",
                                     hibpe_label ~ "dichotomous",
                                     smoken_label ~ "dichotomous"),
                         
                         #Specifying the level of dichotomous the variable that 
                         #  should be displayed
                         value = list(female_label = 'Female',
                                      married_partnered_label = 
                                        "Married/Partnered",
                                      bwc20_label = 
                                        "Backwards count (20): correct",
                                      cactus_label = 
                                        "Item naming (cactus): correct",
                                      scissor_label = 
                                        "Item naming (scissor): correct",
                                      pres_label = "President naming: correct",
                                      stroke_label = "History of stroke",
                                      diabe_label = "History of diabetes",
                                      hearte_label = "History of heart disease",
                                      hibpe_label = "History of hypertension",
                                      smoken_label = "Current smoker"),
                         
                         #Specifying the number of decimal places for 
                         #  categorical vars
                         digits = list(all_continuous() ~ 1, 
                                       female_label ~ c(0, 1),
                                       married_partnered_label ~ c(0, 1), 
                                       bwc20_label ~ c(0, 1), 
                                       cactus_label ~ c(0, 1), 
                                       scissor_label ~ c(0, 1), 
                                       pres_label ~ c(0, 1), 
                                       stroke_label ~ c(0, 1), 
                                       diabe_label ~ c(0, 1), 
                                       hearte_label ~ c(0, 1), 
                                       hibpe_label ~ c(0, 1), 
                                       smoken_label ~ c(0, 1), 
                                       race_label ~ c(0, 1), 
                                       employment_label ~ c(0, 1), 
                                       subjective_cognition_label ~ c(0, 1), 
                                       alcohol_consumption_label ~ c(0, 1), 
                                       impairment_label ~ c(0, 1)), 
                         #order variables
                         include = 
                           c(age, female_label, race_label, edyrs, 
                             employment_label, married_partnered_label, bmi, 
                             stroke_label, diabe_label, hearte_label, 
                             hibpe_label, smoken_label, 
                             alcohol_consumption_label, adl, iadl, immrc, delrc, 
                             ser7, cactus_label, scissor_label, pres_label, 
                             bwc20_label, hrs_cog, subjective_cognition_label, 
                             mmse_norm, animal_naming, wrc_yes, wrc_no, 
                             imm_story, del_story, imm_cp, del_cp, trailsA, 
                             impairment_label)) %>%
  
  #Renaming the header
  gtsummary::modify_header(label = "Variable") %>%
  
  #Moving labels from the bottom to next to each of the variables
  gtsummary::add_stat_label(location = "row") %>%
  
  #setting as a tibble so that it can be output to excel
  gtsummary::as_tibble()

#---- **save output ----
#Going to need to clean up this output in excel
# Need to divide count data by 25 to account for multiply imputed data
writexl::write_xlsx(gtsummary::as_tibble(table), path = paste0(
  path_to_box, "tables/appendix_E/tableE.1_imputed_sample_characteristics.xlsx"))

#---- in-text analyses ----
#---- **% in "Other" category ----
for(wave in c("a")){
  demdx_data_path <- paste0(path_to_box, "data/ADAMS/adams1", wave, "/adams1", 
                            wave, "da/ADAMS1", str_to_upper(wave), "D_R.da")
  demdx_dict_path <- paste0(path_to_box, "data/ADAMS/adams1", wave, "/adams1", 
                            wave, "sta/ADAMS1", str_to_upper(wave), "D_R.dct")
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

ADAMS_demdx %<>% mutate_at("HHIDPN", as.numeric)

ADAMS_Other <- ADAMS %>% left_join(ADAMS_demdx, by = "HHIDPN") %>% 
  filter(Adem_dx_cat == "Other")

#5: Parkinson's (1.2%)
#8: Normal pressure hydrocephalus (1.2%)
#21: Cognitive impairment secondary to vascular disease (15.1%)
#23: Depression (5.8%)
#25: Mental retardation (2.3%)
#26: Alcohol abuse (past) (3.5%)
#27: Alcohol abuse (current) (1.2%)
#28: Stroke (24.4%)
#29: Other neurological conditions (7.0%)
#30: Other medical conditions (38.4%)

table(ADAMS_Other$ADFDX1, useNA = "ifany")/nrow(ADAMS_Other)

#---- **% in "Dementia" category ----
ADAMS_Dementia <- ADAMS %>% left_join(ADAMS_demdx, by = "HHIDPN") %>% 
  filter(Adem_dx_cat == "Dementia")

#1: Probable AD (36.1%)
#2: Possible AD (41.1%)
#3: Probable vascular dementia (8.9%)
#4: Possible vascular dementia (8.2%)
#10: Dementia of undetermined etiology (5.1%)
#15: Alcoholic dementia (0.6%)

table(ADAMS_Dementia$ADFDX1, useNA = "ifany")/nrow(ADAMS_Dementia)

#---- **training sample contingency cell counts ----
contingency_cell_counts <- 
  ADAMS_train %>% unite("cell", c("Black", "Hispanic", "Astroke"), sep = "")

table(contingency_cell_counts$cell, contingency_cell_counts$Adem_dx_cat)
