#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy", "labelled", 
       "gtsummary", "writexl")

options(scipen = 999)

#---- source scripts ----
source(paste0(here::here("functions", "read_da_dct.R")))

#---- Table 3.2 ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_subset_mixed <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_subset_mixed.csv"))

ADAMS_train_IDs <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_train.csv")) %>% 
  dplyr::select("HHIDPN") %>% unlist() %>% unname()

ADAMS_test_IDs <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_test.csv")) %>% 
  dplyr::select("HHIDPN") %>% unlist() %>% unname()

ADAMS <- rbind(ADAMS_subset_mixed %>% filter(HHIDPN %in% c(ADAMS_train_IDs)) %>% 
                 mutate(dataset = "train"), 
               ADAMS_subset_mixed %>% filter(HHIDPN %in% c(ADAMS_test_IDs)) %>% 
                 mutate(dataset = "test")) 

#---- **set factor levels ----
ADAMS %<>% mutate_at("ETHNIC_label", function(x) 
  factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate_at("Adem_dx_cat", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other")))

#---- **label variables ----
ADAMS %<>% 
  labelled::set_variable_labels(AAGE = "Age",
                                ETHNIC_label = "Race/Ethnicity",
                                Abmi = "BMI",
                                Aiadla = "IADLs",
                                Astroke = "History of stroke",
                                ANSER7T = "Serial 7s", 
                                ANIMMCR = "Immediate word recall", 
                                ANDELCOR = "Delayed word recall", 
                                ANMSETOT_norm = "Total MMSE (normalized)", 
                                ANRECYES = "Word recall (yes)", 
                                ANWM1TOT = "Immediate story recall", 
                                proxy_cog = "Average Jorm IQCODE", 
                                Adem_dx_cat = "Adjudicated impairment")

#---- **label categories ----
ADAMS %<>% 
  labelled::set_value_labels(dataset = 
                               c("ADAMS Training\n (Wave A)" = "train", 
                                 "ADAMS Hold Out\n (Wave A)" = "test"),
                             Astroke = c("History of Stroke" = 1,
                                         "No History of Stroke" = 0))

#---- **make table ----
ADAMS_table <- ADAMS %>%
  #set labelled variables as factors 
  purrr::modify_if(labelled::is.labelled, labelled::to_factor) %>%
  gtsummary::tbl_summary(missing = "no", 
                         by = dataset,
                         statistic = list(all_continuous() ~ "{mean} ({sd})"),
                         #If we don't specify these variables types, all levels 
                         #  will be summarized
                         type = list(Astroke ~ "dichotomous", 
                                     Aiadla ~ "continuous", 
                                     ANSER7T ~ "continuous"),
                         
                         #Specifying the level of dichotomous the variable that 
                         #  should be displayed
                         value = list(Astroke = "History of Stroke"),
                         
                         #Specifying the number of decimal places for 
                         #  categorical vars
                         digits = list(all_continuous() ~ 1, 
                                       ETHNIC_label ~ c(0, 1),
                                       Astroke ~ c(0, 1), 
                                       Adem_dx_cat ~ c(0, 1)),
                         
                         #Specifying the exact variables I want in the table
                         include = c(AAGE, ETHNIC_label, Abmi, Astroke, Aiadla, 
                                     ANSER7T, ANIMMCR, ANDELCOR, ANMSETOT_norm, 
                                     ANRECYES, ANWM1TOT, proxy_cog, 
                                     Adem_dx_cat)) %>%
  
  #Renaming the header
  gtsummary::modify_header(label = "Variable") %>%
  
  #Adding in percents
  gtsummary::add_overall() %>%
  
  #Moving labels from the bottom to next to each of the variables
  gtsummary::add_stat_label(location = "row") %>%
  
  # #bolding the variables
  # gtsummary::bold_labels() %>%
  
  #setting as a tibble so that it can be output to excel
  gtsummary::as_tibble()

#---- **save output ----
writexl::write_xlsx(gtsummary::as_tibble(ADAMS_table), path = paste0(
  path_to_box, "tables/chapter_3/table3.2_ADAMS_characteristics.xlsx"))

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
