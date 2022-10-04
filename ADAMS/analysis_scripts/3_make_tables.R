#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy", "labelled", 
       "gtsummary")

options(scipen = 999)

#---- source scripts ----
source(paste0(here::here("functions", "read_da_dct.R")))

#---- Table 3.2 ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_train <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_train.csv")) %>% 
  mutate("dataset" = "ADAMS Training\n(Wave A)")

ADAMS_test <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_test.csv")) %>% 
  mutate("dataset" = "ADAMS Testing\n(Wave A)")

ADAMS <- rbind(ADAMS_train, ADAMS_test)

#---- **select/order vars ----
#remove these variables
ADAMS %<>% dplyr::select(c("AAGE", "ETHNIC_label", "Abmi", "Aiadla", "Astroke",
                           "ANSER7T", "ANIMMCR", "ANDELCOR", "ANMSETOT_norm",  
                           "ANRECYES", "ANWM1TOT", "proxy_cog", "Adem_dx_cat", 
                           "dataset"))

#---- **make table ----
tab1 <- ADAMS %>% 
  tbl_summary(by = dataset, 
              statistic = list(all_continuous() ~ "{mean} ({sd})"))


set_value_labels(high_accult = c("High"=1, "Low"=0)) %>% 
  modify_if(is.labelled, to_factor) %>% 
  modify_header(label = "**Variable**") %>%
  add_overall() %>% 
  add_stat_label(location = "row") %>%
  modify_spanning_header(c("stat_1","stat_2") ~ "**Acculturation Level**") %>%
  modify_spanning_header(starts_with("stat_") ~ "Table 1") %>% 
  bold_labels() 

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
