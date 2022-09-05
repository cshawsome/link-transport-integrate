#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "tidyr")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **clean HRS ----
HRS <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv")) %>% 
  dplyr::select(-one_of("GENDER", "RACE", "HISPANIC", "PIWTYPE", "PJ005M1", 
                        "r13drinkd", "r13drinkn", "r13pstmem", "GENDER_label", 
                        "RACE_label", "RACE_White", "RACE_Black", "RACE_Other", 
                        "HISPANIC_indicator", "ETHNIC_label", "PJ005M1_label", 
                        "PJ005M1_collapsed_label", "drinks_per_week", 
                        "drink_cat", "drink_cat_label"))

#---- **analytic HCAP ----
#need HRS bins to match these
HCAP <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv")) 
