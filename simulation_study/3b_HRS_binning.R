#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "tidyr")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **clean HCAP ----
HCAP <- read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_clean.csv")) %>% 
  dplyr::select(-one_of("H1RMSESCORE", "H1RWLIMM1SCORE", "H1RWLIMM2SCORE", 
                        "H1RWLIMM3SCORE", "H1RBMIMMSCORE", "H1RLMIMMSCORE", 
                        "H1RBMDELSCORE", "H1RLMDELSCORE", "HCAP_SELECT", 
                        "GENDER", "RACE", "HISPANIC", "PIWTYPE", "PJ005M1", 
                        "r13drinkd", "r13drinkn", "r13pstmem", "GENDER_label", 
                        "RACE_label", "RACE_White", "RACE_Black", "RACE_Other", 
                        "HISPANIC_indicator", "ETHNIC_label", "PJ005M1_label", 
                        "PJ005M1_collapsed_label", "r13drinks_per_week", 
                        "r13drink_cat", "r13drink_cat_label"))
