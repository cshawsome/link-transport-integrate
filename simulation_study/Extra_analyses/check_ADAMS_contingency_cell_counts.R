#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **prior imputed clean ----
prior_imputed_clean <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/MI/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

single_ADAMS_analytic <- prior_imputed_clean[[1]]

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(single_ADAMS_analytic))

single_ADAMS_analytic %<>% 
  rename_at(., vars(variable_labels$ADAMS), ~ variable_labels$data_label)

#---- contingency cell IDs ----
#there is a cell with only 2 people in it, sticking with race/ethnicity x stroke instead
W <- c("black", "hispanic", "female", "stroke")

single_ADAMS_analytic %<>% 
  unite("cell_ID", all_of(W), sep = "", remove = FALSE)

table(single_ADAMS_analytic$cell_ID)
