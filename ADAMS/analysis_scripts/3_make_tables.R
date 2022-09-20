#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy", "labelled", 
       "gtsummary")

options(scipen = 999)

#---- read in data ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_train <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_train.csv")) %>% 
  mutate("dataset" = "ADAMS Training\n(Wave A)")

ADAMS_test <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_test.csv")) %>% 
  mutate("dataset" = "ADAMS Testing\n(Wave A)")

ADAMS <- rbind(ADAMS_train, ADAMS_test)

#---- select/order vars ----
#remove these variables
ADAMS %<>% dplyr::select(-one_of(c("HHIDPN", "(Intercept)")))

#variable order


#---- make table ----
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
