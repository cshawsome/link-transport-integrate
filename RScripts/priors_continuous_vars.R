#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "openxlsx", "sjPlot", "here")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv")) %>% 
  mutate_if(is.character, as.factor)

