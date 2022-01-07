#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "MASS")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/"

#---- **ADAMS ----
ADAMS <- rbind(read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_test.csv")), 
               read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_train.csv")))






