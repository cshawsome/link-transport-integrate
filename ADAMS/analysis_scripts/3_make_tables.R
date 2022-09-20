#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy")

options(scipen = 999)

#---- read in data ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
ADAMS_train <- read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_train.csv"))
ADAMS_test <- read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_test.csv"))
