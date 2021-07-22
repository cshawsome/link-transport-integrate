#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr")
install_github("thomasp85/gganimate")
library(gganimate)

#---- read in data ----
#---- **ADAMS ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"))

ADAMS_data <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                              "data/cleaned/ADAMS/ADAMS_subset_mixed.csv"))