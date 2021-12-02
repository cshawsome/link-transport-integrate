#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in ADAMS analytic dataset ----
ADAMS <- read_csv(file = paste0("/Users/CrystalShaw/Box/Dissertation/", 
                        "data/ADAMS/cleaned/ADAMS_subset_mixed.csv"))
