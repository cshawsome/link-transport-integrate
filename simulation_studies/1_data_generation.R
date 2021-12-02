#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in ADAMS analytic dataset ----
ADAMS <- rbind(read_csv(file = paste0("/Users/CrystalShaw/Box/Dissertation/", 
                        "data/ADAMS/cleaned/ADAMS_test.csv")), 
               read_csv(file = paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                      "data/ADAMS/cleaned/ADAMS_train.csv")))


