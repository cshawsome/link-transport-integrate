#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS/ADAMS_train.csv"))

#normal model predictors
norm_preds <- c("AAGE", "Black", "Hispanic", "ANMSETOT", "ANSER7T", "ANIMMCR", 
                "ANRECYES", "ANWM1TOT", "proxy_cog")

#other model predictors
other_preds <- c("AAGE", "ANMSETOT", "ANIMMCR", "ANDELCOR")

#mci model predictors
MCI_preds <- c("ANMSETOT", "ANIMMCR", "Aiadla", "Astroke", "Abmi")





