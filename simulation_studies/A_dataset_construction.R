#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "MASS")

#---- read in ADAMS analytic dataset ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/"
ADAMS <- rbind(read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_test.csv")), 
               read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_train.csv")))

#---- generating synthetic data ----
cont_cols <- c("AAGE", "ANMSETOT_norm", "ANSER7T", "ANIMMCR", "ANRECYES", 
               "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

#---- **normal distribution ----
Sigma <- var(as.matrix(ADAMS %>% dplyr::select(all_of(cont_cols))))

synthetic_data <- mvrnorm(n = nrow(ADAMS), mu = rep(0, length(cont_cols)), 
                          Sigma = Sigma)




