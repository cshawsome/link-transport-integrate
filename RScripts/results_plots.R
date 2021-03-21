#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer")

#---- read in data ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character()))

synthetic_ADAMS <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "analyses/results/ADAMSA/ADAMSA_synthetic.csv"), 
           col_types = cols(HHIDPN = col_character()))
