#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer")

#---- read in data ----
synthetic_ADAMS <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "analyses/results/ADAMSA/ADAMSA_synthetic.csv"), 
           col_types = cols(HHIDPN = col_character()))

ADAMS_columns <- c(colnames(synthetic_ADAMS)[
  which(!colnames(synthetic_ADAMS) %in% 
          c("(Intercept)", "Group", "Black", "Hispanic", "p_Unimpaired", 
            "p_Other", "p_MCI"))], "Adem_dx_cat")

ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character())) %>% 
  #only take those part of modeling (complete cases only for now)
  filter(HHIDPN %in% synthetic_ADAMS$HHIDPN) %>% 
  dplyr::select(all_of(ADAMS_columns))
  
#---- 