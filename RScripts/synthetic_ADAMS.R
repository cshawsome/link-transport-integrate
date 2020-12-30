#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset.csv"), 
                         col_types = cols(HHIDPN = col_character()))
  
#---- **RAND ----
RAND_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/RAND_subset.csv"), 
                        col_types = cols(HHIDPN = col_character()))

#---- join data ----
all_data <- left_join(ADAMS_subset, RAND_subset, by = "HHIDPN")

#---- select variables ----
vars <- c("HHIDPN", "GENDER", "ETHNIC", "AAGE_cat", "EDYRS_cat", "ANMSETOT_cat", 
          paste0("r", c(5, 6, 7), "iadla_cat"))

analytical_sample <- all_data %>% dplyr::select(all_of(vars))

#Variable check
colSums(is.na(analytical_sample))

#---- Mixture Model ----
#---- **Step 1: Take SRS of data ----
samp_size <- floor(0.5*nrow(analytical_sample))
rsamp <- sample_n(analytical_sample, size = samp_size, replace = FALSE)

#---- **Step 2: Set priors ----
#number of group-level latent classes
group_class_n <- 4
#number of sub-latent classes
sub_class_n <- 5





