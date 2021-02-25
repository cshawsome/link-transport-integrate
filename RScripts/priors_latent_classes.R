#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv"))

#---- dem class indicators ----
ADAMS_subset %<>% 
  #Impaired vs. Normal
  mutate("ANormal" = ifelse(Adem_dx_cat == "Normal", 1, 0),
         #Other vs. Dementia/MCI
         "AOther" = ifelse(Adem_dx_cat == "Other", 1, 0), 
         #MCI vs. Dementia
         "AMCI" = ifelse(Adem_dx_cat == "MCI", 1, 0))

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$ANormal, useNA = "ifany")
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AOther, useNA = "ifany")
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$AMCI, useNA = "ifany")
