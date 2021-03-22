#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "openxlsx", "sjPlot", "here")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv")) 

#---- data cleaning: dem dx ----
ADAMS_subset %<>% 
  mutate("Adem_dx_cat_collapse" = 
           case_when(Adem_dx_cat %in% 
                       c("Dementia", "Probable/Possible AD", "Probable Dementia", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia", 
                     TRUE ~ Adem_dx_cat))

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$Adem_dx_cat_collapse,
#       useNA = "ifany")


