#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "openxlsx", "sjPlot", "here")

options(scipen = 999)

#---- read in data ----
#Continuous vars (notation from Schafer 1997)
Z <- c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
       "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv")) %>% 
  dplyr::select()

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

#---- Sigma priors ----
test <- ADAMS_subset %>% filter(Adem_dx_cat_collapse == "Normal")


