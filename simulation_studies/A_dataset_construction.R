#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "MASS")

#---- source scripts ----
source(here("functions", "read_da_dct.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/"

#---- **RAND ----
rand_waves <- seq(5, 9, by = 1) #Corresponding to ADAMS +- 1 wave (for imputation)
rand_variables <- 
  c("hhidpn", 
    #Health and health behaviors (ever/never stroke, ever/never
    # hypertension, ever/never diabetes, ever/never cvd, ever/never cancer, 
    # BMI, IADLs, ADLs, gross motor, fine motor, depressive symptoms, 
    # self-reported change in memory, smokes now, number days drinking per week, 
    # number drinks/day) 
    paste0("r", rand_waves, "stroke"), paste0("r", rand_waves, "hibpe"), 
    paste0("r", rand_waves, "diabe"), paste0("r", rand_waves, "hearte"),
    paste0("r", rand_waves, "cancre"), paste0("r", rand_waves, "bmi"),
    paste0("r", rand_waves, "iadla"), paste0("r", rand_waves, "adla"),
    paste0("r", rand_waves, "grossa"), paste0("r", rand_waves, "finea"),
    paste0("r", rand_waves, "cesd"), paste0("r", rand_waves, "pstmem"), 
    paste0("r", rand_waves, "smoken"), paste0("r", rand_waves, "drinkd"),
    paste0("r", rand_waves, "drinkn"))

RAND <- read_dta(paste0(path_to_box, "Dissertation/data/HRS/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **ADAMS ----
ADAMS <- rbind(read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_test.csv")), 
               read_csv(file = paste0(path_to_box, 
                                      "Dissertation/data/ADAMS/cleaned/", 
                                      "ADAMS_train.csv")))






