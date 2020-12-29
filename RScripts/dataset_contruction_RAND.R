#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "haven", "labelled")

#---- read in RAND data ----
rand_waves <- seq(5, 9, by = 1)
rand_variables <- c("hhidpn", paste0("r", rand_waves, "iadla"))

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- variable check ----
#We need to turn these variables into categories >= 1 
for(wave in rand_waves){
  print(paste0("Wave", wave))
  print(table(RAND[, paste0("r", wave, "iadla")]))
}

#---- clean: functional limitations ----
#---- **iadla ----
for(wave in rand_waves){
  var <- paste0("r", wave, "iadla")
  new_var <- paste0("r", wave, "iadla_cat")
  
  RAND[, new_var] <- as.numeric(unlist(RAND[, var] + 1))
}

#Sanity check
for(wave in rand_waves){
  var <- paste0("r", wave, "iadla")
  new_var <- paste0("r", wave, "iadla_cat")
  
  print(paste0("Wave", wave))
  print(table(as.numeric(unlist(RAND[, var])), 
              as.numeric(unlist(RAND[, new_var])), 
              useNA = "ifany"))
}

#---- drop original variables ----
drop_these <- c(paste0("r", rand_waves, "iadla"))

RAND %<>% dplyr::select(-all_of(drop_these))

#---- save dataset ----
write_csv(RAND, path = paste0("/Users/CrystalShaw/Box/Dissertation/", 
                              "data/cleaned/RAND_subset.csv"))

