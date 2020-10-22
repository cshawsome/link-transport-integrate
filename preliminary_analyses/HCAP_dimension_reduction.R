#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "stringr")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16da/HC16HP_R.da")
dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16sta/HC16HP_R.dct")

HCAP <- read_da_dct(data_path, dict_path, HHIDPN = "TRUE")

#---- variables of interest ----
HCAP_vars <- c("HHIDPN", "MSE", "TICS", "WLIMM1", "WLIMM2", "WLIMM3", "1066", 
               "SCORE")
not_these <- c("TEST", "INTRO", "H1RMSESCORE", "H1RTICSSCORE", "H1R1066SCORE", 
               "SCORESE", "CESDSCORE", "STSCORE")

HCAP_assessment <- HCAP %>% dplyr::select(contains(HCAP_vars)) %>% 
  dplyr::select(-c(contains(not_these), "H1RMSE11T")) #H1RMSE11T1 needs to stay

#---- coding correct/incorrect answers ----
#Need to recode 2 = correct (different level) to 1
recode_correct <- c("H1RMSE17", "H1RMSE21")

#Need to recode 5 = error to 0
recode_incorrect <- 
  colnames(HCAP_assessment)[as.logical(
    str_detect(colnames(HCAP_assessment), "MSE") + 
      str_detect(colnames(HCAP_assessment), "TICS") + 
      str_detect(colnames(HCAP_assessment), "1066"))] 
recode_incorrect <- recode_incorrect[-which(recode_incorrect == "H1RMSE11T1" | 
                                              recode_incorrect == "H1RMSE13")]
for(var in recode_incorrect){
  HCAP_assessment[which(HCAP_assessment[, var] == 5), var] <- 0
}

for(var in recode_correct){
  HCAP_assessment[which(HCAP_assessment[, var] == 2), var] <- 1
}

# #Sanity check
# for(var in colnames(HCAP_assessment)){
#   print(var)
#   print(table(HCAP_assessment[, var], useNA = "ifany"))
# }

#---- code missigness ----
recode_missing_97 <- colnames(HCAP_assessment)[as.logical(
  str_detect(colnames(HCAP_assessment), "MSE") + 
    str_detect(colnames(HCAP_assessment), "H1RCPIMMSCORE") + 
    str_detect(colnames(HCAP_assessment), "H1RCPDELSCORE"))]

recode_missing_8 <- colnames(HCAP_assessment)[as.logical(
  str_detect(colnames(HCAP_assessment), "TICS") + 
    str_detect(colnames(HCAP_assessment), "1066"))]

recode_missing_996 <- c("H1RNSSCORE", "H1RTMASCORE", "H1RTMBSCORE")
  
for(var in recode_missing_97){
  HCAP_assessment[which(HCAP_assessment[, var] >= 97), var] <- NA
}

for(var in recode_missing_8){
  HCAP_assessment[which(HCAP_assessment[, var] >= 8), var] <- NA
}

for(var in recode_missing_996){
  HCAP_assessment[which(HCAP_assessment[, var] >= 996), var] <- NA
}

# #Sanity check
# for(var in colnames(HCAP_assessment)){
#   print(var)
#   print(table(HCAP_assessment[, var], useNA = "ifany"))
# }

#---- sum scores and averages of repeated trials ----
HCAP_assessment %<>% 
  mutate("H1RMSE12SCORE" = rowSums(HCAP_assessment %>% 
                                     dplyr::select(contains("H1RMSE12")), 
                                   na.rm = TRUE), 
         "H1RWLIMMSCORE" = apply(HCAP_assessment %>% 
                                      dplyr::select(contains("H1RWLIMM")), 1, 
                                 max, na.rm = TRUE)) 
HCAP_assessment[is.infinite(HCAP_assessment[, "H1RWLIMMSCORE"]), 
                "H1RWLIMMSCORE"] <- 0

#Sanity check
View(HCAP_assessment %>% dplyr::select(contains("H1RMSE12")))
View(HCAP_assessment %>% dplyr::select(contains("H1RWLIMM")))
table(HCAP_assessment$H1RWLIMMSCORE, useNA = "ifany")

#---- make correlation matrix ----


