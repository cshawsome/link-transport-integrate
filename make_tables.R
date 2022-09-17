#Edit this in the future to just pull in the analytic datasets
#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy")

options(scipen = 999)

#---- ADAMS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

raw_ADAMS <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_subset_mixed.csv"))

ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/chunk_1/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric)) 

#one imputed subset
ADAMS_subset <- ADAMS_imputed_clean[[1]]

#ADAMS training data
ADAMS_train <- read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_train.csv"))
ADAMS_train_IDs <- ADAMS_train$HHIDPN
ADAMS_train_raw <- raw_ADAMS[which(raw_ADAMS$HHIDPN %in% ADAMS_train_IDs), ]

#ADAMS hold-out data
ADAMS_test <- read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_test.csv"))
ADAMS_test_IDs <- ADAMS_test$HHIDPN
ADAMS_test_raw <- raw_ADAMS[which(raw_ADAMS$HHIDPN %in% ADAMS_test_IDs), ]

#---- **age ----
#---- ****train ----
mean(ADAMS_train_raw$AAGE)
sd(ADAMS_train_raw$AAGE)

#---- ****test ----
mean(ADAMS_test_raw$AAGE)
sd(ADAMS_test_raw$AAGE)

#---- **race/ethnicity ----
#---- ****imputed subset ----
mean(ADAMS_subset$White)
mean(ADAMS_subset$Black)
mean(ADAMS_subset$Hispanic)

#---- ****train ----
sum(ADAMS_train_raw$Hispanic)
mean(ADAMS_train_raw$Hispanic)

sum(ADAMS_train_raw$White)
mean(ADAMS_train_raw$White)

sum(ADAMS_train_raw$Black)
mean(ADAMS_train_raw$Black)

#---- ****test ----
sum(ADAMS_test_raw$Hispanic)
mean(ADAMS_test_raw$Hispanic)

sum(ADAMS_test_raw$White)
mean(ADAMS_test_raw$White)

sum(ADAMS_test_raw$Black)
mean(ADAMS_test_raw$Black)

#---- **normalized MMSE ----
#---- ****train ----
mean(ADAMS_train_raw$ANMSETOT_norm)
sd(ADAMS_train_raw$ANMSETOT_norm)

#---- ****test ----
mean(ADAMS_test_raw$ANMSETOT_norm)
sd(ADAMS_test_raw$ANMSETOT_norm)

#---- **immediate word recall ----
#---- ****train ----
mean(ADAMS_train_raw$ANIMMCR)
sd(ADAMS_train_raw$ANIMMCR)

#---- ****test ----
mean(ADAMS_test_raw$ANIMMCR)
sd(ADAMS_test_raw$ANIMMCR)

#---- **serial 7s ----
#---- ****train ----
mean(ADAMS_train_raw$ANSER7T)
sd(ADAMS_train_raw$ANSER7T)

#---- ****test ----
mean(ADAMS_test_raw$ANSER7T)
sd(ADAMS_test_raw$ANSER7T)

#---- **word list recall (yes) ----
#---- ****train ----
mean(ADAMS_train_raw$ANRECYES)
sd(ADAMS_train_raw$ANRECYES)

#---- ****test ----
mean(ADAMS_test_raw$ANRECYES)
sd(ADAMS_test_raw$ANRECYES)

#---- **story recall (immediate) ----
#---- ****train ----
mean(ADAMS_train_raw$ANWM1TOT)
sd(ADAMS_train_raw$ANWM1TOT)

#---- ****test ----
mean(ADAMS_test_raw$ANWM1TOT)
sd(ADAMS_test_raw$ANWM1TOT)

#---- **JORM IQ code ----
#---- ****train ----
mean(ADAMS_train_raw$proxy_cog)
sd(ADAMS_train_raw$proxy_cog)

#---- ****test ----
mean(ADAMS_test_raw$proxy_cog)
sd(ADAMS_test_raw$proxy_cog)

#---- **delayed word recall ----
#---- ****train ----
mean(ADAMS_train_raw$ANDELCOR)
sd(ADAMS_train_raw$ANDELCOR)

#---- ****test ----
mean(ADAMS_test_raw$ANDELCOR)
sd(ADAMS_test_raw$ANDELCOR)

#---- **BMI ----
#---- ****train ----
mean(ADAMS_train_raw$Abmi)
sd(ADAMS_train_raw$Abmi)

#---- ****test ----
mean(ADAMS_test_raw$Abmi)
sd(ADAMS_test_raw$Abmi)

#---- **IADLs ----
#---- ****train ----
mean(ADAMS_train_raw$Aiadla)
sd(ADAMS_train_raw$Aiadla)

#---- ****test ----
mean(ADAMS_test_raw$Aiadla)
sd(ADAMS_test_raw$Aiadla)

#---- **stroke history ----
#---- ****train ----
sum(ADAMS_train_raw$Astroke)
mean(ADAMS_train_raw$Astroke)

#---- ****test ----
sum(ADAMS_test_raw$Astroke)
mean(ADAMS_test_raw$Astroke)

#---- **dementia adjudication ----
#---- ****train ----
table(ADAMS_train_raw$Adem_dx_cat)
table(ADAMS_train_raw$Adem_dx_cat)/nrow(ADAMS_train_raw)

#---- ****test ----
table(ADAMS_test_raw$Adem_dx_cat)
table(ADAMS_test_raw$Adem_dx_cat)/nrow(ADAMS_test_raw)

#---- Synthetic Superpop ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

superpopulations_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpopulations_data <- 
  do.call(rbind, lapply(superpopulations_paths, read_results))

#---- **race/ethnicity ----
mean(superpopulations_data$White)
mean(superpopulations_data$black)
mean(superpopulations_data$hispanic)

#---- Synthetic HRS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_HRS_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HRS_list"))

#filter to normal data for now
dataset_names <- 
  unlist(lapply(synthetic_HRS_list, function(x) unique(x$dataset_name)))

indices <- which(dataset_names %in% 
                   paste0("normal_", c(500, 1000, 2000, 4000, 8000), "_ADAMS"))

#---- **race/ethnicity ----
lapply(synthetic_HRS_list[indices], function(x) mean(x$White))
lapply(synthetic_HRS_list[indices], function(x) mean(x$black))
lapply(synthetic_HRS_list[indices], function(x) mean(x$hispanic))

#---- Synthetic HCAP ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HCAP_list"))

#filter to normal data for now
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

indices <- which(dataset_names %in% 
                   paste0("normal_", c(500, 1000, 2000, 4000, 8000), "_ADAMS"))

#---- **race/ethnicity ----
lapply(synthetic_HCAP_list[indices], function(x) mean(x$White))
lapply(synthetic_HCAP_list[indices], function(x) mean(x$black))
lapply(synthetic_HCAP_list[indices], function(x) mean(x$hispanic))

#---- HRS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

HRS_clean <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv"))

HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

#---- **clean ----
#---- ****married/partnered ----
table(HRS_clean$PCOUPLE, useNA = "ifany")
table(HRS_clean$PCOUPLE, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$PCOUPLE, HRS_clean$White, useNA = "ifany")
table(HRS_clean$PCOUPLE, HRS_clean$White, useNA = "ifany")/sum(HRS_clean$White)

table(HRS_clean$PCOUPLE, HRS_clean$Black, useNA = "ifany")
table(HRS_clean$PCOUPLE, HRS_clean$Black, useNA = "ifany")/sum(HRS_clean$Black)

table(HRS_clean$PCOUPLE, HRS_clean$Hispanic, useNA = "ifany")
table(HRS_clean$PCOUPLE, HRS_clean$Hispanic, useNA = "ifany")/
  sum(HRS_clean$Hispanic)

#---- ****female ----
table(HRS_clean$Female, useNA = "ifany")
table(HRS_clean$Female, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$Female, HRS_clean$White, useNA = "ifany")
table(HRS_clean$Female, HRS_clean$White, useNA = "ifany")/sum(HRS_clean$White)

table(HRS_clean$Female, HRS_clean$Black, useNA = "ifany")
table(HRS_clean$Female, HRS_clean$Black, useNA = "ifany")/sum(HRS_clean$Black)

table(HRS_clean$Female, HRS_clean$Hispanic, useNA = "ifany")
table(HRS_clean$Female, HRS_clean$Hispanic, useNA = "ifany")/
  sum(HRS_clean$Hispanic)

#---- ****age ----
summary(HRS_clean$PAGE)
summary(HRS_clean %>% filter(White == 1) %>% dplyr::select("PAGE"))
summary(HRS_clean %>% filter(Black == 1) %>% dplyr::select("PAGE"))
summary(HRS_clean %>% filter(Hispanic == 1) %>% dplyr::select("PAGE"))

#---- ****race/ethnicity ----
table(HRS_clean$White, useNA = "ifany")
table(HRS_clean$White, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$Black, useNA = "ifany")
table(HRS_clean$Black, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$Hispanic, useNA = "ifany")
table(HRS_clean$Hispanic, useNA = "ifany")/nrow(HRS_clean)

#---- ****employment ----
table(HRS_clean$Working, useNA = "ifany")
table(HRS_clean$Working, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$Working, HRS_clean$White, useNA = "ifany")
table(HRS_clean$Working, HRS_clean$White, useNA = "ifany")/sum(HRS_clean$White)

table(HRS_clean$Working, HRS_clean$Black, useNA = "ifany")
table(HRS_clean$Working, HRS_clean$Black, useNA = "ifany")/sum(HRS_clean$Black)

table(HRS_clean$Working, HRS_clean$Hispanic, useNA = "ifany")
table(HRS_clean$Working, HRS_clean$Hispanic, useNA = "ifany")/
  sum(HRS_clean$Hispanic)

table(HRS_clean$Retired, useNA = "ifany")
table(HRS_clean$Retired, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$Retired, HRS_clean$White, useNA = "ifany")
table(HRS_clean$Retired, HRS_clean$White, useNA = "ifany")/sum(HRS_clean$White)

table(HRS_clean$Retired, HRS_clean$Black, useNA = "ifany")
table(HRS_clean$Retired, HRS_clean$Black, useNA = "ifany")/sum(HRS_clean$Black)

table(HRS_clean$Retired, HRS_clean$Hispanic, useNA = "ifany")
table(HRS_clean$Retired, HRS_clean$Hispanic, useNA = "ifany")/
  sum(HRS_clean$Hispanic)

table(HRS_clean$`Not working`, useNA = "ifany")
table(HRS_clean$`Not working`, useNA = "ifany")/nrow(HRS_clean)

table(HRS_clean$`Not working`, HRS_clean$White, useNA = "ifany")
table(HRS_clean$`Not working`, HRS_clean$White, useNA = "ifany")/
  sum(HRS_clean$White)

table(HRS_clean$`Not working`, HRS_clean$Black, useNA = "ifany")
table(HRS_clean$`Not working`, HRS_clean$Black, useNA = "ifany")/
  sum(HRS_clean$Black)

table(HRS_clean$`Not working`, HRS_clean$Hispanic, useNA = "ifany")
table(HRS_clean$`Not working`, HRS_clean$Hispanic, useNA = "ifany")/
  sum(HRS_clean$Hispanic)

#---- ****missingness ----
colSums(is.na(HRS_clean))/nrow(HRS_clean)
colSums(is.na(HRS_clean %>% filter(White == 1)))/sum(HRS_clean$White)
colSums(is.na(HRS_clean %>% filter(Black == 1)))/sum(HRS_clean$Black)
colSums(is.na(HRS_clean %>% filter(Hispanic == 1)))/sum(HRS_clean$Hispanic)


#---- **CC ----
#---- ****married/partnered ----
table(HRS_analytic$PCOUPLE, useNA = "ifany")
table(HRS_analytic$PCOUPLE, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$PCOUPLE, HRS_analytic$White, useNA = "ifany")
table(HRS_analytic$PCOUPLE, HRS_analytic$White, useNA = "ifany")/
  sum(HRS_analytic$White)

table(HRS_analytic$PCOUPLE, HRS_analytic$Black, useNA = "ifany")
table(HRS_analytic$PCOUPLE, HRS_analytic$Black, useNA = "ifany")/
  sum(HRS_analytic$Black)

table(HRS_analytic$PCOUPLE, HRS_analytic$Hispanic, useNA = "ifany")
table(HRS_analytic$PCOUPLE, HRS_analytic$Hispanic, useNA = "ifany")/
  sum(HRS_analytic$Hispanic)

#---- ****female ----
table(HRS_analytic$Female, useNA = "ifany")
table(HRS_analytic$Female, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$Female, HRS_analytic$White, useNA = "ifany")
table(HRS_analytic$Female, HRS_analytic$White, useNA = "ifany")/
  sum(HRS_analytic$White)

table(HRS_analytic$Female, HRS_analytic$Black, useNA = "ifany")
table(HRS_analytic$Female, HRS_analytic$Black, useNA = "ifany")/
  sum(HRS_analytic$Black)

table(HRS_analytic$Female, HRS_analytic$Hispanic, useNA = "ifany")
table(HRS_analytic$Female, HRS_analytic$Hispanic, useNA = "ifany")/
  sum(HRS_analytic$Hispanic)

#---- ****age ----
summary(HRS_analytic$PAGE)
summary(HRS_analytic %>% filter(White == 1) %>% dplyr::select("PAGE"))
summary(HRS_analytic %>% filter(Black == 1) %>% dplyr::select("PAGE"))
summary(HRS_analytic %>% filter(Hispanic == 1) %>% dplyr::select("PAGE"))

#---- ****race/ethnicity ----
table(HRS_analytic$White, useNA = "ifany")
table(HRS_analytic$White, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$Black, useNA = "ifany")
table(HRS_analytic$Black, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$Hispanic, useNA = "ifany")
table(HRS_analytic$Hispanic, useNA = "ifany")/nrow(HRS_analytic)

#---- ****employment ----
table(HRS_analytic$Working, useNA = "ifany")
table(HRS_analytic$Working, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$Working, HRS_analytic$White, useNA = "ifany")
table(HRS_analytic$Working, HRS_analytic$White, useNA = "ifany")/sum(HRS_analytic$White)

table(HRS_analytic$Working, HRS_analytic$Black, useNA = "ifany")
table(HRS_analytic$Working, HRS_analytic$Black, useNA = "ifany")/sum(HRS_analytic$Black)

table(HRS_analytic$Working, HRS_analytic$Hispanic, useNA = "ifany")
table(HRS_analytic$Working, HRS_analytic$Hispanic, useNA = "ifany")/
  sum(HRS_analytic$Hispanic)

table(HRS_analytic$Retired, useNA = "ifany")
table(HRS_analytic$Retired, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$Retired, HRS_analytic$White, useNA = "ifany")
table(HRS_analytic$Retired, HRS_analytic$White, useNA = "ifany")/sum(HRS_analytic$White)

table(HRS_analytic$Retired, HRS_analytic$Black, useNA = "ifany")
table(HRS_analytic$Retired, HRS_analytic$Black, useNA = "ifany")/sum(HRS_analytic$Black)

table(HRS_analytic$Retired, HRS_analytic$Hispanic, useNA = "ifany")
table(HRS_analytic$Retired, HRS_analytic$Hispanic, useNA = "ifany")/
  sum(HRS_analytic$Hispanic)

table(HRS_analytic$`Not working`, useNA = "ifany")
table(HRS_analytic$`Not working`, useNA = "ifany")/nrow(HRS_analytic)

table(HRS_analytic$`Not working`, HRS_analytic$White, useNA = "ifany")
table(HRS_analytic$`Not working`, HRS_analytic$White, useNA = "ifany")/
  sum(HRS_analytic$White)

table(HRS_analytic$`Not working`, HRS_analytic$Black, useNA = "ifany")
table(HRS_analytic$`Not working`, HRS_analytic$Black, useNA = "ifany")/
  sum(HRS_analytic$Black)

table(HRS_analytic$`Not working`, HRS_analytic$Hispanic, useNA = "ifany")
table(HRS_analytic$`Not working`, HRS_analytic$Hispanic, useNA = "ifany")/
  sum(HRS_analytic$Hispanic)

#---- OLD ----

#---- source scripts ----
#custom read da dct function
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- read in data ----
#---- **HRS ----
#---- ****HRS tracker ----
HRS_tracker_data_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS_tracker/trk2018_3/", 
         "TRK2018TR_R.da")
HRS_tracker_dict_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS_tracker/trk2018_3/", 
         "TRK2018TR_R.dct")

HRS_tracker <- read_da_dct(HRS_tracker_data_path, HRS_tracker_dict_path, 
                           HHIDPN = "TRUE") %>% 
  #Participated in 2016 HRS
  filter(PIWTYPE == 1)  

#---- ****HRS waves ----
#ADAMS waves 5, 6, 7 = 2001, 2002/2003, 2004
#HRS wave 13 = 2016 
#list of variables to read in from HRS data:
# ID, Household Couple Category, Age, Race/Ethnicity, BMI, IADLs, Stroke History, 
# Serial 7s, Immediate Word Recall, Delayed Word Recall

hrs_waves <- c(5, 6, 7, 13)
hrs_vars = c("hhid", "pn", "rahhidpn", "ragender", "raracem", "rahispan", 
             "r13mstat",
             paste0("h", hrs_waves, "cpl"), paste0("r", hrs_waves, "agey_b"),
             paste0("r", hrs_waves, "bmi"), paste0("r", hrs_waves, "iadla"),
             paste0("r", hrs_waves, "stroke"), paste0("r", hrs_waves, "ser7"),
             paste0("r", hrs_waves, "imrc"), paste0("r", hrs_waves, "dlrc"), 
             paste0("r", hrs_waves, "cogtot"))

HRS <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS/", 
                       "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                col_select = all_of(hrs_vars)) 
colnames(HRS)[which(colnames(HRS) == "hhid")] <- "HHID"
colnames(HRS)[which(colnames(HRS) == "pn")] <- "PN"
colnames(HRS)[which(colnames(HRS) == "rahhidpn")] <- "HHIDPN"

HRS_2016 <- HRS %>% 
  #filter to tracker data
  filter(HHIDPN %in% HRS_tracker$HHIDPN) %>% 
  #age 65+
  filter(r13agey_b >= 65) %>%
  #race/ethnicity
  mutate("race_eth" = case_when(rahispan == 1 ~ "Hispanic", 
                                rahispan == 0 & raracem == 1 ~ "White", 
                                rahispan == 0 & raracem == 2 ~ "Black", 
                                rahispan == 0 & raracem == 3 ~ "Other")) %>%
  #missing race/ethnicity data (n = 6)
  filter(!is.na(race_eth))

#---- ******female ----
#1 = male; 2 = female
table(HRS_2016$ragender, useNA = "ifany")
mean(HRS_2016$ragender == 2)

#---- ******hrs cognitive score ----
summary(HRS_2016$r13cogtot, useNA = "ifany")
sd(HRS_2016$r13cogtot, na.rm = TRUE)

#---- ******coupled status ----
table(HRS_2016$h13cpl, useNA = "ifany")
summary(HRS_2016$h13cpl)

#---- ******age ----
hist(HRS_2016$r13agey_b)
summary(HRS_2016$r13agey_b)
sd(HRS_2016$r13agey_b)

#---- ******race/ethnicity ----
table(HRS_2016$race_eth, useNA = "ifany")
table(HRS_2016$race_eth, useNA = "ifany")/nrow(HRS_2016)

#---- ******BMI ----
summary(HRS_2016$r13bmi)
sd(HRS_2016$r13bmi, na.rm = TRUE)

#---- ******IADL ----
hist(HRS_2016$r13iadla)
summary(HRS_2016$r13iadla)
sd(HRS_2016$r13iadla, na.rm = TRUE)

#---- ******stroke history ----
table(HRS_2016$r13stroke, useNA = "ifany")
summary(HRS_2016$r13stroke)

#---- ******serial 7s ----
summary(HRS_2016$r13ser7)
sd(HRS_2016$r13ser7, na.rm = TRUE)

#---- ******immediate word recall ----
summary(HRS_2016$r13imrc)
sd(HRS_2016$r13imrc, na.rm = TRUE)

#---- ******delayed word recall ----
summary(HRS_2016$r13dlrc)
sd(HRS_2016$r13dlrc, na.rm = TRUE)

#---- **ADAMS ----
#---- ****ADAMS tracker ----
ADAMS_tracker_data_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
         "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/",
         "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

#filter to those who completed Wave A assessment
ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") %>% filter(AASSESS == 1) 

#---- ****ADAMS neuropsych ----
for(wave in c("a")){
  neuropsych_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                                 "ADAMS1", str_to_upper(wave), "N_R.da")
  
  neuropsych_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                                 "ADAMS1", str_to_upper(wave), "N_R.dct")
  
  if(wave == "a"){
    ADAMS_neuropsych <- read_da_dct(neuropsych_data_path, neuropsych_dict_path, 
                                    HHIDPN = "TRUE") 
  } else{
    ADAMS_neuropsych %<>% 
      left_join(., read_da_dct(neuropsych_data_path, neuropsych_dict_path, 
                               HHIDPN = "TRUE"), 
                by = "HHIDPN")
  }
}

#---- ****ADAMS informant questionnaire ----
for(wave in c("a")){
  proxy_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                            "ADAMS1", str_to_upper(wave), "G_R.da")
  proxy_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                            "ADAMS1", str_to_upper(wave), "G_R.dct")
  if(wave == "a"){
    ADAMS_proxy <- read_da_dct(proxy_data_path, proxy_dict_path, 
                               HHIDPN = "TRUE") %>% 
      dplyr::select("HHIDPN", paste0(str_to_upper(wave), "GQ", seq(14, 29))) 
  } else{
    ADAMS_proxy %<>% 
      left_join(., read_da_dct(proxy_data_path, proxy_dict_path, 
                               HHIDPN = "TRUE") %>% 
                  dplyr::select("HHIDPN", 
                                paste0(str_to_upper(wave), "GQ", seq(14, 29))))
  }
}

#---- ****ADAMS dem dx ----
for(wave in c("a")){
  demdx_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                            "ADAMS1", str_to_upper(wave), "D_R.da")
  demdx_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                            "ADAMS1", str_to_upper(wave), "D_R.dct")
  if(wave == "a"){
    ADAMS_demdx <- read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
      dplyr::select("HHIDPN", paste0(str_to_upper(wave), "DFDX1")) 
  } else{
    ADAMS_demdx %<>% 
      left_join(., read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
                  dplyr::select("HHIDPN", 
                                paste0(str_to_upper(wave), "DFDX1")))
  }
}

for(wave in c("A")){
  dem_dx_var <- paste0(wave, "DFDX1")
  ADAMS_demdx[, paste0(wave, "dem_dx_cat")] <- 
    case_when(ADAMS_demdx[, dem_dx_var] %in% c(1, 2) ~ "Probable/Possible AD", 
              ADAMS_demdx[, dem_dx_var] %in% c(3, 4) ~ 
                "Probable/Possible Vascular Dementia", 
              ADAMS_demdx[, dem_dx_var] %in% 
                c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 21, 28, 29, 
                  30, 33) ~ "Other",
              ADAMS_demdx[, dem_dx_var] %in% c(18, 32) ~ "Probable Dementia",
              ADAMS_demdx[, dem_dx_var] %in% c(10, 13, 15, 16, 17, 19) ~ 
                "Dementia", 
              ADAMS_demdx[, dem_dx_var] %in% c(20, 22) ~ "MCI", 
              ADAMS_demdx[, dem_dx_var] == 31 ~ "Normal")
}

#Further collapsing categories
for(wave in c("A")){
  ADAMS_demdx[, paste0(wave, "dem_dx_cat")] <- 
    case_when(ADAMS_demdx[, paste0(wave, "dem_dx_cat")] %in% 
                c("Probable/Possible AD", 
                  "Probable/Possible Vascular Dementia", 
                  "Probable Dementia", "Dementia") ~ "Dementia",
              ADAMS_demdx[, paste0(wave, "dem_dx_cat")] == "Other" ~ "Other", 
              ADAMS_demdx[, paste0(wave, "dem_dx_cat")] == "MCI" ~ "MCI", 
              ADAMS_demdx[, paste0(wave, "dem_dx_cat")] == "Normal" ~ 
                "Unimpaired")
}

# #Sanity check
# for(wave in c("A", "B", "C", "D")){
#   dem_dx_var <- paste0(wave, "DFDX1")
#   print(paste0("Wave ", wave))
#   print(table(ADAMS_demdx[, dem_dx_var],
#               ADAMS_demdx[, paste0(wave, "dem_dx_cat")],
#               useNA = "ifany"))
# }

#Remove original variables
ADAMS_demdx %<>% dplyr::select(-c(paste0(c("A"), "DFDX1")))

#---- ****join data ----
ADAMS_A <- left_join(ADAMS_tracker, ADAMS_neuropsych, by = "HHIDPN") %>% 
  left_join(., ADAMS_proxy, by = "HHIDPN") %>% 
  left_join(., ADAMS_demdx, by = "HHIDPN") %>% 
  left_join(., HRS, by = "HHIDPN")

wave_updated_vars <- c("cpl", "bmi", "iadla", "stroke", "cogtot")

for(var in wave_updated_vars){
  if(var == "cpl"){
    ADAMS_A %<>% 
      mutate(!!paste0("A", var) := 
               case_when(AYEAR == 2001 ~ !!sym(paste0("h5", var)), 
                         AYEAR %in% c(2002, 2003) ~ !!sym(paste0("h6", var)), 
                         AYEAR == 2004 ~ !!sym(paste0("h7", var))))
    
  } else{
    ADAMS_A %<>% 
      mutate(!!paste0("A", var) := 
               case_when(AYEAR == 2001 ~ !!sym(paste0("r5", var)), 
                         AYEAR %in% c(2002, 2003) ~ !!sym(paste0("r6", var)), 
                         AYEAR == 2004 ~ !!sym(paste0("r7", var))))
  }
}

ADAMS_A %<>% 
  #race/ethnicity
  mutate("race_eth" = case_when(rahispan == 1 ~ "Hispanic",
                                rahispan == 0 & raracem == 1 ~ "White",
                                rahispan == 0 & raracem == 2 ~ "Black",
                                rahispan == 0 & raracem == 3 ~ "Other")) %>%
  #clean serial 7s data (missing/refused)
  mutate_at(.vars = c("ANSER7T"), function(x) ifelse(x > 5, NA, x)) %>%
  #clean immediate word recall
  mutate_at(.vars = c("ANIMMCR1", "ANIMMCR2", "ANIMMCR3", "ANDELCOR"), 
            #Missing/refused  
            function(x) ifelse(x > 10, NA, x)) %>% 
  #Best of 3 immediate recall trials
  mutate("ANIMMCR" = pmax(ANIMMCR1, ANIMMCR2, ANIMMCR3, na.rm = TRUE)) %>%
  #clean MMSE data
  mutate_at(.vars = "ANMSETOT", function(x) ifelse(x > 30, NA, x)) %>% 
  #normalize MMSE
  mutate("ANMSETOT_norm" = normMMSE(ANMSETOT)) %>%
  #clean word list recognition (missing/refused)
  mutate_at(.vars = c("ANRECYES"), function(x) ifelse(x > 10, NA, x)) %>% 
  #clean immediate story recall (missing/refused)
  mutate_at(.vars = c("ANWM1TOT"), function(x) ifelse(x > 37, NA, x)) %>%
  #proxy cognition 
  mutate("proxy_cog" = ADAMS_A %>% dplyr::select(contains("AGQ")) %>% 
           rowMeans(., na.rm = TRUE))  

#no self-response neuropsych
ADAMS_A %<>% 
  mutate("no_data" = 
           rowSums(!is.na(ADAMS_A %>% 
                            dplyr::select(c("ANSER7T", "ANIMMCR", "ANDELCOR", 
                                            "ANMSETOT_norm", "ANRECYES", 
                                            "ANWM1TOT")))))
#filter data
#no self-response on neuropsych
sum(ADAMS_A$no_data == 0)
ADAMS_A %<>% filter(no_data > 0)

#---- ******female ----
#1 = male; 2 = female
table(ADAMS_A$ragender, useNA = "ifany")
mean(ADAMS_A$ragender == 2)

#---- ******hrs cognitive score ----
summary(ADAMS_A$Acogtot, useNA = "ifany")
sd(ADAMS_A$Acogtot, na.rm = TRUE)

#---- ******coupled status ----
table(ADAMS_A$Acpl, useNA = "ifany")
summary(ADAMS_A$Acpl)

#---- ******age ----
hist(ADAMS_A$AAGE)
summary(ADAMS_A$AAGE)
sd(ADAMS_A$AAGE)

#---- ******race/ethnicity ----
table(ADAMS_A$race_eth, useNA = "ifany")
table(ADAMS_A$race_eth, useNA = "ifany")/nrow(ADAMS_A)

#---- ******BMI ----
summary(ADAMS_A$Abmi)
sd(ADAMS_A$Abmi, na.rm = TRUE)

#---- ******IADL ----
hist(ADAMS_A$Aiadla)
summary(ADAMS_A$Aiadla)
sd(ADAMS_A$Aiadla, na.rm = TRUE)

#---- ******stroke history ----
table(ADAMS_A$Astroke, useNA = "ifany")
summary(ADAMS_A$Astroke)

#---- ******serial 7s ----
summary(ADAMS_A$ANSER7T)
sd(ADAMS_A$ANSER7T, na.rm = TRUE)

#---- ******immediate word recall ----
summary(ADAMS_A$ANIMMCR)
sd(ADAMS_A$ANIMMCR, na.rm = TRUE)

#---- ******delayed word recall ----
summary(ADAMS_A$ANDELCOR)
sd(ADAMS_A$ANDELCOR, na.rm = TRUE)

#---- ******normed MMSE----
summary(ADAMS_A$ANMSETOT_norm)
sd(ADAMS_A$ANMSETOT_norm, na.rm = TRUE)

#---- ******word list recognition (yes)----
summary(ADAMS_A$ANRECYES)
sd(ADAMS_A$ANRECYES, na.rm = TRUE)

#---- ******story recall (immediate)----
summary(ADAMS_A$ANWM1TOT)
sd(ADAMS_A$ANWM1TOT, na.rm = TRUE)

#---- ******proxy cognition----
summary(ADAMS_A$proxy_cog)
sd(ADAMS_A$proxy_cog, na.rm = TRUE)

#---- ******dementia adjudication ----
table(ADAMS_A$Adem_dx_cat)
table(ADAMS_A$Adem_dx_cat)/nrow(ADAMS_A)

#---- **HCAP ----
HCAP_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                         "HCAP/HC16/HC16da/HC16HP_R.da")

HCAP_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                         "HCAP/HC16/HC16sta/HC16HP_R.dct")

#---- **HCAP informant ----
HCAP_informant_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                   "HCAP/HC16/HC16da/HC16HP_I.da")

HCAP_informant_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                   "HCAP/HC16/HC16sta/HC16HP_I.dct")

#---- ****join data ----
HCAP_2016 <- read_da_dct(HCAP_data_path, HCAP_dict_path, HHIDPN = "TRUE") %>%
  #join informant data
  left_join(., read_da_dct(HCAP_informant_data_path, HCAP_informant_dict_path, 
                           HHIDPN = "TRUE"), 
            by = "HHIDPN") %>% 
  #filter to tracker data
  left_join(., HRS, by = "HHIDPN") %>% 
  #round age up to 65 (there's some 64 years olds)
  mutate("age" = ifelse(r13agey_b < 65, 65, r13agey_b)) %>%
  #race/ethnicity
  mutate("race_eth" = case_when(rahispan == 1 ~ "Hispanic", 
                                rahispan == 0 & raracem == 1 ~ "White", 
                                rahispan == 0 & raracem == 2 ~ "Black", 
                                rahispan == 0 & raracem == 3 ~ "Other")) %>% 
  filter(!is.na(race_eth)) %>% 
  #Best of 3 immediate recall trials
  mutate("H1RWLIMMSCORE" = pmax(H1RWLIMM1SCORE, H1RWLIMM2SCORE, H1RWLIMM3SCORE, 
                                na.rm = TRUE)) %>%
  #normalize MMSE
  mutate("H1RMSESCORE_norm" = normMMSE(H1RMSESCORE)) %>%
  #sum score immediate story recall 
  mutate("story_sum_immediate" = H1RBMIMMSCORE + H1RLMIMMSCORE)

#no self-response neuropsych
HCAP_2016 %<>% 
  mutate("no_data" = 
           rowSums(!is.na(HCAP_2016 %>% 
                            dplyr::select(c("r13ser7", "H1RWLIMMSCORE", 
                                            "H1RWLDELSCORE", 
                                            "H1RMSESCORE_norm", "H1RWLRECYSCORE", 
                                            "story_sum_immediate")))))
#filter data
#no self-response on neuropsych
sum(HCAP_2016$no_data == 0)
HCAP_2016 %<>% filter(no_data > 0)

#---- ******female ----
#1 = male; 2 = female
table(HCAP_2016$ragender, useNA = "ifany")
mean(HCAP_2016$ragender == 2)

#---- ******hrs cognitive score ----
summary(HCAP_2016$r13cogtot, useNA = "ifany")
sd(HCAP_2016$r13cogtot, na.rm = TRUE)

#---- ******coupled status ----
table(HCAP_2016$h13cpl, useNA = "ifany")
summary(HCAP_2016$h13cpl)

#---- ******age ----
summary(HCAP_2016$age)
sd(HCAP_2016$age)

#---- ******race/ethnicity ----
table(HCAP_2016$race_eth, useNA = "ifany")
table(HCAP_2016$race_eth, useNA = "ifany")/nrow(HCAP_2016)

#---- ******BMI ----
summary(HCAP_2016$r13bmi)
sd(HCAP_2016$r13bmi, na.rm = TRUE)

#---- ******IADL ----
summary(HCAP_2016$r13iadla)
sd(HRS_2016$r13iadla, na.rm = TRUE)

#---- ******stroke history ----
table(HCAP_2016$r13stroke, useNA = "ifany")
summary(HCAP_2016$r13stroke)

#---- ******serial 7s ----
#had to take from HRS because not sure if item 12 on MMSE is serial7s or WORLD backwards
#**ask Jen Manly
summary(HCAP_2016$r13ser7)
sd(HCAP_2016$r13ser7, na.rm = TRUE)

#---- ******immediate word recall ----
summary(HCAP_2016$H1RWLIMMSCORE)
sd(HCAP_2016$H1RWLIMMSCORE, na.rm = TRUE)

#---- ******delayed word recall ----
summary(HCAP_2016$H1RWLDELSCORE)
sd(HCAP_2016$H1RWLDELSCORE, na.rm = TRUE)

#---- ******normed MMSE----
summary(HCAP_2016$H1RMSESCORE_norm)
sd(HCAP_2016$H1RMSESCORE_norm, na.rm = TRUE)

#---- ******word list recognition (yes)----
summary(HCAP_2016$H1RWLRECYSCORE)
sd(HCAP_2016$H1RWLRECYSCORE, na.rm = TRUE)

#---- ******story recall (immediate)----
summary(HCAP_2016$story_sum_immediate)
sd(HCAP_2016$story_sum_immediate, na.rm = TRUE)

#---- ******proxy cognition----
summary(HCAP_2016$H1IIQSCORE)
sd(HCAP_2016$H1IIQSCORE, na.rm = TRUE)
