#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr", "NormPsy")

options(scipen = 999)

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
hrs_vars = c("hhid", "pn", "rahhidpn", "raracem", "rahispan", "r13mstat",
             paste0("h", hrs_waves, "cpl"), paste0("r", hrs_waves, "agey_b"),
             paste0("r", hrs_waves, "bmi"), paste0("r", hrs_waves, "iadla"),
             paste0("r", hrs_waves, "stroke"), paste0("r", hrs_waves, "ser7"),
             paste0("r", hrs_waves, "imrc"), paste0("r", hrs_waves, "dlrc"))

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

wave_updated_vars <- c("cpl", "bmi", "iadla", "stroke")

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
