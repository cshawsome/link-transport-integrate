#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven", "stringr")

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

#filter to those who completed Wave A assessment
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
hrs_vars = c("rahhidpn", "raracem", "rahispan",
             paste0("h", hrs_waves, "cpl"), paste0("r", hrs_waves, "agey_b"),
             paste0("r", hrs_waves, "bmi"), paste0("r", hrs_waves, "iadla"),
             paste0("r", hrs_waves, "stroke"), paste0("r", hrs_waves, "ser7"),
             paste0("r", hrs_waves, "imrc"), paste0("r", hrs_waves, "dlrc"))

HRS <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS/", 
                       "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                col_select = all_of(hrs_vars)) 
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

wave_updated_vars <- c("cpl", "bmi", "iadla", "stroke", "imrc", "dlrc", "ser7")

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

#race/ethnicity
ADAMS_A %<>% 
  mutate("race_eth" = case_when(rahispan == 1 ~ "Hispanic", 
                                rahispan == 0 & raracem == 1 ~ "White", 
                                rahispan == 0 & raracem == 2 ~ "Black", 
                                rahispan == 0 & raracem == 3 ~ "Other"))

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
summary(ADAMS_A$Aser7)
sd(ADAMS_A$Aser7, na.rm = TRUE)

#---- ******immediate word recall ----
summary(ADAMS_A$Aimrc)
sd(ADAMS_A$Aimrc, na.rm = TRUE)

#---- ******delayed word recall ----
summary(ADAMS_A$Adlrc)
sd(ADAMS_A$Adlrc, na.rm = TRUE)

#---- **HCAP ----


#---- data cleaning ----
HRS_clean <- HRS_data %>%
  #Filter out those younger than 50
  filter(R13AGEY_B >= 50) %>%
  #Remove those missing race/ethnicity data
  filter(!is.na(RARACEM)) %>%
  filter(!is.na(RAHISPAN)) %>%
  #My age categories
  mutate_at("R13AGEY_B", floor) %>%
  mutate("my_age_cat" = case_when(R13AGEY_B %in% seq(50, 54) ~ "50-54",
                                  R13AGEY_B %in% seq(55, 59) ~ "55-59",
                                  R13AGEY_B %in% seq(60, 64) ~ "60-64",
                                  R13AGEY_B %in% seq(65, 69) ~ "65-69",
                                  R13AGEY_B %in% seq(70, 74) ~ "70-74",
                                  R13AGEY_B %in% seq(75, 79) ~ "75-79",
                                  R13AGEY_B %in% seq(80, 84) ~ "80-84",
                                  R13AGEY_B %in% seq(85, 89) ~ "85-89",
                                  R13AGEY_B >= 90 ~ "90+")) %>%
  #My education categories
  mutate("my_edu_cat" = case_when(RAEDYRS == 0 ~ "No school completed",
                                  RAEDYRS %in% seq(1, 8) ~ "1-8",
                                  RAEDYRS %in% seq(9, 11) ~ "Some high school",
                                  RAEDYRS == 12 ~ "High school diploma",
                                  RAEDYRS %in% seq(13, 15) ~ "Some college",
                                  RAEDYRS >= 16 ~ "Bachelor's degree or higher",
                                  TRUE ~ "missing"))


#---- Read in HCAP participant ids (wave 1) ----
# Set path to the data file "*.da"
data_path <- paste0("/Users/CrystalShaw/Box/NIA_F31_April2020/Data/HCAP/",
                    "HC16/HC16da/HC16HP_R.da")

# Set path to the dictionary file "*.dct"
dict_path <- paste0("/Users/CrystalShaw/Box/NIA_F31_April2020/Data/HCAP/",
                    "HC16/HC16sta/HC16HP_R.dct")

# Read the dictionary file
df_dict <- read.table(dict_path, skip = 2, fill = TRUE,
                      stringsAsFactors = FALSE)

#Set column names for dictionary dataframe
colnames(df_dict) <- c("col.num", "col.type", "col.name", "col.width",
                       "col.lbl")

#Remove last row which only contains a closing}
df_dict <- df_dict[-nrow(df_dict), ]

#Extract numeric value from column width field
df_dict$col.width <- as.integer(sapply(df_dict$col.width, gsub,
                                       pattern = "[^0-9\\.]",
                                       replacement = ""))

#Convert column types to format to be used with read_fwf function
df_dict$col.type <-
  sapply(df_dict$col.type,
         function(x) ifelse(x %in% c("int","byte","long"), "i",
                            ifelse(x == "float", "n",
                                   ifelse(x == "double", "d", "c"))))

#Read the data file into a dataframe
HCAP <- read_fwf(file = data_path,
                 fwf_widths(widths = df_dict$col.width,
                            col_names = df_dict$col.name),
                 col_types = paste(df_dict$col.type, collapse = ""))

# Add column labels to headers
attributes(HCAP)$variable.labels <- df_dict$col.lbl

HCAP %<>% as.data.frame() %>%
  unite("HHIDPN", c("HHID", "PN"), sep = "")
HCAP$HHIDPN = str_remove(HCAP$HHIDPN, "^0+")

#Join HCAP with HRS data
HRS_HCAP <- inner_join(HRS_clean, HCAP, by = "HHIDPN")

#---- Filling in F31 Research Strategy Table 2 ----
nrow(HRS_clean)
table(HRS_clean$my_age_cat)
table(HRS_clean$my_age_cat)/nrow(HRS_clean)
# 1 = Male; 2 = Female
table(HRS_clean$RAGENDER)
table(HRS_clean$RAGENDER)/nrow(HRS_clean)
# Race: 1 = White; 2 = Black; 3 = Other
CrossTable(HRS_clean$RARACEM, HRS_clean$RAHISPAN)
table(HRS_clean$my_edu_cat)
table(HRS_clean$my_edu_cat)/nrow(HRS_clean)

nrow(NHATS_clean)
table(NHATS_clean$r6d2intvrage)
table(NHATS_clean$r6d2intvrage)/nrow(NHATS_clean)
# 1 = Male; 2 = Female
table(NHATS_clean$r5dgender)
table(NHATS_clean$r5dgender)/nrow(NHATS_clean)
# 1 = Non-hispanic White; 2 = Non-hispanic Black; 3 = Non-hispanic Other;
# 4 = Hispanic
table(NHATS_clean$rl5dracehisp)
table(NHATS_clean$rl5dracehisp)/nrow(NHATS_clean)
table(NHATS_clean$my_edu_cat)
table(NHATS_clean$my_edu_cat)/nrow(NHATS_clean)

nrow(HRS_HCAP)
table(HRS_HCAP$my_age_cat)
table(HRS_HCAP$my_age_cat)/nrow(HRS_HCAP)
# 1 = Male; 2 = Female
table(HRS_HCAP$RAGENDER)
table(HRS_HCAP$RAGENDER)/nrow(HRS_HCAP)
# Race: 1 = White; 2 = Black; 3 = Other
CrossTable(HRS_HCAP$RARACEM, HRS_HCAP$RAHISPAN)
table(HRS_HCAP$my_edu_cat)
table(HRS_HCAP$my_edu_cat)/nrow(HRS_HCAP)
