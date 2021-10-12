#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "haven")

options(scipen = 999)

#---- source scripts ----
#custom read da dct function
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- read in data ----
#---- **HRS tracker ----
HRS_tracker_data_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS_tracker/trk2018_3/", 
         "TRK2018TR_R.da")
HRS_tracker_dict_path <- 
  paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS_tracker/trk2018_3/", 
         "TRK2018TR_R.dCT")

#filter to those who completed Wave A assessment
HRS_tracker <- read_da_dct(HRS_tracker_data_path, HRS_tracker_dict_path, 
                           HHIDPN = "TRUE") %>% 
  #Participated in 2016 HRS
  filter(PIWTYPE == 1) %>% filter(PAGE >= 65)


#---- **ADAMS ----


#---- **HRS ----
#wave 13 = 2016 wave
#list of variables to read in from HRS data:
# ID, Age, Sex/Gender, Race/Ethnicity

hrs_vars = c("hhidpn", "r13agey_b", "raracem", "rahispan", "") 

HRS <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/HRS/", 
                       "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                col_select = all_of(hrs_vars)) %>% 
  mutate_at("hhidpn", as.character)




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
