#---- Package Loading and Options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "future.apply", "readr", "haven",
       "lubridate", "gmodels", "tidyselect")

options(scipen = 999)

#---- Read in the HRS data (wave 13) ----
#list of variables to read in from HRS data:
# ID, Age, Sex/Gender, Race/Ethnicity, Years of education

vars = c("HHIDPN", "R13AGEY_B", "RAGENDER", "RARACEM", "RAHISPAN", "RAEDYRS")

HRS_data <- read_sas(paste0("/Users/CrystalShaw/Box/NIA_F31_April2020/Data/",
                            "HRS/HRS RAND/randhrs1992_2016v1_SAS_data/",
                            "randhrs1992_2016v1.sas7bdat"),
                     n_max = Inf, col_select = all_of(vars)) %>%
  mutate_at("HHIDPN", as.character)


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

#---- Read in NHATS data (wave 1 and 5 (baseline data), 6) ----
#list of variables to read in from HRS data:
# ID, Age, Sex/Gender, Race/Ethnicity, Years of education
w1_vars = c("spid", "el1higstschl")
w5_vars = c("spid", "el5higstschl")
w6_vars = c("spid", "r6d2intvrage", "r5dgender", "rl5dracehisp")

NHATS_w1 <- read_sas(paste0("/Users/CrystalShaw/Box/NIA_F31_April2020/Data/",
                            "NHATS/Round 1/NHATS_ROUND_1_SP_File.sas7bdat"),
                     n_max = Inf, col_select = all_of(w1_vars)) %>%
  mutate_at("spid", as.character)

NHATS_w5 <- read_sas(paste0("/Users/CrystalShaw/Box/NIA_F31_April2020/Data/",
                            "NHATS/Round 5/NHATS_ROUND_5_SP_File_V2.sas7bdat"),
                     n_max = Inf, col_select = all_of(w5_vars)) %>%
  mutate_at("spid", as.character)

NHATS_w6 <- read_sas(paste0("/Users/CrystalShaw/Box/NIA_F31_April2020/Data/",
                            "NHATS/Round 6/NHATS_ROUND_6_SP_File_V2.sas7bdat"),
                     n_max = Inf, col_select = all_of(w6_vars)) %>%
  mutate_at("spid", as.character)

NHATS_clean <- left_join(NHATS_w6, NHATS_w5, by = "spid") %>%
  left_join(., NHATS_w1, by = "spid") %>%
  #Remove those who didn't take the survey based on age variable
  filter(!is.na(r6d2intvrage) & r6d2intvrage != -1) %>%
  #Remove those who don't know/refuse race/ethnicity question
  filter(!(rl5dracehisp %in% c(5, 6))) %>%
  #Carry forward education data from wave 1
  mutate("elhigstschl" = case_when(el5higstschl == -1 ~ el1higstschl,
                                   TRUE ~ el5higstschl)) %>%
  #Remove those not in the survey based on education variable
  filter(elhigstschl != -1) %>%
  #My education categories
  mutate("my_edu_cat" = case_when(elhigstschl == 1 ~ "No school completed",
                                  elhigstschl == 2 ~ "1-8",
                                  elhigstschl %in% c(3) ~ "Some high school",
                                  elhigstschl == 4 ~ "High school diploma",
                                  elhigstschl %in% c(5, 6, 7) ~ "Some college",
                                  elhigstschl %in% c(8, 9) ~
                                    "Bachelor's degree or higher",
                                  elhigstschl %in% c(-9, -8, -7) ~ "missing"))


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
