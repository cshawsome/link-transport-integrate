#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "broom", "haven", "labelled", "wesanderson")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- color palette ----
green = "#66a182"

#---- import data ----
#---- **ADAMS neuropsych ----
#ADAMS
for(wave in c("a", "b", "c", "d")){
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

#---- **HCAP neuropsych ----
HCAP_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                         "HC16da/HC16HP_R.da")
HCAP_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                         "HC16sta/HC16HP_R.dct")

HCAP <- read_da_dct(HCAP_data_path, HCAP_dict_path, HHIDPN = "TRUE")
cog_test_labels <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                   "cog_test_meaningful_labels.csv"))

#---- **RAND ----
rand_waves <- seq(5, 13, by = 1)
rand_variables <- c("hhidpn", paste0("r", rand_waves, "cogtot"), 
                    paste0("r", rand_waves, "iadla"))

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- **ADAMS sociodemographic ----
ADAMS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                  "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                  "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") 

#---- **HCAP sociodemographics ----
HRS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                "tracker/trk2018_3/TRK2018TR_R.da")
HRS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                "tracker/trk2018_3/TRK2018TR_R.dct")

HRS_tracker <- read_da_dct(HRS_tracker_data_path, HRS_tracker_dict_path, 
                           HHIDPN = "TRUE") 

#---- join data ----
HCAP <- left_join(HCAP, HRS_tracker, by = "HHIDPN")
ADAMS <- left_join(ADAMS_neuropsych, ADAMS_tracker, by = "HHIDPN")

#---- select variables ----
HCAP_vars <- c("HHIDPN", "GENDER", "HISPANIC", "RACE", "SCHLYRS",
               #for age calculation
               "BIRTHYR", "H1RIWYEAR",
               #MMSE
               "H1RMSESCORE")

ADAMS_vars <- c("HHIDPN", "GENDER", "RACE", "HISPANIC", "EDYRS",
                #Interview year, MMSE
                paste0(c("A", "B", "C", "D"), "YEAR"), "NMSETOT", 
                "AAGE")

HCAP_subset <- HCAP %>% dplyr::select(all_of(HCAP_vars))
ADAMS_subset <- ADAMS %>% dplyr::select(contains(ADAMS_vars))

# #Variable check
# colSums(is.na(HCAP_subset))

#Drop people without an interview year
HCAP_subset %<>% filter(!is.na(H1RIWYEAR))
#Set people with missing MMSE data (code = 97) to NA in ADAMS
ADAMS_subset %<>% 
  mutate_at(.vars = paste0(c("A", "B", "C", "D"), "NMSETOT"), 
            function(x) ifelse(x > 30, NA, x))

# #Sanity check
# ADAMS_subset %>% dplyr::select(contains("NMSETOT")) %>% 
#   apply(2, function(x) max(x, na.rm = TRUE))

# table(ADAMS_subset$AYEAR, useNA = "ifany")
# table(ADAMS_subset$BYEAR, useNA = "ifany")
# table(ADAMS_subset$CYEAR, useNA = "ifany")
# table(ADAMS_subset$DYEAR, useNA = "ifany")

#---- format data ----
#---- **ADAMS interview year ----
ADAMS_subset %<>% mutate_at("BYEAR", function(x) ifelse(x == 9997, NA, x))

# #Sanity check
# table(ADAMS_subset$BYEAR, useNA = "ifany")

#---- **age ----
HCAP_subset %<>% mutate("HCAP_age" = H1RIWYEAR - BIRTHYR)

# #Sanity check
# ggplot(data = HCAP_subset) + 
#   geom_histogram(aes(HCAP_age), color = green, fill = green) + 
#   theme_minimal() + xlab("Age") + theme(text = element_text(size = 14))
# ggsave(filename = "HCAP_age.jpeg", plot = last_plot(), device = "jpeg", 
#        path = paste0("/Users/CrystalShaw/Box/Dissertation/", 
#                      "preliminary_analyses/", 
#                      "HCAP_synthetic/figures/"), width = 7, height = 5, 
#        units = "in")
# 
# summary(HCAP_subset$HCAP_age)

# #Sanity check
# colnames(HCAP_subset)

#---- **race/ethnicity ----
# #Variable check
# #0 = "Not Obtained"; 1-3 = diff categories Hispanic; 5 = "Non-Hispanic"
# table(HCAP_subset$HISPANIC) 
# #0 = "Not Obtained"; 1 = "White"; 2 = "Black"; 7 = "Other"
# table(HCAP_subset$RACE)
# table(HCAP_subset$RACE, HCAP_subset$HISPANIC)

#Drop missing race
HCAP_subset %<>% filter(RACE != 0)
ADAMS_subset %<>% filter(RACE != 0)

HCAP_subset %<>% 
  mutate("race_ethnic_cat" = 
           case_when(RACE == 1 & HISPANIC == 5 ~ "Non-hispanic White", 
                     HISPANIC %in% c(1, 2, 3) ~ "Hispanic", 
                     RACE == 2 & (HISPANIC == 5 | HISPANIC == 0) ~ 
                       "Non-hispanic Black", 
                     RACE == 7 & HISPANIC == 5 ~ "Other"))

ADAMS_subset %<>% 
  mutate("race_ethnic_cat" = 
           case_when(RACE == 1 & HISPANIC == 5 ~ "Non-hispanic White", 
                     HISPANIC %in% c(1, 2, 3) ~ "Hispanic", 
                     RACE == 2 & (HISPANIC == 5 | HISPANIC == 0) ~ 
                       "Non-hispanic Black", 
                     RACE == 7 & HISPANIC == 5 ~ "Other"))

# #Sanity check
# table(HCAP_subset$race_ethnic_cat, useNA = "ifany")
# table(HCAP_subset$race_ethnic_cat, useNA = "ifany")/nrow(HCAP_subset)
# 
# table(ADAMS_subset$race_ethnic_cat, useNA = "ifany")
# table(ADAMS_subset$race_ethnic_cat, useNA = "ifany")/nrow(ADAMS_subset)

#Drop original variables
HCAP_subset %<>% dplyr::select(-c("RACE", "HISPANIC"))
ADAMS_subset %<>% dplyr::select(-c("RACE", "HISPANIC"))

# #Sanity check
# colnames(HCAP_subset)
# colnames(ADAMS_subset)

#---- **sex/gender ----
# #Variable check
# #1 = "Male"; 2 = "Female"
# table(HCAP_subset$GENDER)

HCAP_subset %<>% mutate("female" = ifelse(GENDER == 2, 1, 0))
ADAMS_subset %<>% mutate("female" = ifelse(GENDER == 2, 1, 0))

# #Sanity check
# table(HCAP_subset$female, useNA = "ifany")
# table(HCAP_subset$female, useNA = "ifany")/nrow(HCAP_subset)
# table(ADAMS_subset$female, useNA = "ifany")
# table(ADAMS_subset$female, useNA = "ifany")/nrow(ADAMS_subset)

#Drop original variable
HCAP_subset %<>% dplyr::select(-c("GENDER"))
ADAMS_subset %<>% dplyr::select(-c("GENDER"))

# #Sanity check
# colnames(HCAP_subset)
# colnames(ADAMS_subset)

#---- education distributions ----
# #Variable Check
# #0-17 = years of education; 99 = "Not Ascertained"
# table(HCAP_subset$SCHLYRS)
# table(ADAMS_subset$EDYRS, useNA = "ifany")

#Remove missing education
HCAP_subset %<>% filter(SCHLYRS != 99)

ggplot(data = HCAP_subset) + 
  geom_histogram(aes(SCHLYRS), color = green, fill = green) +
  theme_minimal() + xlab("Years of Education") + 
  theme(text = element_text(size = 14))

ggsave(filename = "HCAP_edyrs.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

# #Variable check
# summary(HCAP_subset$SCHLYRS)

#---- neuropsych distributions ----
#HCAP
ggplot(data = HCAP_subset) + 
  geom_histogram(aes(H1RMSESCORE), color = green, fill = green) +
  theme_minimal() + xlab("MMSE Score") + theme(text = element_text(size = 14))

ggsave(filename = "HCAP_MMSE_score.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

#ADAMS
ggplot(data = ADAMS_subset) + 
  geom_histogram(aes(ANMSETOT), color = green, fill = green) +
  theme_minimal() + xlab("MMSE Score") + theme(text = element_text(size = 14))

ggsave(filename = "ADAMS_MMSE_score.jpeg", plot = last_plot(), device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

# #Variable check
# summary(HCAP_subset$H1RMSESCORE)
summary(ADAMS_subset$ANMSETOT)

#---- modeling HCAP data ----
#Formatting race/ethnicity factors
HCAP_subset %<>% mutate_at("race_ethnic_cat", as.factor) %>% 
  mutate_at("race_ethnic_cat", function(x) fct_relevel(x, "Non-hispanic White"))

#---- **closest RAND total cognition measure ----
#join RAND and HCAP
HCAP_subset %<>% left_join(., RAND, by = "HHIDPN")

#closest RAND year
HCAP_subset[, "RANDYEAR"] <- 
  ifelse(HCAP_subset[, "H1RIWYEAR"] %% 2 == 0, 
         HCAP_subset[, "H1RIWYEAR"], 
         HCAP_subset[, "H1RIWYEAR"] - 1)

HCAP_subset[, "RANDWAVE"] <- 
  ((HCAP_subset[, "RANDYEAR"] - 2000)/2) + 5

# #Sanity check
# table(HCAP_subset$H1RIWYEAR, HCAP_subset$RANDYEAR, useNA = "ifany")
# table(HCAP_subset$RANDYEAR, HCAP_subset$RANDWAVE, useNA = "ifany")

#closest RAND measure
for(i in 1:nrow(HCAP_subset)){
  if(!is.na(HCAP_subset[i, "H1RMSESCORE"])){
    RAND_wave <- HCAP_subset[i, "RANDWAVE"]
    RAND_cog_val <- HCAP_subset[i, paste0("r", RAND_wave, "cogtot")]
    if(!is.na(RAND_cog_val)){
      HCAP_subset[i, "RANDcogtot"] <- RAND_cog_val
    } else{
      past_cog_vals <- 
        HCAP_subset[i, paste0("r", seq(5, (RAND_wave - 1)), "cogtot")]
      nearest <- max(which(!is.na(past_cog_vals)))
      if(is.finite(nearest)){
        HCAP_subset[i, "RANDcogtot"] <- past_cog_vals[nearest]
      } 
    }
  }
}

# #Sanity check-- we only have 1,474 people with any RAND cognitive data
# View(HCAP_subset %>% dplyr::select(c("H1RMSESCORE", "RANDcogtot", 
#                                      contains("cogtot"))))
# sum(!is.na(HCAP_subset$RANDcogtot))

#---- HCAP people in model ----
HCAP_RANDcomplete <- HCAP_subset %>% filter(!is.na(RANDcogtot))

hist(HCAP_RANDcomplete$H1RMSESCORE)
summary(HCAP_RANDcomplete$H1RMSESCORE)

#---- **main effects models ----
MMSE_main_effects_no_cog <- lm(H1RMSESCORE ~ HCAP_age + female + 
                                 race_ethnic_cat +
                                 SCHLYRS, data = HCAP_subset)

MMSE_main_effects <- lm(H1RMSESCORE ~ HCAP_age + female + 
                          race_ethnic_cat + RANDcogtot +
                          SCHLYRS, data = HCAP_subset)

summary(MMSE_main_effects)
plot(MMSE_main_effects$residuals)
hist(MMSE_main_effects$fitted.values)
summary(MMSE_main_effects$fitted.values)

summary(MMSE_main_effects_no_cog)
plot(MMSE_main_effects_no_cog$residuals)
hist(MMSE_main_effects_no_cog$fitted.values)
summary(MMSE_main_effects_no_cog$fitted.values)

#---- **2-way interactions ----
MMSE_2way_intx_no_cog <- lm(H1RMSESCORE ~ HCAP_age*female + 
                              HCAP_age*race_ethnic_cat + 
                              female*race_ethnic_cat + female*SCHLYRS + 
                              race_ethnic_cat*SCHLYRS, 
                            data = HCAP_subset)

MMSE_2way_intx <- lm(H1RMSESCORE ~ HCAP_age*female + 
                       HCAP_age*race_ethnic_cat + 
                       female*race_ethnic_cat + female*SCHLYRS + 
                       race_ethnic_cat*SCHLYRS + RANDcogtot*female + 
                       RANDcogtot*race_ethnic_cat + RANDcogtot*SCHLYRS + 
                       RANDcogtot*HCAP_age, 
                     data = HCAP_subset)

summary(MMSE_2way_intx)
tidy(MMSE_2way_intx, conf.int = TRUE)
plot(MMSE_2way_intx$residuals)
hist(MMSE_2way_intx$fitted.values)
summary(MMSE_2way_intx$fitted.values)

summary(MMSE_2way_intx_no_cog)
tidy(MMSE_2way_intx_no_cog, conf.int = TRUE)
plot(MMSE_2way_intx_no_cog$residuals)
hist(MMSE_2way_intx_no_cog$fitted.values)
summary(MMSE_2way_intx_no_cog$fitted.values)

#---- ****performance plots ----
#MMSE Dists
plot_data <- cbind(MMSE_2way_intx$fitted.values, 
                   HCAP_RANDcomplete$H1RMSESCORE) %>% 
  set_colnames(c("fitted", "observed"))
#plot_data[plot_data > 30] <- 30
plot_data %<>% as.data.frame() %>% 
  pivot_longer(everything(), names_to = "Value", values_to = "MMSE")

ggplot(data = plot_data, 
       aes(x = MMSE, color = factor(Value), 
           fill = factor(Value))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Value Type"), 
         color = guide_legend(title = "Value Type")) 

ggsave(filename = "hists_HCAP_obs_pred_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#Scatterplot
plot_data <- cbind(MMSE_2way_intx$fitted.values, 
                   HCAP_RANDcomplete$H1RMSESCORE) %>% 
  set_colnames(c("fitted", "observed"))
#plot_data[plot_data > 30] <- 30
plot_data %<>% as.data.frame()

ggplot(data = plot_data, 
       aes(x = observed, y = fitted)) + 
  geom_point(color = green) + theme_minimal() + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("HCAP MMSE") + 
  ylab("Fitted Value") + 
  theme(text = element_text(size = 14)) 

ggsave(filename = "person-level_HCAP_obs_pred_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#---- **3-way interactions ----
MMSE_3way_intx_no_cog <- lm(H1RMSESCORE ~ HCAP_age*female*race_ethnic_cat + 
                              SCHLYRS*female*race_ethnic_cat, 
                            data = HCAP_subset)

MMSE_3way_intx <- lm(H1RMSESCORE ~ HCAP_age*female*race_ethnic_cat + 
                       SCHLYRS*female*race_ethnic_cat + 
                       RANDcogtot*HCAP_age*female + 
                       RANDcogtot*female*race_ethnic_cat, 
                     data = HCAP_subset)

summary(MMSE_3way_intx_no_cog)
plot(MMSE_3way_intx_no_cog$residuals)
hist(MMSE_3way_intx_no_cog$fitted.values)
summary(MMSE_3way_intx_no_cog$fitted.values)

summary(MMSE_3way_intx)
plot(MMSE_3way_intx$residuals)
hist(MMSE_3way_intx$fitted.values)
summary(MMSE_3way_intx$fitted.values)

#---- mixture distributions ----
#Try to recover the HCAP shape distribution without taking covariates into 
# consideration... yet
samp_size = 10000

high_score_prop = 0.65; high_score_mean = 28; high_score_SD = 2
med_score_prop = 0.30; med_score_mean = 22; med_score_SD = 3.5
low_score_prop = 0.05; low_score_mean = 15; low_score_SD = 5

high_scores <- rnorm(n = floor(high_score_prop*samp_size), 
                     mean = high_score_mean, sd = high_score_SD)
high_scores[high_scores > 30] <- 30
high_scores[high_scores < 0] <- 0

med_scores <- rnorm(n = floor(med_score_prop*samp_size), 
                    mean = med_score_mean, sd = med_score_SD)
med_scores[med_scores > 30] <- 30
med_scores[med_scores < 0] <- 0

low_scores <- rnorm(n = floor(low_score_prop*samp_size), 
                    mean = low_score_mean, sd = low_score_SD)
low_scores[low_scores > 30] <- 30
low_scores[low_scores < 0] <- 0

HCAP_plot_data <- data.frame("values" = HCAP_subset$H1RMSESCORE, 
                             "data_type" = "HCAP")
synthetic_HCAP_plot_data <- 
  data.frame("values" = c(high_scores, med_scores, low_scores), 
             "data_type" = "Synthetic HCAP") 

plot_data <- rbind(HCAP_plot_data, synthetic_HCAP_plot_data)

#Comparison plot
ggplot(data = plot_data, aes(x = values, color = data_type, fill = data_type)) + 
  geom_histogram(alpha = 0.5, position = "identity") + theme_minimal() + 
  xlab("MMSE Total Scores") + theme(text = element_text(size = 14))

ggsave(filename = "synthetic_HCAP_MMSE_score.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5,
       units = "in")

#---- overlap in mixtures ----
plot_data <- 
  data.frame("score" = 
               c(rnorm(n = floor(high_score_prop*samp_size), 
                       mean = high_score_mean, sd = high_score_SD), 
                 rnorm(n = floor(med_score_prop*samp_size), 
                       mean = med_score_mean, sd = med_score_SD), 
                 rnorm(n = floor(low_score_prop*samp_size), 
                       mean = low_score_mean, sd = low_score_SD)), 
             "group" = c(rep("high", floor(high_score_prop*samp_size)), 
                         rep("med", floor(med_score_prop*samp_size)), 
                         rep("low", floor(low_score_prop*samp_size))))
plot_data[which(plot_data$score > 30), "score"] <- 30
plot_data[which(plot_data$score < 0), "score"] <- 0

ggplot(data = plot_data, aes(x = score, color = group, fill = group)) + 
  geom_histogram(alpha = 0.5, position = "identity") + theme_minimal() + 
  xlab("Scores") + theme(text = element_text(size = 14))

ggsave(filename = "mixture_dists.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#---- HCAP MMSE joint distributions ----
#---- **sex x MMSE ----
plot_data <- HCAP_subset %>% dplyr::select("female", "H1RMSESCORE")

ggplot(data = plot_data, aes(x = H1RMSESCORE, color = factor(female), 
                             fill = factor(female))) + 
  geom_histogram(alpha = 0.5, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "female"), 
         color = guide_legend(title = "female"))

ggsave(filename = "HCAP_sex_by_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#---- **race/ethnicity x MMSE ----
plot_data <- HCAP_subset %>% dplyr::select("race_ethnic_cat", "H1RMSESCORE")

ggplot(data = plot_data, aes(x = H1RMSESCORE, color = race_ethnic_cat, 
                             fill = race_ethnic_cat)) + 
  geom_histogram(alpha = 0.5, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Race/Ethnicity"), 
         color = guide_legend(title = "Race/Ethnicity"))

ggsave(filename = "HCAP_race_eth_by_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#---- **age x MMSE ----
plot_data <- HCAP_subset %>% dplyr::select("HCAP_age", "H1RMSESCORE")

ggplot(data = plot_data, aes(x = HCAP_age, y = H1RMSESCORE)) + 
  geom_point(color = green) + theme_minimal() + 
  xlab("Age") + ylab("MMSE") + theme(text = element_text(size = 14)) 

ggsave(filename = "HCAP_age_by_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#---- **age x MMSE ----
plot_data <- HCAP_subset %>% dplyr::select("SCHLYRS", "H1RMSESCORE")

ggplot(data = plot_data, aes(x = SCHLYRS, y = H1RMSESCORE)) + 
  geom_point(color = green) + theme_minimal() + 
  xlab("Years of Education") + ylab("MMSE") + 
  theme(text = element_text(size = 14)) 

ggsave(filename = "HCAP_edu_by_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#---- ADAMS distributions by dem status ----
#---- **read in demdx data ----
for(wave in c("a", "b", "c", "d")){
  demdx_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "da/", 
                            "ADAMS1", str_to_upper(wave), "D_R.da")
  demdx_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1", wave, "/adams1", wave, "sta/", 
                            "ADAMS1", str_to_upper(wave), "D_R.dct")
  if(wave == "a"){
    ADAMS_demdx <- read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
      dplyr::select("HHIDPN", paste0(str_to_upper(wave), "DFDX1")) %>% 
      set_colnames(c("HHIDPN", "dem_dx")) %>% 
      mutate(new_col = 
               case_when(dem_dx %in% c(1, 2) ~ "Probable/Possible AD", 
                         dem_dx %in% c(3, 4) ~ 
                           "Probable/Possible Vascular Dementia", 
                         dem_dx %in% 
                           c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 21, 28, 29, 
                             30, 33) ~ "Other",
                         dem_dx %in% c(18, 32) ~ "Probable Dementia",
                         dem_dx %in% c(10, 13, 15, 16, 17, 19) ~ "Dementia", 
                         dem_dx %in% c(20, 22) ~ "MCI", 
                         dem_dx == 31 ~ "Normal")) %>% 
      dplyr::select(-c("dem_dx")) %>%
      set_colnames(c("HHIDPN", 
                     paste0(str_to_upper(wave), "dem_dx_cat")))
  } else{
    ADAMS_demdx %<>% 
      left_join(., read_da_dct(demdx_data_path, demdx_dict_path, 
                               HHIDPN = "TRUE") %>% 
                  dplyr::select("HHIDPN", 
                                paste0(str_to_upper(wave), "DFDX1")) %>% 
                  set_colnames(c("HHIDPN", "dem_dx")) %>% 
                  mutate(new_col = 
                           case_when(dem_dx %in% c(1, 2) ~ 
                                       "Probable/Possible AD", 
                                     dem_dx %in% c(3, 4) ~ 
                                       "Probable/Possible Vascular Dementia", 
                                     dem_dx %in% 
                                       c(5, 6, 7, 8, 11, 14, 23, 24, 25, 26, 27, 
                                         21, 28, 29, 30, 33) ~ "Other",
                                     dem_dx %in% c(18, 32) ~ 
                                       "Probable Dementia",
                                     dem_dx %in% c(10, 13, 15, 16, 17, 19) ~ 
                                       "Dementia", 
                                     dem_dx %in% c(20, 22) ~ "MCI", 
                                     dem_dx == 31 ~ "Normal")) %>% 
                  dplyr::select(-c("dem_dx")) %>% 
                  set_colnames(c("HHIDPN", 
                                 paste0(str_to_upper(wave), "dem_dx_cat"))), 
                by = "HHIDPN")
  }
}

# #Sanity check-- need to comment out lines 406 and 430 and add "temp" as second 
# # column name in set_colnames for this check
# table(ADAMS_demdx$temp.x, ADAMS_demdx$Adem_dx_cat, useNA = "ifany")
# table(ADAMS_demdx$temp.y, ADAMS_demdx$Bdem_dx_cat, useNA = "ifany")
# table(ADAMS_demdx$temp.x.x, ADAMS_demdx$Cdem_dx_cat, useNA = "ifany")
# table(ADAMS_demdx$temp.y.y, ADAMS_demdx$Ddem_dx_cat, useNA = "ifany")

#---- **join with neurospych ----
ADAMS_subset %<>% left_join(., ADAMS_demdx, by = "HHIDPN") %>% 
  left_join(., RAND, by = "HHIDPN")

#---- **closest RAND year ----
#Wave Year | HRS Core Data | RAND
# 2000 | G | 5
# 2002 | H | 6
# 2004 | J | 7
# 2006 | K | 8
# 2008 | L | 9

#For the even years, take the year itself; for odd years, take the year before
for(wave in c("A", "B", "C", "D")){
  ADAMS_subset[, paste0(wave, "RANDYEAR")] <- 
    ifelse(ADAMS_subset[, paste0(wave, "YEAR")] %% 2 == 0, 
           ADAMS_subset[, paste0(wave, "YEAR")], 
           ADAMS_subset[, paste0(wave, "YEAR")] - 1)
  
  ADAMS_subset[, paste0(wave, "RANDWAVE")] <- 
    ((ADAMS_subset[, paste0(wave, "RANDYEAR")] - 2000)/2) + 5
}

# #Sanity check
# table(ADAMS_subset$AYEAR, ADAMS_subset$ARANDYEAR, useNA = "ifany")
# table(ADAMS_subset$BYEAR, ADAMS_subset$BRANDYEAR, useNA = "ifany")
# table(ADAMS_subset$CYEAR, ADAMS_subset$CRANDYEAR, useNA = "ifany")
# table(ADAMS_subset$DYEAR, ADAMS_subset$DRANDYEAR, useNA = "ifany")
# 
# table(ADAMS_subset$ARANDYEAR, ADAMS_subset$ARANDWAVE, useNA = "ifany")
# table(ADAMS_subset$BRANDYEAR, ADAMS_subset$BRANDWAVE, useNA = "ifany")
# table(ADAMS_subset$CRANDYEAR, ADAMS_subset$CRANDWAVE, useNA = "ifany")
# table(ADAMS_subset$DRANDYEAR, ADAMS_subset$DRANDWAVE, useNA = "ifany")

#---- **closest RAND cog score ----
for(wave in c("A", "B", "C", "D")){
  for(i in 1:nrow(ADAMS_subset)){
    if(!is.na(ADAMS_subset[i, paste0(wave, "NMSETOT")])){
      RAND_wave <- ADAMS_subset[i, paste0(wave, "RANDWAVE")]
      RAND_cog_val <- ADAMS_subset[i, paste0("r", RAND_wave, "cogtot")]
      if(!is.na(RAND_cog_val)){
        ADAMS_subset[i, paste0(wave, "RANDcogtot")] <- RAND_cog_val
      } else{
        if(RAND_wave != 9){
          future_cog_vals <- 
            ADAMS_subset[i, paste0("r", seq(RAND_wave + 1, 9), "cogtot")]
          nearest <- min(which(!is.na(future_cog_vals)))
          if(is.finite(nearest)){
            ADAMS_subset[i, paste0(wave, "RANDcogtot")] <- 
              future_cog_vals[nearest]
          } 
        }
      }
    }
  }
}

#Zscore the cognitive measures
MMSE_Zscore <- ADAMS_subset %>% dplyr::select(contains("NMSETOT")) %>% 
  mutate_all(function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)) 
MMSE_Zscore %<>% set_colnames(paste0(colnames(MMSE_Zscore), "_Zscore"))

cogtot_Zscore <- ADAMS_subset %>% dplyr::select(contains("RANDcogtot")) %>% 
  mutate_all(function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))  
cogtot_Zscore %<>% set_colnames(paste0(colnames(cogtot_Zscore), "_Zscore"))

#cbind Zscores
ADAMS_subset %<>% cbind(MMSE_Zscore) %>% cbind(cogtot_Zscore)

# #Sanity check
# View(ADAMS_subset %>%
#        dplyr::select(contains(c("ARANDWAVE", "ANMSE", "cogtot"))))
# View(ADAMS_subset %>%
#        dplyr::select(contains(c("ARANDWAVE", "ANMSE", "cogtot"))) %>% 
#        filter(is.na(ARANDcogtot)))
# View(ADAMS_subset %>% 
#        dplyr::select(contains(c("BRANDWAVE", "BNMSE", "cogtot"))))
# View(ADAMS_subset %>% 
#        dplyr::select(contains(c("CRANDWAVE", "CNMSE", "cogtot"))))

# head((ADAMS_subset$ANMSETOT - 
#     mean(ADAMS_subset$ANMSETOT, na.rm = TRUE))/
#   sd(ADAMS_subset$ANMSETOT, na.rm = TRUE))
# head((ADAMS_subset$BNMSETOT - 
#         mean(ADAMS_subset$BNMSETOT, na.rm = TRUE))/
#        sd(ADAMS_subset$BNMSETOT, na.rm = TRUE))
# 
# head(MMSE_Zscore)

#---- **plot: MMSE x dem dx ----
plot_data <- ADAMS_subset %>% 
  dplyr::select(contains(c("NMSETOT", "dem_dx_cat"))) %>%
  pivot_longer(everything(), 
               names_to = c("wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  mutate("dem_dx_new_cat" = 
           case_when(dem_dx_cat %in% 
                       c("Dementia", "Probable Dementia", 
                         "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia", 
                     TRUE ~ dem_dx_cat))

# #Variable check
# table(plot_data$dem_dx_cat, useNA = "ifany")
# table(plot_data$dem_dx_cat, plot_data$dem_dx_new_cat, useNA = "ifany")

ggplot(data = plot_data %>% filter(dem_dx_new_cat != "Other"), 
       aes(x = NMSETOT, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx")) + 
  facet_grid(cols = vars(wave))

ggsave(filename = "ADAMS_MMSE_by_dem.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 15, height = 5, 
       units = "in")

ggplot(data = plot_data, 
       aes(x = NMSETOT, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx")) + 
  facet_grid(cols = vars(wave))

ggsave(filename = "ADAMS_MMSE_by_dem_and_other.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 15, height = 5, 
       units = "in")

#---- **for poster ----
ggplot(data = plot_data %>% filter(wave == "A"), 
       aes(x = NMSETOT, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 1, position = "identity") + theme_minimal() + 
  xlab("Total Cognition") + theme(text = element_text(size = 14)) + 
  scale_color_manual(values = rep("black", 4)) + 
  scale_fill_manual(values = rep("black", 4)) + 
  theme(legend.position = "none") 

ggsave(filename = "ADAMS_MMSE_all.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "presentations/admitted_students_day2021"), 
       width = 7, height = 5, units = "in")

ggplot(data = plot_data %>% filter(wave == "A"), 
       aes(x = NMSETOT, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("Total Cognition") + theme(text = element_text(size = 14)) + 
  scale_color_manual(values = c("#ff0000", "#fdab00", "#00a389", "#28bed9")) + 
  scale_fill_manual(values = c("#ff0000", "#fdab00", "#00a389", "#28bed9")) + 
  theme(legend.position = "none") 

ggsave(filename = "ADAMS_MMSE_by_class.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "presentations/admitted_students_day2021"), 
       width = 7, height = 5, units = "in")

#---- **plot: RAND cog score x dem dx ----
plot_data <- ADAMS_subset %>% 
  dplyr::select(contains(c("RANDcogtot", "dem_dx_cat"))) %>%
  pivot_longer(everything(), 
               names_to = c("wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  mutate("dem_dx_new_cat" = 
           case_when(dem_dx_cat %in% 
                       c("Dementia", "Probable Dementia", 
                         "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia", 
                     TRUE ~ dem_dx_cat))

# #Variable check
# table(plot_data$dem_dx_cat, useNA = "ifany")
# table(plot_data$dem_dx_cat, plot_data$dem_dx_new_cat, useNA = "ifany")

ggplot(data = plot_data %>% filter(dem_dx_new_cat != "Other"), 
       aes(x = RANDcogtot, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("HRS Total Cognition Score") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx")) + 
  facet_grid(cols = vars(wave))

ggsave(filename = "HRS_total_cog_by_dem.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 15, height = 5, 
       units = "in")

ggplot(data = plot_data, 
       aes(x = RANDcogtot, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("HRS Total Cognition Score") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx")) + 
  facet_grid(cols = vars(wave))

ggsave(filename = "HRS_total_cog_by_dem_and_other.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 15, height = 5, 
       units = "in")

#---- **plot: MMSE x closest RAND cog score ----
corA <- round(cor(ADAMS_subset$ANMSETOT_Zscore, ADAMS_subset$ARANDcogtot_Zscore, 
                  use = "complete.obs"), 2)
corB <- round(cor(ADAMS_subset$BNMSETOT_Zscore, ADAMS_subset$BRANDcogtot_Zscore, 
                  use = "complete.obs"), 2)
corC <- round(cor(ADAMS_subset$CNMSETOT_Zscore, ADAMS_subset$CRANDcogtot_Zscore, 
                  use = "complete.obs"), 2)
corD <- round(cor(ADAMS_subset$DNMSETOT_Zscore, ADAMS_subset$DRANDcogtot_Zscore, 
                  use = "complete.obs"), 2)

plot_data <- ADAMS_subset %>% 
  dplyr::select(contains(c("NMSETOT_Zscore", "RANDcogtot_Zscore"))) %>% 
  pivot_longer(everything(),
               names_to = c("wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  mutate("wave_corr" = case_when(wave == "A" ~ paste0("A: rho = ", corA), 
                                 wave == "B" ~ paste0("B: rho = ", corB), 
                                 wave == "C" ~ paste0("C: rho = ", corC), 
                                 wave == "D" ~ paste0("D: rho = ", corD)))

ggplot(data = plot_data, 
       aes(x = NMSETOT_Zscore, y = RANDcogtot_Zscore)) + 
  geom_point(color = green) + theme_minimal() + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("ADAMS MMSE Score (Z score)") + 
  ylab("HRS Total Cognition Score (Z score)") + 
  theme(text = element_text(size = 14)) +
  facet_grid(cols = vars(wave_corr))

ggsave(filename = "ADAMS_MMSE_by_total_cog_Zscore.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 15, height = 5, 
       units = "in")

#---- **out-of-sample performance ----
#It looks like the model with 2-way interactions has the best fit in HCAP

#format ADAMS data
ADAMS_test <- ADAMS_subset %>% 
  dplyr::select("HHIDPN", "EDYRS", "AAGE", "ANMSETOT", "race_ethnic_cat", 
                "female", "ARANDcogtot") %>% 
  set_colnames(c("HHIDPN", "SCHLYRS", "HCAP_age", "ANMSETOT", "race_ethnic_cat", 
                 "female", "RANDcogtot")) %>% na.omit()

ADAMS_predicted <- predict.lm(MMSE_2way_intx, newdata = ADAMS_test)

#Results
#MMSE Dists
plot_data <- cbind(ADAMS_predicted, 
                   ADAMS_test$ANMSETOT) %>% 
  set_colnames(c("fitted", "observed"))
#plot_data[plot_data > 30] <- 30
plot_data %<>% as.data.frame() %>% 
  pivot_longer(everything(), names_to = "Value", values_to = "MMSE")

ggplot(data = plot_data, 
       aes(x = MMSE, color = factor(Value), 
           fill = factor(Value))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Value Type"), 
         color = guide_legend(title = "Value Type")) 

ggsave(filename = "hists_ADAMS_obs_pred_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")

#Scatterplot
plot_data <- cbind(ADAMS_predicted, 
                   ADAMS_test$ANMSETOT) %>% 
  set_colnames(c("fitted", "observed"))
#plot_data[plot_data > 30] <- 30
plot_data %<>% as.data.frame()

ggplot(data = plot_data, 
       aes(x = observed, y = fitted)) + 
  geom_point(color = green) + theme_minimal() + 
  geom_abline(intercept = 0, slope = 1) + 
  xlab("ADAMS MMSE") + 
  ylab("Fitted Value") + 
  theme(text = element_text(size = 14)) 

ggsave(filename = "person-level_ADAMS_obs_pred_MMSE.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "preliminary_analyses/",
                     "HCAP_synthetic/figures/"), width = 7, height = 5, 
       units = "in")










