#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "tidyr")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **clean HCAP ----
HCAP <- read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_clean.csv")) %>% 
  dplyr::select(-one_of("H1RMSESCORE", "H1RWLIMM1SCORE", "H1RWLIMM2SCORE", 
                        "H1RWLIMM3SCORE", "H1RBMIMMSCORE", "H1RLMIMMSCORE", 
                        "H1RBMDELSCORE", "H1RLMDELSCORE", "HCAP_SELECT", 
                        "GENDER", "RACE", "HISPANIC", "PIWTYPE", "PJ005M1", 
                        "r13drinkd", "r13drinkn", "r13pstmem", "GENDER_label", 
                        "RACE_label", "RACE_White", "RACE_Black", "RACE_Other", 
                        "HISPANIC_indicator", "ETHNIC_label", "PJ005M1_label", 
                        "PJ005M1_collapsed_label", "r13drinks_per_week", 
                        "r13drink_cat", "r13drink_cat_label"))

#---- **imputation matrix ----
hotdeck_vars_mat <- 
  read_csv(paste0(path_to_box, "analyses/HCAP/hotdeck_impute_mat.csv")) %>% 
  column_to_rownames("var_names")

#---- **summarize missingness ----
#double check that all of these are in the rownames of the imputation matrix
#do 10 imputations because max missingness is 9.6%
colMeans(is.na(HCAP))[which(colMeans(is.na(HCAP)) > 0)]

#---- source functions ----
source(here::here("HCAP", "functions", "hotdeck_function.R"))

#---- impute: derive variable bins ----
HCAP %<>% 
  mutate(
    #---- **age ----
    "age_cat" = cut(age, breaks = c(70, 85, max(HCAP$age)), 
                    include.lowest = TRUE, right = FALSE), 
    #---- **education ----
    "edyrs_cat" = cut(edyrs, 
                      breaks = c(0, 11, 12, max(HCAP$edyrs)),
                      include.lowest = TRUE, right = TRUE),
    #---- **bmi ----
    "bmi_cat" = cut_number(bmi, n = 5),
    #---- **immediate word recall ----
    "immrc_cat" = cut(immrc, breaks = c(0, 6, 8, max(HCAP$immrc, na.rm = TRUE)), 
                      include.lowest = TRUE, right = FALSE),
    #---- **delayed word recall ----
    "delrc_cat" = cut(delrc, breaks = c(0, 5, 7, max(HCAP$delrc, na.rm = TRUE)), 
                      include.lowest = TRUE, right = FALSE),
    #---- **serial 7s ----
    "ser7_cat" = cut(ser7, breaks = c(0, 4, max(HCAP$ser7, na.rm = TRUE)), 
                     include.lowest = TRUE, right = TRUE), 
    #---- **hrs cognition ----
    "hrs_cog_HCAP_cat" = cut_number(hrs_cog, n = 5),
    "hrs_cog_superpop_cat" = cut_number(hrs_cog, n = 4),
    #---- **mmse_norm ----
    "mmse_norm_cat" = cut_number(mmse_norm, n = 4),
    #---- **animal naming ----
    "animal_naming_cat" = cut_number(animal_naming, n = 5),
    #---- **word recall yes ----
    "wrc_yes_cat" = cut(wrc_yes, 
                        breaks = c(0, 8, max(HCAP$wrc_yes, na.rm = TRUE)), 
                        include.lowest = TRUE, right = TRUE),
    #---- **word recall no ----
    "wrc_no_cat" = cut(wrc_no, 
                       breaks = c(0, 9, max(HCAP$wrc_no, na.rm = TRUE)), 
                       include.lowest = TRUE, right = TRUE), 
    #---- **immediate story recall ----
    "imm_story_cat" = cut_number(imm_story, n = 4), 
    #---- **delayed story recall ----
    "del_story_cat" = cut_number(del_story, n = 5),
    #---- **immediate constructional praxis ----
    "imm_cp_cat" = cut(imm_cp, 
                       breaks = c(0, 7, 10, max(HCAP$imm_cp, na.rm = TRUE)), 
                       include.lowest = TRUE, right = TRUE),
    #---- **delayed constructional praxis ----
    "del_cp_cat" = cut(del_cp, 
                       breaks = c(0, 4, 6, 8, max(HCAP$imm_cp, na.rm = TRUE)), 
                       include.lowest = TRUE, right = TRUE),
    #---- **trails A ----
    "trailsA_cat" = cut_number(trailsA, n = 5))

# #Sanity check bins
# check_vars <- c("age_cat", "edyrs_cat", "immrc_cat", "delrc_cat", "ser7_cat",
#                 "hrs_cog_HCAP_cat", "hrs_cog_superpop_cat", "mmse_norm_cat",
#                 "wrc_yes_cat", "wrc_no_cat", "animal_naming_cat",
#                 "imm_story_cat", "del_story_cat", "imm_cp_cat", "del_cp_cat",
#                 "trailsA_cat")
# 
# for(var in check_vars){
#   print(table(HCAP[, var]))
# }

#---- impute: MI hotdecking ----
set.seed(20221116)
m = 10
HCAP_impute_list <- list()

for(i in 1:m){
  HCAP_impute_list[[i]] <- 
    hotdeck(dataset_to_impute = HCAP, hotdeck_dataset = HCAP, 
            imputation_mat = hotdeck_vars_mat, 
            binary_vars = c("moderate_drinking", "heavy_drinking", 
                            "not_working", "retired", "smoken", "hibpe", 
                            "diabe", "hearte", "stroke", "subj_cog_better", 
                            "subj_cog_worse", "bwc20", "scissor", "cactus", 
                            "pres"), 
            dataset_specific_vars = c("hrs_cog"), 
            dataset_specific_label = "HCAP")
}

# #Sanity check-- do this on one element of the list at a time
# colMeans(is.na(HCAP))[which(colMeans(is.na(HCAP)) > 0)]
# #pool size
# apply(HCAP[, str_detect(colnames(HCAP), "pool")], 2, 
#       function(x) min(x, na.rm = TRUE))

#---- clean dummy vars ----
clean_dummy_vars <- function(data){
  #---- clean: subjective cog decline vars ----
  data %<>% mutate(subj_cog_count = subj_cog_better + subj_cog_worse)
  
  # #check counts
  # table(data$subj_cog_count)
  
  #---- **subjective cog same ----
  data %<>% mutate("subj_cog_same" = ifelse(subj_cog_count == 0, 1, 0))
  
  #---- **subjective cog better/worse ----
  fix_these <- which(data$subj_cog_count > 1)
  
  for(index in fix_these){
    data[index, sample(c("subj_cog_worse", "subj_cog_better"), size = 1)] <- 0
  }
  
  #Sanity check-- should only have sums equal to 1
  data %<>%
    mutate(subj_cog_count = subj_cog_better + subj_cog_worse + subj_cog_same)

  #table(data$subj_cog_count)
  
  #---- clean: drinking vars ----
  data %<>% mutate(drinking_count = moderate_drinking + heavy_drinking)
  
  # #check counts
  # table(data$drinking_count)
  
  #---- **no drinking ----
  data %<>% mutate("no_drinking" = ifelse(drinking_count == 0, 1, 0))
  
  #---- **moderate/heavy drinking ----
  fix_these <- which(data$drinking_count > 1)
  
  for(index in fix_these){
    data[index, sample(c("moderate_drinking", "heavy_drinking"), size = 1)] <- 0
  }
  
  #Sanity check-- should only have sums equal to 1
  data %<>%
    mutate(drinking_count = no_drinking + moderate_drinking + heavy_drinking)

  #table(data$drinking_count)
  
  #---- clean: employment vars ----
  data %<>% mutate(employment_count = not_working + retired)
  
  # #check counts
  # table(data$employment_count)
  
  #---- **working ----
  data %<>% mutate("working" = ifelse(employment_count == 0, 1, 0))
  
  #---- **not working/retired ----
  fix_these <- which(data$employment_count > 1)
  
  for(index in fix_these){
    data[index, sample(c("not_working", "retired"), size = 1)] <- 0
  }
  
  #Sanity check-- should only have sums equal to 1
  data %<>%
    mutate(employment_count = working + not_working + retired)

  #table(data$employment_count)
  
  return(data)
}

HCAP_impute_list %<>% lapply(., clean_dummy_vars)

# #Sanity check-- should all be 1
# for(var in c("subj_cog_count", "drinking_count", "employment_count")){
#   print(lapply(HCAP_impute_list, function(x) table(x[, var])))
# }

#---- standardize continuous vars ----
standardize_vars <- c("animal_naming", "delrc", "wrc_yes", "wrc_no", "imm_cp", 
                      "del_cp", "trailsA", "age", "edyrs", "bmi", "ser7", 
                      "hrs_cog", "adl", "iadl", "mmse_norm", "immrc", "imm_story", 
                      "del_story")

Z_score <- function(data, vars){
  subset <- data %>% dplyr::select(all_of(vars)) %>% 
    mutate_all(scale) %>% set_colnames(paste0(all_of(vars), "_Z"))
  
  data %<>% cbind(., subset)
  
  return(data)
}

HCAP_impute_list <- lapply(HCAP_impute_list, Z_score, standardize_vars)


#---- save dataset ----
HCAP_impute_list %>% 
  saveRDS(paste0(path_to_box, "analyses/HCAP/HCAP_MI_hotdeck"))
