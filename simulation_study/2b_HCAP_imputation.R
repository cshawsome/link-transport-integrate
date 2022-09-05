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
  read_csv(paste0(path_to_box, "data/HCAP/hotdeck_impute_mat.csv")) %>% 
  column_to_rownames("var_names")

#---- **summarize missingness ----
#double check that all of these are in the rownames of the imputation matrix
colMeans(is.na(HCAP))[which(colMeans(is.na(HCAP)) > 0)]

#---- source functions ----
source(here::here("simulation_study", "functions", "hotdeck_function.R"))

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
    "hrs_cog_cat" = cut_number(hrs_cog, n = 5), 
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
    "trailsA_cat" = cut_number(trailsA, n = 5), 
    #---- **imputed memory ----
    "memimp16_cat" = cut_number(memimp16, n = 4))

#Sanity check bins
check_vars <- c("age_cat", "edyrs_cat", "immrc_cat", "delrc_cat", "ser7_cat", 
                "hrs_cog_cat", "mmse_norm_cat", "wrc_yes_cat", "wrc_no_cat", 
                "animal_naming_cat", "imm_story_cat", "del_story_cat", 
                "imm_cp_cat", "del_cp_cat", "trailsA_cat", "memimp16_cat")

for(var in check_vars){
  print(table(HCAP[, var]))
}

#---- impute: hotdecking ----
set.seed(20220904)
HCAP %<>% 
  hotdeck(dataset_to_impute = ., hotdeck_dataset = HCAP, 
          imputation_mat = hotdeck_vars_mat, 
          binary_vars = c("subj_cog_better", "subj_cog_worse", "bwc20", "scissor", 
                          "cactus", "pres"))

#Sanity check
colMeans(is.na(HCAP))[which(colMeans(is.na(HCAP)) > 0)]

#---- clean: subjective cog decline vars ----
HCAP %<>% mutate(subj_cog_count = subj_cog_better + subj_cog_worse)

#check counts
table(HCAP$subj_cog_count)

#---- **subjective cog same ----
HCAP %<>% mutate("subj_cog_same" = ifelse(subj_cog_count == 0, 1, 0))

#---- **subjective cog better/worse ----
fix_these <- which(HCAP$subj_cog_count > 1)

for(index in fix_these){
  HCAP[index, sample(c("subj_cog_worse", "subj_cog_better"), size = 1)] <- 0
}

# #Sanity check-- should only have sums equal to 1
# HCAP %<>% 
#   mutate(subj_cog_count = subj_cog_better + subj_cog_worse + subj_cog_same)
# 
# table(HCAP$subj_cog_count)

#---- last bin check ----
# #sum should be 2235
# check_vars <- c("age_cat", "edyrs_cat", "immrc_cat", "delrc_cat", "ser7_cat", 
#                 "adl_cat", "iadl_cat", "wrc_yes_cat", "wrc_no_cat", "imm_cp_cat", 
#                 "del_cp_cat")
# 
# for(var in check_vars){
#   print(sum(table(HCAP[, var])))
# }

#---- save dataset ----
HCAP %>% dplyr::select(-one_of("subj_cog_count")) %>% 
  write_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv"))
