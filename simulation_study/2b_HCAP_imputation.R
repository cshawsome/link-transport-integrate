#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "tidyr", "miceFast", "ggforce")

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

#---- **summarize missingness ----
#double check that all of these are in the rownames of the imputation matrix
needs_imputing <- names(colMeans(is.na(HCAP))[which(colMeans(is.na(HCAP)) > 0)])

#---- source functions ----
source(here::here("functions", "fast_impute.R"))

#---- define imputation var types ----
not_predictors <- c("HHIDPN", "White", "Other", "Working", "r13no_drinking",
                    "subj_cog_same")

#---- predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(HCAP)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(HCAP))

#---- **cannot predict themselves ----
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- **non-predictors ----
predict[, not_predictors] <- 0

#---- imputation ----
#About 6 seconds
set.seed(20220202)
start <- Sys.time()
fast_impute(predictor_matrix = predict, data = HCAP, 
            path_for_output = paste0(path_to_box, "data/HCAP/cleaned/"),
            method = "PMM", m = 1, maxit = 15, chunk = 1)
end <- Sys.time() - start

#---- read in results ----
HCAP_imputed <- 
  readRDS(paste0(path_to_box, "data/HCAP/cleaned/MI/chunk_1/MI_datasets")) %>% 
  as.data.frame()

# #Sanity check
# colMeans(is.na(HCAP_imputed))[which(colMeans(is.na(HCAP_imputed)) > 0)]

#---- OLD ----


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
HCAP %>% dplyr::select(-one_of("subj_cog_count", "contingency_cell")) %>% 
  write_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv"))
