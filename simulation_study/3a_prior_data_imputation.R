#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "miceFast", "ggforce", "future.apply")

options(scipen = 999)

#---- source scripts ----
source(here::here("functions", "fast_impute.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **ADAMS analytic ----
ADAMS_analytic <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))

HCAP_analytic <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))

#---- ADAMS imputation ----
#---- **pare down list of vars based on ADAMS variable selection ----
remove <- c("avg_proxy_cog_Better", "avg_proxy_cog_Same", "avg_proxy_cog_Worse", 
            "proxy_Spouse", "proxy_Child", "proxy_Other_Relative", 
            "proxy_Other")

ADAMS_subset <- ADAMS_analytic %>% dplyr::select(-all_of(remove))

#---- **check missings ----
needs_imputing <- names(which(colSums(is.na(ADAMS_subset)) != 0))

#---- **define imputation var types ----
hrs_waves <- seq(4, 7)
cog_waves <- seq(5, 7)

sociodem_vars <- c("AAGE", "EDYRS", "Female", "Black", "Hispanic", 
                   "Married/partnered", "Working", "Retired", "Not working")

ADAMS_vars <- c("SELFCOG", "ANMSETOT", "ANSER7T", "ANSCISOR", "ANCACTUS", 
                "ANPRES", "ANAFTOT", "ANDELCOR", "ANRECYES", "ANRECNO", 
                "ANWM1TOT", "ANWM2TOT", "ANCPTOT", "ANRCPTOT", "ANTMASEC", 
                "ANBWC20", "ANIMMCR", "ANSMEM2_Better", "ANSMEM2_Same", 
                "ANSMEM2_Worse", "Dementia", "Other", "MCI")

HRS_vars <- c(paste0("r", cog_waves, "mpart"), paste0("r", hrs_waves, "height"), 
              paste0("r", hrs_waves, "weight"), paste0("r", hrs_waves, "smoken"), 
              paste0("r", hrs_waves, "drinkd"), paste0("r", hrs_waves, "drinkn"), 
              paste0("r", hrs_waves, "hibpe"), paste0("r", hrs_waves, "diabe"), 
              paste0("r", hrs_waves, "hearte"), paste0("r", hrs_waves, "stroke"),
              paste0("r", hrs_waves, "adla"), paste0("r", hrs_waves, "iadla"), 
              paste0("r", cog_waves, "imrc"), paste0("r", cog_waves, "dlrc"), 
              paste0("r", cog_waves, "ser7"), paste0("r", cog_waves, "bwc20"), 
              paste0("r", cog_waves, "scis"), paste0("r", cog_waves, "cact"), 
              paste0("r", cog_waves, "pres"), paste0("r", cog_waves, "cogtot"), 
              paste0("r", cog_waves, "pstmem_Better"), 
              paste0("r", cog_waves, "pstmem_Same"), 
              paste0("r", cog_waves, "pstmem_Worse"))

not_predictors <- c("HHIDPN", "AYEAR", "White", "ANMSETOT_norm", 
                    "Adem_dx_label_collapsed", "Unimpaired", "Astroke", "Ahibpe", 
                    "Adiabe", "Ahearte", "Abmi", "Aiadla", "Aadla", "Asmoken", 
                    "Adrinkd", "Adrinkn")

#Sanity check
union_var_names <- c(sociodem_vars, ADAMS_vars, HRS_vars, not_predictors)
ncol(ADAMS_subset) == length(union_var_names)

colnames(ADAMS_subset)[which(!colnames(ADAMS_subset) %in% union_var_names)]
union_var_names[which(!union_var_names %in% colnames(ADAMS_subset))]

#---- **predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(ADAMS_subset)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(ADAMS_subset))

#---- ****cannot predict themselves ----
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- ****non-predictors ----
predict[, not_predictors] <- 0

#---- ****ADAMS-ADAMS prediction ----
#use ADAMS variables to predict each other: remove HRS predictors for now
predict[needs_imputing[which(needs_imputing %in% ADAMS_vars)], HRS_vars] <- 0

#---- ****HRS-HRS prediction ----
#use HRS variables to predict each other: remove ADAMS predictors 
predict[needs_imputing[which(needs_imputing %in% HRS_vars)], ADAMS_vars] <- 0

#---- ****ADAMS extra predictors ----
predict["ANBWC20", paste0("r", cog_waves, "bwc20")] <- 1
predict["ANCACTUS", paste0("r", cog_waves, "cact")] <- 1
predict["SELFCOG", paste0("r", cog_waves, "cogtot")] <- 1
predict["ANDELCOR", paste0("r", cog_waves, "dlrc")] <- 1
predict["ANIMMCR", paste0("r", cog_waves, "imrc")] <- 1
predict["ANPRES", paste0("r", cog_waves, "pres")] <- 1
predict["ANSCISOR", paste0("r", cog_waves, "scis")] <- 1
predict["ANSER7T", paste0("r", cog_waves, "ser7")] <- 1
predict[paste0("ANSMEM2_", c("Better", "Same", "Worse")), 
        paste0("r", cog_waves, "pstmem_Better")] <- 1
predict[paste0("ANSMEM2_", c("Better", "Same", "Worse")), 
        paste0("r", cog_waves, "pstmem_Same")] <- 1
predict[paste0("ANSMEM2_", c("Better", "Same", "Worse")), 
        paste0("r", cog_waves, "pstmem_Worse")] <- 1

#---- ****derived vars (do not impute) ----
remove <- c("ANMSETOT_norm", 
            paste0("A", c("stroke", "hibpe", "diabe", "hearte", "bmi", "iadla", 
                          "adla", "smoken", "drinkd", "drinkn"))) 

predict <- predict[!rownames(predict) %in% remove, ]

#---- ****sanity check ----
#non-predictors should all be 0
colSums(predict[, not_predictors])

#---- **impute data ----
per_run_num_impute = 100
total_num_impute = 10000
chunks <- seq(1, total_num_impute/per_run_num_impute) %>% as.matrix()

#About 13.5 hours
start <- Sys.time()
plan(multisession, workers = (availableCores() - 2))
set.seed(20220329)
future_apply(chunks, 1, function(x) 
  fast_impute(predictor_matrix = predict, data = ADAMS_subset, 
              path_for_output = paste0(path_to_box, "data/ADAMS/prior_data/"),
              method = "PMM", m = per_run_num_impute, maxit = 20, chunk = x), 
  future.seed = TRUE)
end <- Sys.time() - start
plan(sequential)

