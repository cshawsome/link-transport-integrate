#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "miceFast", "ggforce", "stringr", "magrittr")

#---- source scripts ----
source(here::here("functions", "fast_impute.R"))

#---- read data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
ADAMS_analytic <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))

#---- check missings ----
needs_imputing <- names(which(colSums(is.na(ADAMS_analytic)) != 0))

#---- define imputation var types ----
hrs_waves <- seq(4, 7)
cog_waves <- seq(5, 7)

sociodem_vars <- c("AAGE", "EDYRS", "Female", "Black", "Hispanic", 
                   "Married/partnered", "Working", "Retired", "Not working")

ADAMS_vars <- c("SELFCOG", "ANMSETOT", "ANSER7T", "ANSCISOR", "ANCACTUS", 
                "ANPRES", "ANAFTOT", "ANBNTTOT", "ANDELCOR", "ANRECYES", 
                "ANRECNO", "ANWM1TOT", "ANWM2TOT", "ANCPTOT", "ANRCPTOT", 
                "ANTMASEC", "ANBWC20", "ANBWC86", "ANIMMCR", "ANSMEM2_Better", 
                "ANSMEM2_Same", "ANSMEM2_Worse", "avg_proxy_cog_Better", 
                "avg_proxy_cog_Same", "avg_proxy_cog_Worse", "Dementia", "Other", 
                "MCI", "proxy_Spouse", "proxy_Child", "proxy_Other_Relative", 
                "proxy_Other")

HRS_vars <- c(paste0("r", cog_waves, "mpart"), paste0("r", hrs_waves, "height"), 
              paste0("r", hrs_waves, "weight"), paste0("r", hrs_waves, "smoken"), 
              paste0("r", hrs_waves, "drinkd"), paste0("r", hrs_waves, "drinkn"), 
              paste0("r", hrs_waves, "hibpe"), paste0("r", hrs_waves, "diabe"), 
              paste0("r", hrs_waves, "hearte"), paste0("r", hrs_waves, "stroke"), 
              paste0("r", hrs_waves, "adla"), paste0("r", hrs_waves, "iadla"), 
              paste0("r", cog_waves, "imrc"), paste0("r", cog_waves, "dlrc"), 
              paste0("r", cog_waves, "ser7"), paste0("r", cog_waves, "bwc20"), 
              paste0("r", seq(5, 6), "bwc86"), paste0("r", cog_waves, "scis"), 
              paste0("r", cog_waves, "cact"), paste0("r", cog_waves, "pres"), 
              paste0("r", cog_waves, "cogtot"), 
              paste0("r", cog_waves, "pstmem_Better"), 
              paste0("r", cog_waves, "pstmem_Same"), 
              paste0("r", cog_waves, "pstmem_Worse"))

not_predictors <- c("HHIDPN", "AYEAR", "White", "ANMSETOT_norm", 
                    "Adem_dx_label_collapsed", "Unimpaired", "Astroke", "Ahibpe", 
                    "Adiabe", "Ahearte", "Abmi", "Aiadla", "Aadla", "Asmoken", 
                    "Adrinkd", "Adrinkn")

# #Sanity check
# union_var_names <- c(sociodem_vars, ADAMS_vars, HRS_vars, not_predictors)
# ncol(ADAMS_analytic) == length(union_var_names)
# 
# colnames(ADAMS_analytic)[which(!colnames(ADAMS_analytic) %in% union_var_names)]
# union_var_names[which(!union_var_names %in% colnames(ADAMS_analytic))]

#---- predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(ADAMS_analytic)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(ADAMS_analytic))

#---- **cannot predict themselves ----
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- **non-predictors ----
predict[, not_predictors] <- 0

#---- **ADAMS-ADAMS prediction ----
#use ADAMS variables to predict each other: remove HRS predictors for now
predict[needs_imputing[which(needs_imputing %in% ADAMS_vars)], HRS_vars] <- 0

#---- **HRS-HRS prediction ----
#use HRS variables to predict each other: remove ADAMS predictors 
predict[needs_imputing[which(needs_imputing %in% HRS_vars)], ADAMS_vars] <- 0

#---- **ADAMS extra predictors ----
predict["ANBWC20", paste0("r", cog_waves, "bwc20")] <- 1
predict["ANBWC86", paste0("r", seq(5, 6), "bwc86")] <- 1
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

#---- **derived vars (do not impute) ----
remove <- c("ANMSETOT_norm", 
            paste0("A", c("stroke", "hibpe", "diabe", "hearte", "bmi", "iadla", 
                          "adla", "smoken", "drinkd", "drinkn"))) 

predict <- predict[!rownames(predict) %in% remove, ]

#---- **sanity check ----
#non-predictors should all be 0
colSums(predict[, not_predictors])

#---- imputation ----
start <- Sys.time()
fast_impute(predictor_matrix = predict, data = ADAMS_analytic, 
            path_for_output = paste0(path_to_box, "data/ADAMS/cleaned/"),
            method = "PMM", m = 25, maxit = 15)
end <- Sys.time() - start
