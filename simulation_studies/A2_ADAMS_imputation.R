#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "mice")

#---- read data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
ADAMS_analytic <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))

#---- check missings ----
needs_imputing <- names(which(colSums(is.na(ADAMS_analytic)) != 0))

#---- define imputation var types ----
ADAMS_vars <- c("SELFCOG", "AAGE", "EDYRS", "ANMSETOT", "ANSER7T", "ANSCISOR", 
                "ANCACTUS", "ANPRES", "ANAFTOT", "ANBNTTOT", "ANDELCOR", 
                "ANRECYES", )

not_predictors <- c("HHIDPN", "AYEAR")

#---- predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(ADAMS_analytic)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(ADAMS_analytic))

#cannot predict themselves
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- OLD ----
  time_updated_vars <- c("married_partnered", "not_married_partnered", 
                         "widowed", "drinking_cat", "memrye_impute", 
                         "stroke_impute", "hearte_impute", "lunge_impute", 
                         "cancre_impute", "hibpe_impute", "diabe_impute", 
                         "cesd", "BMI")

time_invariant_vars <- c("ed_cat", "black", "hispanic", "other", 
                         "female", "survtime", "death2018", "smoker", 
                         "observed")

blocks <- c(apply(expand.grid("r", seq(4, 9), time_updated_vars), 1, 
                  paste, collapse = ""))


#---- **time-updated var models ----
for(var in time_updated_vars){
  #can't predict itself
  predict[paste0("r", seq(4, 9), var), paste0("r", seq(4, 9), var)] <- 
    (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
  
  #can't predict in the same wave (all missing)
  predictors <- time_updated_vars[which(time_updated_vars != var)]
  for(predictor in predictors){
    predict[paste0("r", seq(4, 9), var), 
            paste0("r", seq(4, 9), predictor)] <- 
      (diag(x = 1, nrow = 6, ncol = 6) == 0)*1
  }
  
  #use time-updated age
  predict[paste0("r", seq(4, 9), var), 
          paste0("r", seq(4, 9), "age_y_int")] <- 
    diag(x = 1, nrow = 6, ncol = 6)
}

#Don't use these as predictors
predict[, c("HHIDPN", "intercept", 
            paste0("r", seq(3, 9), "conde_impute"), "white", "r3cesd", 
            paste0("r", seq(3, 9), "shlt"), "age_death_y", 
            "r4cesd_elevated", "r9cesd_elevated", "total_elevated_cesd",
            "prop_elevated_cesd", "avg_cesd", "avg_cesd_elevated", 
            "observed", paste0("r", seq(4, 9), "cesd_death2018"), 
            paste0("r", seq(3, 8), "cesd_conde_impute"))] <- 0

# #Sanity check
# colSums(predict)
# 
#---- imputation ----
data_imputed <- 
  mice::mice(data = ADAMS_analytic, m = 2, maxit = 2, method = "pmm",
             predictorMatrix = predict, where = is.na(ADAMS_analytic), 
             blocks = as.list(rownames(predict)), seed = 20220126)

