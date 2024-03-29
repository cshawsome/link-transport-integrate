#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "miceFast", "ggforce", "stringr", "magrittr")

#---- source scripts ----
source(here::here("functions", "fast_impute.R"))

#---- read data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

#---- check missings ----
needs_imputing <- names(which(colSums(is.na(HRS_analytic)) != 0))
needs_imputing <- needs_imputing[-1]

#---- define imputation var types ----
not_predictors <- c("HHIDPN", "HCAP_SELECT")

#---- predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(HRS_analytic)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(HRS_analytic))

#---- **cannot predict themselves ----
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- **non-predictors ----
predict[, not_predictors] <- 0

#---- **sanity check ----
#non-predictors should all be 0
colSums(predict[, not_predictors])

#---- imputation ----
set.seed(20220211)
start <- Sys.time()
fast_impute(predictor_matrix = predict, data = HRS_analytic, 
            path_for_output = paste0(path_to_box, "data/HRS/cleaned/"),
            method = "PMM", m = 2, maxit = 15)
end <- Sys.time() - start

#---- read in results ----
HRS_imputed <- 
  readRDS(paste0(path_to_box, "data/HRS/cleaned/MI/MI_datasets"))

#---- OLD ----

#---- derive post-imputation variables ----
hrs_waves <- seq(4, 7)

derive_vars <- function(data, waves){
  for(wave in waves){
    data %<>% 
      #---- **bmi ----
    mutate(!!paste0("r", wave, "bmi_derived") := 
             !!sym(paste0("r", wave, "weight"))/
             (!!sym(paste0("r", wave, "height"))*
                !!sym(paste0("r", wave, "height")))) %>% 
      #---- **drinks per week ----
    mutate(!!paste0("r", wave, "drinks_per_week") := 
             !!sym(paste0("r", wave, "drinkd"))*
             !!sym(paste0("r", wave, "drinkn")))
  }
  return(data)
}

ADAMS_imputed <- lapply(ADAMS_imputed, derive_vars, hrs_waves)
ADAMS_imputed <- lapply(ADAMS_imputed, function(x) x %<>% 
                          mutate("ANMSETOT_norm" = normMMSE(ANMSETOT)))

# #Sanity check
# test <- ADAMS_imputed[[1]]
# View(test[, c("r4weight", "r4height", "r4bmi_derived")])
# head(test[, "r4weight"]/(test[, "r4height"]^2))
# View(test[, c("r4drinkd", "r4drinkn", "r4drinks_per_week")])
# tail(colnames(test))

#---- representative health variables ----
health_vars <- c("adla", "bmi_derived", "diabe", "drinks_per_week", "hearte", 
                 "hibpe", "iadla", "smoken", "stroke")

rep_vars <- function(data, vars){
  for(var in vars){
    data %<>% 
      mutate(!!paste0("A", var) := 
               case_when(AYEAR %in% c(2001, 2002) ~ !!sym(paste0("r5", var)), 
                         AYEAR %in% c(2003, 2004) ~ !!sym(paste0("r6", var))))
  }
  return(data)
}

ADAMS_imputed <- lapply(ADAMS_imputed, rep_vars, health_vars)

# #Sanity check
# test <- ADAMS_imputed[[1]]
# View(test[, c("AYEAR", "r5adla", "r6adla", "Aadla")])
# tail(colnames(test)) #only expect bmi_derived and drinks_per_week at end

#---- **derive drinking category ----
drinking_stat <- function(data){
  data %<>% 
    mutate("Adrink_cat" = 
             case_when(Adrinks_per_week == 0 ~ 1,
                       Female == 0 & 
                         (Adrinks_per_week >= 1 & Adrinks_per_week < 14) ~ 2,
                       Female == 1 &
                         (Adrinks_per_week >= 1 & Adrinks_per_week < 7) ~ 2,
                       Female == 0 &
                         (Adrinks_per_week >= 14 | Adrinkn >= 4) ~ 3,
                       Female == 1 &
                         (Adrinks_per_week >= 7 | Adrinkn >= 3) ~ 3)) %>% 
    mutate("Adrink_cat_label" = 
             case_when(Adrink_cat == 1 ~ "No Drinking", 
                       Adrink_cat == 2 ~ "Moderate Drinking", 
                       Adrink_cat == 3 ~ "Heavy Drinking")) %>% 
    mutate("Ano_drinking" = ifelse(Adrink_cat_label == "No Drinking", 1, 0), 
           "Amoderate_drinking" = 
             ifelse(Adrink_cat_label == "Moderate Drinking", 1, 0), 
           "Aheavy_drinking" = 
             ifelse(Adrink_cat_label == "Heavy Drinking", 1, 0))
  
  return(data)
}

ADAMS_imputed <- lapply(ADAMS_imputed, drinking_stat)

# #Sanity check
# test <- ADAMS_imputed[[1]]
# table(test$Adrink_cat_label, test$Ano_drinking, useNA = "ifany")
# table(test$Adrink_cat_label, test$Amoderate_drinking, useNA = "ifany")
# table(test$Adrink_cat_label, test$Aheavy_drinking, useNA = "ifany")
# View(test %>% 
#        dplyr::select("Female", "Adrinks_per_week", "Adrinkn", "Adrink_cat_label"))
# tail(colnames(test)) 

#---- clean: multi-cat vars ----
#working status, subjective cognitive decline, proxy cognition,  
# #make sure none are missing
# test <- ADAMS_imputed[[1]]
# colSums(is.na(test %>% 
#                 dplyr::select("Working", "Not working", "Retired", 
#                               "avg_proxy_cog_Better", "avg_proxy_cog_Same", 
#                               "avg_proxy_cog_Worse", "ANSMEM2_Better", 
#                               "ANSMEM2_Same", "ANSMEM2_Worse")))

#---- check if any datasets have more than one dummy indicator ----
Working <- c("Working", "Not working", "Retired")
Subjective <- c("ANSMEM2_Better", "ANSMEM2_Same", "ANSMEM2_Worse")
Proxy <- c("avg_proxy_cog_Better", "avg_proxy_cog_Same", 
           "avg_proxy_cog_Worse")

dummy_var_check <- function(data, working_vars, subjective_vars, proxy_vars){
  indicator <- vector(length = 3) 
  names(indicator) = c("Working", "Subjective", "Proxy")
  indicator["Working"] <- data %>% dplyr::select(all_of(working_vars)) %>% 
    rowSums() %>% table() %>% as.data.frame() %>% nrow()
  indicator["Subjective"] <- data %>% dplyr::select(all_of(subjective_vars)) %>% 
    rowSums() %>% table() %>% as.data.frame() %>% nrow()
  indicator["Proxy"] <- data %>% dplyr::select(all_of(proxy_vars)) %>% 
    rowSums() %>% table() %>% as.data.frame() %>% nrow()
  
  return(indicator)
}

#The subjective and proxy subjective cognitive decline variable has some 
# duplicates
dummy_check_results <- 
  lapply(ADAMS_imputed, dummy_var_check, Working, Subjective, Proxy) %>% 
  do.call(rbind, .)

#might move this out to its own script
select_cat <- function(data, vars){
  max_cats <- length(vars)
  data %<>% mutate("count" = rowSums(data[, vars]))
  data[which(data$count %in% c(0, max_cats)), sample(vars, size = 1)] <- 1
  
  probs <- which(!data$count %in% c(0, 1, max_cats))
  if(length(probs) != 0){
    for(row in 1:length(probs)){
      which_cats <- vars[which(data[row, vars] == 1)]
      data[row, vars] <- 0
      data[row, sample(which_cats, size = 1)] <- 1
    }
  }
  
  return(data %>% dplyr::select(-one_of("count")))
}

ADAMS_imputed <- lapply(ADAMS_imputed, select_cat, Subjective)
ADAMS_imputed <- lapply(ADAMS_imputed, select_cat, Proxy)

# #Sanity check
# dummy_recheck_results <- 
#   lapply(ADAMS_imputed, dummy_var_check, Working, Subjective, Proxy) %>% 
#   do.call(rbind, .)

#---- define predictor variable types ----
sociodemographics <- c("Female", "Married/partnered", "Not working", "Retired", 
                       "AAGE", "EDYRS", "Black", "Hispanic")
neuropsych <- c("ANAFTOT", "ANBNTTOT", "ANCPTOT", "ANRCPTOT", "ANMSETOT_norm", 
                "ANRECNO", "ANRECYES", "ANTMASEC", "ANWM1TOT", "ANWM2TOT")
gen_cog <- c("SELFCOG", "ANBWC20", "ANBWC86", "ANCACTUS", "ANDELCOR", "ANIMMCR",
             "ANPRES", "ANSCISOR", "ANSER7T", "ANSMEM2_Better", "ANSMEM2_Worse", 
             "avg_proxy_cog_Better", "avg_proxy_cog_Worse")
functional <- c("Aadla", "Aiadla")
health <- c("Adiabe", "Ahearte", "Ahibpe", "Asmoken", "Astroke", 
            "Amoderate_drinking", "Aheavy_drinking", "Abmi_derived")

outcome <- c("Unimpaired", "MCI", "Dementia", "Other")

#ID vars + vars that won't be used in models, but want to keep for summary stats 
keep <- c("HHIDPN", "Working", "White", "ANSMEM2_Same", "avg_proxy_cog_Same", 
          "Ano_drinking") 

#---- keep vars in analytic datasets ----
analysis_vars <- c(keep, sociodemographics, neuropsych, gen_cog, functional, 
                   health, outcome)

ADAMS_imputed_clean <- 
  lapply(ADAMS_imputed, function(x) x %<>% dplyr::select(all_of(analysis_vars)))

# #Sanity check
# lapply(ADAMS_imputed_clean, ncol)

#---- standardize continuous vars ----
standardize_vars <- c("AAGE", "EDYRS", "ANMSETOT_norm", "ANAFTOT", "ANBNTTOT", 
                      "ANCPTOT", "ANRCPTOT","ANRECNO", "ANRECYES", "ANTMASEC", 
                      "ANWM1TOT", "ANWM2TOT", "SELFCOG", "ANBWC20", "ANBWC86", 
                      "ANDELCOR", "ANIMMCR", "ANSER7T", "Aadla", "Aiadla", 
                      "Abmi_derived")

Z_score <- function(data, vars){
  subset <- data %>% dplyr::select(all_of(vars)) %>% 
    mutate_all(scale) %>% set_colnames(paste0(all_of(vars), "_Z"))
  
  data %<>% cbind(., subset)
  
  return(data)
}

ADAMS_imputed_clean <- lapply(ADAMS_imputed_clean, Z_score, standardize_vars)

# #Sanity check
# sapply(ADAMS_imputed_clean, ncol)
# test <- ADAMS_imputed_clean[[1]]
# tail(colnames(test))

#---- **save clean datasets ----
saveRDS(ADAMS_imputed_clean, 
        file = paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned"))
