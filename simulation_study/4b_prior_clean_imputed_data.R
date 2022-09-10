#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "NormPsy")

options(scipen = 999)

#---- read in results ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

data_folders <-  
  list.dirs(path = paste0(path_to_box, "data/ADAMS/prior_data/MI"), 
            full.names = TRUE, recursive = FALSE)

for(directory in data_folders){
  if(!exists("ADAMS_imputed")){
    ADAMS_imputed <- readRDS(paste0(directory, "/MI_datasets"))
  } else{
    ADAMS_imputed %<>% append(., readRDS(paste0(directory, "/MI_datasets")))
  }
}

#variable labels
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))
  
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
health_vars <- c("adla", "iadla", "bmi_derived", "diabe", "drinks_per_week", 
                 "hearte", "hibpe", "smoken", "stroke")

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

dummy_var_check <- function(data, working_vars, subjective_vars){
  indicator <- vector(length = 2) 
  names(indicator) = c("Working", "Subjective")
  indicator["Working"] <- data %>% dplyr::select(all_of(working_vars)) %>% 
    rowSums() %>% table() %>% as.data.frame() %>% nrow()
  indicator["Subjective"] <- data %>% dplyr::select(all_of(subjective_vars)) %>% 
    rowSums() %>% table() %>% as.data.frame() %>% nrow()
  
  return(indicator)
}

dummy_check_results <- 
  lapply(ADAMS_imputed, dummy_var_check, Working, Subjective) %>% 
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

# #Sanity check
# dummy_recheck_results <- 
#   lapply(ADAMS_imputed, dummy_var_check, Working, Subjective, Proxy) %>% 
#   do.call(rbind, .)

#---- define predictor variable types ----
sociodemographics <- c("Female", "Married/partnered", "Not working", "Retired", 
                       "AAGE", "EDYRS", "Black", "Hispanic")
neuropsych <- c("ANAFTOT", "ANCPTOT", "ANRCPTOT", "ANMSETOT_norm", "ANRECYES", 
                "ANRECNO", "ANTMASEC", "ANWM1TOT", "ANWM2TOT")
gen_cog <- c("SELFCOG", "ANBWC20", "ANCACTUS", "ANDELCOR", "ANIMMCR","ANPRES", 
             "ANSCISOR", "ANSER7T", "ANSMEM2_Better", "ANSMEM2_Worse")
functional <- c("Aadla", "Aiadla")
health <- c("Adiabe", "Ahearte", "Ahibpe", "Asmoken", "Astroke", 
            "Amoderate_drinking", "Aheavy_drinking", "Abmi_derived")

outcome <- c("Unimpaired", "MCI", "Dementia", "Other")

#ID vars + vars that won't be used in models, but want to keep for summary stats 
keep <- c("HHIDPN", "Working", "White", "ANSMEM2_Same") 

#---- keep vars in analytic datasets ----
analysis_vars <- c(keep, sociodemographics, neuropsych, gen_cog, functional, 
                   health, outcome)

ADAMS_imputed_clean <- 
  lapply(ADAMS_imputed, function(x) x %<>% dplyr::select(all_of(analysis_vars)))

# #Sanity check
# lapply(ADAMS_imputed_clean, ncol)

#---- standardize continuous vars ----
standardize_vars <- c("AAGE", "EDYRS", "ANMSETOT_norm", "ANAFTOT", "ANCPTOT", 
                      "ANRCPTOT", "ANRECYES", "ANRECNO", "ANTMASEC", "ANWM1TOT", 
                      "ANWM2TOT", "SELFCOG", "ANDELCOR", "ANIMMCR", "ANSER7T", 
                      "Aadla", "Aiadla", "Abmi_derived")

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

#---- **name variables ----
variable_labels_ADAMS <- variable_labels %>% 
  filter(ADAMS %in% colnames(ADAMS_imputed_clean[[1]]))

ADAMS_imputed_clean <- 
  lapply(ADAMS_imputed_clean, 
         function(x) rename_at(x, vars(variable_labels_ADAMS$ADAMS), ~ 
                                 variable_labels_ADAMS$data_label)) 

#---- **save clean datasets ----
saveRDS(ADAMS_imputed_clean, 
        file = paste0(path_to_box, "data/prior_data/MI/MI_datasets_cleaned"))
