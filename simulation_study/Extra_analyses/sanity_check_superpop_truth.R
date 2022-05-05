#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "broom", "miceFast", "ggforce")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_analytic <-  
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv")) %>% 
  mutate("age_Z" = scale(AAGE))

#---- ADAMS analysis ----
glm(Dementia ~ age_Z + Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ age_Z + Female + Black + Hispanic + ANMSETOT, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

#---- try one imputation of superpop ----
#---- **source functions ----
source(here::here("functions", "fast_impute.R"))

#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
HRS_clean <- read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv"))

#---- **check missings ----
needs_imputing <- names(which(colSums(is.na(HRS_clean)) != 0))

#---- **define imputation var types ----
not_predictors <- c("HHIDPN", "White")

#---- **predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(HRS_clean)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(HRS_clean))

#---- **cannot predict themselves ----
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- **non-predictors ----
predict[, not_predictors] <- 0

#---- **sanity check ----
#non-predictors should all be 0
colSums(predict[, not_predictors])

#---- **imputation ----
set.seed(20220202)
start <- Sys.time()
fast_impute(predictor_matrix = predict, data = HRS_clean, 
            path_for_output = paste0(path_to_box, "data/HRS/cleaned/"),
            method = "PMM", m = 1, maxit = 15, chunk = 1)
end <- Sys.time() - start

#---- read in results ----
HRS_imputed <- readRDS(paste0(path_to_box, 
                              "data/HRS/cleaned/MI/chunk_1/MI_datasets"))



