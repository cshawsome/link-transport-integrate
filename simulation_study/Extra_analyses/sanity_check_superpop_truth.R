#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "broom", "miceFast", "ggforce", 
       "locfit", "LaplacesDemon", "MCMCpack")

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

glm(Dementia ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ age_Z + Female + Black + Hispanic + ANMSETOT, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

#---- **sex/gender + race/ethnicity only ----
glm(Unimpaired ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(MCI ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Other ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

#---- **sex/gender + race/ethnicity only | impairment class ----
glm(Dementia ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(MCI ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic %>% filter(Dementia == 0)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Other ~ Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic %>% 
      filter(Dementia == 0 & MCI == 0)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))


#---- try one imputation of superpop with selected vars betas ----
#---- **source functions ----
source(here::here("functions", "fast_impute.R"))
source(here("simulation_study", "functions", 
            "generate_synthetic_continuous_function.R"))

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

#---- **read in results ----
HRS_imputed <- readRDS(paste0(path_to_box, 
                              "data/HRS/cleaned/MI/chunk_1/MI_datasets")) %>% 
  do.call(rbind, .)

#---- **generate synthetic pop ----
#---- ****selected vars betas ----
selected_vars_betas <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "selected_vars_model_coefficients.csv"))

#---- ****variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(HRS %in% colnames(HRS_imputed)) 

#---- **continuous distribution parameters ----
normal_parameter_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/", 
                 "continuous_distribution_parameters/", 
                 "normal_parameters"))

#---- **format data ----
HRS_imputed %<>% 
  rename_at(vars(variable_labels$HRS), ~ variable_labels$data_label) %>% 
  rename("Intercept" = "intercept")

#---- **synthetic data ----
set.seed(20220303)

start <- Sys.time()
for(dist_name in c("normal")){
  #---- ****compare with ADAMS ----
  generate_synthetic_continuous(HRS_imputed, sample_size = 1000000, 
                                unimpaired_prop = 0.35, mci_prop = 0.10, 
                                dementia_prop = 0.35, dist = dist_name, 
                                parameters = normal_parameter_list, 
                                selected_vars_estimates = selected_vars_betas,
                                scenario_name = "ADAMS",
                                path_to_results = 
                                  paste0(path_to_box, "analyses/", 
                                         "simulation_study/", 
                                         "test_superpopulations/"))
}

end <- Sys.time() - start

#---- **read in data ----
superpop <- read_csv(paste0(path_to_box, "analyses/simulation_study/", 
                            "test_superpopulations/normal_1000000_ADAMS.csv"))

#---- **test overall model ----
glm(Dementia ~ age_Z + female + black + hispanic, 
    family = "poisson", data = superpop) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

#---- **test stratified model ----
glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% filter(female == 1)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% filter(female == 1 & age_Z > 2)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% 
      filter(female == 1 & age_Z > 1 & age_Z < 2)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% 
      filter(female == 1 & age_Z < -1)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

#---- try one imputation of superpop with fixed vars betas ----
#---- **source functions ----
source(here("simulation_study", "functions", 
            "generate_synthetic_continuous_function.R"))

#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

HRS_imputed <- readRDS(paste0(path_to_box, 
                              "data/HRS/cleaned/MI/chunk_1/MI_datasets")) %>% 
  do.call(rbind, .)

#---- **generate synthetic pop ----
#---- ****selected vars betas ----
fixed_betas <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "fixed_model_coefficients.csv"))

#---- ****variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(HRS %in% colnames(HRS_imputed)) 

#---- **continuous distribution parameters ----
normal_parameter_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/", 
                 "continuous_distribution_parameters/", 
                 "normal_parameters"))

#---- **format data ----
HRS_imputed %<>% 
  rename_at(vars(variable_labels$HRS), ~ variable_labels$data_label) %>% 
  rename("Intercept" = "intercept")

#---- **synthetic data ----
set.seed(20220303)

start <- Sys.time()
for(dist_name in c("normal")){
  #---- ****compare with ADAMS ----
  generate_synthetic_continuous(HRS_imputed, sample_size = 1000000, 
                                dementia_prop = 0.35, mci_prop = 0.10, 
                                other_prop = 0.20, dist = dist_name, 
                                parameters = normal_parameter_list, 
                                selected_vars_estimates = fixed_betas,
                                scenario_name = "ADAMS_fixed_betas",
                                path_to_results = 
                                  paste0(path_to_box, "analyses/", 
                                         "simulation_study/", 
                                         "test_superpopulations/"))
}

end <- Sys.time() - start

#---- **read in data ----
superpop <- read_csv(paste0(path_to_box, "analyses/simulation_study/", 
                            "test_superpopulations/", 
                            "normal_1000000_ADAMS_fixed_betas.csv"))

#---- **test overall model ----
glm(Dementia ~ age_Z + female + black + hispanic + stroke, 
    family = "poisson", data = superpop) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

#---- **test stratified model ----
glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% filter(female == 1)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% filter(female == 1 & age_Z > 2)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% 
      filter(female == 1 & age_Z > 1 & age_Z < 2)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ black + hispanic, 
    family = "poisson", data = superpop %>% 
      filter(female == 1 & age_Z < -1)) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))




