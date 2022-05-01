#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "LaplacesDemon", "boot")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

#---- **selected vars betas ----
selected_vars_betas <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "selected_vars_model_coefficients.csv"))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(HRS %in% colnames(HRS_analytic)) 

#---- **continuous distribution parameters ----
normal_parameter_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/", 
                 "continuous_distribution_parameters/", 
                 "normal_parameters"))

#---- source functions ----
source(here("simulation_study", "functions", 
            "generate_synthetic_continuous_function.R"))

#---- format data ----
HRS_analytic %<>% mutate("Intercept" = 1) %>% 
  rename_at(vars(variable_labels$HRS), ~ variable_labels$data_label)

#---- synthetic data ----
set.seed(20220303)

start <- Sys.time()
for(dist_name in c("normal", "lognormal", "bathtub")){
  #---- ****compare with ADAMS ----
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000000, 
                                unimpaired_prop = 0.35, mci_prop = 0.10, 
                                dementia_prop = 0.35, dist = dist_name, 
                                parameters = normal_parameter_list, 
                                selected_vars_estimates = selected_vars_betas,
                                scenario_name = "ADAMS",
                                path_to_results = 
                                  paste0(path_to_box, "analyses/", 
                                         "simulation_study/superpopulations/"))
  #---- ****mostly unimpaired ----
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000000, 
                                unimpaired_prop = 0.50, mci_prop = 0.20, 
                                dementia_prop = 0.20, dist = dist_name, 
                                parameters = normal_parameter_list, 
                                selected_vars_estimates = selected_vars_betas,
                                path_to_results = 
                                  paste0(path_to_box, "analyses/", 
                                         "simulation_study/superpopulations/")) 
  #---- ****mostly dementia ----
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000000, 
                                unimpaired_prop = 0.20, mci_prop = 0.20, 
                                dementia_prop = 0.50, dist = dist_name, 
                                parameters = normal_parameter_list, 
                                selected_vars_estimates = selected_vars_betas,
                                path_to_results = 
                                  paste0(path_to_box, "analyses/", 
                                         "simulation_study/superpopulations/")) 
}
end <- Sys.time() - start
