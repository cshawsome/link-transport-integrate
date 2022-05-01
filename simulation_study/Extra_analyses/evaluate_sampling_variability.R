#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "LaplacesDemon")

#---- source scripts ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv")) 

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(HRS %in% colnames(HRS_analytic)) 

#---- **selected vars estimates ----
selected_vars_betas <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "selected_vars_model_coefficients.csv"))

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
for(n in c(500, 1000, 2000, 4000, 8000)){
  for(dist_name in c("normal", "lognormal", "bathtub")){
    #---- ****compare with ADAMS ----
    generate_synthetic_continuous(HRS_analytic, sample_size = n, 
                                  unimpaired_prop = 0.35, mci_prop = 0.10, 
                                  dementia_prop = 0.35, dist = dist_name, 
                                  parameters = normal_parameter_list, 
                                  scenario_name = "ADAMS",
                                  path_to_results = 
                                    paste0(path_to_box, "analyses/", 
                                           "simulation_study/synthetic_data/"))
    #---- ****mostly unimpaired ----
    generate_synthetic_continuous(HRS_analytic, sample_size = n, 
                                  unimpaired_prop = 0.50, mci_prop = 0.20, 
                                  dementia_prop = 0.20, dist = dist_name, 
                                  parameters = normal_parameter_list, 
                                  path_to_results = 
                                    paste0(path_to_box, "analyses/", 
                                           "simulation_study/synthetic_data/")) 
    #---- ****mostly dementia ----
    generate_synthetic_continuous(HRS_analytic, sample_size = n, 
                                  unimpaired_prop = 0.20, mci_prop = 0.20, 
                                  dementia_prop = 0.50, dist = dist_name, 
                                  parameters = normal_parameter_list, 
                                  path_to_results = 
                                    paste0(path_to_box, "analyses/", 
                                           "simulation_study/synthetic_data/")) 
  }
}

#---- read in results ----
#---- **data paths ----
synthetic_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/synthetic_data/HRS"), 
             full.names = TRUE, pattern = "*.csv")

synthetic_data <- 
  do.call(rbind, lapply(synthetic_data_paths, read_results))

#----