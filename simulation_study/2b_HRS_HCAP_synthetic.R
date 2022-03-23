#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "LaplacesDemon", "MBSP")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv")) 

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
#---- **normal ----
#---- ****compare with ADAMS ----
generate_synthetic_continuous(HRS_analytic, sample_size = 1000, 
                              unimpaired_prop = 0.35, mci_prop = 0.10, 
                              dementia_prop = 0.35, dist = "normal", 
                              parameters = normal_parameter_list, 
                              path_to_results = 
                                paste0(path_to_box, "analyses/", 
                                       "simulation_study/synthetic_data/", 
                                       "ADAMS_props/"))

#---- ****500 mostly unimpaired ----
generate_synthetic_continuous(HRS_analytic, sample_size = 500, 
                              unimpaired_prop = 0.50, mci_prop = 0.20, 
                              dementia_prop = 0.20, dist = "normal", 
                              parameters = normal_parameter_list, 
                              path_to_results = 
                                paste0(path_to_box, "analyses/", 
                                       "simulation_study/synthetic_data/")) 

#---- ****500 mostly dementia ----
generate_synthetic_continuous(HRS_analytic, sample_size = 500, 
                              unimpaired_prop = 0.20, mci_prop = 0.20, 
                              dementia_prop = 0.50, dist = "normal", 
                              parameters = normal_parameter_list, 
                              path_to_results = 
                                paste0(path_to_box, "analyses/", 
                                       "simulation_study/synthetic_data/")) 

#---- ****1000 mostly unimpaired ----
generate_synthetic_continuous(HRS_analytic, sample_size = 1000, 
                              unimpaired_prop = 0.50, mci_prop = 0.20, 
                              dementia_prop = 0.20, dist = "normal", 
                              parameters = normal_parameter_list, 
                              path_to_results = 
                                paste0(path_to_box, "analyses/", 
                                       "simulation_study/synthetic_data/")) 

#---- ****1000 mostly dementia ----
generate_synthetic_continuous(HRS_analytic, sample_size = 1000, 
                              unimpaired_prop = 0.20, mci_prop = 0.20, 
                              dementia_prop = 0.50, dist = "normal", 
                              parameters = normal_parameter_list, 
                              path_to_results = 
                                paste0(path_to_box, "analyses/", 
                                       "simulation_study/synthetic_data/")) 



