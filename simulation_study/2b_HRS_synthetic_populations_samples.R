#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "LaplacesDemon", "boot", "MCMCpack")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

#---- **fixed betas ----
fixed_betas <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "fixed_model_coefficients.csv"))

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

#---- synthetic superpopulations ----
set.seed(20220303)

#About 2.2 hours for all superpops
start <- Sys.time()
for(dist_name in c("normal", "lognormal", "bathtub")){
  #---- ****compare with ADAMS ----
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000000, 
                                dementia_prop = 0.35, mci_prop = 0.10, 
                                other_prop = 0.20, dist = dist_name, 
                                parameters = normal_parameter_list, 
                                selected_vars_estimates = fixed_betas,
                                scenario_name = "ADAMS",
                                path_to_results = 
                                  paste0(path_to_box, "analyses/", 
                                         "simulation_study/superpopulations/"))
  #---- ****mostly unimpaired ----
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000000,
                                dementia_prop = 0.20, mci_prop = 0.20,
                                other_prop = 0.10, dist = dist_name,
                                parameters = normal_parameter_list,
                                selected_vars_estimates = fixed_betas,
                                path_to_results =
                                  paste0(path_to_box, "analyses/",
                                         "simulation_study/superpopulations/"))
  #---- ****mostly dementia ----
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000000,
                                dementia_prop = 0.50, mci_prop = 0.20,
                                other_prop = 0.10, dist = dist_name,
                                parameters = normal_parameter_list,
                                selected_vars_estimates = fixed_betas,
                                path_to_results =
                                  paste0(path_to_box, "analyses/",
                                         "simulation_study/superpopulations/"))
}
end <- Sys.time() - start

#---- synthetic HRS ----
#---- **source functions ----
source(here::here("functions", "read_results.R"))

#---- **read in superpop data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **data paths ----
superpop_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpop_data_list <- lapply(superpop_data_paths, read_results)

#---- **create one set of synthetic HRS ----
create_HRS_datasets <- function(superpop, n){
  sample_n(superpop, size = n) %>% 
    separate(dataset_name, sep = "_", into = c("dist", "size", "prior")) %>% 
    mutate_at("size", as.numeric) %>% mutate(size = n) %>% 
    unite("dataset_name", c("dist", "size", "prior"), sep = "_")
}

set.seed(20220507)

for(n in c(500, 1000, 2000, 4000, 8000)){
  if(!exists("synthetic_HRS_list")){
    synthetic_HRS_list <- 
      lapply(superpop_data_list, function(x) create_HRS_datasets(x, n))
  } else{
    synthetic_HRS_list <- 
      append(synthetic_HRS_list, lapply(superpop_data_list, function(x) 
        create_HRS_datasets(x, n)))
  }
}

#---- **save data ----
saveRDS(synthetic_HRS_list, 
        file = paste0(path_to_box, 
                      "analyses/simulation_study/synthetic_HRS_list"))

#---- synthetic HCAP ----
#---- **read in synthetic HRS data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_HRS_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HRS_list"))

#---- **create one set of synthetic HCAP ----
set.seed(20220507)

synthetic_HCAP_list <- 
  lapply(synthetic_HRS_list, 
         function(x) 
           x %>% group_by(married_partnered) %>% slice_sample(prop = 0.5) %>% 
           mutate("(Intercept)" = 1) %>% ungroup())

#---- **save data ----
saveRDS(synthetic_HCAP_list, 
        file = paste0(path_to_box, 
                      "analyses/simulation_study/synthetic_HCAP_list"))


