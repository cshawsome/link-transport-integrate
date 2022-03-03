#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "LaplacesDemon", "MBSP")

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
source(here("simulation_study", "functions", "generate_synthetic_continuous.R"))

#---- format data ----
HRS_analytic %<>% mutate("Intercept" = 1) %>% 
  rename_at(vars(variable_labels$HRS), ~ variable_labels$data_label)

#---- synthetic data ----
set.seed(20220303)
#---- **multivariate normal ----
synthetic_normal_1000 <- 
  generate_synthetic_continuous(HRS_analytic, sample_size = 1000, 
                                unimpaired_prop = 0.35, mci_prop = 0.10, 
                                dementia_prop = 0.35, dist = "normal", 
                                parameters = normal_parameter_list)

write_csv(as.data.frame(synthetic_normal_1000), 
          paste0(path_to_box, "analyses/simulation_study/", "synthetic_data/", 
                 "synthetic_normal_1000.csv"))

