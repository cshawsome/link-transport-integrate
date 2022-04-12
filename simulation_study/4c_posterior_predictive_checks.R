#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "transformr", "moments", 
       "qdapRegex")

#---- source functions ----
source(here::here("functions", "read_results.R"))
source(here::here("simulation_study", "functions", 
                  "posterior_predictive_checks_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **data paths ----
synthetic_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/synthetic_data/HRS"), 
             full.names = TRUE, pattern = "*.csv")

#subset for now
synthetic_data_paths <- 
  synthetic_data_paths[
    str_detect(synthetic_data_paths, 
               "synthetic_normal_500_ADAMS|synthetic_normal_1000_ADAMS")]


synthetic_data_list <- lapply(synthetic_data_paths, read_results)

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(here::here("color_palette.csv")) 

#---- define vars ----
#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- colnames(synthetic_data_list[[1]])[str_detect(
  colnames(synthetic_data_list[[1]]), "_Z")]

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- **test 1 dataset ----
posterior_predictive_checks(dataset_to_copy = synthetic_data_list[[1]] %>% 
                              group_by(married_partnered) %>% 
                              slice_sample(prop = 0.5) %>% 
                              mutate("(Intercept)" = 1) %>% ungroup() , 
                            categorical_covariates = W, 
                            continuous_covariates = Z, contrasts_matrix = A,
                            cell_ID_key, variable_labels, color_palette,
                            path_to_analyses_folder = 
                              paste0(path_to_box, 
                                     "analyses/simulation_study/HCAP_HRS_", 
                                     unique(synthetic_data_list[[1]][, "dataset_name"]), "/"), 
                            path_to_figures_folder = 
                              paste0(path_to_box,
                                     "figures/simulation_study/HCAP_HRS_", 
                                     unique(synthetic_data_list[[1]][, "dataset_name"]), "/"))

#---- posterior predictive checks ----
set.seed(20220329)
start <- Sys.time()
lapply(synthetic_data_list, function(x)
  posterior_predictive_checks(dataset_to_copy = x %>% 
                                group_by(married_partnered) %>% 
                                slice_sample(prop = 0.5) %>% 
                                mutate("(Intercept)" = 1) %>% ungroup() , 
                              categorical_covariates = W, 
                              continuous_covariates = Z, contrasts_matrix = A,
                              cell_ID_key, variable_labels, num_samples = 10, 
                              num_chains = 1, color_palette,
                              path_to_analyses_folder = 
                                paste0(path_to_box, 
                                       "analyses/simulation_study/HCAP_HRS_", 
                                       unique(x[, "dataset_name"]), "/"), 
                              path_to_figures_folder = 
                                paste0(path_to_box,
                                       "figures/simulation_study/HCAP_HRS_", 
                                       unique(x[, "dataset_name"]), "/")))

end <- Sys.time() - start

