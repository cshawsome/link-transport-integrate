#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "LaplacesDemon", "locfit", "wesanderson", "vroom", "mvnfast", 
       "future.apply")

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("simulation_study", "functions", 
                  "generate_synthetic_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_box, "data/tuning/nu_0_matrix.csv")) 
#scaling for inverse wishart as variance of Beta
kappa_0_mat <- 
  read_csv(paste0(path_to_box, "data/tuning/kappa_0_matrix.csv"))

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/", 
                  "model_coefficients.csv"))
#---- **categorical vars (notation from Schafer 1997) ----
W <- c("black", "hispanic", "stroke")

#---- **continuous vars (notation from Schafer 1997) ----
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

#---- generate datasets in serial ----
#About 1.5 hours to generate data for all datasets in serial
#---- **generate datasets: no calibration ----
#---- ****read in data ----
synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "data/HCAP/synthetic_HCAP_list"))

#---- ****rename datasets based on calibration scenario ----
calibration_scenario = "no_calibration"

synthetic_HCAP_list <- 
  lapply(synthetic_HCAP_list, function(x)
    x %<>% mutate("dataset_name_stem" = unlist(unique(x[, "dataset_name"]))))

synthetic_HCAP_list <- 
  lapply(synthetic_HCAP_list, function(x)
    x %<>% mutate("dataset_name" = 
                    paste0(unlist(unique(x[, "dataset_name_stem"])), "_", 
                           calibration_scenario)))

#---- ****dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

#---- ****specify indices ----
indices <-
  which(dataset_names %in% paste0("HRS_", c(500), "_", calibration_scenario))

set.seed(20220329)
start <- Sys.time()

lapply(synthetic_HCAP_list[indices], function(x)
  generate_synthetic(warm_up = 100, run_number = 1, 
                     starting_props = c(0.25, 0.25, 0.25, 0.25),
                     dataset_to_copy = x, calibration_sample = FALSE, 
                     calibration_prop = NA, calibration_sample_name = NA, 
                     path_to_raw_prior_sample = NA,
                     path_to_data = paste0(path_to_box,"data/"), 
                     path_to_analyses_folder = 
                       paste0(path_to_box, "analyses/simulation_study/HCAP_", 
                              unique(x[, "dataset_name_stem"]), "/"), 
                     path_to_figures_folder = 
                       paste0(path_to_box, "figures/simulation_study/HCAP_", 
                              unique(x[, "dataset_name_stem"]), "/"), 
                     categorical_vars = W, continuous_vars = Z, 
                     id_var = "HHIDPN", variable_labels = variable_labels, 
                     cell_ID_key = cell_ID_key, color_palette = color_palette, 
                     contrasts_matrix = A, kappa_0_mat = kappa_0_mat, 
                     nu_0_mat = nu_0_mat, num_synthetic = 1000, 
                     data_only = FALSE))

end <- Sys.time() - start

#---- **generate datasets: HCAP_50 calibration ----
#---- ****read in data ----
synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "data/HCAP/synthetic_HCAP_list"))

#---- ****rename datasets based on calibration scenario ----
calibration_scenario = "HCAP_50"

synthetic_HCAP_list <- 
  lapply(synthetic_HCAP_list, function(x)
    x %<>% mutate("dataset_name_stem" = unlist(unique(x[, "dataset_name"]))))

synthetic_HCAP_list <- 
  lapply(synthetic_HCAP_list, function(x)
    x %<>% mutate("dataset_name" = 
                    paste0(unlist(unique(x[, "dataset_name_stem"])), "_", 
                           calibration_scenario)))

#---- ****dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

#---- ****specify indices ----
indices <-
  which(dataset_names %in% paste0("HRS_", c(500), "_", calibration_scenario))

set.seed(20220329)
start <- Sys.time()

lapply(synthetic_HCAP_list[indices], function(x)
  generate_synthetic(warm_up = 100, run_number = 1, 
                     starting_props = c(0.25, 0.25, 0.25, 0.25),
                     dataset_to_copy = x, calibration_sample = TRUE, 
                     calibration_prop = 
                       as.numeric(str_remove(calibration_scenario, 
                                             "HCAP_"))/100,
                     calibration_sample_name = calibration_scenario, 
                     path_to_raw_prior_sample = 
                       paste0(path_to_box, "data/prior_data/MI/", 
                              "MI_datasets_cleaned"),
                     path_to_data = paste0(path_to_box,"data/"), 
                     path_to_analyses_folder = 
                       paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
                              unique(x[, "dataset_name_stem"]), "/"), 
                     path_to_figures_folder = 
                       paste0(path_to_box, "figures/simulation_study/HCAP_HRS_", 
                              unique(x[, "dataset_name_stem"]), "/"), 
                     categorical_vars = W, continuous_vars = Z, 
                     id_var = "HHIDPN", variable_labels = variable_labels, 
                     cell_ID_key = cell_ID_key, color_palette = color_palette, 
                     contrasts_matrix = A, kappa_0_mat = kappa_0_mat, 
                     nu_0_mat = nu_0_mat, num_synthetic = 1000, 
                     data_only = FALSE))

end <- Sys.time() - start

# #---- generate datasets in parallel ----
# #About 1.5 hours to generate data for all datasets in serial
# 
# #---- **specify indices ----
# indices <- 
#   which(dataset_names %in% paste0("normal_", c(500), "_ADAMS"))
# 
# set.seed(20220329)
# start <- Sys.time()
# plan(multisession, workers = (availableCores() - 2))
# 
# future_lapply(synthetic_HCAP_list[indices], function(x)
#   generate_synthetic(warm_up = 100, run_number = 1, 
#                      starting_props = c(0.25, 0.25, 0.25, 0.25),
#                      dataset_to_copy = x, calibration_sample = TRUE, 
#                      calibration_prop = 0.5, calibration_sample_name = "HCAP_50", 
#                      path_to_raw_prior_sample = 
#                        paste0(path_to_box, "analyses/simulation_study/", 
#                               "prior_data/MI/MI_datasets_cleaned"),
#                      path_to_data = path_to_box, 
#                      path_to_analyses_folder = 
#                        paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
#                               unique(x[, "dataset_name"]), "/"), 
#                      path_to_figures_folder = 
#                        paste0(path_to_box, "figures/simulation_study/HCAP_HRS_", 
#                               unique(x[, "dataset_name"]), "/"), 
#                      categorical_vars = W, continuous_vars = Z, 
#                      id_var = "HHIDPN", variable_labels = variable_labels, 
#                      cell_ID_key = cell_ID_key, color_palette = color_palette, 
#                      contrasts_matrix = A, kappa_0_mat = kappa_0_mat, 
#                      nu_0_mat = nu_0_mat, num_synthetic = 1000, 
#                      data_only = FALSE), future.seed = TRUE)
# 
# end <- Sys.time() - start
# plan(sequential)
