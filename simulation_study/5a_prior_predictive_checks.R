#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr", "moments", "qdapRegex", "future.apply", "mvnfast", 
       "LaplacesDemon", "vroom")
install_github("thomasp85/gganimate")
#library(gganimate)

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("simulation_study", "functions", 
                  "prior_predictive_checks_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **impairment class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/model_coefficients.csv"))

#---- **categorical ----
#notation from Schafer 1997
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_at("cell_ID", as.character)

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_box, "data/tuning/nu_0_matrix.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_box, "data/tuning/kappa_0_matrix.csv"))

#---- data formatting ----
#---- **user input ----
#calibration scenario options: "no_calibration", "calibration_50_SRS", 
# "calibration_50_design", "calibration_100"
calibration_scenario = "calibration_50_SRS" 

#HCAP sample prop options: 0.25, 0.50
HCAP_sample_prop = 0.50

#---- **read in data ----
synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "data/HCAP/synthetic_HCAP_list"))

#---- **rename datasets based on calibration scenario ----
synthetic_HCAP_list <- 
  lapply(synthetic_HCAP_list, function(x)
    x %<>% mutate("dataset_name_stem" = unlist(unique(x[, "dataset_name"]))))

synthetic_HCAP_list <- 
  lapply(synthetic_HCAP_list, function(x)
    x %<>% mutate("dataset_name" = 
                    paste0(unlist(unique(x[, "dataset_name_stem"])), "_", 
                           calibration_scenario)))

#---- **dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

#---- **specify indices ----
indices <-
  which(dataset_names %in% 
          paste0("HRS_", c(2000, 4000, 8000), "_sample_", HCAP_sample_prop*100, 
                 "_", calibration_scenario))

#---- serial checks ----
set.seed(20220329)
start <- Sys.time()
lapply(synthetic_HCAP_list[indices], function(x)
  prior_predictive_checks(dataset_to_copy = x, calibration_sample = 
                            !(calibration_scenario == "no_calibration"), 
                          calibration_prop = 
                            suppressWarnings(
                              parse_number(calibration_scenario)/100), 
                          calibration_sample_name = calibration_scenario,
                          path_to_data = path_to_box, 
                          path_to_output_folder =
                            paste0(path_to_box,
                                   "figures/chapter_4/simulation_study/HCAP_",
                                   unique(x[, "dataset_name_stem"]),
                                   "/prior_predictive_checks/"), 
                          continuous_check_test = TRUE,
                          continuous_check = c("Unimpaired", "MCI", "Dementia", 
                                               "Other"),
                          categorical_vars = W, continuous_vars = Z,
                          variable_labels = variable_labels, 
                          color_palette = color_palette, contrasts_matrix = A,
                          weights_matrix = cell_ID_key, 
                          kappa_0_mat = kappa_0_mat, nu_0_mat = nu_0_mat,
                          num_synthetic = 1000))
end <- Sys.time() - start

#---- parallel checks ----
set.seed(20220329)
start <- Sys.time()
future_lapply(synthetic_HCAP_list[indices], function(x)
  prior_predictive_checks(dataset_to_copy = x, calibration_sample = 
                            !(calibration_scenario == "no_calibration"), 
                          calibration_prop = 
                            suppressWarnings(
                              parse_number(calibration_scenario)/100), 
                          calibration_sample_name = calibration_scenario,
                          path_to_data = path_to_box, 
                          path_to_output_folder =
                            paste0(path_to_box,
                                   "figures/chapter_4/simulation_study/HCAP_",
                                   unique(x[, "dataset_name_stem"]),
                                   "/prior_predictive_checks/"), 
                          continuous_check_test = TRUE,
                          continuous_check = c("Unimpaired", "MCI", "Dementia", 
                                               "Other"),
                          categorical_vars = W, continuous_vars = Z,
                          variable_labels = variable_labels, 
                          color_palette = color_palette, contrasts_matrix = A,
                          weights_matrix = cell_ID_key, 
                          kappa_0_mat = kappa_0_mat, nu_0_mat = nu_0_mat,
                          num_synthetic = 1000), future.seed = TRUE)
end <- Sys.time() - start
plan(sequential)