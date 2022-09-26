#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "wesanderson", "RColorBrewer", "transformr", "moments", 
       "qdapRegex", "future.apply")

#---- source functions ----
source(here::here("functions", "read_results.R"))
source(here::here("simulation_study", "functions", 
                  "posterior_predictive_checks_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))

#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/", 
                  "model_coefficients.csv"))

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- define vars ----
#---- **categorical vars (notation from Schafer 1997) ----
W <- c("black", "hispanic", "stroke")

#---- **continuous vars (notation from Schafer 1997) ----
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- data formatting ----
#---- **user input ----
#calibration scenario options: "no_calibration", "HCAP_50", "HCAP_100"
calibration_scenario = "no_calibration" 

#HCAP sample prop options: 0.25, 0.50
HCAP_sample_prop = 0.25

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
start <- Sys.time()

lapply(synthetic_HCAP_list[indices], function(x)
  posterior_predictive_checks(dataset_to_copy = x, 
                              calibration_sample = 
                                !(calibration_scenario == "no_calibration"),
                              calibration_prop = suppressWarnings(
                                as.numeric(str_remove(calibration_scenario, 
                                                      "HCAP_"))/100), 
                              calibration_sample_name = calibration_scenario,
                              categorical_covariates = W, 
                              continuous_covariates = Z, contrasts_matrix = A,
                              cell_ID_key = cell_ID_key, 
                              variable_labels = variable_labels, 
                              color_palette = color_palette,
                              path_to_analyses_folder = 
                                paste0(path_to_box, 
                                       "analyses/simulation_study/HCAP_", 
                                       unique(x[, "dataset_name_stem"]), "/"), 
                              path_to_figures_folder = 
                                paste0(path_to_box,
                                       "figures/chapter_4/simulation_study/HCAP_", 
                                       unique(x[, "dataset_name_stem"]), "/")))

end <- Sys.time() - start

#---- parallel checks ----
start <- Sys.time()

future_lapply(synthetic_HCAP_list[indices], function(x)
  posterior_predictive_checks(dataset_to_copy = x, 
                              calibration_sample = 
                                !(calibration_scenario == "no_calibration"),
                              calibration_prop = suppressWarnings(
                                as.numeric(str_remove(calibration_scenario, 
                                                      "HCAP_"))/100), 
                              calibration_sample_name = calibration_scenario,
                              categorical_covariates = W, 
                              continuous_covariates = Z, contrasts_matrix = A,
                              cell_ID_key = cell_ID_key, 
                              variable_labels = variable_labels, 
                              color_palette = color_palette,
                              path_to_analyses_folder = 
                                paste0(path_to_box, 
                                       "analyses/simulation_study/HCAP_", 
                                       unique(x[, "dataset_name_stem"]), "/"), 
                              path_to_figures_folder = 
                                paste0(path_to_box,
                                       "figures/chapter_4/simulation_study/HCAP_", 
                                       unique(x[, "dataset_name_stem"]), "/")), 
  future.seed = TRUE)

end <- Sys.time() - start
plan(sequential)