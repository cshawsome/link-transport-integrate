#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr", "moments", "qdapRegex", "future.apply", "mvnfast", 
       "LaplacesDemon", "vroom")
install_github("thomasp85/gganimate")
library(gganimate)

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("simulation_study", "functions", 
                  "prior_predictive_checks_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **impairment class color palette ----
color_palette <- read_csv(here::here("color_palette.csv")) 

#---- **synthetic HCAP ----
synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HCAP_list"))

#---- dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv"))

#---- **categorical ----
#notation from Schafer 1997
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_box, "analyses/nu_0_matrix.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_box, "analyses/kappa_0_matrix.csv"))

#---- run checks in parallel ----
#1.7 days for all checks to run in serial
#1.3 hours for 5 runs in parallel
set.seed(20220329)
start <- Sys.time()
plan(multisession, workers = (availableCores() - 2))

#---- **specify indices ----
indices <-
  which(dataset_names %in%
          paste0("normal_", c(500, 1000, 2000, 4000, 8000), "_ADAMS"))

#---- **run parallel checks ----
future_lapply(synthetic_HCAP_list[indices], function(x)
  prior_predictive_checks(dataset_to_copy = x, calibration_sample = TRUE, 
                          calibration_prop = 0.50, 
                          calibration_sample_name = "HCAP_50",
                          path_to_raw_prior_sample = 
                            paste0(path_to_box, "analyses/simulation_study/", 
                                   "prior_data/MI/MI_datasets_cleaned"), 
                          path_to_data = path_to_box, 
                          path_to_output_folder =
                            paste0(path_to_box,
                                   "figures/simulation_study/HCAP_HRS_",
                                   unique(x[, "dataset_name"]),
                                   "/prior_predictive_checks/"), 
                          continuous_check_test = TRUE,
                          continuous_check = 
                            c("Unimpaired", "MCI", "Dementia", "Other"),
                          categorical_vars = W, continuous_vars = Z,
                          variable_labels, color_palette, contrasts_matrix = A,
                          kappa_0_mat = kappa_0_mat, nu_0_mat = nu_0_mat,
                          num_synthetic = 1000), future.seed = TRUE)
end <- Sys.time() - start
plan(sequential)