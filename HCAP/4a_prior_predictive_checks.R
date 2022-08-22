#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr", "moments", "qdapRegex", "future.apply", "mvnfast", 
       "LaplacesDemon", "vroom")
# install_github("thomasp85/gganimate")
library(gganimate)

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("HCAP", "functions", "prior_predictive_checks_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HCAP data ----
HCAP_analytic <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic.csv")) 

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(HCAP %in% colnames(HCAP_analytic)) 

#label data
HCAP_analytic %<>% 
  rename_at(vars(variable_labels$HCAP), ~ variable_labels$data_label) %>% 
  mutate("(Intercept)" = 1)

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

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_box, "data/tuning/nu_0_matrix_HCAP.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- 
  read_csv(paste0(path_to_box, "data/tuning/kappa_0_matrix_HCAP.csv"))

#---- run checks in serial ----
set.seed(20220819)
start <- Sys.time()

#---- **run checks ----
prior_predictive_checks(dataset_to_copy = HCAP_analytic, 
                        calibration_sample = FALSE, 
                        calibration_prop = NA, calibration_sample_name = NA,
                        path_to_raw_prior_sample = NA, 
                        path_to_data = path_to_box, 
                        path_to_output_folder =
                          paste0(path_to_box,
                                 "figures/HCAP/prior_predictive_checks/"), 
                        continuous_check_test = TRUE,
                        continuous_check = c("Unimpaired", "MCI", "Dementia", 
                                             "Other"),
                        categorical_vars = W, continuous_vars = Z,
                        variable_labels = variable_labels, 
                        color_palette = color_palette, contrasts_matrix = A,
                        kappa_0_mat = kappa_0_mat, nu_0_mat = nu_0_mat,
                        num_synthetic = 1000)
end <- Sys.time() - start