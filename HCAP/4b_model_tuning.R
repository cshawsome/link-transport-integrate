#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "LaplacesDemon", "locfit", "wesanderson", "vroom", "mvnfast", 
       "future.apply")

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("HCAP", "functions", "generate_synthetic_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HCAP data ----
HCAP_imputed <- 
  readRDS(paste0(path_to_box, "analyses/HCAP/HCAP_MI_hotdeck")) %>% 
  lapply(., function(x) x %<>% mutate("(Intercept)" = 1))

#just do checks on one dataset
HCAP_analytic <- HCAP_imputed[[1]] 

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairment class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- define vars ----
#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
all_vars <- colnames(HCAP_analytic)
Z <- all_vars[str_detect(all_vars, "_Z")]

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_box, "data/tuning/nu_0_matrix_HCAP.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- 
  read_csv(paste0(path_to_box, "data/tuning/kappa_0_matrix_HCAP.csv"))

#---- generate datasets ----
set.seed(20221116)
start <- Sys.time()

generate_synthetic(warm_up = 100, run_number = 1, 
                   starting_props = c(0.25, 0.25, 0.25, 0.25),
                   dataset_to_copy = HCAP_analytic, calibration_sample = FALSE, 
                   calibration_prop = NA, calibration_sample_name = NA, 
                   path_to_raw_prior_sample = NA,
                   path_to_data = paste0(path_to_box,"data/"), 
                   path_to_analyses_folder = 
                     paste0(path_to_box, "analyses/HCAP/"), 
                   path_to_figures_folder = 
                     paste0(path_to_box, "figures/chapter_6/"), 
                   categorical_vars = W, continuous_vars = Z, 
                   id_var = "HHIDPN", variable_labels = variable_labels, 
                   cell_ID_key = cell_ID_key, color_palette = color_palette, 
                   contrasts_matrix = A, kappa_0_mat = kappa_0_mat, 
                   nu_0_mat = nu_0_mat, num_synthetic = 1000, data_only = FALSE)

end <- Sys.time() - start
