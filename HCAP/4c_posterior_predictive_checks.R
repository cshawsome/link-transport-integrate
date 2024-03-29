#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "wesanderson", "RColorBrewer", "transformr", "moments", 
       "qdapRegex", "future.apply")

#---- source functions ----
source(here::here("functions", "read_results.R"))
source(here::here("HCAP", "functions", "posterior_predictive_checks_function.R"))

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

#---- **impairement class color palette ----
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
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- run checks ----
set.seed(20221116)
start <- Sys.time()

posterior_predictive_checks(dataset_to_copy = HCAP_analytic, 
                            calibration_sample = FALSE,
                            calibration_prop = NA, calibration_sample_name = NA,
                            categorical_covariates = W, 
                            continuous_covariates = Z, contrasts_matrix = A,
                            cell_ID_key = cell_ID_key, 
                            variable_labels = variable_labels, 
                            color_palette = color_palette,
                            path_to_analyses_folder = 
                              paste0(path_to_box, "analyses/HCAP/"), 
                            path_to_figures_folder = 
                              paste0(path_to_box,"figures/chapter_6/"))

end <- Sys.time() - start
