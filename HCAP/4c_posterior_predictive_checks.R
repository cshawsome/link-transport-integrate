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

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- define vars ----
#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
all_vars <- colnames(HCAP_analytic)
Z <- all_vars[str_detect(all_vars, "_Z")]

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- run checks ----
set.seed(20220329)
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
                              paste0(path_to_box,"figures/HCAP/"))

end <- Sys.time() - start
