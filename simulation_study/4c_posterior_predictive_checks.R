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

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(here::here("color_palette.csv")) 

#---- **synthetic HCAP ----
synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HCAP_list"))

#---- dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

#---- define vars ----
#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
all_vars <- colnames(synthetic_HCAP_list[[1]])
Z <- all_vars[str_detect(all_vars, "_Z")]

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- run checks in parallel ----
set.seed(20220329)
start <- Sys.time()
plan(multisession, workers = (availableCores() - 2))

#---- **specify indices ----
indices <- which(dataset_names %in% 
                   paste0("normal_", c(500), "_ADAMS"))

future_lapply(synthetic_HCAP_list[indices], function(x)
  posterior_predictive_checks(dataset_to_copy = x, calibration_sample = FALSE,
                              calibration_prop = 0.50, 
                              calibration_sample_name = "HCAP_50",
                              categorical_covariates = W, 
                              continuous_covariates = Z, contrasts_matrix = A,
                              cell_ID_key = cell_ID_key, 
                              variable_labels = variable_labels, 
                              color_palette = color_palette,
                              path_to_analyses_folder = 
                                paste0(path_to_box, 
                                       "analyses/simulation_study/HCAP_HRS_", 
                                       unique(x[, "dataset_name"]), "/"), 
                              path_to_figures_folder = 
                                paste0(path_to_box,
                                       "figures/simulation_study/HCAP_HRS_", 
                                       unique(x[, "dataset_name"]), "/")), 
  future.seed = TRUE)

end <- Sys.time() - start
plan(sequential)
