#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- source functions ----
source(here::here("HCAP", "functions", "analysis_function.R"))
source(here::here("HCAP", "functions", "generate_synthetic_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HCAP data ----
HCAP_imputed <- 
  readRDS(paste0(path_to_box, "analyses/HCAP/HCAP_MI_hotdeck")) %>% 
  lapply(., function(x) x %<>% mutate("(Intercept)" = 1))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairment class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/", 
                  "model_coefficients.csv"))

#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

#---- run analysis ----
set.seed(20221125)

MI_results <- lapply(HCAP_imputed, function(x) 
  analysis_function(warm_up = 100, starting_props = c(0.25, 0.25, 0.25, 0.25),
                    dataset_to_copy = x,
                    orig_means = 
                      x %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
                      colMeans() %>% t() %>% as.data.frame(),
                    orig_sds = 
                      x %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
                      apply(., 2, sd) %>% t() %>% as.data.frame(),
                    calibration_sample = FALSE, calibration_prop = NA,
                    calibration_sample_name = NA,
                    path_to_data = paste0(path_to_box,"data/"),
                    categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN",
                    variable_labels = variable_labels, cell_ID_key = cell_ID_key,
                    color_palette = color_palette, contrasts_matrix = A,
                    kappa_0_mat = 
                      read_csv(paste0(path_to_box, 
                                      "data/tuning/kappa_0_matrix_HCAP.csv")),
                    nu_0_mat = 
                      read_csv(paste0(path_to_box, 
                                      "data/tuning/nu_0_matrix_HCAP.csv")),
                    num_synthetic = 1000))

#---- plot results ----