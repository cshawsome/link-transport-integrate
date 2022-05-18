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

#---- **data paths ----
superpop_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpop_data_list <- lapply(superpop_data_paths, read_results)

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(here::here("color_palette.csv")) 

#---- define vars ----
#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- colnames(synthetic_data_list[[1]])[str_detect(
  colnames(synthetic_data_list[[1]]), "_Z")]

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- create one set of synthetic HRS ----
create_HRS_datasets <- function(superpop, n){
  sample_n(superpop, size = n) %>% 
    separate(dataset_name, sep = "_", into = c("dist", "size", "prior")) %>% 
    mutate_at("size", as.numeric) %>% mutate(size = n) %>% 
    unite("dataset_name", c("dist", "size", "prior"), sep = "_")
}

set.seed(20220507)

for(n in c(500, 1000, 2000, 4000, 8000)){
  if(!exists("synthetic_data_list")){
    synthetic_data_list <- 
      lapply(superpop_data_list, function(x) create_HRS_datasets(x, n))
  } else{
    synthetic_data_list <- 
      append(synthetic_data_list, lapply(superpop_data_list, function(x) 
        create_HRS_datasets(x, n)))
  }
}

# #Sanity check
# lapply(synthetic_data_list, function(x) unique(x$dataset_name))

#---- **dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_data_list, function(x) unique(x$dataset_name)))

#---- run checks in parallel ----
set.seed(20220329)
start <- Sys.time()
plan(multisession, workers = (availableCores() - 2))

#---- **specify indices ----
indices <- which(dataset_names %in% 
                   paste0("normal_", c(500, 1000), "_ADAMS"))

future_lapply(synthetic_data_list[indices], function(x)
  posterior_predictive_checks(dataset_to_copy = x %>% 
                                group_by(married_partnered) %>% 
                                slice_sample(prop = 0.5) %>% 
                                mutate("(Intercept)" = 1) %>% ungroup(), 
                              categorical_covariates = W, 
                              continuous_covariates = Z, 
                              contrasts_matrix = A,
                              cell_ID_key, variable_labels, color_palette,
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
