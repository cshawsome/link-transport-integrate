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

#---- **data paths ----
superpop_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpop_data_list <- lapply(superpop_data_paths, read_results)

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **impairment class color palette ----
color_palette <- read_csv(here::here("color_palette.csv")) 

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv"))

#---- **detect categorical ----
all_cat_vars <- selected_vars[!str_detect(selected_vars$data_label, "_Z"), 
                              "data_label"] %>% unlist()
#remove "Intercept"
all_cat_vars <- as.vector(all_cat_vars[-1])

#categorical vars that continuous vars depend on (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

# #Sanity check-- need to add one for intercept
# nrow(selected_vars) == sum(length(all_cat_vars), length(Z), 1)

#---- specifying priors ----
#---- **latent classes ----
for(group in c("unimpaired", "mci", "other")){
  assign(paste0(group, "_betas"), 
         vroom(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                      "latent_class_", group, "_betas.csv"), delim = ","))
  assign(paste0(group, "_cov"), 
         readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                        "latent_class_", group, "_cov")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                 "imputation_cell_props")) 

#--- ****beta and sigma ----
priors_beta <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                              "prior_data/priors_beta")) 
prior_V_inv <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                              "prior_data/priors_V_inv"))  
prior_Sigma <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                              "prior_data/priors_Sigma")) 

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0 <- read_csv(paste0(path_to_box, "analyses/nu_0.csv")) %>% 
  column_to_rownames("dataset_name") %>% t()

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_box, "analyses/kappa_0_matrix.csv"))

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

#---- run checks ----
#1.7 days for all checks to run in serial
set.seed(20220329)
start <- Sys.time()
plan(multisession, workers = (availableCores() - 2))

future_lapply(synthetic_data_list[seq(31, 43, by = 3)], function(x)
  prior_predictive_checks(unimpaired_preds, other_preds, mci_preds, 
                          all_categorical_vars = all_cat_vars,
                          selected_categorical_vars = W, continuous_vars = Z, 
                          variable_labels, color_palette, 
                          dataset_to_copy = x %>% 
                            group_by(married_partnered) %>% 
                            slice_sample(prop = 0.5) %>% 
                            mutate("(Intercept)" = 1) %>% ungroup(), 
                          num_synthetic = 1000, unimpaired_betas, 
                          unimpaired_cov, other_betas, other_cov, 
                          mci_betas, mci_cov, alpha_0_dist, prior_Sigma, 
                          prior_V_inv, prior_beta, 
                          nu_0 = nu_0[, unique(x$dataset_name)], kappa_0_mat, 
                          contrasts_matrix = A, 
                          path_to_folder = 
                            paste0(path_to_box,
                                   "figures/simulation_study/HCAP_HRS_", 
                                   unique(x[, "dataset_name"]), 
                                   "/prior_predictive_checks/")), 
  future.seed = TRUE)
end <- Sys.time() - start
plan(sequential)
