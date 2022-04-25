#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "vroom", "mvnfast")

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("simulation_study", "functions", 
                  "generate_synthetic_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **data paths ----
synthetic_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/synthetic_data/HRS"), 
             full.names = TRUE, pattern = "*.csv")

#Focus on full analysis with normal, ADAMS props data
synthetic_data_paths <- 
  synthetic_data_paths[unlist(lapply(synthetic_data_paths, 
                                     function(x) str_detect(x, "ADAMS") & 
                                       str_detect(x, "_normal")))]

synthetic_data_list <- lapply(synthetic_data_paths, read_results)

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
nu_0 <- 65
#scaling for inverse wishart as variance of Beta
kappa_0 <- c(0.85, 0.85, 0.85, 0.85) %>% 
  set_names(c("Unimpaired", "MCI", "Dementia", "Other"))

#---- code testing ----
#---- **test 1: serial; one at a time ----
#n = 500 dataset (1.05 mins laptop)
#n = 1000 dataset (XX mins desktop, 1.13 min laptop)
#n = 2000 dataset (XX mins desktop, XX min laptop)
#n = 4000 dataset (XX mins desktop, XX min laptop)
#n = 8000 dataset (XX mins desktop, XX min laptop)
#
#Just data generation
#n = 500 dataset 
#n = 1000 dataset (0.8 mins laptop)
#n = 2000 dataset 
#n = 4000 dataset
#n = 8000 dataset
set.seed(20220329)
start <- Sys.time()
generate_synthetic(warm_up = 100, run_number = 1, 
                   starting_props = c(0.25, 0.25, 0.25, 0.25),
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, 
                   id_var = "HHIDPN", variable_labels, 
                   dataset_to_copy = synthetic_data_list[[1]] %>% 
                     group_by(married_partnered) %>% 
                     slice_sample(prop = 0.5) %>% 
                     mutate("(Intercept)" = 1) %>% ungroup(), cell_ID_key, 
                   color_palette, num_synthetic = 1000, unimpaired_betas, 
                   unimpaired_cov, other_betas, other_cov, mci_betas, mci_cov, 
                   alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                   kappa_0, contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
                            unique(synthetic_data_list[[1]][, "dataset_name"]), 
                            "/"), 
                   path_to_figures_folder = 
                     paste0(path_to_box,
                            "figures/simulation_study/HCAP_HRS_", 
                            unique(synthetic_data_list[[1]][, "dataset_name"]), 
                            "/"), data_only = FALSE)
end <- Sys.time() - start

#---- ****code timing ----
library(profvis)
profvis::profvis(
  generate_synthetic(warm_up = 100, run_number = 1, 
                     starting_props = c(0.25, 0.25, 0.25, 0.25),
                     unimpaired_preds, other_preds, mci_preds, 
                     categorical_vars = W, continuous_vars = Z, 
                     id_var = "HHIDPN", variable_labels, 
                     dataset_to_copy = synthetic_data_list[[2]] %>% 
                       group_by(married_partnered) %>% 
                       slice_sample(prop = 0.5) %>% 
                       mutate("(Intercept)" = 1) %>% ungroup(), cell_ID_key, 
                     color_palette, num_synthetic = 1000, unimpaired_betas, 
                     unimpaired_cov, other_betas, other_cov, mci_betas, mci_cov, 
                     alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                     kappa_0, contrasts_matrix = A,
                     path_to_analyses_folder = 
                       paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
                              unique(synthetic_data_list[[2]][, "dataset_name"]), 
                              "/"), 
                     path_to_figures_folder = 
                       paste0(path_to_box,
                              "figures/simulation_study/HCAP_HRS_", 
                              unique(synthetic_data_list[[2]][, "dataset_name"]), 
                              "/"), 
                     data_only = FALSE))

#---- **test 2: parallel ----
set.seed(20220329)
start <- Sys.time()
plan(multisession, workers = (availableCores() - 2))

future_lapply(synthetic_data_list[1:2], function(x)
  generate_synthetic(warm_up = 100, run_number = 1, 
                     starting_props = c(0.25, 0.25, 0.25, 0.25),
                     unimpaired_preds, other_preds, mci_preds, 
                     categorical_vars = W, continuous_vars = Z, 
                     id_var = "HHIDPN", variable_labels, 
                     dataset_to_copy = x %>% 
                       group_by(married_partnered) %>% 
                       slice_sample(prop = 0.5) %>% 
                       mutate("(Intercept)" = 1) %>% ungroup(), cell_ID_key, 
                     color_palette, num_synthetic = 1000, unimpaired_betas, 
                     unimpaired_cov, other_betas, other_cov, mci_betas, mci_cov, 
                     alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                     kappa_0, contrasts_matrix = A,
                     path_to_analyses_folder = 
                       paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
                              unique(x[, "dataset_name"]), "/"), 
                     path_to_figures_folder = 
                       paste0(path_to_box,
                              "figures/simulation_study/HCAP_HRS_", 
                              unique(x[, "dataset_name"]), "/")), 
  future.seed = TRUE)
end <- Sys.time() - start
plan(sequential)

#---- generate synthetic data ----
set.seed(20220329)
start <- Sys.time()

lapply(synthetic_data_list[1:2], function(x)
  generate_synthetic(warm_up = 100, run_number = 1, 
                     starting_props = c(0.25, 0.25, 0.25, 0.25),
                     unimpaired_preds, other_preds, mci_preds, 
                     categorical_vars = W, continuous_vars = Z, 
                     id_var = "HHIDPN", variable_labels, 
                     dataset_to_copy = x %>% 
                       group_by(married_partnered) %>% 
                       slice_sample(prop = 0.5) %>% 
                       mutate("(Intercept)" = 1) %>% ungroup(), cell_ID_key, 
                     color_palette, num_synthetic = 1000, unimpaired_betas, 
                     unimpaired_cov, other_betas, other_cov, mci_betas, mci_cov, 
                     alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                     kappa_0, contrasts_matrix = A,
                     path_to_analyses_folder = 
                       paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
                              unique(x[, "dataset_name"]), "/"), 
                     path_to_figures_folder = 
                       paste0(path_to_box,
                              "figures/simulation_study/HCAP_HRS_", 
                              unique(x[, "dataset_name"]), "/")))
end <- Sys.time() - start
