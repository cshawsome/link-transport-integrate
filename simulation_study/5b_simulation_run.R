#This script is formatted to run on UCLA Hoffman, so all file paths reference the
# directory setup on that computing cluster

#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MCMCpack", "locfit", 
       "vroom", "mvnfast", "mice")

#---- source functions ----
source(here::here("functions", "read_results.R"))
source(here::here("simulation_study", "functions", 
                  "generate_synthetic_function.R"))
source(here::here("simulation_study", "functions", "simulation_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **data paths ----
superpop_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpop_data_list <- lapply(superpop_data_paths, read_results)

#---- **truth table ----
truth <- read_csv(paste0(path_to_box, "analyses/simulation_study/truth.csv"))

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(here::here("color_palette.csv"))

#---- **all sim scenarios matrix ----
all_sim_scenarios <- read_csv(paste0(path_to_box, "analyses/simulation_study/", 
                                     "sim_study_scenarios.csv"))

#---- define vars ----
#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- colnames(superpop_data_list[[1]])[str_detect(
  colnames(superpop_data_list[[1]]), "_Z")]

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
nu_0_vec <- read_csv(paste0(path_to_box, "analyses/nu_0.csv")) %>% 
  column_to_rownames("dataset_name") %>% t()

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_box, "analyses/kappa_0_matrix.csv"))

#---- run sim ----
set.seed(20220512)

replicate(2,
          simulation_function(warm_up = 100, starting_props = rep(0.25, 4), 
                              unimpaired_preds, other_preds, mci_preds, 
                              categorical_vars = W, continuous_vars = Z, 
                              id_var = "HHIDPN", variable_labels, scenario = 1,
                              superpops_list = superpop_data_list,
                              all_scenarios_list = all_sim_scenarios, 
                              cell_ID_key, color_palette, 
                              num_synthetic = 1000, unimpaired_betas, 
                              unimpaired_cov, other_betas, other_cov, mci_betas, 
                              mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                              prior_beta, nu_0_vec, kappa_0_mat, 
                              contrasts_matrix = A, truth,
                              path_to_results = 
                                paste0(path_to_box, 
                                       "analyses/simulation_study/results/")))




