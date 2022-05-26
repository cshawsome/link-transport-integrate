#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "LaplacesDemon", "locfit", "wesanderson", "vroom", "mvnfast")

#---- source functions ----
source(here::here("functions", "read_results.R"))

source(here::here("simulation_study", "functions", 
                  "generate_synthetic_function.R"))

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
A <- read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_box, "analyses/nu_0_matrix.csv")) 
#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_box, "analyses/kappa_0_matrix.csv"))

#---- generate sets for tuning ----
#About 1.5 hours to generate data for all datasets

#---- **specify indices ----
indices <- which(dataset_names %in% 
                   paste0("normal_", c(500, 1000, 2000), "_ADAMS"))

set.seed(20220329)
start <- Sys.time()

lapply(synthetic_data_list[indices], 
       function(x)
         generate_synthetic(warm_up = 100, run_number = 1, 
                            starting_props = c(0.25, 0.25, 0.25, 0.25),
                            unimpaired_preds, other_preds, mci_preds, 
                            categorical_vars = W, continuous_vars = Z, 
                            id_var = "HHIDPN", variable_labels, 
                            dataset_to_copy = x, cell_ID_key, color_palette, 
                            num_synthetic = 1000, unimpaired_betas, 
                            unimpaired_cov, other_betas, other_cov, mci_betas, 
                            mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                            prior_beta, nu_0_mat, kappa_0_mat, 
                            contrasts_matrix = A,
                            path_to_analyses_folder = 
                              paste0(path_to_box, "analyses/simulation_study/", 
                                     "HCAP_HRS_", 
                                     unique(x[, "dataset_name"]), "/"), 
                            path_to_figures_folder = 
                              paste0(path_to_box,
                                     "figures/simulation_study/HCAP_HRS_", 
                                     unique(x[, "dataset_name"]), "/"), 
                            data_only = FALSE))

end <- Sys.time() - start

