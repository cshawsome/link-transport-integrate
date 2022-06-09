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

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **impairment class color palette ----
color_palette <- read_csv(here::here("color_palette.csv")) 

#---- **synthetic HCAP ----
synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HCAP_list"))

#---- dataset names ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv"))

#---- **categorical ----
#notation from Schafer 1997
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

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
nu_0_mat <- read_csv(paste0(path_to_box, "analyses/nu_0_matrix.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_box, "analyses/kappa_0_matrix.csv"))

#---- run single check ----
#About 1 hour per scenario
set.seed(20220329)
start <- Sys.time()

lapply(synthetic_HCAP_list[which(dataset_names == "normal_500_ADAMS")],
       function(x)
         prior_predictive_checks(unimpaired_preds, other_preds, mci_preds,
                                 categorical_vars = W, continuous_vars = Z,
                                 continuous_check_test = FALSE,
                                 continuous_check = 
                                   c("Unimpaired", "MCI", "Dementia", "Other"),
                                 variable_labels, color_palette,
                                 dataset_to_copy = x,
                                 num_synthetic = 1000, unimpaired_betas,
                                 unimpaired_cov, other_betas, other_cov,
                                 mci_betas, mci_cov, alpha_0_dist, prior_Sigma,
                                 prior_V_inv, prior_beta, nu_0_mat, kappa_0_mat,
                                 contrasts_matrix = A,
                                 path_to_folder =
                                   paste0(path_to_box,
                                          "figures/simulation_study/HCAP_HRS_",
                                          unique(x[, "dataset_name"]),
                                          "/prior_predictive_checks/")))
end <- Sys.time() - start

# #---- run checks in parallel ----
# #1.7 days for all checks to run in serial
# #1.3 hours for 5 runs in parallel
# set.seed(20220329)
# start <- Sys.time()
# plan(multisession, workers = (availableCores() - 2))
# 
# #---- **specify indices ----
# indices <- 
#   which(dataset_names %in% 
#           paste0("normal_", c(500, 1000, 2000, 4000, 8000), "_ADAMS"))
# 
# #---- **run parallel checks ----
# future_lapply(synthetic_HCAP_list[indices], function(x)
#   prior_predictive_checks(unimpaired_preds, other_preds, mci_preds, 
#                           categorical_vars = W, continuous_vars = Z, 
#                           variable_labels, color_palette, 
#                           dataset_to_copy = x , 
#                           num_synthetic = 1000, unimpaired_betas, 
#                           unimpaired_cov, other_betas, other_cov, 
#                           mci_betas, mci_cov, alpha_0_dist, prior_Sigma, 
#                           prior_V_inv, prior_beta, 
#                           nu_0_mat, kappa_0_mat, 
#                           contrasts_matrix = A, 
#                           path_to_folder = 
#                             paste0(path_to_box,
#                                    "figures/simulation_study/HCAP_HRS_", 
#                                    unique(x[, "dataset_name"]), 
#                                    "/prior_predictive_checks/")), 
#   future.seed = TRUE)
# end <- Sys.time() - start
# plan(sequential)
