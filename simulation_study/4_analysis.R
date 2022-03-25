#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr", "moments", "qdapRegex")
install_github("thomasp85/gganimate")
library(gganimate)

#---- source functions ----
source(here::here("simulation_study", "functions", 
                  "prior_predictive_checks_function.R"))

source(here::here("simulation_study", "functions", 
                  "generate_synthetic_function.R"))

source(here::here("simulation_study", "functions", 
                  "posterior_predictive_checks_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **dataset to copy ----
#dataset we're trying to create synthetic versions of
dataset_to_copy <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/synthetic_data/HCAP/", 
                  "HCAP_synthetic_normal_250_unimpaired.csv")) %>% 
  mutate("(Intercept)" = 1)

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
Z <- colnames(dataset_to_copy)[str_detect(colnames(dataset_to_copy), "_Z")]

#---- prior predictive checks ----
#---- **specifying priors ----
#---- ****latent classes ----
for(group in c("unimpaired", "mci", "other")){
  assign(paste0(group, "_betas"), 
         read_csv(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                         "latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         read_csv(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                         "latent_class_", group, "_cov.csv")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                  "imputation_cell_props.csv")) 

#--- ****beta and sigma ----
priors_beta <- read_csv(paste0(path_to_box, "analyses/simulation_study/",
                               "prior_data/priors_beta.csv")) 
prior_V_inv <- read_csv(paste0(path_to_box, "analyses/simulation_study/",
                               "prior_data/priors_V_inv.csv"))  
prior_Sigma <- read_csv(paste0(path_to_box, "analyses/simulation_study/",
                               "prior_data/priors_Sigma.csv")) 

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0 <- 65
#scaling for inverse wishart as variance of Beta
kappa_0 <- c(0.85, 0.85, 0.85, 0.85) %>% 
  set_names(c("Unimpaired", "MCI", "Dementia", "Other"))

#---- **run checks ----
prior_predictive_checks(unimpaired_preds, other_preds, mci_preds, 
                        categorical_vars = W, continuous_vars = Z, 
                        id_var = "HHIDPN", variable_labels = variable_labels,
                        dementia_var = "Adem_dx_cat", 
                        dataset_to_copy, num_synthetic = 10, 
                        unimpaired_betas, unimpaired_cov, 
                        other_betas, other_cov, mci_betas, mci_cov, 
                        alpha_0_dist, prior_Sigma, prior_V_inv, 
                        prior_beta, nu_0, kappa_0, 
                        contrasts_matrix = A, 
                        path_to_folder = 
                          paste0(path_to_box,
                                 "figures/simulation_study/", 
                                 "HCAP_normal_250_unimpaired/", 
                                 "prior_predictive_checks/"))

#---- generate synthetic data ----
generate_synthetic(warm_up = 2, run_number = 1, 
                   starting_props = c(0.25, 0.25, 0.25, 0.25),
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN", 
                   variable_labels, dataset_to_copy, cell_ID_key, color_palette,
                   num_synthetic = 10, unimpaired_betas, unimpaired_cov, 
                   other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
                   prior_Sigma, prior_V_inv, prior_beta, nu_0, kappa_0, 
                   contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0(path_to_box, "analyses/simulation_study/", 
                            "HCAP_normal_250_unimpaired/"), 
                   path_to_figures_folder = 
                     paste0(path_to_box, "figures/simulation_study/", 
                            "HCAP_normal_250_unimpaired/"))

#---- posterior predictive checks ----
posterior_predictive_checks(dataset_to_copy, continuous_covariates = Z, 
                            num_samples = 1000, num_chains = 1, 
                            path_to_analyses_folder = 
                              paste0("/Users/CrystalShaw/Box/",
                                     "Dissertation/analyses/ADAMS_test/"), 
                            path_to_figures_folder = 
                              paste0("/Users/CrystalShaw/Box/", 
                                     "Dissertation/figures/ADAMS_test/"))


