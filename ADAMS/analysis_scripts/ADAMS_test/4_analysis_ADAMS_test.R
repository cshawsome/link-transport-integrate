#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr", "moments", "vroom", "qdapRegex")
install_github("thomasp85/gganimate")
library(gganimate)

#---- source functions ----
source(here::here("ADAMS", "analysis_scripts", "functions", 
                  "ADAMS_prior_predictive_checks_props_function.R"))
source(here::here("ADAMS", "analysis_scripts", "functions", 
                  "generate_synthetic_function.R"))
source(here::here("ADAMS", "analysis_scripts", "functions", 
                  "ADAMS_posterior_predictive_checks.R"))

#---- read in data ----
#dataset we're trying to copy
dataset_to_copy <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/ADAMS/cleaned/ADAMS_test.csv"))

#complete dataset
ADAMS_data <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                           "data/ADAMS/cleaned/ADAMS_subset_mixed.csv"))

#---- define relevant vars ----
#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- cbind(c("AAGE", "ANMSETOT_norm", "ANSER7T", "ANIMMCR", "ANRECYES", 
             "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Abmi"), 
           c("Age", "Total MMSE", "Serial 7s", "Immediate Word Recall", 
             "Wordlist Recall (Yes)", "Story Recall I", "Proxy Cognition (Avg)", 
             "Delayed Word Recall", "IADLs", "BMI")) %>% 
  set_colnames(c("var", "label"))

#---- prior predictive checks ----
#---- **specifying priors ----
#---- ****latent classes ----
#based on analysis in priors_latent_classes.R
for(group in c("unimpaired", "mci", "other")){
  assign(paste0(group, "_betas"), 
         vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                      "ADAMS_test/latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                      "ADAMS_test/latent_class_", group, "_cov.csv")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- 
  vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
               "ADAMS_test/bootstrap_cell_props.csv")) 

#--- ****beta and sigma ----
priors_beta <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                            "priors/ADAMS_test/priors_beta.csv")) 
prior_V_inv <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "priors/ADAMS_test/priors_V_inv.csv")) 
prior_Sigma <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "priors/ADAMS_test/priors_Sigma.csv")) 

#---- **contrasts matrix ----
A = do.call(cbind, list(
  #intercept
  rep(1, 6),
  #race/ethnicity main effect: Black
  rep(c(1, 0, 0), 2),
  #race/ethnicity main effect: Hispanic
  rep(c(0, 1, 0), 2),
  #stroke main effect
  rep(c(0, 1), each = 3)))

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0 <- 65
#scaling for inverse wishart as variance of Beta
kappa_0 <- c(0.85, 0.85, 0.85, 0.85)

#---- **original means and variances ----
ADAMS_data %<>% 
  dplyr::select("HHIDPN", all_of(W), all_of(Z[, "var"]), "Adem_dx_cat") %>% 
  filter(HHIDPN %in% dataset_to_copy$HHIDPN)

ADAMS_means <- colMeans(ADAMS_data %>% dplyr::select(all_of(Z[, "var"])))
ADAMS_sds <- apply(ADAMS_data %>% dplyr::select(all_of(Z[, "var"])), 2, sd)

#---- **run checks ----
ADAMS_prior_predictive_checks_props(unimpaired_preds, other_preds, mci_preds, 
                                    categorical_vars = W, continuous_vars = Z, 
                                    id_var = "HHIDPN", 
                                    dementia_var = "Adem_dx_cat", 
                                    dataset_to_copy, num_synthetic = 1000, 
                                    unimpaired_betas, unimpaired_cov, 
                                    other_betas, other_cov, mci_betas, mci_cov, 
                                    alpha_0_dist, prior_Sigma, prior_V_inv, 
                                    prior_beta, nu_0, kappa_0, 
                                    contrasts_matrix = A, 
                                    orig_means = ADAMS_means, 
                                    orig_sds = ADAMS_sds, 
                                    path_to_folder = 
                                      paste0("/Users/CrystalShaw/Box/",
                                             "Dissertation/figures/ADAMS_test/", 
                                             "prior_predictive_checks/"))

#---- synthetic datasets ----
#starting_props are for (normal, other, mci, dementia)
generate_synthetic_props(warm_up = 500, run_number = 1, 
                         #warm start
                         starting_props = c(0.40, 0.13, 0.17, 0.30), 
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars = W, continuous_vars = Z, 
                         id_var = "HHIDPN", dementia_var = "Adem_dx_cat", 
                         dataset_to_copy, num_synthetic = 1000, 
                         unimpaired_betas, unimpaired_cov, other_betas, 
                         other_cov, mci_betas, mci_cov, alpha_0_dist, 
                         prior_Sigma, prior_V_inv, prior_beta, nu_0, kappa_0, 
                         contrasts_matrix = A,
                         path_to_analyses_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                  "analyses/ADAMS_test/"), 
                         path_to_figures_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                  "figures/ADAMS_test/"))

generate_synthetic_props(warm_up = 500, run_number = 2, 
                         starting_props = c(0.25, 0.25, 0.25, 0.25),
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars = W, continuous_vars = Z, 
                         id_var = "HHIDPN", dementia_var = "Adem_dx_cat", 
                         dataset_to_copy, num_synthetic = 1000, 
                         unimpaired_betas, unimpaired_cov, other_betas, 
                         other_cov, mci_betas, mci_cov, alpha_0_dist, 
                         prior_Sigma, prior_V_inv, prior_beta, nu_0, kappa_0, 
                         contrasts_matrix = A,
                         path_to_analyses_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                  "analyses/", "ADAMS_test/"), 
                         path_to_figures_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                                  "ADAMS_test/"))

generate_synthetic_props(warm_up = 500, run_number = 3, 
                         starting_props = c(0.10, 0.20, 0.30, 0.40),
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars = W, continuous_vars = Z, 
                         id_var = "HHIDPN", dementia_var = "Adem_dx_cat", 
                         dataset_to_copy, num_synthetic = 1000, unimpaired_betas, 
                         unimpaired_cov, other_betas, other_cov, mci_betas, 
                         mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                         prior_beta, nu_0, kappa_0, contrasts_matrix = A,
                         path_to_analyses_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                                  "ADAMS_test/"), 
                         path_to_figures_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                                  "ADAMS_test/"))

generate_synthetic_props(warm_up = 500, run_number = 4, 
                         starting_props = c(0.10, 0.30, 0.40, 0.20),
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars = W, continuous_vars = Z, 
                         id_var = "HHIDPN", dementia_var = "Adem_dx_cat", 
                         dataset_to_copy, num_synthetic = 1000, unimpaired_betas, 
                         unimpaired_cov, other_betas, other_cov, mci_betas, 
                         mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                         prior_beta, nu_0, kappa_0, contrasts_matrix = A,
                         path_to_analyses_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                                  "ADAMS_test/"), 
                         path_to_figures_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                                  "ADAMS_test/"))

generate_synthetic_props(warm_up = 500, run_number = 5, 
                         starting_props = c(0.05, 0.15, 0.25, 0.55),
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars = W, continuous_vars = Z, 
                         id_var = "HHIDPN", dementia_var = "Adem_dx_cat", 
                         dataset_to_copy, num_synthetic = 1000, unimpaired_betas, 
                         unimpaired_cov, other_betas, other_cov, mci_betas, 
                         mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                         prior_beta, nu_0, kappa_0, contrasts_matrix = A,
                         path_to_analyses_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                                  "ADAMS_test/"), 
                         path_to_figures_folder = 
                           paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                                  "ADAMS_test/"))

#---- posterior predictive checks ----
ADAMS_posterior_predictive_checks(dataset_to_copy, continuous_covariates = Z, 
                                  orig_means = ADAMS_means, 
                                  orig_sds = ADAMS_sds, 
                                  num_samples = 1000, num_chains = 5, 
                                  dementia_var = "Adem_dx_cat", 
                                  path_to_analyses_folder = 
                                    paste0("/Users/CrystalShaw/Box/",
                                           "Dissertation/analyses/ADAMS_test/"), 
                                  path_to_figures_folder = 
                                    paste0("/Users/CrystalShaw/Box/", 
                                           "Dissertation/figures/ADAMS_test/"))
