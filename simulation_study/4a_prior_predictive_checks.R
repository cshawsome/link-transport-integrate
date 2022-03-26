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

#---- **impairment class color palette ----
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
                        id_var = "HHIDPN", variable_labels, color_palette, 
                        dataset_to_copy, num_synthetic = 10, unimpaired_betas, 
                        unimpaired_cov, other_betas, other_cov, mci_betas, 
                        mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                        prior_beta, nu_0, kappa_0, contrasts_matrix = A, 
                        path_to_folder = 
                          paste0(path_to_box,
                                 "figures/simulation_study/", 
                                 "HCAP_normal_250_unimpaired/", 
                                 "prior_predictive_checks/"))
