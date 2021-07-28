#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr")
install_github("thomasp85/gganimate")
library(gganimate)

#---- read in data ----
#dataset we're trying to copy
dataset_to_copy <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                   "data/ADAMS/cleaned/ADAMS_train.csv"))

#complete dataset
ADAMS_data <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
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
for(group in c("normal", "mci", "other")){
  assign(paste0(group, "_betas"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "ADAMS_train/latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "ADAMS_train/latent_class_", group, "_cov.csv")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                  "ADAMS_train/bootstrap_cell_counts.csv")) 

#--- ****beta and sigma ----
priors_beta <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/ADAMS_train/priors_beta.csv")) 
prior_V_inv <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                               "priors/ADAMS_train/priors_V_inv.csv")) 
prior_Sigma <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                               "priors/ADAMS_train/priors_Sigma.csv")) 

#---- ****hyperparameters ----
#DOF for inverse wishart
nu_0 <- 65
#scaling for inverse wishart as variance of Beta
kappa_0 <- c(0.85, 0.85, 0.85, 0.85)


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

#---- **original means and variances ----
ADAMS_data %<>% 
  dplyr::select("HHIDPN", all_of(W), all_of(Z[, "var"]), "Adem_dx_cat") %>% 
  filter(HHIDPN %in% ADAMS_train$HHIDPN)

ADAMS_means <- colMeans(ADAMS_data %>% dplyr::select(all_of(Z[, "var"])))
ADAMS_sds <- apply(ADAMS_data %>% dplyr::select(all_of(Z[, "var"])), 2, sd)

#---- **run checks ----
path_to_folder <- paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "prior_predictive_checks")
