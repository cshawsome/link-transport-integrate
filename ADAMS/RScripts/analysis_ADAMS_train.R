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
#---- **ADAMS ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"))

ADAMS_data <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                              "data/cleaned/ADAMS/ADAMS_subset_mixed.csv"))

#---- prior predictive checks ----
#---- **specifying priors ----
#---- ****latent classes ----
for(group in c("normal", "mci", "other")){
  assign(paste0(group, "_betas"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "latent_class_", group, "_cov.csv")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                  "bootstrap_cell_counts.csv")) 

#--- ****beta and sigma ----
priors_beta <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/priors_beta.csv")) 
prior_V_inv <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                               "priors/priors_V_inv.csv")) 
prior_Sigma <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                               "priors/priors_Sigma.csv")) 

#---- ****hyperparameters ----
#DOF for inverse wishart
nu_0 <- 65
#scaling for inverse wishart as variance of Beta
kappa_0 <- c(0.85, 0.85, 0.85, 0.85)