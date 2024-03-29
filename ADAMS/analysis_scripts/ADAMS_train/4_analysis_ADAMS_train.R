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
                  "ADAMS_prior_predictive_checks_counts_function.R"))
source(here::here("ADAMS", "analysis_scripts", "functions", 
                  "generate_synthetic_function.R"))
source(here::here("ADAMS", "analysis_scripts", "functions", 
                  "ADAMS_posterior_predictive_checks.R"))

#---- read in data ----
#dataset we're trying to copy
dataset_to_copy <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/ADAMS/cleaned/ADAMS_train.csv"))

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
                      "ADAMS_train/latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                      "ADAMS_train/latent_class_", group, "_cov.csv")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- 
  vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
               "ADAMS_train/bootstrap_cell_counts.csv")) 

#--- ****beta and sigma ----
priors_beta <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                            "priors/ADAMS_train/priors_beta.csv")) 
prior_V_inv <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "priors/ADAMS_train/priors_V_inv.csv")) 
prior_Sigma <- vroom(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "priors/ADAMS_train/priors_Sigma.csv")) 

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
ADAMS_prior_predictive_checks_counts(unimpaired_preds, other_preds, mci_preds, 
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
                                              "Dissertation/figures/", 
                                              "ADAMS_train/", 
                                              "prior_predictive_checks/"))

#---- synthetic datasets ----
#starting_props are for (normal, other, mci, dementia)
#---- **visualize starting points ----
groups <- c("Unimpaired", "Other", "MCI", "Dementia") 
plot_data <- 
  cbind(expand.grid(paste0("Chain ", seq(1, 5)), groups),
        expand.grid(seq(1, 5), seq(1, 4))) %>% 
  set_colnames(c("Chain", "Group", "y", "x")) %>% 
  mutate("props" = NA)

plot_data[which(plot_data$Chain == "Chain 1" & plot_data$Group %in% groups), 
          "props"] <- c(0.40, 0.20, 0.10, 0.30)
plot_data[which(plot_data$Chain == "Chain 2" & plot_data$Group %in% groups), 
          "props"] <- c(0.25, 0.25, 0.25, 0.25)
plot_data[which(plot_data$Chain == "Chain 3" & plot_data$Group %in% groups), 
          "props"] <- c(0.10, 0.20, 0.30, 0.40)
plot_data[which(plot_data$Chain == "Chain 4" & plot_data$Group %in% groups), 
          "props"] <- c(0.10, 0.30, 0.40, 0.20)
plot_data[which(plot_data$Chain == "Chain 5" & plot_data$Group %in% groups), 
          "props"] <- c(0.05, 0.15, 0.25, 0.55)

ggplot(plot_data, aes(fill = Group, y = props, x = Chain)) + 
  geom_bar(position = "fill", stat = "identity") + theme_minimal() + 
  scale_fill_manual(values = c("#00a389", "#28bed9", "#fdab00", "#ff0000")) + 
  ylab("Proportion") + xlab(" ")

ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "ADAMS_train/diagnostics/chain_starts.jpeg"), 
       width = 5, height = 3, units = "in")

#---- **generate datasets ----
generate_synthetic(warm_up = 500, run_number = 1, 
                   #warm start
                   starting_props = c(0.40, 0.20, 0.10, 0.30), 
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN", 
                   dementia_var = "Adem_dx_cat", dataset_to_copy, 
                   num_synthetic = 1000, unimpaired_betas, unimpaired_cov, 
                   other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
                   count = "yes", prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                   kappa_0, contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "analyses/ADAMS_train/"), 
                   path_to_figures_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "figures/ADAMS_train/"))

generate_synthetic(warm_up = 500, run_number = 2, 
                   starting_props = c(0.25, 0.25, 0.25, 0.25),
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN", 
                   dementia_var = "Adem_dx_cat", dataset_to_copy, 
                   num_synthetic = 1000, unimpaired_betas, unimpaired_cov, 
                   other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
                   count = "yes", prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                   kappa_0, contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "analyses/ADAMS_train/"), 
                   path_to_figures_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "figures/ADAMS_train/"))

generate_synthetic(warm_up = 500, run_number = 3, 
                   starting_props = c(0.10, 0.20, 0.30, 0.40),
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN", 
                   dementia_var = "Adem_dx_cat", dataset_to_copy, 
                   num_synthetic = 1000, unimpaired_betas, unimpaired_cov, 
                   other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
                   count = "yes",prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                   kappa_0, contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "analyses/ADAMS_train/"), 
                   path_to_figures_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "figures/ADAMS_train/"))

generate_synthetic(warm_up = 500, run_number = 4, 
                   starting_props = c(0.10, 0.30, 0.40, 0.20),
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN", 
                   dementia_var = "Adem_dx_cat", dataset_to_copy, 
                   num_synthetic = 1000, unimpaired_betas, unimpaired_cov, 
                   other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
                   count = "yes", prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                   kappa_0, contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "analyses/ADAMS_train/"), 
                   path_to_figures_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "figures/ADAMS_train/"))

generate_synthetic(warm_up = 500, run_number = 5, 
                   starting_props = c(0.05, 0.15, 0.25, 0.55),
                   unimpaired_preds, other_preds, mci_preds, 
                   categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN", 
                   dementia_var = "Adem_dx_cat", dataset_to_copy, 
                   num_synthetic = 1000, unimpaired_betas, unimpaired_cov, 
                   other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
                   count = "yes", prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                   kappa_0, contrasts_matrix = A,
                   path_to_analyses_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "analyses/ADAMS_train/"), 
                   path_to_figures_folder = 
                     paste0("/Users/CrystalShaw/Box/Dissertation/", 
                            "figures/ADAMS_train/"))

#---- posterior predictive checks ----
ADAMS_posterior_predictive_checks(dataset_to_copy, continuous_covariates = Z, 
                                  orig_means = ADAMS_means, 
                                  orig_sds = ADAMS_sds, 
                                  num_samples = 1000, num_chains = 5, 
                                  dementia_var = "Adem_dx_cat", 
                                  path_to_analyses_folder = 
                                    paste0("/Users/CrystalShaw/Box/",
                                           "Dissertation/analyses/ADAMS_train/"), 
                                  path_to_figures_folder = 
                                    paste0("/Users/CrystalShaw/Box/", 
                                           "Dissertation/figures/ADAMS_train/"))


