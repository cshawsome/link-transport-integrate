#This script is formatted to run on UCLA Hoffman, so all file paths reference the
# directory setup on that computing cluster

#---- hoffman setup ----
## Read in the arguments listed in the:
## R CMD BATCH --no-save --no-restore '--args mechanism="MNAR" method="JMVN" mask_percent="10%"' test_args2.R
## expression:
args=(commandArgs(TRUE))

## args is now a list of character vectors

## Check to see if arguments are passed and set default values if not,
## then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  scenario_num = 1
  num_replicates = 1
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
## Now print values just to make sure:
print(scenario_num)
print(num_replicates)

#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "MCMCpack", "locfit", "vroom", 
       "mvnfast", "mice", "LaplacesDemon")

#---- source functions ----
path_to_RScripts <- "/u/home/c/cshaw343/link_transport_integrate/RScripts/"

source(paste0(path_to_RScripts, "read_results.R"))
source(paste0(path_to_RScripts, "generate_synthetic_function.R"))
source(paste0(path_to_RScripts, "simulation_function.R"))

#---- read in data ----
path_to_data <- "/u/home/c/cshaw343/link_transport_integrate/data/"

#---- **data paths ----
superpop_data_paths <- 
  list.files(path = paste0(path_to_data, "superpopulations"), full.names = TRUE, 
             pattern = "*.csv")

superpop_data_list <- lapply(superpop_data_paths, read_results)

#---- **truth table ----
truth <- read_csv(paste0(path_to_data, "truth.csv"))

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_data, "variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_data, "cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairement class color palette ----
color_palette <- read_csv(paste0(path_to_data, "color_palette.csv"))

#---- **all sim scenarios matrix ----
all_sim_scenarios <- read_csv(paste0(path_to_data, "sim_study_scenarios.csv"))

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
         vroom(paste0(path_to_data, "latent_class_", group, "_betas.csv"), 
               delim = ","))
  assign(paste0(group, "_cov"), 
         readRDS(paste0(path_to_data, "latent_class_", group, "_cov")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- ****contingency cells ----
alpha_0_dist <- readRDS(paste0(path_to_data, "imputation_cell_props")) 

#--- ****beta and sigma ----
priors_beta <- readRDS(paste0(path_to_data, "priors_beta")) 
prior_V_inv <- readRDS(paste0(path_to_data, "priors_V_inv"))  
prior_Sigma <- readRDS(paste0(path_to_data, "priors_Sigma"))

#---- **contrasts matrix ----
A = read_csv(paste0(path_to_data, "contrasts_matrix.csv")) %>% as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_data, "nu_0_matrix.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_data, "kappa_0_matrix.csv"))

#---- set seed ----
seed <- as.numeric(all_sim_scenarios[scenario_num, "seed"])
set.seed(seed)

#---- run sim ----
replicate(num_replicates,
          simulation_function(warm_up = 100, starting_props = rep(0.25, 4), 
                              unimpaired_preds, other_preds, mci_preds, 
                              categorical_vars = W, continuous_vars = Z, 
                              id_var = "HHIDPN", variable_labels, 
                              scenario = scenario_num,
                              superpops_list = superpop_data_list,
                              all_scenarios_list = all_sim_scenarios, 
                              cell_ID_key, color_palette, 
                              num_synthetic = 1000, unimpaired_betas, 
                              unimpaired_cov, other_betas, other_cov, mci_betas, 
                              mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                              prior_beta, nu_0_mat, kappa_0_mat, 
                              contrasts_matrix = A, truth, seed,
                              path_to_results = 
                                paste0("/u/home/c/cshaw343/", 
                                       "link_transport_integrate/results/")))




