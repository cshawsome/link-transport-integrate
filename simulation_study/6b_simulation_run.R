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
library("tidyverse") 
library("DirichletReg")
library("magrittr")
library("MCMCpack") 
library("locfit")
library("vroom")
library("mvnfast")
library("mice")
library("LaplacesDemon")

#---- source functions ----
path_to_RScripts <- "/u/home/c/cshaw343/link_transport_integrate/RScripts/"

source(paste0(path_to_RScripts, "read_results.R"))
source(paste0(path_to_RScripts, "generate_synthetic_function.R"))
source(paste0(path_to_RScripts, "simulation_function.R"))

#---- read in data ----
path_to_data <- "/u/home/c/cshaw343/link_transport_integrate/data/"

#---- **data paths ----
superpop <- 
  read_results(paste0(path_to_data, "superpopulations/superpop_1000000.csv"))

#---- **truth table ----
truth <- read_csv(paste0(path_to_data, 
                         "superpopulations/agesex_standardized_prevs.csv"))

#---- **superpop means and sds ----
means <- 
  read_csv(paste0(path_to_box, "data/superpopulations/superpop_means.csv"))

sds <- 
  read_csv(paste0(path_to_box, "data/superpopulations/superpop_sds.csv"))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_data, "variable_crosswalk.csv")) 

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
Z <- colnames(superpop)[str_detect(colnames(superpop), "_Z")]

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_data, "contrasts_matrix.csv")) %>% as.matrix()

#---- **hyperparameters (tune these) ----
#DOF for inverse wishart
nu_0_mat <- read_csv(paste0(path_to_data, "tuning/nu_0_matrix.csv"))

#scaling for inverse wishart as variance of Beta
kappa_0_mat <- read_csv(paste0(path_to_data, "tuning/kappa_0_matrix.csv"))

#---- set seed ----
seed <- as.numeric(all_sim_scenarios[scenario_num, "seed"])
set.seed(seed)

#---- run sim ----
replicate(num_replicates,
          simulation_function(warm_up = 100, starting_props = rep(0.25, 4), 
                              categorical_vars = W, continuous_vars = Z, 
                              id_var = "HHIDPN", 
                              variable_labels = variable_labels, 
                              scenario = scenario_num,
                              superpopulation = superpop, orig_means = means, 
                              orig_sds = sds, 
                              all_scenarios_list = all_sim_scenarios, 
                              cell_ID_key = cell_ID_key, 
                              color_palette = color_palette, 
                              num_synthetic = 1000, contrasts_matrix = A,
                              kappa_0_mat = kappa_0_mat, nu_0_mat = nu_0_mat,
                              truth = truth, seed = seed, 
                              path_to_raw_prior_sample = 
                                paste0(path_to_data, 
                                       "prior_data/MI/MI_datasets_cleaned"), 
                              path_to_data = path_to_data, 
                              path_to_results = 
                                paste0("/u/home/c/cshaw343/", 
                                       "link_transport_integrate/results/")))
