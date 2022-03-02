#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "LaplacesDemon")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv")) 

#---- **continuous distribution parameters ----
normal_parameter_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/", 
                 "continuous_distribution_parameters/", 
                 "normal_parameters"))

#---- source functions ----
source(here("simulation_study", "functions", "generate_synthetic_continuous.R"))

#---- synthetic data ----
#---- **multivariate normal ----


for(ID in cell_IDs){
  HRS_analytic
}

ID = "000"
start_1 <- Sys.time()
test_1 <- matrixNormal::rmatnorm(s = 1, 
                               M = as.matrix(normal_parameter_list[[ID]]$M), 
                               U = as.matrix(normal_parameter_list[[ID]]$U),
                               V = as.matrix(normal_parameter_list[[ID]]$V))
end_1 <- Sys.time() - start_1

start_2 <- Sys.time()
test_2 <- LaplacesDemon::rmatrixnorm(M = as.matrix(normal_parameter_list[[ID]]$M), 
                                     U = as.matrix(normal_parameter_list[[ID]]$U),
                                     V = as.matrix(normal_parameter_list[[ID]]$V))
end_2 <- Sys.time() - start_2



#---- predicted impairment ----
