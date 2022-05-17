#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- all scenarios ----
dists <- c("normal", "lognormal", "bathtub")
samples_sizes <- c(500, 1000, 2000, 4000, 8000)
prior_props <- c("ADAMS", "unimpaired", "dementia")
seeds <- c()

sim_scenarios <- 
  as.data.frame(expand.grid(dists, samples_sizes, prior_props)) %>% 
  set_colnames(c("distribution", "sample_size", "prior_prop"))

#---- **save output ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

write_csv(sim_scenarios, paste0(path_to_box, "analyses/simulation_study/", 
                                "sim_study_scenarios.csv"))

