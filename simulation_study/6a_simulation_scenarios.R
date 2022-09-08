#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- all scenarios ----
dists <- c("normal", "lognormal", "bathtub")
samples_sizes <- c(500, 1000, 2000, 4000, 8000)
prior_props <- c("ADAMS", "unimpaired", "dementia")
calibration <- c("none", "HCAP_50")
batch_size <- 20
total_runs <- 1000

max_batch_count <- total_runs/batch_size
seeds <- seq(1, max_batch_count*batch_size, by = batch_size)

sim_scenarios <- 
  as.data.frame(expand.grid(dists, samples_sizes, prior_props, calibration, 
                            seeds)) %>% 
  set_colnames(c("distribution", "sample_size", "prior_prop", "calibration", 
                 "seed")) %>%
  arrange(prior_prop, calibration, distribution, sample_size, seed) 

sim_scenarios %<>% mutate("job" = seq(1, nrow(sim_scenarios)))

#---- **save output ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

write_csv(sim_scenarios, paste0(path_to_box, "analyses/simulation_study/", 
                                "sim_study_scenarios.csv"))

