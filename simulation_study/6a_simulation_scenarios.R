#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- all scenarios ----
HRS_samples_sizes <- c(2000, 4000, 8000)
HCAP_sampling_prop <- c(25, 50)
calibration_scenario <- c("ADAMS_prior", "calibration_50_SRS", 
                          "calibration_50_design")
batch_size <- 15
total_runs <- 1000

max_batch_count <- total_runs/batch_size
seeds <- seq(1, max_batch_count*batch_size, by = batch_size)

sim_scenarios <- 
  as.data.frame(expand.grid(HRS_samples_sizes, HCAP_sampling_prop, 
                            calibration_scenario, seeds)) %>% 
  set_colnames(c("sample_size", "HCAP_prop", "calibration", "seed")) %>%
  arrange(calibration, sample_size, seed) 

sim_scenarios %<>% mutate("job" = seq(1, nrow(sim_scenarios)))

#---- **save output ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

write_csv(sim_scenarios, 
          paste0(path_to_box, "data/", "sim_study_scenarios.csv"))

