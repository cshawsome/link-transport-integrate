#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- all scenarios ----
samples_sizes <- c(500, 1000, 2000, 4000, 8000)
calibration <- c("none", "HCAP_50")
batch_size <- 20
total_runs <- 1000

max_batch_count <- total_runs/batch_size
seeds <- seq(1, max_batch_count*batch_size, by = batch_size)

sim_scenarios <- 
  as.data.frame(expand.grid(samples_sizes, calibration, seeds)) %>% 
  set_colnames(c("sample_size", "calibration", "seed")) %>%
  arrange(calibration, sample_size, seed) 

sim_scenarios %<>% mutate("job" = seq(1, nrow(sim_scenarios)))

#---- **save output ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

write_csv(sim_scenarios, 
          paste0(path_to_box, "data/", "sim_study_scenarios.csv"))

