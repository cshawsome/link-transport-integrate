#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- all scenarios ----
HRS_samples_sizes <- c(2000, 4000, 8000)
HCAP_sampling_prop <- c(25, 50)
calibration_scenario <- c("ADAMS_prior", 
                          paste0("calibration_", c(20, 35, 50), "_SRS"), 
                          paste0("calibration_", c(20, 35, 50), "_SRS_race"))
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

#---- missing runs ----
missing_runs_function <- 
  function(missing_scenario_number, n_missing_runs, HRS_sample_size, 
           HCAP_sampling_prop, calibration_scenario, 
           sim_scenarios_dataframe = NA){
    
    batch_size = 15
    total_runs <- n_missing_runs*1/((1000 - n_missing_runs)/1000)
    
    max_batch_count <- total_runs/batch_size
    seeds <- 
      seq(1, max_batch_count*batch_size, by = batch_size)*
      missing_scenario_number*1000
    
    additional_missing_runs <- 
      as.data.frame(expand.grid(HRS_sample_size, HCAP_sampling_prop, 
                                calibration_scenario, seeds)) %>% 
      set_colnames(c("sample_size", "HCAP_prop", "calibration", "seed")) %>%
      arrange(calibration, sample_size, seed) 
    
    if(sum(is.na(sim_scenarios_dataframe)) > 0){
      additional_missing_runs %<>% 
        mutate("job" = seq(1, nrow(additional_missing_runs)))
      
      return(additional_missing_runs)
      
    } else{
      job_start <- max(sim_scenarios_missing_runs$job) + 1
      
      additional_missing_runs %<>% 
        mutate("job" = seq(job_start, 
                           job_start + nrow(additional_missing_runs) - 1))
      
      return(rbind(sim_scenarios_dataframe, additional_missing_runs))
    }
  }

#---- **1: n = 2000; HCAP_prop = 25; calibration_scenario = calibration_35_SRS ----
sim_scenarios_missing_runs <- 
  missing_runs_function(missing_scenario_number = 1, n_missing_runs = 838, 
                        HRS_sample_size = 2000, HCAP_sampling_prop = 25, 
                        calibration_scenario = 
                          c(paste0("calibration_", c(35), "_SRS")), 
                        sim_scenarios_dataframe = NA)

#---- **save output ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

write_csv(sim_scenarios_missing_runs, 
          paste0(path_to_box, "data/", "sim_study_scenarios_missing_runs.csv"))
