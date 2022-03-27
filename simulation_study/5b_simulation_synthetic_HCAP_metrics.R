#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- source scripts ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

simulated_data_directories <- 
  list.dirs(path = paste0(path_to_box, 
                          "analyses/simulation_study/HCAP_normal_250_unimpaired/", 
                          "synthetic_data"), full.names = TRUE, 
            recursive = FALSE)

simulated_data_directories <- 
  simulated_data_directories[str_detect(simulated_data_directories, 
                                        pattern = "sim_")]

for(path in simulated_data_directories){
  simulated_data_paths <- 
    list.files(path = path, full.names = TRUE, pattern = "*.csv")
  
  if(!exists("simulated_data")){
    simulated_data <- do.call(rbind, lapply(simulated_data_paths, read_results))
  } else{
    simulated_data %<>% 
      rbind(., do.call(rbind, lapply(simulated_data_paths, read_results)))
  }
}

#---- impairment class coverage ----
#---- **true counts ----
test <- 



