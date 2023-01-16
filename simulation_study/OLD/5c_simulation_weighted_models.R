#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "broom", "mice")

#---- source scripts ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **synthetic true HRS----
simulation_truth <-  
  read_csv(paste0(path_to_box, "analyses/simulation_study/synthetic_data/HRS/", 
                  "HRS_synthetic_normal_500_unimpaired.csv")) 

#---- **synthetic HCAP ----
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
    simulated_data <- lapply(simulated_data_paths, read_results)
  } else{
    simulated_data %<>% append(., lapply(simulated_data_paths, read_results))
  }
}

#---- ****add weights ----
#the only variable that contributed to selection was marital status
# we selected half of married people and half of unmarried people, thus the 
# weight for every observation should be 2

simulated_data <- 
  lapply(simulated_data, function(x) x %<>% mutate("weight" = 2))

#---- ****predicted dementia status ----
simulated_data <- 
  lapply(simulated_data, function(x) x %<>% 
           mutate("predicted_Dementia" = ifelse(Group == "Dementia", 1, 0)))

#---- models ----
models <- lapply(simulated_data, 
                 function(dataset) glm(predicted_Dementia ~ 
                                         age_Z + female + black + hispanic, 
                                       family = "poisson", data = dataset))

#---- **pooled ----
pooled_model <- summary(mice::pool(models[1:10]))


