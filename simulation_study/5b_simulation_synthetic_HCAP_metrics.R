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
count_data <- simulated_data %>% filter(sim_name == "sim_1" & 
                                          dataset_number == "synthetic_1") %>% 
  dplyr::select(c("Unimpaired", "MCI", "Dementia", "Other")) %>% colSums() %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  set_colnames(c("Group", "true_counts")) %>% 
  mutate("sample_size" = sum(true_counts))

#---- **simulated counts ----
test <- simulated_data %>% group_by(dataset_name, sim_name, dataset_number) %>% 
  count(Group) %>% group_by



