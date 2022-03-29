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

#---- coverage: impairment class ----
#---- **true counts ----
true_counts <- simulated_data %>% filter(sim_name == "sim_1" & 
                                           dataset_number == "synthetic_1") %>% 
  dplyr::select(c("Unimpaired", "MCI", "Dementia", "Other")) %>% colSums() %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  set_colnames(c("Group", "true_counts")) %>% 
  mutate("sample_size" = sum(true_counts))

#---- **simulated counts ----
simulated_counts <- simulated_data %>% 
  group_by(dataset_name, sim_name, dataset_number) %>% 
  count(Group) %>% ungroup() %>% group_by(dataset_name, sim_name, Group) %>% 
  summarise_at("n", list("LB" = function(x) quantile(x, 0.25), 
                         "UB" = function(x) quantile(x, 0.975))) %>% 
  left_join(., count_data) %>% 
  mutate("truth_capture" = 
           ifelse(LB <= true_counts & UB >= true_counts, 1, 0)) %>% 
  group_by(sample_size, Group) %>% summarize_at("truth_capture", mean) %>% 
  rename("Normal" = "truth_capture")

#---- **formatted table ----
coverage_table <- left_join(true_counts, simulated_counts)

write_csv(coverage_table, 
          paste0(path_to_box, "analyses/simulation_study/results/", 
                 "HCAP_metrics_coverage_table.csv"))

#---- bias: impairment class ----


