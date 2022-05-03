#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "LaplacesDemon")

#---- source scripts ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **superpopulations ----
superpopulations_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpopulations_data <- 
  do.call(rbind, lapply(superpopulations_paths, read_results)) 

#---- **HRS synthetic datasets ----
HRS_synthetic_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/OLD/synthetic_data/HRS"), 
             full.names = TRUE, pattern = "*.csv")

HRS_synthetic_data <- 
  do.call(rbind, lapply(HRS_synthetic_paths, read_results)) 

#---- draws from superpopulation ----
for(sample_size in c(500, 1000, 2000, 4000, 8000)){
  if(!exists("HRS_draws")){
    HRS_draws <- superpopulations_data %>% group_by(dataset_name) %>% 
      sample_n(size = sample_size, replace = TRUE) %>% 
      mutate_at("dataset_name", 
                function(x) str_replace(x, "1000000", paste0(sample_size)))
  } else{
    HRS_draws %<>% 
      rbind(., superpopulations_data %>% group_by(dataset_name) %>% 
              sample_n(size = sample_size, replace = TRUE) %>% 
              mutate_at("dataset_name", 
                        function(x) str_replace(x, "1000000", 
                                                paste0(sample_size))))
  }
}

#---- age variability plots ----
#---- **superpopulations ----
ggplot(data = superpopulations_data, aes(x = age_Z)) + 
  geom_histogram() + theme_bw() + 
  facet_wrap(vars(dataset_name), ncol = 3, scales = "free")

ggsave(paste0(path_to_box, "figures/simulation_study/extra_analyses/", 
              "superpopulations_age_dist.jpeg"))

#---- **superpopulation draws ----
ggplot(data = HRS_draws, aes(x = age_Z)) + 
  geom_histogram() + theme_bw() + 
  facet_wrap(vars(dataset_name), ncol = 3, scales = "free")

ggsave(paste0(path_to_box, "figures/simulation_study/extra_analyses/", 
              "HRS_draws_age_dist.jpeg"), width = 11, height = 14, 
       units = "in")

#---- **synthetic HRS ----
ggplot(data = HRS_synthetic_data, aes(x = age_Z)) + 
  geom_histogram() + theme_bw() + 
  facet_wrap(vars(dataset_name), ncol = 3, scales = "free")

ggsave(paste0(path_to_box, "figures/simulation_study/extra_analyses/", 
              "HRS_synthetic_age_dist.jpeg"), width = 11, height = 14, 
       units = "in")

#---- continuous vars summary stats ----
#These are very similar
HRS_synthetic_data_sds <- HRS_synthetic_data %>% 
  dplyr::select(contains("Z")) %>% summarize_all(sd)

HRS_draws_sds <- HRS_draws %>% ungroup() %>% 
  dplyr::select(names(HRS_synthetic_data_sds)) %>% summarize_all(sd) %>% 
  unlist()
  
