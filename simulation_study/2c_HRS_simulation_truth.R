#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "broom")

#---- source scripts ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

superpopulations_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpopulations_data <- 
  do.call(rbind, lapply(superpopulations_paths, read_results))

#---- analytic models ----
for(dataset in unique(superpopulations_data$dataset_name)){
  subset <- superpopulations_data %>% filter(dataset_name == dataset)
  
  if(!exists("truth")){
    truth <- glm(Dementia ~ age_Z + female + black + hispanic, 
                 family = "poisson", data = subset) %>% 
      tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
      mutate("dataset_name" = dataset, 
             "RR" = exp(estimate))
  } else{
    truth <- rbind(truth,  glm(Dementia ~ age_Z + female + black + hispanic, 
                 family = "poisson", data = subset) %>% 
      tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
      mutate("dataset_name" = dataset, 
             "RR" = exp(estimate)))
  }
}

write_csv(truth, paste0(path_to_box, "analyses/simulation_study/truth.csv"))

