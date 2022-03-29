#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "broom")

#---- source scripts ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/synthetic_data/HRS"), 
             full.names = TRUE, pattern = "*.csv")

synthetic_data <- do.call(rbind, lapply(synthetic_data_paths, read_results))

#---- analytic models ----
for(dataset in unique(synthetic_data$dataset_name)){
  subset <- synthetic_data %>% filter(dataset_name == dataset)
  
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



