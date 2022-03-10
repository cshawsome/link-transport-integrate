#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "broom")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_data <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/synthetic_data/", 
                  "synthetic_normal_500_unimpaired.csv"))

#---- analytic models ----
truth <- glm(Dementia ~ age_Z + female + edyrs_Z + black + hispanic, 
             family = "poisson", data = synthetic_data) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = TRUE)

write_csv(truth, paste0(path_to_box, "analyses/simulation_study/truth/", 
                        "truth_normal_500_unimpaired.csv"))



