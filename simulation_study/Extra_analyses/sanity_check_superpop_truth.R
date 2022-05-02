#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

ADAMS_analytic <-  
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv")) %>% 
  mutate("age_Z" = scale(AAGE))


#---- ADAMS analysis ----
glm(Dementia ~ age_Z + Female + Black + Hispanic, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))

glm(Dementia ~ age_Z + Female + Black + Hispanic + ANMSETOT, 
    family = "poisson", data = ADAMS_analytic) %>% 
  tidy(., conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE) %>% 
  mutate("RR" = exp(estimate))
