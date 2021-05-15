#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS/ADAMS_train.csv"))

#---- predictors ----
#normal model predictors
norm_preds <- c("AAGE", "Black", "Hispanic", "ANMSETOT", "ANSER7T", "ANIMMCR", 
                "ANRECYES", "ANWM1TOT", "proxy_cog")

#other model predictors
other_preds <- c("AAGE", "ANMSETOT", "ANIMMCR", "ANDELCOR")

#mci model predictors
MCI_preds <- c("ANMSETOT", "ANIMMCR", "Aiadla", "Astroke", "Abmi")

#---- model ----
bootstrap_models <- function(){
  subsample <- sample_frac(ADAMS_subset, size = 0.5, replace = TRUE)
  
  normal_model <- 
    glm(formula(paste("ANormal ~ ", paste(norm_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample)
  
  other_model <- 
    glm(formula(paste("AOther ~ ", paste(other_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(ANormal == 0))
  
  mci_model <- 
    glm(formula(paste("AMCI ~ ", paste(MCI_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(ANormal == 0 & AOther == 0)) 
  
  return(list("normal_betas" = coefficients(normal_model),
              "other_betas" = coefficients(other_model), 
              "mci_betas" = coefficients(mci_model), 
              "normal_cov" = as.vector(vcov(normal_model)), 
              "other_cov" = as.vector(vcov(other_model)), 
              "mci_cov" = as.vector(vcov(mci_model))))
}

bootstrap_runs <- replicate(2, bootstrap_models(), simplify = FALSE)

#---- format output ----




