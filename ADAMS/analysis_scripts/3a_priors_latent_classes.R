#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/ADAMS/cleaned/ADAMS_train.csv"))

#---- predictors ----
#unimpaired model predictors
unimpaired_preds <- c("AAGE", "Black", "Hispanic", "ANMSETOT_norm", "ANSER7T", 
                  "ANIMMCR", "ANRECYES", "ANWM1TOT", "proxy_cog")

#other model predictors
other_preds <- c("AAGE", "ANMSETOT_norm", "ANIMMCR", "ANDELCOR")

#mci model predictors
mci_preds <- c("ANMSETOT_norm", "ANIMMCR", "Aiadla", "Astroke", "Abmi")

#---- model ----
bootstrap_models <- function(prop){
  subsample <- sample_frac(ADAMS_subset, size = prop, replace = TRUE)
  
  unimpaired_model <- 
    glm(formula(paste("AUnimpaired ~ ", paste(unimpaired_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample)
  
  other_model <- 
    glm(formula(paste("AOther ~ ", paste(other_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(AUnimpaired == 0))
  
  mci_model <- 
    glm(formula(paste("AMCI ~ ", paste(mci_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(AUnimpaired == 0 & AOther == 0)) 
  
  return(list("unimpaired_betas" = coefficients(unimpaired_model),
              "other_betas" = coefficients(other_model), 
              "mci_betas" = coefficients(mci_model), 
              "unimpaired_cov" = as.vector(vcov(unimpaired_model)), 
              "other_cov" = as.vector(vcov(other_model)), 
              "mci_cov" = as.vector(vcov(mci_model))))
}

bootstrap_runs <- replicate(10000, bootstrap_models(prop = 1), 
                            simplify = FALSE)

#---- check distributions ----
for(group in c("unimpaired", "other", "mci")){
  data <- lapply(bootstrap_runs, "[[", paste0(group, "_betas")) %>% 
    do.call(rbind, .) %>% t() %>% as.data.frame() 
  
  for(var in rownames(data)){
    show(hist(as.numeric(data[var, ]), main = paste0(group, " ", var)))
  }
}

#---- format output ----
for(est in c("betas", "cov")){
  for(group in c("unimpaired", "other", "mci")){
    data <- lapply(bootstrap_runs, "[[", paste0(group, "_", est)) %>% 
      do.call(rbind, .) %>% t() %>% as.data.frame() 
    
    if(est == "betas"){
      data %<>% mutate("preds" = c("(Intercept)", get(paste0(group, "_preds"))))
    }
    
    data %>% 
      write_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                       "latent_class_", group, "_", est, ".csv"))
  }
}




