#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "locfit")

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character()))

#---- **models for priors ----
normal_prior <- readRDS(here::here("priors", "normal_model_25.rds"))
other_prior <- readRDS(here::here("priors", "other_model_25.rds"))
MCI_prior <- readRDS(here::here("priors", "MCI_model_25.rds"))

normal_preds <- names(coefficients(normal_prior))
normal_preds[which(normal_preds == "ETHNIC_labelBlack")] <- "Black"
normal_preds[which(normal_preds == "ETHNIC_labelHispanic")] <- "Hispanic"
other_preds <- names(coefficients(other_prior))
MCI_preds <- names(coefficients(MCI_prior))

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- c("AAGE", "ETHNIC_label", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", 
          "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Astroke", "Abmi")

analytical_sample <- ADAMS_subset %>% dplyr::select("HHIDPN", all_of(vars)) %>% 
  mutate("Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0)) %>% 
  #use complete data for now
  na.omit() %>% 
  #add intercept
  mutate("(Intercept)" = 1) %>%
  #pre-allocate columns
  mutate("Group" = 0, "p_Unimpaired" = 0, "p_Other" = 0, "p_MCI" = 0)

# #Sanity check
# table(analytical_sample$ETHNIC_label, analytical_sample$Black, useNA = "ifany")
# table(analytical_sample$ETHNIC_label, analytical_sample$Hispanic,
#       useNA = "ifany")
# 
# colSums(is.na(analytical_sample))

#---- all-way contingency table ----
#We have one small cell-- Hispanics who have had a stroke
cross_class_label <- table(analytical_sample$ETHNIC_label, 
                           analytical_sample$Astroke) %>% as.data.frame()

#---- Bayes Stuff ----
#---- **parameters ----
#number of runs
B = 2

#---- **chain storage ----
normal_beta_chain <- 
  matrix(nrow = length(coefficients(normal_prior)) , ncol = B)
other_beta_chain <- 
  matrix(nrow = length(coefficients(other_prior)) , ncol = B)
MCI_beta_chain <- 
  matrix(nrow = length(coefficients(MCI_prior)) , ncol = B)

latent_class_chain <- matrix(nrow = 4, ncol = B)

#---- **priors ----

#---- **sampling ----
for(b in 1:B){
  #---- ****latent class betas ----
  normal_beta_chain[, b] <- mvrnorm(n = 1, mu = coefficients(normal_prior), 
                                    Sigma = vcov(normal_prior))
  other_beta_chain[, b] <- mvrnorm(n = 1, mu = coefficients(other_prior), 
                                   Sigma = vcov(other_prior))
  MCI_beta_chain[, b] <- mvrnorm(n = 1, mu = coefficients(MCI_prior), 
                                 Sigma = vcov(MCI_prior))
  
  #---- ****group: Unimpaired vs. Impaired ----
  analytical_sample[, "p_Unimpaired"] <- 
    expit(as.matrix(analytical_sample[, normal_preds]) %*% 
            as.matrix(normal_beta_chain[, b]))
  
  analytical_sample[, "Group"] <- 
    rbernoulli(n = nrow(analytical_sample), 
               p = analytical_sample[, "p_Unimpaired"])*1
  
  #---- ****group: Other vs. MCI/Dementia ----
  subset_index <- which(analytical_sample$Group == 0)
  
  analytical_sample[subset_index, "p_Other"] <- 
    expit(as.matrix(analytical_sample[subset_index, other_preds]) %*% 
            as.matrix(other_beta_chain[, b]))
  
  analytical_sample[subset_index, "Group"] <- 
    rbernoulli(n = length(subset_index), 
               p = analytical_sample[subset_index, "p_Other"])*2
  
  #---- ****group: MCI vs. Dementia ----
  
  
  #---- ****group: summary ----
}



