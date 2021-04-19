#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "MCMCpack")

options(scipen = 999)

#---- read in data ----
#Categorical vars (notation from Schafer 1997)
W <- c("ETHNIC_label", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
       "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

group <- c("Adem_dx_cat")

ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv")) %>% 
  dplyr::select(c("HHIDPN", all_of(group), all_of(W), all_of(Z))) %>% 
  na.omit() %>% 
  #don't standardize this
  mutate_at("Astroke", as.character) %>%
  #Z-score continuous
  mutate_if(is.numeric, scale) %>%
  #transform to correct type
  mutate_at("Astroke", as.numeric) %>%
  mutate("Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0),
         #Add intercept
         "(Intercept)" = 1) %>% 
  mutate("group_class" = 
           case_when(Adem_dx_cat %in% 
                       c("Dementia", "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia",
                     Adem_dx_cat == "Normal" ~ "Unimpaired",
                     TRUE ~ Adem_dx_cat))

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$group_class, useNA = "ifany")

#---- arrange data ----
ADAMS_subset %<>% arrange(Astroke, desc(Black), desc(Hispanic))

#---- A (contrasts) ----
#categorical vars contrasts matrix
A = do.call(cbind, list(
  #intercept
  rep(1, 6),
  #race/ethnicity main effect: Black
  rep(c(1, 0, 0), 2),
  #race/ethnicity main effect: Hispanic
  rep(c(0, 1, 0), 2),
  #stroke main effect
  rep(c(0, 1), each = 3)))

#---- pre-allocate ----
priors_beta <- matrix(nrow = (10*4) , ncol = 4)
i = 1

for(group in c("Unimpaired", "Other", "MCI", "Dementia")){
  #---- filter data ----
  subset <- ADAMS_subset %>% filter(group_class == group)
  
  #---- U (contingency cell) ----
  contingency_table <- table(subset$ETHNIC_label, 
                             subset$Astroke) %>% as.data.frame()
  
  U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
  
  for(j in 1:nrow(contingency_table)){
    if(contingency_table[j, "Freq"] == 0){next}
    if(j == 1){
      index = 1
    } else{
      index = sum(contingency_table[1:(j - 1), "Freq"]) + 1
    }
    U[index:(index - 1 + contingency_table[j, "Freq"]), j] <- 1
  }
  
  UtU <- diag(contingency_table[, "Freq"])
  
  #---- beta hat ----
  continuous_covariates <- subset %>% dplyr::select(all_of(Z)) %>% as.matrix
  
  V <- solve(t(A) %*% UtU %*% A)
  
  beta_hat <- V %*% t(A) %*% t(U) %*% continuous_covariates
  priors_beta[, i] <- as.vector(beta_hat)
  
  if(group == "Unimpaired"){
    #---- Sigma prior (based on normal group) ----
    eps_hat <- continuous_covariates - (U %*% A %*% beta_hat)
    Sigma_hat <- (t(eps_hat) %*% eps_hat)/(nrow(ADAMS_subset) - ncol(beta_hat))
  }
  i = i + 1
}

#---- save ----
saveRDS(Sigma_hat, file = here::here("priors", "Sigma.rds"))
saveRDS(priors_beta, file = here::here("priors", "beta.rds"))
  

