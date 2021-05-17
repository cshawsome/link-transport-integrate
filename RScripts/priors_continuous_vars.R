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
                  "data/cleaned/ADAMS/ADAMS_train.csv")) 

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
cells <- 
  as.data.frame(table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke)) %>% 
  unite("cell", c("Var1", "Var2"), sep = ":")

B = 10000
priors_beta <- as.data.frame(matrix(nrow = (10*4*4) , ncol = B)) %>% 
  set_colnames(seq(1, B)) %>% 
  mutate("group" = rep(unique(ADAMS_subset$Adem_dx_cat), each = 40)) 
priors_V_inv <- as.data.frame(matrix(nrow = (4*4*4), ncol = B)) %>% 
  set_colnames(seq(1, B)) %>% 
  mutate("group" = rep(unique(ADAMS_subset$Adem_dx_cat), each = 16)) 
priors_Sigma <- as.data.frame(matrix(nrow = (10*10*4), ncol = B)) %>% 
  set_colnames(seq(1, B)) %>% 
  mutate("group" = rep(unique(ADAMS_subset$Adem_dx_cat), each = 100)) 

for(b in 1:B){
  sample <- sample_n(ADAMS_subset, size = nrow(ADAMS_subset), replace = TRUE)
  
  for(group in unique(ADAMS_subset$Adem_dx_cat)){
    #---- filter data ----
    subset <- sample %>% filter(Adem_dx_cat == group) %>% 
      arrange(Astroke, desc(Black), desc(Hispanic))
    
    #---- U (contingency cell) ----
    contingency_table <- 
      as.data.frame(table(subset$ETHNIC_label, subset$Astroke)) 
    
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
    
    V_inv <- t(A) %*% UtU %*% A
    V <- solve(V_inv)
    priors_V_inv[which(priors_V_inv$group == group), b] <- as.vector(V_inv)
    
    beta_hat <- V %*% t(A) %*% t(U) %*% continuous_covariates
    priors_beta[which(priors_beta$group == group), b] <- as.vector(beta_hat)
    
    residual <- continuous_covariates - U %*% A %*% beta_hat
    priors_Sigma[which(priors_Sigma$group == group), b] <- 
      as.vector(t(residual) %*% residual)
  }
}

#---- save ----
write_csv(priors_beta, paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                              "priors/priors_beta.csv"))
write_csv(priors_V_inv, paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/priors_V_inv.csv"))
write_csv(priors_Sigma, paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/priors_Sigma.csv"))


