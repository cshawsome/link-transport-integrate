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
Z <- c("AAGE", "ANMSETOT_norm", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
       "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

group <- c("Adem_dx_cat")

ADAMS_train <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS/ADAMS_train.csv")) %>% 
  mutate("group_number" = case_when(Adem_dx_cat == "Normal" ~ 1, 
                                    Adem_dx_cat == "Other" ~ 2, 
                                    Adem_dx_cat == "MCI" ~ 3, 
                                    Adem_dx_cat == "Dementia" ~ 4))

alpha_0_dist <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                  "bootstrap_cell_counts.csv")) 

#---- arrange data ----
ADAMS_train %<>% arrange(Astroke, desc(Black), desc(Hispanic))

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
  as.data.frame(table(ADAMS_train$ETHNIC_label, ADAMS_train$Astroke)) %>% 
  unite("cell", c("Var1", "Var2"), sep = ":")

B = 10000
priors_beta <- as.data.frame(matrix(nrow = (10*4*4) , ncol = B)) %>% 
  set_colnames(seq(1, B)) %>% 
  mutate("group" = rep(unique(ADAMS_train$Adem_dx_cat), each = 40), 
         "group_number" = rep(c(4, 3, 2, 1), each = 40)) 
priors_V_inv <- as.data.frame(matrix(nrow = (4*4*4), ncol = B)) %>% 
  set_colnames(seq(1, B)) %>% 
  mutate("group" = rep(unique(ADAMS_train$Adem_dx_cat), each = 16), 
         "group_number" = rep(c(4, 3, 2, 1), each = 16)) 
priors_Sigma <- as.data.frame(matrix(nrow = (10*10*4), ncol = B)) %>% 
  set_colnames(seq(1, B)) %>% 
  mutate("group" = rep(unique(ADAMS_train$Adem_dx_cat), each = 100), 
         "group_number" = rep(c(4, 3, 2, 1), each = 100)) 

for(b in 1:B){
  sample <- sample_n(ADAMS_train, size = nrow(ADAMS_train), replace = TRUE)
  
  for(group in seq(1, 4)){
    #---- filter data ----
    subset <- sample %>% filter(group_number == group) %>% 
      arrange(Astroke, desc(Black), desc(Hispanic))
    
    #---- U (contingency cell) ----
    contingency_table_temp <- 
      as.data.frame(table(subset$ETHNIC_label, subset$Astroke)) %>% 
      unite("cell", c("Var1", "Var2"), sep = ":")
    
    if(nrow(contingency_table_temp) < 6){
      contingency_table <- tibble("cells" = cells$cell, "Freq" = 0)
      contingency_table[which(contingency_table_temp$cell %in% 
                                contingency_table$cells), "Freq"] <- 
        contingency_table_temp$Freq
    } else{
      contingency_table <- contingency_table_temp
    }
    
    U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
    
    for(j in 1:nrow(contingency_table)){
      if(contingency_table[j, "Freq"] == 0){next}
      if(j == 1){
        index = 1
      } else{
        index = sum(contingency_table[1:(j - 1), "Freq"]) + 1
      }
      U[index:(index - 1 + as.numeric(contingency_table[j, "Freq"])), j] <- 1
    }
    
    UtU <- diag(unlist(contingency_table[, "Freq"]))
    
    #---- **draw new UtU if needed ----
    while(det(t(A) %*% UtU %*% A) < 1e-9){
      random_draw <- sample(seq(1, 10000), size = 1)
      new_counts <- alpha_0_dist[, c(random_draw, ncol(alpha_0_dist))] %>% 
        filter(group_number == group)
      
      UtU <- diag(unlist(new_counts[, 1]))
    }
    
    #---- beta hat ----
    continuous_covariates <- subset %>% dplyr::select(all_of(Z)) %>% as.matrix
    
    V_inv <- t(A) %*% UtU %*% A
    V <- solve(V_inv)
    priors_V_inv[which(priors_V_inv$group_number == group), b] <- 
      as.vector(V_inv)
    
    beta_hat <- V %*% t(A) %*% t(U) %*% continuous_covariates
    priors_beta[which(priors_beta$group_number == group), b] <- 
      as.vector(beta_hat)
    
    #---- Sigma hat ----
    residual <- continuous_covariates - U %*% A %*% beta_hat
    priors_Sigma[which(priors_Sigma$group_number == group), b] <- 
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


