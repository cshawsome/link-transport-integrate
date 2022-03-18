#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "MCMCpack")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **prior imputed clean ----
prior_imputed_clean <- 
  readRDS(paste0(path_to_box, 
                 "data/ADAMS/prior_data/MI/MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(prior_imputed_clean[[1]]))

#---- relabel columns ----
prior_imputed_clean <- 
  lapply(prior_imputed_clean, 
         function(x) rename_at(x, vars(variable_labels$ADAMS), ~ 
                                 variable_labels$data_label)) 

#---- **variable selection results ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv")) %>% 
  dplyr::select("data_label") %>% unlist()

#---- ****classify vars ----
#notation from Schafer 1997
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars, "Z")]

#---- **prior: cell counts ----
alpha_0_dist <- 
  read_csv(paste0(path_to_box, "data/ADAMS/prior_data/", 
                  "imputation_cell_props.csv")) 

#---- define A (contrasts) ----
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

cells <- A %>% as.data.frame() %>% unite("cells", -1, sep = "")

#---- estimates ----
test <- prior_imputed_clean[[1]]

estimate_cont_priors <- function(data, W, Z){
  for(group in c("Unimpaired", "MCI", "Other", "Dementia")){
    #---- **filter data ----
    subset <- data %>% filter(!!sym(group) == 1)
    
    #---- U (contingency cell) ----
    contingency_table_temp <- subset %>% 
      unite("cell_ID", all_of(W), sep = "") %>% dplyr::select(cell_ID) %>% 
      table() %>% as.data.frame() %>% set_colnames(c("cell", "Freq"))
  
    if(nrow(contingency_table_temp) < 6){
      contingency_table <- data.frame("cell" = cells$cell) %>% 
        left_join(., contingency_table_temp)
      
      contingency_table[is.na(contingency_table)] <- 0
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
      
      UtU <- diag(round(unlist(new_counts[, 1])*nrow(subset)))
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



#---- OLD ----

group <- c("Adem_dx_cat")

ADAMS_train <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/ADAMS/cleaned/ADAMS_train.csv")) %>% 
  mutate("group_number" = case_when(Adem_dx_cat == "Unimpaired" ~ 1, 
                                    Adem_dx_cat == "Other" ~ 2, 
                                    Adem_dx_cat == "MCI" ~ 3, 
                                    Adem_dx_cat == "Dementia" ~ 4))



#---- arrange data ----
ADAMS_train %<>% arrange(Astroke, desc(Black), desc(Hispanic))



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
  
  
}

#---- save ----
write_csv(priors_beta, paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                              "priors/ADAMS_test/priors_beta.csv"))
write_csv(priors_V_inv, paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/ADAMS_test/priors_V_inv.csv"))
write_csv(priors_Sigma, paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/ADAMS_test/priors_Sigma.csv"))
