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
  readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/MI/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(prior_imputed_clean[[1]]))

#---- ****relabel columns ----
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
  readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                 "imputation_cell_props")) 

#---- **contrasts matrix ----
#categorical vars contrasts matrix
A = read_csv(paste0(path_to_box, "analyses/contrasts_matrix.csv")) %>% 
  as.matrix()

cells <- A %>% as.data.frame() %>% unite("cells", -1, sep = "") %>% 
  dplyr::select(-"Intercept") %>% table() %>% as.data.frame() %>% 
  dplyr::select(-"Freq") %>% set_colnames("cells")

#---- estimates ----
estimate_cont_priors <- function(data, W, Z, A){
  #---- **pre-allocate ----
  priors_V_inv <- 
    as.data.frame(matrix(nrow = (ncol(A)*4), ncol = ncol(A))) %>%
    mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                         each = ncol(A))) 
  
  priors_beta <- 
    as.data.frame(matrix(nrow = (ncol(A)*4) , ncol = length(Z))) %>% 
    mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                         each = ncol(A)))
  
  priors_Sigma <- 
    as.data.frame(matrix(nrow = (length(Z)*4), ncol = length(Z))) %>% 
    mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                         each = length(Z))) 
  
  for(class in c("Unimpaired", "MCI", "Other", "Dementia")){
    #---- **filter data ----
    subset <- data %>% filter(!!sym(class) == 1) 
    
    #---- **U (contingency cell) ----
    contingency_table_temp <- subset %>% 
      unite("cell_ID", all_of(W), sep = "") %>% dplyr::select(cell_ID) %>% 
      table() %>% as.data.frame() %>% set_colnames(c("cells", "Freq")) 
    
    if(nrow(contingency_table_temp) < 6){
      contingency_table <- data.frame("cells" = cells$cell) %>% 
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
    
    #---- NEEDS TO BE REFACTORED ----
    # #---- **draw new UtU if needed ----
    # while(det(t(A) %*% UtU %*% A) < 1e-9){
    #   max_index <- 
    #     colnames(alpha_0_dist)[str_detect(colnames(alpha_0_dist), "[0-9]+")] %>% 
    #     as.numeric() %>% max()
    #   random_draw <- sample(seq(1, max_index), size = 1)
    #   new_props <- alpha_0_dist[, c(as.character(random_draw), "group")] %>% 
    #     filter(group == class)
    #   
    #   UtU <- diag(round(unlist(new_counts[, 1])*nrow(subset)))
    # }
    
    #---- **beta hat ----
    continuous_covariates <- subset %>% dplyr::select(all_of(Z)) %>% as.matrix
    
    V_inv <- t(A) %*% UtU %*% A
    V <- solve(V_inv)
    priors_V_inv[which(priors_V_inv$group == class), 
                 seq(1, ncol(V_inv))] <- V_inv
    
    beta_hat <- V %*% t(A) %*% t(U) %*% continuous_covariates
    priors_beta[which(priors_beta$group == class), seq(1, ncol(beta_hat))] <- 
      beta_hat
    
    #---- **Sigma hat ----
    residual <- continuous_covariates - U %*% A %*% beta_hat
    priors_Sigma[which(priors_Sigma$group == class), seq(1, ncol(residual))] <- 
      t(residual) %*% residual
  }
  
  #---- **return values ----
  return(list("priors_V_inv" = priors_V_inv %<>% group_by(group) %>% 
                group_split() %>% 
                set_names(c("Dementia", "MCI", "Other", "Unimpaired")), 
              "priors_beta" = priors_beta %<>% group_by(group) %>% 
                group_split() %>% 
                set_names(c("Dementia", "MCI", "Other", "Unimpaired")), 
              "priors_Sigma" = priors_Sigma %<>% group_by(group) %>% 
                group_split() %>% 
                set_names(c("Dementia", "MCI", "Other", "Unimpaired"))))
}

#---- **estimate values ----
all_priors <- lapply(prior_imputed_clean, estimate_cont_priors, W, Z, A)

#---- format output ----
for(est in c("V_inv", "beta", "Sigma")){
  
  data <- lapply(all_priors, "[[", paste0("priors_", est)) 
  
  saveRDS(data, paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                       "priors_", est))
}

#---- OLD ----

#---- **pre-allocate ----
priors_V_inv <- 
  as.data.frame(matrix(nrow = (ncol(A)*ncol(A)*4), 
                       ncol = length(prior_imputed_clean))) %>% 
  set_colnames(seq(1:length(prior_imputed_clean))) %>% 
  mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                       each = ncol(A)*ncol(A))) 

priors_beta <- 
  as.data.frame(matrix(nrow = (length(Z)*ncol(A)*4) , 
                       ncol = length(prior_imputed_clean))) %>% 
  set_colnames(seq(1:length(prior_imputed_clean))) %>% 
  mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                       each = length(Z)*ncol(A)))

priors_Sigma <- 
  as.data.frame(matrix(nrow = (length(Z)*length(Z)*4), 
                       ncol = length(prior_imputed_clean))) %>% 
  set_colnames(seq(1:length(prior_imputed_clean))) %>% 
  mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                       each = length(Z)*length(Z))) 

for(index in 1:length(all_priors)){
  list <- all_priors[[index]]
  priors_V_inv[which(priors_V_inv$group == list[[1]]$group), index] <- 
    list[[1]]$V_inv
  priors_beta[which(priors_beta$group == list[[2]]$group), index] <- 
    list[[2]]$beta
  priors_Sigma[which(priors_Sigma$group == list[[3]]$group), index] <- 
    list[[3]]$Sigma
}

#---- **save results ----
write_csv(priors_V_inv, paste0(path_to_box, "analyses/simulation_study/", 
                               "prior_data/priors_V_inv.csv"))
write_csv(priors_beta, paste0(path_to_box, "analyses/simulation_study/", 
                              "prior_data/priors_beta.csv"))
write_csv(priors_Sigma, paste0(path_to_box, "analyses/simulation_study/", 
                               "prior_data/priors_Sigma.csv"))
