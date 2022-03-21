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
  rep(c(0, 1, 0), 2),
  #race/ethnicity main effect: Hispanic
  rep(c(0, 0, 1), 2),
  #stroke main effect
  rep(c(0, 1), each = 3)))

cells <- A %>% as.data.frame() %>% unite("cells", -1, sep = "") %>% 
  dplyr::select(-"V1") %>% table() %>% as.data.frame() %>% 
  dplyr::select(-"Freq") %>% set_colnames("cells")

#---- estimates ----
test <- prior_imputed_clean[[1]]

estimate_cont_priors <- function(data, W, Z, A){
  #---- **pre-allocate ----
  priors_V_inv <- 
    as.data.frame(matrix(nrow = (ncol(A)*ncol(A)*4), ncol = 1)) %>% 
    set_colnames("V_inv") %>% 
    mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                         each = ncol(A)*ncol(A))) 
  
  priors_beta <- 
    as.data.frame(matrix(nrow = (length(Z)*ncol(A)*4) , ncol = 1)) %>% 
    set_colnames("beta") %>% 
    mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                         each = length(Z)*ncol(A)))
  
  priors_Sigma <- 
    as.data.frame(matrix(nrow = (length(Z)*length(Z)*4), ncol = 1)) %>% 
    set_colnames("Sigma") %>% 
    mutate("group" = rep(c("Unimpaired", "MCI", "Dementia", "Other"), 
                         each = length(Z)*length(Z))) 
    
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
    
    #---- **draw new UtU if needed ----
    while(det(t(A) %*% UtU %*% A) < 1e-9){
      random_draw <- sample(seq(1, (ncol(alpha_0_dist) - 2)), size = 1)
      new_props <- alpha_0_dist[, c(random_draw, "group")] %>% 
        filter(group == class)
      
      UtU <- diag(round(unlist(new_counts[, 1])*nrow(subset)))
    }
    
    #---- **beta hat ----
    continuous_covariates <- subset %>% dplyr::select(all_of(Z)) %>% as.matrix
    
    V_inv <- t(A) %*% UtU %*% A
    V <- solve(V_inv)
    priors_V_inv[which(priors_V_inv$group == class), "V_inv"] <- 
      as.vector(V_inv)
    
    beta_hat <- V %*% t(A) %*% t(U) %*% continuous_covariates
    priors_beta[which(priors_beta$group == class), "beta"] <- 
      as.vector(beta_hat)
    
    #---- **Sigma hat ----
    residual <- continuous_covariates - U %*% A %*% beta_hat
    priors_Sigma[which(priors_Sigma$group == class), "Sigma"] <- 
      as.vector(t(residual) %*% residual)
  }
  
  #---- **return values ----
  return(list(priors_V_inv, priors_beta, priors_Sigma))
}

#---- **estimate values ----
all_priors <- lapply(prior_imputed_clean, estimate_cont_priors, W, Z, A)

for(list in 1:length(all_priors[[1]])){
  name <- names(all_priors[[1]][[list]])[1]
  lapply(all_priors, function(x) x[[list]]) %>% 
    reduce(left_join, by = "group") %>% 
    set_colnames(c("group", seq(1, length(prior_imputed_clean)))) %>% 
    write_csv(paste0(path_to_box, "data/ADAMS/prior_data/priors_", name, ".csv"))
}
