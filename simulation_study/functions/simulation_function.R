simulation_function <- 
  function(warm_up, starting_props, unimpaired_preds, other_preds, 
           mci_preds, categorical_vars, continuous_vars, id_var, variable_labels, 
           dataset, cell_ID_key, color_palette, num_synthetic, unimpaired_betas, 
           unimpaired_cov, other_betas, other_cov, mci_betas, mci_cov, 
           alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, kappa_0, 
           contrasts_matrix, seed, path_to_results){
    
    #---- pre-allocated results ----
    result_names <- 
      c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other", "seed", 
        "time")
    
    results <- matrix(ncol = length(result_names), nrow = 1) %>% 
      set_colnames(all_of(result_names))
    
    #---- seed setting ----
    set.seed(seed)
    results[, "seed"] <- seed
    
    #---- start time ----
    start <- Sys.time()
    
    #---- create synthetic HCAP ----
    dataset_to_copy <- dataset %>% group_by(married_partnered) %>% 
      slice_sample(prop = 0.5) %>% 
      mutate("(Intercept)" = 1) %>% ungroup()
    
    synthetic_HCAP <- 
      generate_synthetic(warm_up, run_number = NA, starting_props,
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars, continuous_vars, id_var, 
                         variable_labels, dataset_to_copy , cell_ID_key, 
                         color_palette, num_synthetic, unimpaired_betas, 
                         unimpaired_cov, other_betas, other_cov, mci_betas, 
                         mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                         prior_beta, nu_0, kappa_0, contrasts_matrix,
                         path_to_analyses_folder = NA, 
                         path_to_figures_folder = NA, data_only = TRUE)
    
    #---- end time ----
    results[, "time"] <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    
    #---- write results ----
    write_csv(as.data.frame(results), 
              file = paste0(path_to_results, dataset$dataset_name[1], 
                            ".csv"), append = TRUE, col_names = TRUE)
  }

#---- testing ----
warm_up = 100 
starting_props = rep(0.25, 4) 
unimpaired_preds 
other_preds
mci_preds 
categorical_vars = W
continuous_vars = Z
id_var = "HHIDPN"
variable_labels 
dataset = synthetic_data_list[[1]]
cell_ID_key
color_palette 
num_synthetic = 1000
unimpaired_betas 
unimpaired_cov
other_betas
other_cov
mci_betas
mci_cov
alpha_0_dist
prior_Sigma
prior_V_inv
prior_beta
nu_0
kappa_0 
contrasts_matrix = A
seed = 1
path_to_results <- paste0(path_to_box, "analyses/simulation_study/results/")
