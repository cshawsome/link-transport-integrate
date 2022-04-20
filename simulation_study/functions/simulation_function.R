simulation_function <- function(){
  generate_synthetic(warm_up = 100, run_number = 1, 
                     starting_props = c(0.25, 0.25, 0.25, 0.25),
                     unimpaired_preds, other_preds, mci_preds, 
                     categorical_vars = W, continuous_vars = Z, 
                     id_var = "HHIDPN", variable_labels, 
                     dataset_to_copy = synthetic_data_list[[2]] %>% 
                       group_by(married_partnered) %>% 
                       slice_sample(prop = 0.5) %>% 
                       mutate("(Intercept)" = 1) %>% ungroup(), cell_ID_key, 
                     color_palette, num_synthetic = 1000, unimpaired_betas, 
                     unimpaired_cov, other_betas, other_cov, mci_betas, mci_cov, 
                     alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, 
                     kappa_0, contrasts_matrix = A,
                     path_to_analyses_folder = 
                       paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_", 
                              unique(synthetic_data_list[[2]][, "dataset_name"]), 
                              "/"), 
                     path_to_figures_folder = 
                       paste0(path_to_box,
                              "figures/simulation_study/HCAP_HRS_", 
                              unique(synthetic_data_list[[2]][, "dataset_name"]), 
                              "/"), data_only = FALSE)
  return()
}