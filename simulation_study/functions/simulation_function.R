simulation_function <- 
  function(warm_up, starting_props, unimpaired_preds, other_preds, 
           mci_preds, categorical_vars, continuous_vars, id_var, variable_labels, 
           scenario, superpops_list, all_scenarios_list, cell_ID_key, 
           color_palette, num_synthetic, unimpaired_betas, unimpaired_cov, 
           other_betas, other_cov, mci_betas, mci_cov, alpha_0_dist, 
           prior_Sigma, prior_V_inv, prior_beta, nu_0_vec, kappa_0_mat, 
           contrasts_matrix, truth, seed, path_to_results){
    
    #---- pre-allocated results ----
    result_names <- 
      c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other", 
        "mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other",
        "bias_Unimpaired", "bias_MCI", "bias_Dementia", "bias_Other",
        "LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other",
        "UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other",
        "Unimpaired_coverage", "MCI_coverage", "Dementia_coverage", 
        "Other_coverage", "black_beta", "black_se", "black_LCI", "black_UCI",
        "hispanic_beta", "hispanic_se", "hispanic_LCI", "hispanic_UCI",
        "black_coverage", "hispanic_coverage", "time", "seed")
    
    results <- matrix(ncol = length(result_names), nrow = 1) %>% 
      set_colnames(all_of(result_names))
    
    #---- superpop names ----
    superpop_names <- 
      unlist(lapply(superpops_list, function(x) unique(x[, "dataset_name"])))
    
    #---- scenario name ----
    scenario_name <- 
      unlist(
        unite(all_sim_scenarios[scenario, 
                                c("distribution", "sample_size", "prior_prop")], 
              col = "name", sep = "_"))
    
    superpop_index <- which(superpop_names == 
                              paste0(all_sim_scenarios[scenario, "distribution"], 
                                     "_1000000_", 
                                     all_sim_scenarios[scenario, "prior_prop"]))
    
    superpop_name <- superpops_list[[superpop_index]]$dataset_name[1]
    
    #---- start time ----
    start <- Sys.time()
    
    dataset_to_copy <- 
      #---- create synthetic HRS ----
    superpops_list[[superpop_index]] %>% 
      slice_sample(n = unlist(all_sim_scenarios[scenario, "sample_size"]), 
                   replace = FALSE) %>%
      #---- create synthetic HCAP ----
    group_by(married_partnered) %>% 
      slice_sample(prop = 0.5, replace = FALSE) %>% 
      mutate("(Intercept)" = 1) %>% ungroup() %>% 
      mutate("dataset_name" = scenario_name)
    
    #---- **true impairment class counts ----
    results[, 
            c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other")] <- 
      colSums(dataset_to_copy[, c("Unimpaired", "MCI", "Dementia", "Other")])
    
    #---- generate synthetic data ----
    synthetic_HCAP <- 
      generate_synthetic(warm_up, run_number = NA, starting_props,
                         unimpaired_preds, other_preds, mci_preds, 
                         categorical_vars, continuous_vars, id_var, 
                         variable_labels, dataset_to_copy , cell_ID_key, 
                         color_palette, num_synthetic, unimpaired_betas, 
                         unimpaired_cov, other_betas, other_cov, mci_betas, 
                         mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
                         prior_beta, nu_0 = nu_0_vec[, scenario_name], 
                         kappa_0_mat, contrasts_matrix,
                         path_to_analyses_folder = NA, 
                         path_to_figures_folder = NA, data_only = TRUE)
    
    #---- impairment class counts ----
    counts <- 
      lapply(synthetic_HCAP, function(x) table(x[, "Group"])) %>% 
      do.call(rbind, .)
    
    #---- **mean counts ----
    mean_counts <- round(colMeans(counts))
    results[, paste0("mean_", colnames(counts))] <- mean_counts
    
    #---- **CI ----
    results[, paste0("LCI_", colnames(counts))] <- 
      apply(counts, 2, function(x) quantile(x, 0.25))
    
    results[, paste0("UCI_", colnames(counts))] <- 
      apply(counts, 2, function(x) quantile(x, 0.975))
    
    #---- **bias ----
    results[, paste0("bias_", colnames(counts))] <- 
      mean_counts - results[, paste0("true_", colnames(counts))]
    
    #---- **coverage ----
    for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
      results[, paste0(class, "_coverage")] <- 
        (results[, paste0("true_", class)] >= results[, paste0("LCI_", class)])*
        (results[, paste0("true_", class)] <= results[, paste0("UCI_", class)])
    }
    
    #---- models ----
    #---- **add weights ----
    #the only variable that contributed to selection was marital status
    # we selected half of married people and half of unmarried people, thus the 
    # weight for every observation should be 2
    
    synthetic_HCAP %<>% lapply(function(x) x %<>% mutate("weight" = 2))
    
    #---- **predicted dementia status ----
    synthetic_HCAP %<>% 
      lapply(function(x) x %<>% 
               mutate("predicted_Dementia" = ifelse(Group == "Dementia", 1, 0)))
    
    #---- models ----
    models <- 
      lapply(synthetic_HCAP, 
             function(dataset) glm(predicted_Dementia ~ 
                                     age_Z + female + black + hispanic, 
                                   family = "poisson", data = dataset)) 
    
    #---- **pooled ----
    pooled_model <- summary(mice::pool(models))
    
    #---- ****truth table ----
    truth_table <- 
      truth[which(truth$dataset_name == superpop_name), ] 
    
    for(race_eth in c("black", "hispanic")){
      #---- ****estimates ----
      results[, paste0(race_eth, "_beta")] <- 
        pooled_model[which(pooled_model$term == race_eth), "estimate"]
      
      results[, paste0(race_eth, "_se")] <- 
        pooled_model[which(pooled_model$term == race_eth), "std.error"]
      
      results[, paste0(race_eth, "_LCI")] <- 
        results[, paste0(race_eth, "_beta")] - 
        1.96*results[, paste0(race_eth, "_se")]
      
      results[, paste0(race_eth, "_UCI")] <- 
        results[, paste0(race_eth, "_beta")] + 
        1.96*results[, paste0(race_eth, "_se")]
      
      #---- ****coverage ----
      results[, paste0(race_eth, "_coverage")] <- 
        (truth_table[which(truth_table$term == race_eth), "estimate"] >= 
           results[, paste0(race_eth, "_LCI")])*
        (truth_table[which(truth_table$term == race_eth), "estimate"] <= 
           results[, paste0(race_eth, "_UCI")])
    }
    
    #---- end time ----
    results[, "time"] <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    
    #---- seed ----
    results[, "seed"] <- seed
    
    #---- write results ----
    file_path <- 
      paste0(path_to_results, scenario_name, ".csv")
    
    if(file.exists(file_path)){
      write_csv(as.data.frame(results), file = file_path, append = TRUE)
    } else{
      write_csv(as.data.frame(results), file = file_path, col_names = TRUE)
    }
  }

# #---- testing ----
# warm_up = 100
# starting_props = rep(0.25, 4)
# unimpaired_preds
# other_preds
# mci_preds
# categorical_vars = W
# continuous_vars = Z
# id_var = "HHIDPN"
# variable_labels
# scenario = 1
# superpops_list = superpop_data_list
# all_scenarios_list = all_sim_scenarios
# cell_ID_key
# color_palette
# num_synthetic = 1000
# unimpaired_betas
# unimpaired_cov
# other_betas
# other_cov
# mci_betas
# mci_cov
# alpha_0_dist
# prior_Sigma
# prior_V_inv
# prior_beta
# nu_0_vec
# kappa_0
# contrasts_matrix = A
# truth
# path_to_results <- paste0(path_to_box, "analyses/simulation_study/results/")
# 
# set.seed(20220512)
# 
# replicate(2, 
#           simulation_function(warm_up, starting_props, unimpaired_preds, 
#                               other_preds, mci_preds, categorical_vars, 
#                               continuous_vars, id_var, variable_labels, 
#                               scenario = 1, 
#                               superpops_list = superpop_data_list, 
#                               all_sim_scenarios, 
#                               cell_ID_key, color_palette, num_synthetic, 
#                               unimpaired_betas, unimpaired_cov, other_betas, 
#                               other_cov, mci_betas, mci_cov, alpha_0_dist, 
#                               prior_Sigma, prior_V_inv, prior_beta, nu_0_vec, 
#                               kappa_0_mat, contrasts_matrix, truth, 
#                               path_to_results))
# 
