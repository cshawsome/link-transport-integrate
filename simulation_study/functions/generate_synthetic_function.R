generate_synthetic <- 
  function(warm_up, run_number, starting_props, unimpaired_preds, other_preds, 
           mci_preds, categorical_vars, continuous_vars, id_var, variable_labels, 
           dataset_to_copy, cell_ID_key, color_palette, num_synthetic, 
           unimpaired_betas, unimpaired_cov, other_betas, other_cov, mci_betas, 
           mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, prior_beta, nu_0, 
           kappa_0, contrasts_matrix, path_to_analyses_folder, 
           path_to_figures_folder){
    
    #---- check subfolders for results ----
    if(!dir.exists(paste0(path_to_analyses_folder, "synthetic_data/run_", 
                          run_number))){
      dir.create(paste0(path_to_analyses_folder, "synthetic_data/run_", 
                        run_number), recursive = TRUE)
      
    }
    
    if(!dir.exists(paste0(path_to_figures_folder, "diagnostics/run_", 
                          run_number))){
      dir.create(paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
                 recursive = TRUE)
    }
    
    if(!dir.exists(paste0(path_to_analyses_folder, "diagnostics_data/run_", 
                          run_number))){
      dir.create(paste0(path_to_analyses_folder, "diagnostics_data/run_", 
                        run_number), recursive = TRUE)
    }
    
    #---- sampling counts ----
    warm_up = warm_up
    synthetic_sets = num_synthetic
    B = warm_up + synthetic_sets
    
    #---- count contingency cells ----
    cross_class_label <- dataset_to_copy %>% 
      dplyr::select(all_of(categorical_vars)) %>% 
      unite("cell_ID", everything(), sep = "") %>% table() %>% 
      as.data.frame() %>% set_colnames(c("cell_ID", "count")) %>% 
      left_join(cell_ID_key)
    
    #---- chain storage ----
    model_gamma_chain <- 
      matrix(nrow = sum(length(unimpaired_preds), length(other_preds), 
                        length(mci_preds)), ncol = B) %>% as.data.frame() %>%
      mutate("model" = c(rep("unimpaired", length(unimpaired_preds)), 
                         rep("other", length(other_preds)), 
                         rep("mci", length(mci_preds))), 
             "pred" = c(unimpaired_preds, mci_preds, other_preds))
    
    latent_class_chain <- matrix(nrow = 4, ncol = B) %>% 
      set_rownames(c("Unimpaired", "Other", "MCI", "Dementia"))
    
    pi_chain <- matrix(nrow = nrow(cross_class_label), ncol = 4*B) %>% 
      set_colnames(gsub(" ", "", 
                        apply(expand.grid(
                          c("Unimpaired", "MCI", "Dementia", "Other"), 
                          seq(1, B)), 1, paste, collapse = ":"))) %>% 
      set_rownames(cross_class_label$cell_ID)
    
    Sigma_chain <- matrix(nrow = length(Z), ncol = 4*B) %>%
      set_colnames(gsub(" ", "", 
                        apply(expand.grid(
                          c("Unimpaired", "MCI", "Dementia", "Other"), 
                          seq(1, B)), 1, paste, collapse = ":"))) %>% 
      set_rownames(unlist(variable_labels[variable_labels$data_label %in% Z, 
                                          "figure_label"]))
    
    mu_chain <-
      matrix(nrow = length(Z), ncol = 4*nrow(cross_class_label)*B) %>%
      set_colnames(gsub(" ", "", 
                        apply(expand.grid(
                          c("Unimpaired", "MCI", "Dementia", "Other"),
                          seq(1:nrow(cross_class_label)), seq(1, B)), 1, paste, 
                          collapse = ":"))) %>% 
      set_rownames(unlist(variable_labels[variable_labels$data_label %in% Z, 
                                          "figure_label"]))
    
    #---- start sampling ----
    for(b in 1:B){
      if(b == 1){
        #---- **init group membership ----
        dataset_to_copy[, "group_num"] <- 
          sample(seq(1, 4), size = nrow(dataset_to_copy) , replace = TRUE, 
                 prob = starting_props)
      } else{
        #---- **latent class gammas ----
        for(model in c("unimpaired", "other", "mci")){
          max_index <- 
            colnames(priors_beta)[str_detect(
              colnames(priors_beta), "[0-9]+")] %>% as.numeric() %>% max()
          
          random_draw <- sample(seq(1, max_index), size = 1)
          
          prior_betas <- 
            as.vector(get(paste0(model, "_betas"))[, as.character(random_draw)])
          prior_cov <- 
            matrix(unlist(get(paste0(model, "_cov"))[, as.character(random_draw)]), 
                   nrow = nrow(prior_betas))
          
          model_gamma_chain[which(model_gamma_chain$model == model), b] <- 
            mvrnorm(n = 1, mu = unlist(prior_betas), Sigma = prior_cov)
        }
        
        #---- ****group membership ----
        group_num = 1
        dataset_to_copy[, "group_num"] <- 0
        for(model in c("unimpaired", "other", "mci")){
          subset_index <- which(dataset_to_copy$group_num == 0)
          
          dataset_to_copy[subset_index, paste0("p_", model)] <- 
            expit(as.matrix(dataset_to_copy[subset_index, 
                                            get(paste0(model, "_preds"))]) %*% 
                    as.matrix(model_gamma_chain[which(model_gamma_chain$model == 
                                                        model), b]))
          
          dataset_to_copy[subset_index, "group_num"] <- 
            rbernoulli(n = length(subset_index), 
                       p = dataset_to_copy[subset_index, 
                                           paste0("p_", model)])*group_num
          
          group_num = group_num + 1
        }
      }
      
      dataset_to_copy[, "Group"] <- 
        case_when(dataset_to_copy$group_num == 1 ~ "Unimpaired", 
                  dataset_to_copy$group_num == 2 ~ "Other", 
                  dataset_to_copy$group_num == 3 ~ "MCI", 
                  dataset_to_copy$group_num %in% c(0, 4) ~ "Dementia")
      
      #---- ****group: summary ----
      summary <- table(dataset_to_copy$Group)/sum(table(dataset_to_copy$Group)) 
      if(length(summary) < 4){
        missing <- which(!seq(1, 4) %in% names(summary))
        new_summary <- vector(length = 4)
        new_summary[missing] <- 0
        new_summary[-missing] <- summary
        latent_class_chain[, b] <- new_summary 
      } else{
        latent_class_chain[, b] <- summary 
      }
      
      for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
        subset <- dataset_to_copy %>% filter(Group == class) 
        if(nrow(subset) == 0){
          next
        } else{
          max_index <- 
            colnames(alpha_0_dist)[str_detect(
              colnames(alpha_0_dist), "[0-9]+")] %>% as.numeric() %>% max()
          
          random_draw <- sample(seq(1, max_index), size = 1)
          posterior_first_count <- 
            alpha_0_dist[which(alpha_0_dist$group == class), 
                         as.character(random_draw)]*nrow(subset)
          
          posterior_count <- posterior_first_count +
            subset %>% dplyr::select(all_of(W)) %>% 
            unite("cell_ID", sep = "") %>% table() %>% as.data.frame() %>% 
            dplyr::select("Freq") %>% unlist()
          
          #---- **p(contingency table cell) ----
          pi_chain[, paste0(class, ":", b)] <- 
            MCMCpack::rdirichlet(1, alpha = as.numeric(unlist(posterior_count)))
          
          #---- ****contingency table count ----
          contingency_table <- 
            rmultinom(n = 1, size = nrow(subset), 
                      prob = pi_chain[, paste0(class, ":", b)])
          
          UtU <- diag(contingency_table[, 1])
          
          #---- ****draw new UtU if needed ----
          while(det(t(A) %*% UtU %*% A) < 1e-9){
            
            random_draw <- sample(seq(1, max_index), size = 1)
            
            posterior_first_count <- 
              alpha_0_dist[which(alpha_0_dist$group == class), 
                           as.character(random_draw)]*nrow(subset)
            
            posterior_count <- posterior_first_count +
              subset %>% dplyr::select(all_of(W)) %>% 
              unite("cell_ID", sep = "") %>% table() %>% as.data.frame() %>% 
              dplyr::select("Freq") %>% unlist()
            
            pi_chain[, paste0(class, ":", b)] <- 
              MCMCpack::rdirichlet(1, 
                                   alpha = as.numeric(unlist(posterior_count)))
            
            contingency_table <- 
              rmultinom(n = 1, size = nrow(subset), 
                        prob = pi_chain[, paste0(class, ":", b)])
            
            UtU <- diag(contingency_table[, 1])
          }
          
          #---- ****make U matrix ----
          U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
          
          for(j in 1:nrow(contingency_table)){
            if(contingency_table[j, ] == 0){next}
            if(j == 1){
              index = 1
            } else{
              index = sum(contingency_table[1:(j - 1), ]) + 1
            }
            U[index:(index - 1 + contingency_table[j, ]), j] <- 1
          }
          
          #---- **Mm ----
          continuous_covariates <- subset %>% dplyr::select(all_of(Z)) %>% 
            as.matrix
          
          V_inv <- t(A) %*% UtU %*% A 
          random_draw <- sample(seq(1, max_index), size = 1)
          V_0_inv <- 
            matrix(unlist(prior_V_inv[which(prior_V_inv$group == class), 
                                      as.character(random_draw)]), 
                   nrow = nrow(V_inv), ncol = ncol(V_inv))
          beta_0 <- 
            matrix(unlist(priors_beta[which(priors_beta$group == class), 
                                      as.character(random_draw)]), 
                   nrow = nrow(V_inv),  
                   ncol = ncol(continuous_covariates))
          
          M <- solve(V_inv + kappa_0[class]*V_0_inv)
          m <-  t(A) %*% t(U) %*% continuous_covariates - 
            kappa_0[class]*V_0_inv %*% beta_0
          
          Mm <- M %*% m
          
          #---- ****draw Sigma | Y ----
          ZtZ <- t(continuous_covariates) %*% continuous_covariates
          third_term <- kappa_0[class]*t(beta_0) %*% V_0_inv %*% beta_0
          
          random_draw <- sample(seq(1, max_index), size = 1)
          Sigma_prior <- 
            matrix(unlist(prior_Sigma[which(prior_Sigma$group == class), 
                                      as.character(random_draw)]), 
                   nrow = ncol(continuous_covariates))
          sig_Y <- riwish(v = (nu_0 + nrow(subset)), 
                          S = Sigma_prior + ZtZ + third_term)
          
          Sigma_chain[, paste0(class, ":", b)] <- diag(sig_Y)
          
          #---- ****draw beta | Sigma, Y ----
          beta_Sigma_Y <- matrix.normal(Mm, M, sig_Y/kappa_0[class])
          
          #---- ****compute mu ----
          mu_chain[, paste0(class, ":", 
                            seq(1, nrow(cross_class_label)), ":", b)] <- 
            t(A %*% matrix(beta_Sigma_Y, nrow = ncol(A), ncol = length(Z), 
                           byrow = FALSE))
          
          #---- ****draw data ----
          #reformat contingency table
          contingency_table %<>% set_colnames("Count") %>% 
            cbind(contrasts_matrix[, -1])
          
          for(j in 1:nrow(contingency_table)){
            if(contingency_table[j, "Count"] == 0){next}
            if(j == 1){
              index = 1
            } else{
              index = sum(contingency_table[1:(j - 1), "Count"]) + 1
            }
            #Z (continuous data)
            if(contingency_table[j, "Count"] == 1){
              subset[index:(index - 1 + contingency_table[j, "Count"]), 
                     colnames(sig_Y)] <- 
                t(as.matrix(mvrnorm(n = contingency_table[j, "Count"],
                                    mu = mu_chain[, paste0(class, ":", j, ":", b)], 
                                    Sigma = sig_Y)))
            } else{
              subset[index:(index - 1 + contingency_table[j, "Count"]), 
                     colnames(sig_Y)] <- 
                mvrnorm(n = contingency_table[j, "Count"],
                        mu = mu_chain[, paste0(class, ":", j, ":", b)], 
                        Sigma = sig_Y)
            }
            
            #W (categorical data)
            subset[index:(index - 1 + contingency_table[j, "Count"]), 
                   colnames(contingency_table)[-1]] <- 
              matrix(rep(contingency_table[j, colnames(contingency_table)[-1]], 
                         contingency_table[j, "Count"]), 
                     ncol = 3, byrow = TRUE)
          }
          
          #---- ****replace synthetic data ----
          dataset_to_copy[which(dataset_to_copy$HHIDPN %in% subset$HHIDPN),
                          c(W, Z)] <- subset[, c(W, Z)]
        }
      }
      
      #---- ****save synthetic sample ----
      if(b > warm_up){
        # write_csv(dataset_to_copy, 
        #           file = paste0(path_to_analyses_folder, "synthetic_data/run_", 
        #                         run_number, "/synthetic_", b - warm_up, 
        #                         ".csv"))
        #test a different way of saving the datasets
        if(!exists("dataset_list")){
          dataset_list <- list()
        }
        dataset_list[[b - warm_up]] <- dataset_to_copy
      }
    }
    
    saveRDS(dataset_list, 
            file = paste0(path_to_analyses_folder, "synthetic_data/run_",
                          run_number, "/synthetic_dataset_list"))
    
    #---- **dx plots ----
    #---- ****gamma chains ----
    gamma_plot_data <- model_gamma_chain %>% as.data.frame() %>% 
      set_colnames(c(seq(1:B), "model", "pred")) %>%
      pivot_longer(cols = paste0(seq(1:B)), 
                   names_to = "Run", values_to = "gamma") %>% 
      filter(pred != "(Intercept)") %>% 
      filter(!is.na(gamma)) %>% 
      left_join(., variable_labels[, c("data_label", "figure_label")], 
                by = c("pred" = "data_label"))
    
    gamma_chain_plot <- 
      ggplot(data = gamma_plot_data, 
             aes(x = reorder(Run, sort(as.numeric(Run))), y = gamma, 
                 colour = figure_label)) + geom_line(aes(group = figure_label)) + 
      facet_grid(rows = vars(factor(model, 
                                    levels = c("unimpaired", "mci", "other"))), 
                 scales = "free") + 
      geom_vline(xintercept = warm_up, size = 1) + theme_bw() + xlab("Run") + 
      scale_x_discrete(breaks = seq(0, B, by = 100)) + 
      scale_color_manual(values = 
                           rev(colorRampPalette(wes_palette("Darjeeling1"))(
                             length(unique(model_gamma_chain$pred)))), 
                         name = "Predictors") + 
      theme(legend.position = "bottom")
    
    ggsave(filename = "gamma_chain.jpeg", plot = gamma_chain_plot, 
           path = paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
           height = 7, width = 12, units = "in", device = "jpeg")
    
    #---- ****latent class chain ----
    latent_class_data <- t(latent_class_chain) %>% as.data.frame() %>%
      mutate("run" = seq(1:B)) %>% 
      pivot_longer(-c("run"), names_to = c("Group"), values_to = "prob") %>% 
      arrange(desc(prob)) %>%
      mutate_at("Group", as.factor) %>% left_join(color_palette)
    latent_class_data$Group <- 
      fct_relevel(latent_class_data$Group, 
                  paste0(unique(latent_class_data$Group))) 
    
    latent_class_chain_plot <- 
      ggplot(data = latent_class_data, 
             aes(x = run, y = prob, colour = Group)) +       
      geom_line(aes(group = Group)) + 
      geom_vline(xintercept = warm_up, size = 1) + 
      theme_minimal() + xlab("Run") + ylab("Proportion of Sample") +  
      scale_color_manual(values = unique(latent_class_data$Color)) + 
      scale_x_continuous(breaks = seq(0, B, by = 100)) + 
      theme(legend.position = "bottom")
    
    ggsave(filename = "latent_class_chain.jpeg", plot = latent_class_chain_plot, 
           path = paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
           height = 7, width = 10, units = "in", device = "jpeg")
    
    #---- ****pi chain ----
    pi_chain_data <- pi_chain %>% as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(-c("Cell"), names_to = c("Group", "Run"), names_sep = ":", 
                   values_to = "probability") %>% arrange(desc(probability)) %>%
      left_join(., cell_ID_key, by = c("Cell" = "cell_ID")) %>%
      mutate_at("Run", as.numeric) %>% mutate_if(is.character, as.factor) %>% 
      rename("cell_ID" = "Cell", "Cell" = "cell_name")
    
    pi_chain_data$Cell <- 
      fct_relevel(pi_chain_data$Cell, paste0(unique(pi_chain_data$Cell)))
    
    pi_chain_plot <- ggplot(data = pi_chain_data, 
                            aes(x = Run, y = probability, colour = Cell)) +       
      geom_line(aes(group = Cell)) + 
      geom_vline(xintercept = warm_up, size = 1) + 
      xlab("Run") + ylab("Probability of cell membership") +  
      scale_color_manual(values = 
                           rev(colorRampPalette(wes_palette("Darjeeling1"))(
                             nrow(cell_ID_key)))) + 
      scale_x_continuous(breaks = seq(0, B, by = 100)) +
      facet_grid(rows = vars(factor(Group, 
                                    levels = c("Unimpaired", "MCI", "Dementia", 
                                               "Other")))) + theme_bw() + 
      theme(legend.position = "bottom")
    
    ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
           path = paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
           device = "jpeg")
    
    #---- ****Sigma chain ----
    Sigma_chain_data <- Sigma_chain %>% as.data.frame() %>% 
      rownames_to_column("Z") %>% 
      pivot_longer(-c("Z"), names_to = c("Group", "Run"), names_sep = ":", 
                   values_to = "variance") %>% mutate_at("Run", as.numeric) %>%
      mutate_if(is.character, as.factor) 
    
    Sigma_chain_plot <- ggplot(data = Sigma_chain_data, 
                               aes(x = Run, y = variance, colour = Z)) +       
      geom_line(aes(group = Z)) + geom_vline(xintercept = warm_up, size = 1) +
      xlab("Run") + ylab("Variance") +  
      scale_color_manual(values = 
                           rev(colorRampPalette(wes_palette("Darjeeling1"))(
                             nrow(Sigma_chain)))) + 
      scale_x_continuous(breaks = seq(0, B, by = 100)) + 
      facet_grid(rows = vars(factor(Group, 
                                    levels = c("Unimpaired", "MCI", "Dementia", 
                                               "Other")))) + theme_bw() + 
      theme(legend.position = "bottom")
    
    ggsave(filename = "Sigma_chain.jpeg", plot = Sigma_chain_plot, 
           path = paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
           height = 7, width = 12, units = "in", device = "jpeg")
    
    #---- ****mu chain ----
    mu_chain_data <- mu_chain %>% as.data.frame() %>% 
      rownames_to_column("Z") %>% 
      pivot_longer(-c("Z"), names_to = c("Group", "Cell", "Run"), 
                   names_sep = ":", values_to = "mu") %>% 
      left_join(color_palette) %>% 
      left_join(cell_ID_key, by = c("Cell" = "cell_order")) %>%
      mutate_at("Run", as.numeric) %>% mutate_if(is.character, as.factor) %>% 
      mutate_at("Color", as.character)
    
    mu_chain_data$Group <- 
      fct_relevel(mu_chain_data$Group, 
                  c("Unimpaired", "MCI", "Dementia", "Other"))
    
    
    mu_chain_plot <- ggplot(data = mu_chain_data, 
                            aes(x = Run, y = mu, colour = Group)) +       
      geom_line(aes(group = Group), alpha = 0.75) + 
      xlab("Run") + ylab("mu") + geom_vline(xintercept = warm_up, size = 1) + 
      scale_color_manual(values = unique(mu_chain_data$Color)) +
      scale_x_continuous(breaks = seq(0, B, by = 100)) + 
      facet_grid(rows = vars(factor(Z)), cols = vars(factor(cell_name)), 
                 scales = "free") + theme_bw() + 
      theme(legend.position = "bottom")
    
    ggsave(filename = "mu_chain.jpeg", plot = mu_chain_plot, 
           path = paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
           height = 15, width = 12, units = "in", device = "jpeg")
    
    #---- save datasets ----
    write_csv(gamma_plot_data, 
              file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "run_", run_number, "/gamma_plot_data.csv"))
    
    write_csv(latent_class_data, 
              file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "run_", run_number, "/latent_class_data.csv"))
    
    write_csv(pi_chain_data, 
              file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "run_", run_number, "/pi_chain_data.csv"))
    
    write_csv(Sigma_chain_data, 
              file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "run_", run_number, "/Sigma_chain_data.csv"))
    
    write_csv(mu_chain_data, 
              file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "run_", run_number, "/mu_chain_data.csv"))
  }

# #---- test function ----
# warm_up = 2 
# run_number = 1 
# starting_props = c(0.25, 0.25, 0.25, 0.25)
# unimpaired_preds = unimpaired_preds
# other_preds = other_preds
# mci_preds = mci_preds
# categorical_vars = W 
# continuous_vars = Z 
# id_var = "HHIDPN" 
# variable_labels = variable_labels
# dataset_to_copy = dataset_to_copy
# cell_ID_key = cell_ID_key
# color_palette = color_palette
# num_synthetic = 10 
# unimpaired_betas = unimpaired_betas
# unimpaired_cov = unimpaired_cov
# other_betas = other_betas
# other_cov = other_cov
# mci_betas = mci_betas
# mci_cov = mci_cov
# alpha_0_dist = alpha_0_dist 
# prior_Sigma = prior_Sigma
# prior_V_inv = prior_V_inv
# prior_beta = priors_beta
# nu_0 = nu_0
# kappa_0 = kappa_0 
# contrasts_matrix = A
# path_to_analyses_folder = 
#   paste0(path_to_box, "analyses/simulation_study/HCAP_normal_250_unimpaired/") 
# path_to_figures_folder = 
#   paste0(path_to_box, "figures/simulation_study/HCAP_normal_250_unimpaired/") 










