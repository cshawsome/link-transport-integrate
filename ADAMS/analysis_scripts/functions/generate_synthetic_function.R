generate_synthetic <- 
  function(warm_up, run_number, starting_props, 
           unimpaired_preds, other_preds, mci_preds, categorical_vars, 
           continuous_vars, id_var, dementia_var, dataset_to_copy, 
           num_synthetic, unimpaired_betas, unimpaired_cov, other_betas, 
           other_cov, mci_betas, mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
           prior_beta, nu_0, kappa_0, contrasts_matrix,
           path_to_analyses_folder, path_to_figures_folder){
    #---- generate subfolders for results ----
    dir.create(paste0(path_to_analyses_folder, "synthetic_data/run_", 
                      run_number), recursive = TRUE)
    dir.create(paste0(path_to_figures_folder, "diagnostics/run_", run_number), 
               recursive = TRUE)
    dir.create(paste0(path_to_analyses_folder, "diagnostics_data/run_", 
                      run_number), recursive = TRUE)
    
    #---- sampling counts ----
    warm_up = warm_up
    synthetic_sets = num_synthetic
    B = warm_up + synthetic_sets
    
    #---- count contingency cells ----
    cross_class_label <- 
      table(dataset_to_copy$ETHNIC_label, dataset_to_copy$Astroke) %>% 
      as.data.frame() %>% 
      mutate("Stroke" = ifelse(Var2 == 0, "No Stroke", "Stroke")) %>% 
      unite("Cell Label", c("Var1", "Stroke"), sep = " | ", remove = FALSE)
    
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
      set_colnames(apply(expand.grid(seq(1, 4), seq(1:B)), 1, paste, 
                         collapse = ":")) %>% 
      set_rownames(cross_class_label$`Cell Label`)
    
    Sigma_chain <- matrix(nrow = nrow(Z), ncol = 4*B) %>%
      set_colnames(apply(expand.grid(seq(1, 4), seq(1:B)), 1, paste,
                         collapse = ":")) %>%
      set_rownames(Z[, "label"])
    
    mu_chain <-
      matrix(nrow = nrow(Z), ncol = 4*nrow(cross_class_label)*B) %>%
      set_colnames(apply(
        expand.grid(seq(1, 4), seq(1:nrow(cross_class_label)), seq(1:B)), 1, 
        paste,collapse = ":")) %>% set_rownames(Z[, "label"])
    
    #---- start sampling ----
    for(b in 1:B){
      if(b == 1){
        #---- ****init group membership ----
        dataset_to_copy[, "Group"] <- 
          sample(seq(1, 4), size = nrow(dataset_to_copy) , replace = TRUE, 
                 prob = starting_props)
      } else{
        #---- ****latent class gammas ----
        for(model in c("unimpaired", "other", "mci")){
          random_draw <- sample(seq(1, 10000), size = 1)
          
          prior_betas <- as.vector(get(paste0(model, "_betas"))[, random_draw])
          prior_cov <- matrix(unlist(get(paste0(model, "_cov"))[, random_draw]), 
                              nrow = nrow(prior_betas))
          
          model_gamma_chain[which(model_gamma_chain$model == model), b] <- 
            mvrnorm(n = 1, mu = unlist(prior_betas), Sigma = prior_cov)
        }
        
        #---- ****group membership ----
        group = 1
        dataset_to_copy[, "Group"] <- 0
        for(model in c("unimpaired", "other", "mci")){
          subset_index <- which(dataset_to_copy$Group == 0)
          
          dataset_to_copy[subset_index, paste0("p_", model)] <- 
            expit(as.matrix(dataset_to_copy[subset_index, 
                                             get(paste0(model, "_preds"))]) %*% 
                    as.matrix(model_gamma_chain[which(model_gamma_chain$model == 
                                                        model), b]))
          
          dataset_to_copy[subset_index, "Group"] <- 
            rbernoulli(n = length(subset_index), 
                       p = dataset_to_copy[subset_index, 
                                            paste0("p_", model)])*group
          
          group = group + 1
        }
        
        dataset_to_copy[which(dataset_to_copy$Group == 0), "Group"] <- 4
      }
      
      #---- ****group: summary ----
      latent_class_chain[, b] <- 
        table(dataset_to_copy$Group)/sum(table(dataset_to_copy$Group))
      
      for(i in 1:4){
        subset <- dataset_to_copy %>% filter(Group == i) 
        random_draw <- sample(seq(1, 10000), size = 1)
        posterior_counts <- 
          alpha_0_dist[which(alpha_0_dist$group_number == i), random_draw] + 
          table(subset$ETHNIC_label, subset$Astroke) %>% as.data.frame() %>% 
          dplyr::select("Freq") %>% unlist()
        
        #---- ****p(contingency table cell) ----
        pi_chain[, paste0(i, ":", b)] <- 
          rdirichlet(1, alpha = as.numeric(unlist(posterior_counts)))
        
        #---- ****contingency table count ----
        contingency_table <- rmultinom(n = 1, size = nrow(subset), 
                                       prob = pi_chain[, paste0(i, ":", b)])
        UtU <- diag(contingency_table[, 1])
        
        #---- ****draw new UtU if needed ----
        while(det(t(A) %*% UtU %*% A) < 1e-9){
          random_draw <- sample(seq(1, 10000), size = 1)
          new_counts <- alpha_0_dist[, c(random_draw, ncol(alpha_0_dist))] %>% 
            filter(group_number == i) + 
            table(subset$ETHNIC_label, subset$Astroke) %>% as.data.frame() %>% 
            dplyr::select("Freq") %>% unlist()
          
          UtU <- diag(unlist(new_counts[, 1]))
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
        
        #---- ****Mm ----
        continuous_covariates <- subset %>% 
          dplyr::select(all_of(Z[, "var"])) %>% as.matrix
        
        V_inv <- t(A) %*% UtU %*% A 
        random_draw <- sample(seq(1, 10000), size = 1)
        V_0_inv <- 
          matrix(unlist(prior_V_inv[which(prior_V_inv$group_number == i), 
                                    random_draw]), 
                 nrow = nrow(V_inv), ncol = ncol(V_inv))
        beta_0 <- matrix(unlist(priors_beta[which(priors_beta$group_number == i), 
                                            random_draw]), 
                         nrow = nrow(V_inv),  ncol = ncol(continuous_covariates))
        
        M <- solve(V_inv + kappa_0[i]*V_0_inv)
        m <-  t(A) %*% t(U) %*% continuous_covariates - 
          kappa_0[i]*V_0_inv %*% beta_0
        
        Mm <- M %*% m
        
        #---- ****draw Sigma | Y ----
        ZtZ <- t(continuous_covariates) %*% continuous_covariates
        third_term <- kappa_0[i]*t(beta_0) %*% V_0_inv %*% beta_0
        
        random_draw <- sample(seq(1, 10000), size = 1)
        Sigma_prior <- 
          matrix(unlist(prior_Sigma[which(prior_Sigma$group_number == i), 
                                    random_draw]), 
                 nrow = ncol(continuous_covariates))
        sig_Y <- riwish(v = (nu_0 + nrow(subset)), 
                        S = Sigma_prior + ZtZ + third_term)
        
        Sigma_chain[, paste0(i, ":", b)] <- diag(sig_Y)
        
        #---- ****draw beta | Sigma, Y ----
        beta_Sigma_Y <- matrix.normal(Mm, M, sig_Y/kappa_0[i])
        
        #---- ****compute mu ----
        mu_chain[, paste0(i, ":", seq(1, nrow(cross_class_label)), ":", b)] <- 
          t(A %*% matrix(beta_Sigma_Y, nrow = ncol(A), ncol = nrow(Z), 
                         byrow = FALSE))
        
        #---- ****draw data ----
        #reformat contingency table
        contingency_table %<>% cbind(do.call(cbind, list(
          #Black              #Hispanic           #Stroke
          rep(c(1, 0, 0), 2), rep(c(0, 1, 0), 2), c(rep(0, 3), rep(1, 3))))) %>% 
          set_colnames(c("Count", W))
        
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
                                  mu = mu_chain[, paste0(i, ":", j, ":", b)], 
                                  Sigma = sig_Y)))
          } else{
            subset[index:(index - 1 + contingency_table[j, "Count"]), 
                   colnames(sig_Y)] <- 
              mvrnorm(n = contingency_table[j, "Count"],
                      mu = mu_chain[, paste0(i, ":", j, ":", b)], 
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
                         c(W, Z[, "var"])] <- subset[, c(W, Z[, "var"])]
      }
      #---- ****post-processing ----
      #---- ******race/ethnicity ----
      dataset_to_copy %<>% 
        mutate("ETHNIC_label" = case_when(Black == 1 ~ "Black", 
                                          Hispanic == 1 ~ "Hispanic", 
                                          TRUE ~ "White"))
      #---- ****save synthetic sample ----
      if(b > warm_up){
        write_csv(dataset_to_copy, 
                  file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                                "results/ADAMSA/standard_normal/run_", 
                                run_number, "/ADAMSA_synthetic_", b - warm_up, 
                                ".csv"))
      }
    }
    
    #---- **dx plots ----
    #---- ****color palettes ----
    extended_pallette10 <- colorRampPalette(wes_palette("Darjeeling1"))(10)
    extended_pallette14 <- colorRampPalette(wes_palette("Darjeeling1"))(14)
    extended_pallette6 <- colorRampPalette(wes_palette("Darjeeling1"))(6)
    
    #---- ****gamma chains ----
    gamma_plot_data <- model_gamma_chain %>% as.data.frame() %>% 
      set_colnames(c(seq(1:B), "model", "pred")) %>%
      pivot_longer(cols = paste0(seq(1:B)), 
                   names_to = "Run", values_to = "gamma") %>% 
      filter(pred != "(Intercept)")
    
    gamma_chain_plot <- 
      ggplot(data = gamma_plot_data, 
             aes(x = reorder(Run, sort(as.numeric(Run))), y = gamma, 
                 colour = pred)) + geom_line(aes(group = pred)) + 
      facet_grid(rows = vars(factor(model, 
                                    levels = c("unimpaired", "mci", "other"))), 
                 scales = "free") + 
      geom_vline(xintercept = warm_up, size = 1) + theme_bw() + xlab("Run") + 
      scale_x_discrete(breaks = seq(0, B, by = 100)) + 
      scale_color_manual(values = rev(extended_pallette14))
    
    ggsave(filename = "gamma_chain.jpeg", plot = gamma_chain_plot, 
           path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "diagnostics/standard_normal/run_", run), 
           width = 7, height = 6, units = "in", device = "jpeg")
    
    #---- **latent class chain ----
    latent_class_data <- t(latent_class_chain) %>% as.data.frame() %>%
      mutate("run" = seq(1:B)) %>% 
      pivot_longer(-c("run"), names_to = c("Group"), values_to = "prob") %>% 
      arrange(desc(prob)) %>%
      mutate_at("Group", as.factor)
    latent_class_data$Group <- 
      fct_relevel(latent_class_data$Group, 
                  paste0(unique(latent_class_data$Group)))
    
    latent_class_chain_plot <- 
      ggplot(data = latent_class_data, 
             aes(x = run, y = prob, colour = Group)) +       
      geom_line(aes(group = Group)) + 
      geom_vline(xintercept = warm_up, size = 1) + 
      theme_minimal() + xlab("Run") + ylab("Proportion of Sample") +  
      scale_color_manual(values = c(wes_palette("Darjeeling1")[1], 
                                    wes_palette("Darjeeling1")[3], 
                                    wes_palette("Darjeeling1")[2], 
                                    wes_palette("Darjeeling1")[5])) + 
      scale_x_continuous(breaks = seq(0, B, by = 100)) 
    
    ggsave(filename = "latent_class_chain.jpeg", plot = latent_class_chain_plot, 
           path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "diagnostics/standard_normal/run_", run), 
           width = 7, height = 3, units = "in", device = "jpeg")
    
    #---- **pi chain ----
    pi_chain_data <- pi_chain %>% as.data.frame() %>% 
      rownames_to_column("Cell") %>% 
      pivot_longer(-c("Cell"), names_to = c("Group", "Run"), names_sep = ":", 
                   values_to = "probability") %>% 
      mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                       Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                       Group == 4 ~ "Dementia")) %>% 
      mutate_at("Run", as.numeric) %>%
      mutate_if(is.character, as.factor) 
    
    pi_chain_plot <- ggplot(data = pi_chain_data, 
                            aes(x = Run, y = probability, colour = Cell)) +       
      geom_line(aes(group = Cell)) + 
      geom_vline(xintercept = warm_up, size = 1) + 
      xlab("Run") + ylab("Probability of cell membership") +  
      scale_color_manual(values = extended_pallette6) + 
      scale_x_continuous(breaks = seq(0, B, by = 100)) +
      facet_grid(rows = vars(factor(Group_label, 
                                    levels = c("Unimpaired", "MCI", "Dementia", 
                                               "Other")))) + theme_bw() 
    
    ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
           path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "diagnostics/standard_normal/run_", run), 
           width = 7, height = 5, units = "in", device = "jpeg")
    
    #---- **Sigma chain ----
    Sigma_chain_data <- Sigma_chain %>% as.data.frame() %>% 
      rownames_to_column("Z") %>% 
      pivot_longer(-c("Z"), names_to = c("Group", "Run"), names_sep = ":", 
                   values_to = "variance") %>% 
      mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                       Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                       Group == 4 ~ "Dementia")) %>%
      mutate_at("Run", as.numeric) %>%
      mutate_if(is.character, as.factor) 
    
    Sigma_chain_plot <- ggplot(data = Sigma_chain_data, 
                               aes(x = Run, y = variance, colour = Z)) +       
      geom_line(aes(group = Z)) + geom_vline(xintercept = warm_up, size = 1) +
      xlab("Run") + ylab("Variance") +  
      scale_color_manual(values = rev(extended_pallette10)) + 
      scale_x_continuous(breaks = seq(0, B, by = 100)) + 
      facet_grid(rows = vars(factor(Group_label, 
                                    levels = c("Unimpaired", "MCI", "Dementia", 
                                               "Other")))) + theme_bw() 
    
    ggsave(filename = "Sigma_chain.jpeg", plot = Sigma_chain_plot, 
           path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "diagnostics/standard_normal/run_", run), 
           width = 7, height = 5, units = "in", device = "jpeg")
    
    #---- **mu chain ----
    mu_chain_data <- mu_chain %>% as.data.frame() %>% 
      rownames_to_column("Z") %>% 
      pivot_longer(-c("Z"), names_to = c("Group", "Cell", "Run"), 
                   names_sep = ":", values_to = "mu") %>% 
      mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                       Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                       Group == 4 ~ "Dementia")) %>% 
      mutate_at("Run", as.numeric) %>%
      mutate_if(is.character, as.factor) 
    
    mu_chain_plot <- ggplot(data = mu_chain_data, 
                            aes(x = Run, y = mu, colour = Group_label)) +       
      geom_line(aes(group = Group_label), alpha = 0.75) + 
      xlab("Run") + ylab("mu") + geom_vline(xintercept = warm_up, size = 1) + 
      scale_color_manual(values = c(wes_palette("Darjeeling1")[1], 
                                    wes_palette("Darjeeling1")[3], 
                                    wes_palette("Darjeeling1")[5], 
                                    wes_palette("Darjeeling1")[2])) +
      scale_x_continuous(breaks = seq(0, B, by = 100)) + 
      facet_grid(rows = vars(factor(Z)), cols = vars(factor(Cell)), 
                 scales = "free") + theme_bw() 
    
    ggsave(filename = "mu_chain.jpeg", plot = mu_chain_plot, 
           path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "diagnostics/standard_normal/run_", run), 
           width = 14, height = 10, units = "in", device = "jpeg")
    
    #---- save datasets ----
    write_csv(gamma_plot_data, 
              file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/diagnostics_data/", 
                            "run_", run, "/ADAMSA_gamma_plot_data.csv"))
    
    write_csv(latent_class_data, 
              file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/diagnostics_data/", 
                            "run_", run,  "/ADAMSA_latent_class_data.csv"))
    
    write_csv(pi_chain_data, 
              file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/diagnostics_data/", 
                            "run_", run,  "/ADAMSA_pi_chain_data.csv"))
    
    write_csv(Sigma_chain_data, 
              file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/diagnostics_data/", 
                            "run_", run, "/ADAMSA_Sigma_chain_data.csv"))
    
    write_csv(mu_chain_data, 
              file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/diagnostics_data/", 
                            "run_", run, "/ADAMSA_mu_chain_data.csv"))
  }










