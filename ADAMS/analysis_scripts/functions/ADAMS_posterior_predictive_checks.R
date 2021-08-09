ADAMS_posterior_predictive_checks <- 
  function(dataset_to_copy, continuous_covariates, orig_means, orig_sds, 
           num_samples, num_chains, dementia_var, path_to_analyses_folder, 
           path_to_figures_folder){
    
    #---- create directories for results ----
    for(chain in 1:num_chains){
      for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
        dir.create(paste0(path_to_figures_folder, 
                          "posterior_predictive_checks/", 
                          "run_", chain, "/cell_counts/group_specific/", 
                          tolower(group)), recursive = TRUE)
      }
    }
    
    for(chain in 1:num_chains){
      for(metric in c("median", "skew")){
        dir.create(paste0(path_to_figures_folder, 
                          "posterior_predictive_checks/", 
                          "run_", chain, "/continuous_vars/", metric, 
                          "/overall"), recursive = TRUE)
        for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
          dir.create(paste0(path_to_figures_folder, 
                            "posterior_predictive_checks/", 
                            "run_", chain, "/continuous_vars/", metric, "/",
                            tolower(group)), recursive = TRUE)
        }
      }
      dir.create(paste0(path_to_figures_folder, "posterior_predictive_checks/", 
                        "run_", chain, "/impairment_classes"), recursive = TRUE)
    }
    
    #---- read in synthetic data ----
    for(sample in 1:num_samples){
      for(chain in 1:num_chains){
        if(sample == 1 & chain == 1){
          synthetic_sample <- 
            vroom(paste0(path_to_analyses_folder, "synthetic_data/run_", 
                         chain, "/ADAMSA_synthetic_", sample, ".csv")) %>% 
            mutate("sample" = sample, "chain" = chain)
        } else{
          synthetic_sample %<>% 
            rbind(., 
                  vroom(paste0(path_to_analyses_folder, 
                               "synthetic_data/run_", chain, 
                               "/ADAMSA_synthetic_", sample, ".csv")) %>% 
                    mutate("sample" = sample, "chain" = chain))
        }
      }
    }
    
    #---- chain convergence ----
    for(chain_num in 1:num_chains){
      if(chain_num == 1){
        group_membership <- 
          vroom(paste0(path_to_analyses_folder, "diagnostics_data/", 
                       "run_", chain_num, 
                       "/ADAMSA_latent_class_data.csv")) %>% 
          mutate("chain" = chain_num)
      } else{
        group_membership %<>% 
          rbind(., 
                vroom(paste0(path_to_analyses_folder, "diagnostics_data/", 
                             "run_", chain_num, 
                             "/ADAMSA_latent_class_data.csv")) %>% 
                  mutate("chain" = chain_num))
      }
    }
    
    #---- **plot ----
    chain_convergence <- 
      ggplot(data = group_membership, 
             aes(x = run, y = prob, colour = as.factor(chain))) +       
      geom_line(aes(group = as.factor(chain)), alpha = 0.75) + xlab("Run") + 
      ylab("Proportion") + guides(color = guide_legend(title = "Chain")) + 
      scale_color_manual(values = wes_palette("Darjeeling2")) +
      facet_wrap(facets = 
                   vars(factor(Group, 
                               levels = c("Unimpaired", "MCI", "Dementia", 
                                          "Other"))), nrow = 2, ncol = 2) + 
      theme_bw() 
    
    ggsave(filename = "group_membership_chains.jpeg", plot = chain_convergence, 
           path = paste0(path_to_figures_folder, "diagnostics/"), 
           width = 10, height = 5, units = "in", device = "jpeg")
    
    #---- pre-allocate results matrix ----
    #---- **synthetic counts ----
    synthetic_counts <- 
      matrix(0, nrow = 6*4*num_chains, ncol = (num_samples + 3)) %>% 
      as.data.frame() %>% 
      set_colnames(c(seq(1, num_samples), "group", "cell", "chain"))
    
    synthetic_counts[, "group"] <- 
      rep(rep(unlist(unique(dataset_to_copy[, dementia_var])), each = 6), 
          num_chains)
    
    cells <- 
      as.data.frame(table(dataset_to_copy$ETHNIC_label, 
                          dataset_to_copy$Astroke)) %>% 
      unite("cell", c("Var1", "Var2"), sep = ":")
    
    synthetic_counts[, "cell"] <- rep(rep(cells$cell, 4), num_chains)
    
    synthetic_counts[, "chain"] <- rep(seq(1, num_chains), each = 6*4) 
    
    #---- categorical checks ----
    #---- **race/ethnicity x stroke ----
    #true counts
    for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
      subset <- dataset_to_copy %>% filter(!!as.symbol(dementia_var) == group)
      assign(paste0(group, "_data_counts"), 
             as.data.frame(table(subset$ETHNIC_label, 
                                 subset$Astroke)) %>% 
               mutate("percent" = round((Freq/sum(Freq))*100, 1)) %>% 
               unite("cell", c("Var1", "Var2"), sep = ":")) 
    }
    
    for(chain_num in unique(synthetic_sample$chain)){
      #counts from synthetic datasets
      for(num in 1:num_samples){
        subsample <- synthetic_sample %>% 
          filter(chain == chain_num & sample == num) 
        
        for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
          subset <- subsample %>% filter(!!as.symbol(dementia_var) == group) 
          
          counts <- 
            as.data.frame(table(subset$ETHNIC_label, subset$Astroke)) %>% 
            unite("cell", c("Var1", "Var2"), sep = ":")
          
          synthetic_counts[
            which(synthetic_counts$group == group & 
                    synthetic_counts$chain == chain_num & 
                    synthetic_counts$cell %in% counts$cell), num] <- counts$Freq
        }
      }
    }
    
    synthetic_count_plot_data <- synthetic_counts %>% mutate("truth" = 0) 
    
    for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
      true_counts <- get(paste0(group, "_data_counts"))
      synthetic_count_plot_data[
        which(synthetic_count_plot_data$group == group & 
                synthetic_count_plot_data$cell %in% true_counts$cell), 
        "truth"] <- true_counts$Freq
    }
    
    synthetic_count_plot_data %<>% 
      mutate("cat" = rep(rep(c("Black + No Stroke", "Hispanic + No Stroke", 
                               "White + No Stroke", "Black + Stroke", 
                               "Hispanic + Stroke", "White + Stroke"), 4), 
                         num_chains)) %>% 
      pivot_longer(-c("group", "cell", "chain", "truth", "cat")) %>% 
      mutate("color" = case_when(group == "Unimpaired" ~ "#00a389", 
                                 group == "Other" ~ "#28bed9", 
                                 group == "MCI" ~ "#fdab00", 
                                 group == "Dementia" ~ "#ff0000"))
    
    #---- **plot ----
    for(chain in 1:num_chains){
      for(dem_group in unique(synthetic_count_plot_data$group)){
        for(category in unique(synthetic_count_plot_data$cat)){
          subset <- synthetic_count_plot_data %>% 
            filter(group == dem_group & cat == category & chain == chain)
          ggplot(data = subset , aes(x = value)) + 
            geom_histogram(fill = "black", color = "black") + theme_minimal() + 
            xlab("Count") + ggtitle(category) + 
            geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                       size = 2)
          
          ggsave(filename = paste0(path_to_figures_folder, 
                                   "posterior_predictive_checks/run_", chain, 
                                   "/cell_counts/group_specific/", 
                                   tolower(dem_group), "/", tolower(dem_group), 
                                   "_", category, "_count.jpeg"), 
                 width = 5, height = 3, units = "in")
        } 
      }
    }
    
    #---- continuous checks ----
    #---- **density ----
    #---- **median ----
    #---- ****overall ----
    synthetic_continuous <- 
      matrix(0, nrow = nrow(Z)*num_chains, ncol = (num_samples + 2)) %>% 
      as.data.frame() %>% set_colnames(c(seq(1, num_samples), "var", "chain"))
    synthetic_continuous[, "var"] <- rep(Z[, "var"], num_chains)
    synthetic_continuous[, "chain"] <- 
      rep(seq(1, num_chains), each = length(Z[, "var"]))
    
    #medians from synthetic datasets
    for(var in Z[, "var"]){
      subset <- synthetic_sample %>% 
        dplyr::select(!!as.symbol(var), "sample", "chain") 
      subset[, var] <- subset[, var]*orig_sds[var] + orig_means[var]
      synthetic_continuous[which(synthetic_continuous$var == var & 
                                   synthetic_continuous$chain %in% 
                                   seq(1, num_chains)), 
                           seq(1, num_samples)] <- 
        matrix(subset %>% group_by(chain, sample) %>% 
                 summarize_all(median) %>% ungroup() %>% 
                 dplyr::select(all_of(var)) %>% unlist(), 
               nrow = num_chains, ncol = num_samples, byrow = TRUE)
    }
    
    synthetic_continuous %<>% 
      mutate("truth" = 0, 
             "label" = rep(Z[, "label"], num_chains))
    
    for(var in Z[, "var"]){
      synthetic_continuous[which(synthetic_continuous$var == var), "truth"] <- 
        median(unlist(dataset_to_copy[, var]*orig_sds[var] + orig_means[var]))
    }
    
    synthetic_continuous %<>% pivot_longer(seq(1:num_samples))
    
    #---- ****plot ----
    for(chain_num in seq(1:num_chains)){
      for(var_name in Z[, "label"]){
        subset <- synthetic_continuous %>% 
          filter(label == var_name & chain == chain_num)
        ggplot(data = subset , aes(x = value)) + 
          geom_histogram(fill = "black", color = "black") + theme_minimal() + 
          xlab("Median") + ggtitle(var_name) + 
          geom_vline(xintercept = subset$truth, color = "#f2caaa" , size = 2)
        
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/run_", chain_num, 
                                 "/continuous_vars/median/overall/", var_name, 
                                 ".jpeg"), width = 5, height = 3, units = "in")
      }
    }
    
    #---- ****by class ----
    synthetic_continuous <- 
      matrix(0, nrow = nrow(Z)*4*num_chains, ncol = (num_samples + 3)) %>% 
      as.data.frame() %>% 
      set_colnames(c(seq(1, num_samples), "group", "var", "chain"))
    synthetic_continuous[, "group"] <- rep(seq(1, 4), each = nrow(Z))
    synthetic_continuous[, "var"] <- rep(Z[, "var"], 4)
    synthetic_continuous[, "chain"] <- rep(seq(1, num_chains), each = nrow(Z)*4)
    
    #medians from synthetic datasets
    for(group in 1:4){
      subsample <- synthetic_sample %>% filter(Group == group)
      
      for(var in Z[, "var"]){
        subset <- subsample %>% 
          dplyr::select(!!as.symbol(var), "sample", "chain") 
        subset[, var] <- subset[, var]*orig_sds[var] + orig_means[var]
        synthetic_continuous[which(synthetic_continuous$group == group & 
                                     synthetic_continuous$var == var & 
                                     synthetic_continuous$chain %in% 
                                     seq(1, num_chains)), 
                             seq(1, num_samples)] <- 
          matrix(subset %>% group_by(chain, sample) %>%
                   summarize_all(median) %>% ungroup() %>%
                   dplyr::select(all_of(var)) %>% unlist(),
                 nrow = num_chains, ncol = num_samples, byrow = TRUE)
      }
    }
    
    synthetic_continuous %<>% 
      mutate("truth" = 0, 
             "label" = rep(Z[, "label"], 4*num_chains), 
             "group_label" = case_when(group == 1 ~ "Unimpaired", 
                                       group == 2 ~ "Other", 
                                       group == 3 ~ "MCI", 
                                       group == 4 ~ "Dementia")) %>% 
      mutate("color" = case_when(group_label == "Unimpaired" ~ "#00a389", 
                                 group_label == "Other" ~ "#28bed9", 
                                 group_label == "MCI" ~ "#fdab00", 
                                 group_label == "Dementia" ~ "#ff0000"))
    
    for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
      subsample <- dataset_to_copy %>% 
        filter(!!as.symbol(dementia_var) == group)
      
      for(var in Z[, "var"]){
        synthetic_continuous[which(synthetic_continuous$group_label == group & 
                                     synthetic_continuous$var == var), 
                             "truth"] <- 
          median(unlist(subsample[, var])*orig_sds[var] + orig_means[var])
      }
    }
    
    synthetic_continuous %<>% pivot_longer(seq(1:num_samples))
    
    #---- ****plot ----
    for(chain_num in 1:num_chains){
      for(dem_group in unique(synthetic_continuous$group_label)){
        for(var_name in Z[, "label"]){
          subset <- synthetic_continuous %>% 
            filter(group_label == dem_group & label == var_name & 
                     chain == chain_num)
          ggplot(data = subset , aes(x = value)) + 
            geom_histogram(fill = "black", color = "black") + theme_minimal() + 
            xlab("Median") + ggtitle(var_name) + 
            geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                       size = 2)
          
          ggsave(filename = paste0(path_to_figures_folder, 
                                   "posterior_predictive_checks/run_", 
                                   chain_num, "/continuous_vars/median/", 
                                   tolower(dem_group), "/", tolower(dem_group), 
                                   "_", var_name, ".jpeg"), 
                 width = 5, height = 3, units = "in")
        } 
      }
    }
    
    #---- **skew ----
    #---- ****overall ----
    synthetic_continuous <- 
      matrix(0, nrow = nrow(Z)*num_chains, ncol = (num_samples + 2)) %>% 
      as.data.frame() %>% set_colnames(c(seq(1, num_samples), "var", "chain"))
    synthetic_continuous[, "var"] <- rep(Z[, "var"], num_chains)
    synthetic_continuous[, "chain"] <- 
      rep(seq(1, num_chains), each = length(Z[, "var"]))
    
    #skew from synthetic datasets
    for(var in Z[, "var"]){
      subset <- synthetic_sample %>% 
        dplyr::select(!!as.symbol(var), "sample", "chain") 
      subset[, var] <- subset[, var]*orig_sds[var] + orig_means[var]
      synthetic_continuous[which(synthetic_continuous$var == var & 
                                   synthetic_continuous$chain %in% 
                                   seq(1, num_chains)), 
                           seq(1, num_samples)] <- 
        matrix(subset %>% group_by(chain, sample) %>% 
                 summarize_all(skewness) %>% ungroup() %>% 
                 dplyr::select(all_of(var)) %>% unlist(), 
               nrow = num_chains, ncol = num_samples, byrow = TRUE)
    }
    
    synthetic_continuous %<>% 
      mutate("truth" = 0, 
             "label" = rep(Z[, "label"], num_chains))
    
    for(var in Z[, "var"]){
      synthetic_continuous[which(synthetic_continuous$var == var), "truth"] <- 
        skewness(unlist(dataset_to_copy[, var]*orig_sds[var] + orig_means[var]))
    }
    
    synthetic_continuous %<>% pivot_longer(seq(1:num_samples))
    
    #---- ****plot ----
    for(chain_num in seq(1:num_chains)){
      for(var_name in Z[, "label"]){
        subset <- synthetic_continuous %>% 
          filter(label == var_name & chain == chain_num)
        ggplot(data = subset , aes(x = value)) + 
          geom_histogram(fill = "black", color = "black") + theme_minimal() + 
          xlab("Skew") + ggtitle(var_name) + 
          geom_vline(xintercept = subset$truth, color = "#f2caaa" , size = 2)
        
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/run_", chain_num, 
                                 "/continuous_vars/skew/overall/", var_name, 
                                 ".jpeg"), width = 5, height = 3, units = "in")
      }
    }
    
    #---- ****by class ----
    synthetic_continuous <- 
      matrix(0, nrow = nrow(Z)*4*num_chains, ncol = (num_samples + 3)) %>% 
      as.data.frame() %>% 
      set_colnames(c(seq(1, num_samples), "group", "var", "chain"))
    synthetic_continuous[, "group"] <- rep(seq(1, 4), each = nrow(Z))
    synthetic_continuous[, "var"] <- rep(Z[, "var"], 4)
    synthetic_continuous[, "chain"] <- rep(seq(1, num_chains), each = nrow(Z)*4)
    
    #skewness from synthetic datasets
    for(group in 1:4){
      subsample <- synthetic_sample %>% filter(Group == group)
      
      for(var in Z[, "var"]){
        subset <- subsample %>% 
          dplyr::select(!!as.symbol(var), "sample", "chain") 
        subset[, var] <- subset[, var]*orig_sds[var] + orig_means[var]
        synthetic_continuous[which(synthetic_continuous$group == group & 
                                     synthetic_continuous$var == var & 
                                     synthetic_continuous$chain %in% 
                                     seq(1, num_chains)), 
                             seq(1, num_samples)] <- 
          matrix(subset %>% group_by(chain, sample) %>%
                   summarize_all(skewness) %>% ungroup() %>%
                   dplyr::select(all_of(var)) %>% unlist(),
                 nrow = num_chains, ncol = num_samples, byrow = TRUE)
      }
    }
    
    synthetic_continuous %<>% 
      mutate("truth" = 0, 
             "label" = rep(Z[, "label"], 4*num_chains), 
             "group_label" = case_when(group == 1 ~ "Unimpaired", 
                                       group == 2 ~ "Other", 
                                       group == 3 ~ "MCI", 
                                       group == 4 ~ "Dementia")) %>% 
      mutate("color" = case_when(group_label == "Unimpaired" ~ "#00a389", 
                                 group_label == "Other" ~ "#28bed9", 
                                 group_label == "MCI" ~ "#fdab00", 
                                 group_label == "Dementia" ~ "#ff0000"))
    
    for(group in unlist(unique(dataset_to_copy[, dementia_var]))){
      subsample <- dataset_to_copy %>% 
        filter(!!as.symbol(dementia_var) == group)
      
      for(var in Z[, "var"]){
        synthetic_continuous[which(synthetic_continuous$group_label == group & 
                                     synthetic_continuous$var == var), 
                             "truth"] <- 
          skewness(unlist(subsample[, var])*orig_sds[var] + orig_means[var])
      }
    }
    
    synthetic_continuous %<>% pivot_longer(seq(1:num_samples))
    
    #---- ****plot ----
    for(chain_num in 1:num_chains){
      for(dem_group in unique(synthetic_continuous$group_label)){
        for(var_name in Z[, "label"]){
          subset <- synthetic_continuous %>% 
            filter(group_label == dem_group & label == var_name & 
                     chain == chain_num)
          ggplot(data = subset , aes(x = value)) + 
            geom_histogram(fill = "black", color = "black") + theme_minimal() + 
            xlab("Skew") + ggtitle(var_name) + 
            geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                       size = 2)
          
          ggsave(filename = paste0(path_to_figures_folder, 
                                   "posterior_predictive_checks/run_", 
                                   chain_num, "/continuous_vars/skew/", 
                                   tolower(dem_group), "/", tolower(dem_group), 
                                   "_", var_name, ".jpeg"), 
                 width = 5, height = 3, units = "in")
        } 
      }
    }
    
    #---- impairment classification ----
    #truth
    dementia_plot_data <- 
      as.data.frame(table(dataset_to_copy[, dementia_var])) %>% 
      mutate("prop" = Freq/sum(Freq))
    
    #synthetic
    dem_sub <- synthetic_sample[, c("Group", "sample", "chain")] %>% 
      mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                       Group == 2 ~ "Other", 
                                       Group == 3 ~ "MCI", 
                                       TRUE ~ "Dementia"))
    
    synthetic_dementia_plot_data <- 
      dem_sub %>% dplyr::count(chain, sample, Group_label) %>%
      group_by(chain, sample) %>%
      mutate(prop = n/sum(n)) %>% 
      mutate_at("sample", as.numeric) %>% 
      mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                                 Group_label == "Other" ~ "#28bed9", 
                                 Group_label == "MCI" ~ "#fdab00", 
                                 TRUE ~ "#ff0000"))
    
    #---- **combined plot ----
    combined_plot_data <- synthetic_dementia_plot_data %>% ungroup() %>% 
      group_by(chain, Group_label) %>% 
      summarize_at("n", list("mean" = mean, 
                             "lower" = ~ quantile(.x, probs = 0.025), 
                             "upper" = ~ quantile(.x, probs = 0.975))) %>% 
      mutate("truth" = dementia_plot_data$Freq) %>% 
      mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                                 Group_label == "Other" ~ "#28bed9", 
                                 Group_label == "MCI" ~ "#fdab00", 
                                 Group_label == "Dementia" ~ "#ff0000")) %>% 
      mutate("chain" = paste0("Chain ", chain))
    
    
    combined_plot_data$Group_label <- 
      factor(combined_plot_data$Group_label, 
             levels = c("Unimpaired", "Other", "MCI", "Dementia"))
    
    ggplot(data = combined_plot_data, 
           aes(x = Group_label, y = mean, color = Group_label)) + 
      geom_point(aes(size = 1)) + theme_bw() + 
      geom_point(aes(x = Group_label, y = truth, color = Group_label, size = 1), 
                 shape = 18) + 
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.10) + 
      xlab("") + ylab("Mean Count") + theme(legend.position = "none") + 
      scale_color_manual(values = rev(c(combined_plot_data$color))) + 
      facet_wrap(facets = as.factor(combined_plot_data$chain), 
                 nrow = 2, ncol = 3) +
      ggtitle(paste0("95% Credible intervals from ", num_samples, 
                     " synthetic datasets"))
    
    ggsave(filename = paste0(path_to_figures_folder, 
                             "posterior_predictive_checks/", 
                             "impairment_classes_combined_all_runs.jpeg"), 
           height = 5, width = 10, units = "in")
    
    #---- ** individual plots ----
    for(chain_num in 1:num_chains){
      ggplot(data = combined_plot_data %>% 
               filter(chain == paste0("Chain ", chain_num)), 
             aes(x = Group_label, y = mean, color = Group_label)) + 
        geom_point(aes(size = 1)) + theme_minimal() + 
        geom_point(aes(x = Group_label, y = truth, color = Group_label, 
                       size = 1), shape = 18) + 
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.10) + 
        xlab("") + ylab("Mean Count") + theme(legend.position = "none") + 
        scale_color_manual(values = rev(c(combined_plot_data$color))) + 
        ggtitle(paste0("95% Credible intervals from ", num_samples, 
                       " synthetic datasets"))
      
      ggsave(filename = paste0(path_to_figures_folder, 
                               "posterior_predictive_checks/run_", chain_num,  
                               "/impairment_classes.jpeg"), 
             height = 5, width = 10, units = "in")
      
    }
  }

