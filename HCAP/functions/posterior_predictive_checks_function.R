posterior_predictive_checks <- 
  function(dataset_to_copy, calibration_sample = FALSE, calibration_prop, 
           calibration_sample_name, categorical_covariates, 
           continuous_covariates, contrasts_matrix, cell_ID_key, variable_labels, 
           color_palette, path_to_analyses_folder, path_to_figures_folder){
    
    if(!calibration_sample){
      #---- read in synthetic data ----
      file_paths <- 
        list.dirs(path = paste0(path_to_analyses_folder, "synthetic_data/", 
                                "ADAMS_prior"), 
                  full.names = TRUE, recursive = FALSE)
      
      for(chain in 1:length(file_paths)){
        data_file <- list.files(file_paths[chain], full.names = TRUE)
        synthetic_sample <- readRDS(data_file)
        
        #---- **formatting ----
        #add chain number
        synthetic_sample <-
          lapply(synthetic_sample, function(x) x %<>% mutate("chain" = chain))
        
        #add sample number
        sample_nums <- as.list(seq(1, length(synthetic_sample)))
        synthetic_sample <- Map(cbind, synthetic_sample, sample = sample_nums)
      }
      
      synthetic_sample %<>% do.call(rbind, .)
    } else{
      #---- read in synthetic data ----
      file_paths <- 
        list.dirs(path = paste0(path_to_analyses_folder, "synthetic_data/", 
                                calibration_sample_name), 
                  full.names = TRUE, recursive = FALSE)
      
      for(chain in 1:length(file_paths)){
        data_file <- list.files(file_paths[chain], full.names = TRUE)
        synthetic_sample <- readRDS(data_file)
        
        #---- **formatting ----
        #add chain number
        synthetic_sample <-
          lapply(synthetic_sample, function(x) x %<>% mutate("chain" = chain))
        
        #add sample number
        sample_nums <- as.list(seq(1, length(synthetic_sample)))
        synthetic_sample <- Map(cbind, synthetic_sample, sample = sample_nums)
      }
      
      synthetic_sample %<>% do.call(rbind, .)
    }
    
    #---- number of chains and samples ----
    num_chains <- max(synthetic_sample$chain)
    num_samples <- max(synthetic_sample$sample)
    
    #---- create directories for results ----
    for(chain in 1:num_chains){
      if(!dir.exists(paste0(path_to_figures_folder, 
                            "posterior_predictive_checks/", 
                            "ADAMS_prior/run_", chain, 
                            "/cell_counts/"))){
        dir.create(paste0(path_to_figures_folder, 
                          "posterior_predictive_checks/", 
                          "ADAMS_prior/run_", chain, 
                          "/cell_counts/"), recursive = TRUE)
      }
    }
    
    for(chain in 1:num_chains){
      for(metric in c("median", "skew")){
        if(!dir.exists(paste0(path_to_figures_folder, 
                              "posterior_predictive_checks/", 
                              "ADAMS_prior/run_", chain, 
                              "/continuous_vars/", metric, "/"))){
          dir.create(paste0(path_to_figures_folder, 
                            "posterior_predictive_checks/", 
                            "ADAMS_prior/run_", chain, 
                            "/continuous_vars/", metric, "/"), 
                     recursive = TRUE)
        }
      }
    }
    
    #---- chain convergence ----
    if(!calibration_sample){
      #---- **read in data ----
      file_paths <- 
        list.dirs(path = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "ADAMS_prior"), full.names = TRUE, 
                  recursive = FALSE)
      
      group_membership <- 
        lapply(file_paths, function(x) 
          read_csv(paste0(x, "/latent_class_data.csv")))
    } 
    
    #---- **formatting ----
    #add chain number
    chain_nums <- as.list(seq(1, length(file_paths)))
    group_membership <- Map(cbind, group_membership, chain = chain_nums)
    
    group_membership %<>% do.call(rbind, .)
    
    #---- **plot ----
    chain_convergence <- 
      ggplot(data = group_membership, 
             aes(x = run, y = prob, colour = as.factor(chain))) +       
      geom_line(aes(group = as.factor(chain)), alpha = 0.60) + xlab("Run") + 
      ylab("Proportion") + guides(color = guide_legend(title = "Chain")) + 
      scale_color_manual(values = rev(wes_palette("Darjeeling2"))) +
      theme_bw() +
      facet_wrap(facets = 
                   vars(factor(Group, 
                               levels = c("Unimpaired", "MCI", "Dementia", 
                                          "Other"))), nrow = 2, ncol = 2) + 
      theme(text = element_text(size = 12),
            strip.text = element_text(size = 12)) + 
      theme(legend.position = "bottom")  
    
    if(!calibration_sample){
      ggsave(filename = "group_membership_chains.jpeg", plot = chain_convergence, 
             path = paste0(path_to_figures_folder, "posterior_predictive_checks/", 
                           "ADAMS_prior/"), 
             width = 10, height = 5, units = "in", device = "jpeg")  
    } else{
      ggsave(filename = "group_membership_chains.jpeg", plot = chain_convergence, 
             path = paste0(path_to_figures_folder, "posterior_predictive_checks/", 
                           calibration_sample_name, "/"), 
             width = 10, height = 5, units = "in", device = "jpeg")  
    }
    
    #---- categorical checks ----
    #---- **race/ethnicity x stroke overall ----
    #---- ****pre-allocate results matrix ----
    #---- ****synthetic counts ----
    synthetic_counts <- 
      matrix(0, nrow = nrow(contrasts_matrix)*num_chains, 
             ncol = (num_samples + 2)) %>% as.data.frame() %>% 
      set_colnames(c(seq(1, num_samples), "cell", "chain"))
    
    cells <- contrasts_matrix[, -1] %>% as.data.frame() %>% 
      unite("cell_ID", sep = "") %>% left_join(cell_ID_key, by = "cell_ID")
    
    synthetic_counts[, "cell"] <- rep(cells$cell_name, num_chains)
    
    synthetic_counts[, "chain"] <- rep(seq(1, num_chains), 
                                       each = nrow(contrasts_matrix)) 
    
    for(chain_num in unique(synthetic_sample$chain)){
      #counts from synthetic datasets
      for(num in 1:num_samples){
        subsample <- synthetic_sample %>% 
          filter(chain == chain_num & sample == num) 
        
        counts <- subsample %>% 
          dplyr::select(all_of(categorical_covariates)) %>% 
          unite("cell_ID", sep = "") %>% table() %>% as.data.frame() %>% 
          set_colnames(c("cell_ID", "Freq")) %>% 
          left_join(., cell_ID_key, by = "cell_ID")
        
        counts[which(is.na(counts$Freq)), "Freq"] <- 0
        
        synthetic_counts[
          which(synthetic_counts$chain == chain_num & 
                  synthetic_counts$cell %in% counts$cell_name), num] <- 
          counts$Freq
      }
    }
    
    synthetic_count_plot_data <- synthetic_counts %>% mutate("truth" = 0) 
    
    #true counts
    counts <- dataset_to_copy %>% 
      dplyr::select(all_of(categorical_covariates)) %>% 
      unite("cell_ID", sep = "") %>% table() %>% as.data.frame() %>% 
      set_colnames(c("cell_ID", "Freq")) %>% 
      left_join(cell_ID_key, ., by = "cell_ID")  
    
    counts[which(is.na(counts$Freq)), "Freq"] <- 0
    
    synthetic_count_plot_data[
      which(synthetic_count_plot_data$cell %in% counts$cell_name), "truth"] <- 
      counts$Freq
    
    synthetic_count_plot_data %<>% 
      pivot_longer(-c("cell", "chain", "truth")) 
    
    #---- **plot ----
    for(chain in 1:num_chains){
      subset <- synthetic_count_plot_data %>% filter(chain == chain)
      ggplot(data = subset , aes(x = value)) + 
        geom_histogram(fill = "black", color = "black") + theme_bw() + 
        xlab("Contingency Cell Count") + ylab("") + 
        facet_wrap(facets = vars(cell), ncol = 2, scales = "free") +
        geom_vline(aes(xintercept = truth), color = "#f2caaa", size = 1) +
        theme(text = element_text(size = 6), 
              strip.text = element_text(size = 6))  
      
      if(!calibration_sample){
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 "ADAMS_prior/run_", chain, 
                                 "/cell_counts/overall_count.jpeg"), 
               width = 3, height = 4, units = "in")
      } else{
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 calibration_sample_name, "/run_", chain, 
                                 "/cell_counts/overall_count.jpeg"), 
               width = 3, height = 4, units = "in")
      }
    }
    
    #---- continuous checks ----
    #do these on the original scale
    continuous_covariates <- str_remove(continuous_covariates, "_Z")
    
    #---- **density ----
    #---- **median ----
    #---- ****overall ----
    synthetic_continuous <- 
      matrix(0, nrow = length(continuous_covariates)*num_chains, 
             ncol = (num_samples + 2)) %>% 
      as.data.frame() %>% set_colnames(c(seq(1, num_samples), "var", "chain"))
    synthetic_continuous[, "var"] <- rep(continuous_covariates, num_chains)
    synthetic_continuous[, "chain"] <- 
      rep(seq(1, num_chains), each = length(continuous_covariates))
    
    #medians from synthetic datasets
    for(var in continuous_covariates){
      subset <- synthetic_sample %>% 
        dplyr::select(!!as.symbol(var), "sample", "chain") 
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
             "label" = 
               rep(unlist(variable_labels[
                 which(unique(variable_labels$data_label) %in% 
                         continuous_covariates), 
                 "figure_label"]), num_chains))
    
    for(var in continuous_covariates){
      synthetic_continuous[which(synthetic_continuous$var == var), "truth"] <- 
        median(unlist(dataset_to_copy[, var]))  
    }
    
    synthetic_continuous %<>% pivot_longer(seq(1:num_samples))
    
    #---- ****plot ----
    for(chain_num in seq(1:num_chains)){
      subset <- synthetic_continuous %>% filter(chain == chain_num)
      
      ggplot(data = subset , aes(x = value)) + 
        geom_histogram(fill = "black", color = "black") + theme_bw() + 
        xlab("Median") + 
        facet_wrap(facets = vars(label), ncol = 2, scales = "free") +
        geom_vline(aes(xintercept = truth), color = "#f2caaa", size = 2)
      
      if(!calibration_sample){
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 "ADAMS_prior/run_", chain_num, 
                                 "/continuous_vars/median/overall.jpeg"), 
               width = 8, height = 10, units = "in")
      } else{
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 calibration_sample_name, "/run_", chain_num, 
                                 "/continuous_vars/median/overall.jpeg"), 
               width = 8, height = 10, units = "in")
      }
    }
    
    #---- **skew ----
    #---- ****overall ----
    synthetic_continuous <- 
      matrix(0, nrow = length(continuous_covariates)*num_chains, 
             ncol = (num_samples + 2)) %>% 
      as.data.frame() %>% set_colnames(c(seq(1, num_samples), "var", "chain"))
    synthetic_continuous[, "var"] <- rep(continuous_covariates, num_chains)
    synthetic_continuous[, "chain"] <- 
      rep(seq(1, num_chains), each = length(continuous_covariates))
    
    #skew from synthetic datasets
    for(var in continuous_covariates){
      subset <- synthetic_sample %>% 
        dplyr::select(!!as.symbol(var), "sample", "chain") 
      
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
             "label" = rep(unlist(variable_labels[
               which(unique(variable_labels$data_label) %in% 
                       continuous_covariates), 
               "figure_label"]), num_chains))
    
    for(var in continuous_covariates){
      synthetic_continuous[which(synthetic_continuous$var == var), "truth"] <- 
        skewness(unlist(dataset_to_copy[, var]))  
    }
    
    synthetic_continuous %<>% pivot_longer(seq(1:num_samples))
    
    #---- ****plot ----
    for(chain_num in seq(1:num_chains)){
      subset <- synthetic_continuous %>% filter(chain == chain_num)
      
      ggplot(data = subset , aes(x = value)) + 
        geom_histogram(fill = "black", color = "black") + theme_bw() + 
        xlab("Skew") + 
        facet_wrap(facets = vars(label), ncol = 2, scales = "free") +
        geom_vline(aes(xintercept = truth), color = "#f2caaa" , size = 2)
      
      if(!calibration_sample){
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 "ADAMS_prior/run_", chain_num, 
                                 "/continuous_vars/skew/overall.jpeg"), 
               width = 8, height = 10, units = "in")  
      } else{
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 calibration_sample_name, "/run_", chain_num, 
                                 "/continuous_vars/skew/overall.jpeg"), 
               width = 8, height = 10, units = "in")
      }
    }
    
    #---- impairment classification ----
    #synthetic
    synthetic_dementia_plot_data <- 
      synthetic_sample[, c("Group", "sample", "chain")]  %>% 
      mutate("Unimpaired" = ifelse(Group == "Unimpaired", 1, 0), 
             "MCI" = ifelse(Group == "MCI", 1, 0), 
             "Dementia" = ifelse(Group == "Dementia", 1, 0), 
             "Other" = ifelse(Group == "Other", 1, 0)) %>%
      dplyr::count(chain, sample, Unimpaired, MCI, Dementia, Other) %>% 
      pivot_longer(c("Unimpaired", "MCI", "Dementia", "Other"), 
                   names_to = "Group") %>% filter(value == 1) %>%
      group_by(chain, sample) %>%
      mutate_at("sample", as.numeric) 
    
    #---- **combined plot ----
    combined_plot_data <- synthetic_dementia_plot_data %>% ungroup() %>% 
      group_by(chain, Group) %>% 
      summarize_at("n", list("mean" = mean, 
                             "lower" = ~ quantile(.x, probs = 0.025), 
                             "upper" = ~ quantile(.x, probs = 0.975))) %>% 
      left_join(color_palette, by = "Group") %>%
      mutate("chain" = paste0("Chain ", chain))
    
    combined_plot_data$Group <- 
      factor(combined_plot_data$Group, 
             levels = c("Unimpaired", "MCI", "Dementia", "Other"))
    
    #---- ****count ----
    ggplot(data = combined_plot_data, 
           aes(x = Group, y = mean)) + theme_bw() + 
      geom_errorbar(aes(ymin = lower, ymax = upper, color = Group), 
                    width = 0.10) + 
      geom_point(aes(size = 1, color = Group)) +
      xlab("") + ylab("Mean Count") + theme(legend.position = "none") + 
      scale_color_manual(
        values = combined_plot_data$Color[order(combined_plot_data$Group)]) + 
      facet_wrap(facets = as.factor(combined_plot_data$chain), 
                 nrow = 2, ncol = 3) +
      ggtitle(paste0("95% Credible intervals from ", num_samples, 
                     " synthetic datasets"))
    
    if(!calibration_sample){
      ggsave(filename = paste0(path_to_figures_folder, 
                               "posterior_predictive_checks/", 
                               "ADAMS_prior/", 
                               "impairment_classes_combined_all_runs_count.jpeg"), 
             height = 5, width = 10, units = "in")
    } else{
      ggsave(filename = paste0(path_to_figures_folder, 
                               "posterior_predictive_checks/",
                               calibration_sample_name,
                               "/impairment_classes_combined_all_runs_count.jpeg"), 
             height = 5, width = 10, units = "in")
    }
    
    #---- ****prop ----
    ggplot(data = combined_plot_data %>% 
             mutate_at(c("mean", "lower", "upper"), 
                       function(x) x/nrow(dataset_to_copy)), 
           aes(x = Group, y = mean)) + theme_bw() + 
      geom_errorbar(aes(ymin = lower, ymax = upper, color = Group), 
                    width = 0.10) + 
      geom_point(aes(size = 1, color = Group)) +
      xlab("") + ylab("Mean Proportion") + theme(legend.position = "none") + 
      scale_color_manual(
        values = combined_plot_data$Color[order(combined_plot_data$Group)]) + 
      facet_wrap(facets = as.factor(combined_plot_data$chain), 
                 nrow = 2, ncol = 3) +
      ggtitle(paste0("95% Credible intervals from ", num_samples, 
                     " synthetic datasets"))
    
    if(!calibration_sample){
      ggsave(filename = paste0(path_to_figures_folder, 
                               "posterior_predictive_checks/", 
                               "ADAMS_prior/", 
                               "impairment_classes_combined_all_runs_prop.jpeg"), 
             height = 5, width = 10, units = "in")
    } else{
      ggsave(filename = paste0(path_to_figures_folder, 
                               "posterior_predictive_checks/",
                               calibration_sample_name,
                               "/impairment_classes_combined_all_runs_prop.jpeg"), 
             height = 5, width = 10, units = "in")
    }
    
    #---- **individual plots ----
    for(chain_num in 1:num_chains){
      #---- ****count ----
      ggplot(data = combined_plot_data %>% 
               filter(chain == paste0("Chain ", chain_num)), 
             aes(x = Group, y = mean)) + 
        geom_errorbar(aes(ymin = lower, ymax = upper, color = Group), 
                      width = 0.10) + 
        geom_point(aes(size = 1, color = Group)) + 
        theme_minimal() + 
        xlab("") + ylab("Count") + 
        theme(text = element_text(size = 10), legend.title = element_blank(), 
              legend.position = "none", 
              plot.title = element_text(size = 10)) +
        scale_color_manual(
          values = combined_plot_data$Color[order(combined_plot_data$Group)]) + 
        ggtitle(paste0("95% Credible intervals from ", num_samples, 
                       " synthetic datasets")) + 
        guides(color = "none") 
      
      if(!calibration_sample){
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 "ADAMS_prior/run_", chain_num,  
                                 "/impairment_classes_count.jpeg"), 
               height = 4, width = 5.5, units = "in")
      } else{
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 calibration_sample_name, "/run_", chain_num,  
                                 "/impairment_classes_count.jpeg"), 
               height = 4, width = 5.5, units = "in")
      }
      
      #---- ****prop ----
      ggplot(data = combined_plot_data %>% 
               filter(chain == paste0("Chain ", chain_num)) %>% 
               mutate_at(c("mean", "lower", "upper"), 
                         function(x) x/nrow(dataset_to_copy)), 
             aes(x = Group, y = mean)) + 
        geom_errorbar(aes(ymin = lower, ymax = upper, color = Group), 
                      width = 0.10) + 
        geom_point(aes(size = 1, color = Group)) + 
        theme_minimal() + 
        xlab("") + ylab("Count") + 
        theme(text = element_text(size = 10), legend.title = element_blank(), 
              legend.position = "none", 
              plot.title = element_text(size = 10)) +
        scale_color_manual(
          values = combined_plot_data$Color[order(combined_plot_data$Group)]) + 
        ggtitle(paste0("95% Credible intervals from ", num_samples, 
                       " synthetic datasets")) + 
        guides(color = "none") 
      
      if(!calibration_sample){
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 "ADAMS_prior/run_", chain_num,  
                                 "/impairment_classes_prop.jpeg"), 
               height = 4, width = 5.5, units = "in")
      } else{
        ggsave(filename = paste0(path_to_figures_folder, 
                                 "posterior_predictive_checks/", 
                                 calibration_sample_name, "/run_", chain_num,  
                                 "/impairment_classes_prop.jpeg"), 
               height = 4, width = 5.5, units = "in")
      }
    }
  }

# #---- test function ----
# dataset_to_copy = HCAP_analytic
# calibration_sample = FALSE
# calibration_prop = NA
# calibration_sample_name = NA
# categorical_covariates = W
# continuous_covariates = Z
# contrasts_matrix = A
# cell_ID_key = cell_ID_key
# variable_labels = variable_labels
# color_palette = color_palette
# path_to_analyses_folder =
#   paste0(path_to_box, "analyses/HCAP/")
# path_to_figures_folder =
#   paste0(path_to_box,"figures/chapter_6/")
