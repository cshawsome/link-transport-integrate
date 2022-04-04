prior_predictive_checks <- 
  function(unimpaired_preds, other_preds, mci_preds, categorical_vars, 
           continuous_vars, variable_labels, color_palette, dataset_to_copy, 
           num_synthetic, unimpaired_betas, unimpaired_cov, other_betas, 
           other_cov, mci_betas, mci_cov, alpha_0_dist, prior_Sigma, 
           prior_V_inv, prior_beta, nu_0, kappa_0, contrasts_matrix, 
           path_to_folder){
    
    #---- create folders for results ----
    if(!dir.exists(paste0(path_to_folder, "impairment_classes/"))){
      dir.create(paste0(path_to_folder, "impairment_classes/"), 
                 recursive = TRUE) 
    }
    
    for(class in c("unimpaired", "mci", "dementia", "other")){
      if(!dir.exists(paste0(path_to_folder, "continuous_vars/", 
                            tolower(class)))){
        dir.create(paste0(path_to_folder, "continuous_vars/", 
                          tolower(class)), recursive = TRUE)
      }
      
    }
    
    #---- select variables ----
    vars <- unique(c(unimpaired_preds, other_preds, mci_preds, 
                     "Unimpaired", "MCI", "Dementia", "Other"))
    
    synthetic_sample <- dataset_to_copy %>% dplyr::select(all_of(vars)) %>% 
      #pre-allocate columns
      mutate("group_num" = 0, "p_unimpaired" = 0, "p_other" = 0, "p_mci" = 0)
    
    generate_data <- function(color_palette){
      #---- latent class ----
      group_num = 1
      for(model in c("unimpaired", "other", "mci")){
        subset_index <- which(synthetic_sample$group_num == 0)
        random_draw <- sample(seq(1, ncol(get(paste0(model, "_cov")))), 
                              size = 1)
        
        prior_betas <- as.vector(get(paste0(model, "_betas"))[, random_draw])
        prior_cov <- matrix(unlist(get(paste0(model, "_cov"))[, random_draw]), 
                            nrow = nrow(prior_betas))
        
        betas <- mvrnorm(n = 1, mu = t(prior_betas), Sigma = prior_cov)
        
        synthetic_sample[subset_index, paste0("p_", model)] <- 
          expit(as.matrix(synthetic_sample[subset_index, 
                                           get(paste0(model, "_preds"))]) %*% 
                  as.matrix(betas))
        
        synthetic_sample[subset_index, "group_num"] <- 
          rbernoulli(n = length(subset_index), 
                     p = synthetic_sample[subset_index, 
                                          paste0("p_", model)])*group_num
        
        group_num = group_num + 1
      }
      
      synthetic_sample[, "Group"] <- 
        case_when(synthetic_sample$group_num == 1 ~ "Unimpaired", 
                  synthetic_sample$group_num == 2 ~ "Other", 
                  synthetic_sample$group_num == 3 ~ "MCI", 
                  synthetic_sample$group_num == 0 ~ "Dementia")
      
      #pre-allocate: ncol = num impairement groups * num contingency cells
      # these are contingency-cell specific means for continuous variables by
      # impairment class
      mu <- matrix(0, ncol = 4*nrow(contrasts_matrix), 
                   nrow = length(continuous_vars)) %>%
        set_colnames(apply(
          expand.grid(c("Unimpaired", "MCI", "Dementia", "Other"), 
                      seq(1, nrow(contrasts_matrix))), 1, paste0,
          collapse = ":"))
      
      #figure out max sampling index-- could use any prior distribution
      max_index <- 
        colnames(alpha_0_dist)[str_detect(colnames(alpha_0_dist), "[0-9]+")] %>% 
        as.numeric() %>% max()
      
      for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
        #---- **index for random draws ----
        random_draw <- sample(seq(1, max_index), size = 1)
        
        #---- **contingency cells ----
        subset <- synthetic_sample %>% filter(Group == class)
        prior_counts <- 
          alpha_0_dist[, c(as.character(random_draw), "group")] %>% 
          filter(group == class)
        
        #---- **p(contingency table cell) ----
        pi <- 
          MCMCpack::rdirichlet(1, alpha = as.numeric(unlist(prior_counts[, 1]))*
                                 nrow(subset))
        
        #---- **contingency table count ----
        contingency_table <- rmultinom(n = 1, size = nrow(subset), prob = pi)
        UtU <- diag(contingency_table[, 1])
        
        #---- **draw new UtU if needed ----
        while(det(t(contrasts_matrix) %*% UtU %*% contrasts_matrix) < 1e-9){
          #---- ****new random draw index ----
          random_draw <- sample(seq(1, max_index), size = 1)
          new_counts <- 
            alpha_0_dist[, c(as.character(random_draw), "group")] %>% 
            filter(group == class)
          
          UtU <- diag(unlist(new_counts[, 1])*nrow(subset))
        }
        
        #---- **make U matrix ----
        U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
        
        for(j in 1:nrow(contingency_table)){
          if(contingency_table[j, 1] == 0){next}
          if(j == 1){
            index = 1
          } else{
            index = sum(contingency_table[1:(j - 1), 1]) + 1
          }
          U[index:(index - 1 + contingency_table[j, 1]), j] <- 1
        }
        
        #---- **draw Sigma_0----
        Sigma_prior <- prior_Sigma[, c(as.character(random_draw), "group")] %>% 
          filter(group == class)
        sig_Y <- riwish(v = (nu_0), S = matrix(unlist(Sigma_prior[, 1]), 
                                               nrow = length(continuous_vars)))
        
        #---- **beta_0 ----
        V_0_inv <- prior_V_inv[, c(as.character(random_draw), "group")] %>% 
          filter(group == class)
        beta_0 <- priors_beta[, c(as.character(random_draw), "group")] %>% 
          filter(group == class)
        
        #as matrices
        V_0_inv <- matrix(unlist(V_0_inv[, 1]), nrow = ncol(contrasts_matrix))
        beta_0 <- matrix(unlist(beta_0[, 1]), nrow = nrow(V_0_inv))
        
        #---- **draw beta | Sigma----
        beta_Sigma_Y <- matrix.normal(beta_0, solve(V_0_inv), 
                                      sig_Y/kappa_0[class])
        
        #---- **compute mu ----
        mu[, paste0(class, ":", seq(1, 6))] <-
          t(contrasts_matrix %*% 
              matrix(beta_Sigma_Y, nrow = ncol(contrasts_matrix), 
                     ncol = length(continuous_vars), byrow = FALSE))
        
        #---- **draw data ----
        #reformat contingency table
        table <- contingency_table %>% as.data.frame() %>%
          cbind(., contrasts_matrix[, -1]) %>% rename("Count" = "V1")
        
        for(j in 1:nrow(table)){
          if(table[j, "Count"] == 0){next}
          if(j == 1){
            index = 1
          } else{
            index = sum(table[1:(j - 1), "Count"]) + 1
          }
          #continuous data
          if(table[j, "Count"] == 1){
            subset[index:(index - 1 + table[j, "Count"]), continuous_vars] <-
              t(as.matrix(mvrnorm(n = table[j, "Count"],
                                  mu = mu[, paste0(class, ":", j)], 
                                  Sigma = sig_Y)))
          } else{
            subset[index:(index - 1 + table[j, "Count"]), continuous_vars] <-
              mvrnorm(n = table[j, "Count"],
                      mu = mu[, paste0(class, ":", j)], Sigma = sig_Y)
          }
        }
        assign(paste0("Z_", class), subset[, all_of(continuous_vars)])
      }
      
      #---- **return ----
      return(list("Group" = synthetic_sample$Group, 
                  "Z_unimpaired" = Z_Unimpaired %>% 
                    mutate("Group" = "Unimpaired") %>% left_join(color_palette), 
                  "Z_other" = Z_Other %>% 
                    mutate("Group" = "Other") %>% left_join(color_palette), 
                  "Z_mci" = Z_MCI %>% 
                    mutate("Group" = "MCI") %>% left_join(color_palette), 
                  "Z_dementia" = Z_Dementia %>% 
                    mutate("Group" = "Dementia") %>% left_join(color_palette)))
    }
    
    #---- multiruns ----
    #start <- Sys.time()
    runs = num_synthetic
    synthetic <- replicate(num_synthetic, generate_data(color_palette), 
                           simplify = FALSE) 
    #stop <- Sys.time() - start
    
    #---- plots ----
    #---- **dem class ----
    to_copy_dementia_plot_data <- 
      as.data.frame(table(
        dataset_to_copy[, c("Unimpaired", "MCI", "Dementia", "Other")])) %>% 
      pivot_longer(-"Freq") %>% filter(value == 1 & Freq != 0)
    
    dem_sub <- lapply(synthetic, "[[", "Group") %>% do.call(cbind, .) %>% 
      set_colnames(seq(1, runs)) %>% as.data.frame() %>%
      pivot_longer(everything()) 
    
    synthetic_dementia_plot_data <- 
      dem_sub %>% dplyr::count(name, value) %>%
      group_by(name) %>%
      mutate_at("name", as.numeric) %>% 
      left_join(., color_palette, by = c("value" = "Group"))
    
    for(class in unique(synthetic_dementia_plot_data$value)){
      subset <- synthetic_dementia_plot_data %>% filter(value == class)
      
      #original data only
      ggplot(data = subset) + 
        geom_histogram(aes(x = n),
                       color = "white", fill = "white") +
        theme_minimal() + ggtitle(class) + xlab("Count") + ylab("Frequency") +
        geom_vline(xintercept = to_copy_dementia_plot_data[[
          which(to_copy_dementia_plot_data$name == class), "Freq"]], size = 2, 
          color = unique(subset$Color))
      
      ggsave(filename = paste0(path_to_folder, "impairment_classes/", class, 
                               "_line_only.jpeg"), 
             height = 3, width = 5, units = "in")
      
      #Prior predictive check
      ggplot(data = subset) + 
        geom_histogram(aes(x = n), 
                       color = unique(subset$Color), 
                       fill = unique(subset$Color)) + 
        theme_minimal() + ggtitle(class) + xlab("Count") + ylab("Frequency") +
        geom_vline(xintercept = to_copy_dementia_plot_data[[
          which(to_copy_dementia_plot_data$name == class), "Freq"]], size = 2)
      
      ggsave(filename = paste0(path_to_folder, "impairment_classes/", class, 
                               ".jpeg"), 
             height = 3, width = 5, units = "in")
    }
    
    #---- **class-specific continuous ----
    for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
      true_data <- dataset_to_copy %>% 
        dplyr::select(c(all_of(continuous_vars), all_of(class))) %>% 
        filter(!!as.symbol(class) == 1) %>% mutate("Color" = "black")
      
      continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
      
      for(i in 1:length(continuous_list)){
        continuous_list[[i]] <- continuous_list[[i]] %>% 
          dplyr::select(-c("Group")) %>%
          mutate("run" = i, "Type" = "Synthetic") %>% 
          rbind(., true_data %>% dplyr::select(-!!as.symbol(class)) %>% 
                  mutate("run" = i, "Type" = "Observed"))
      }
      
      continuous_list %<>% do.call(rbind, .) %>% as.data.frame() 
      
      for(var in continuous_vars){
        data <- continuous_list[, c(var, "run", "Type", "Color")]
        
        continuous_plot <- 
          ggplot(data = data, aes(color = Type, fill = Type)) + 
          geom_density(aes(x = data[, 1]), alpha = 0.5) + 
          theme_minimal() + 
          xlab(variable_labels[variable_labels$data_label == var, 
                               "figure_label"]) + 
          scale_color_manual(values = rev(unique(data$Color))) + 
          scale_fill_manual(values = rev(unique(data$Color))) + 
          transition_states(data$run, transition_length = 1, state_length = 1) +
          labs(title = "Synthetic {round(frame_time)}") + transition_time(run) +
          ease_aes('linear')
        
        animate(continuous_plot, fps = 2, height = 4, width = 5, units = "in", 
                res = 150, renderer = gifski_renderer())
        
        anim_save(filename = paste0(path_to_folder, "continuous_vars/", 
                                    tolower(class), "/", var, ".gif"),
                  animation = last_animation(),
                  renderer = gifski_renderer())
      }
    }
  }

#---- test function ----
unimpaired_preds = unimpaired_preds
other_preds = other_preds
mci_preds = mci_preds
categorical_vars = W
continuous_vars = Z
variable_labels = variable_labels
color_palette = color_palette
dataset_to_copy = synthetic_data_list[[1]] %>% group_by(married_partnered) %>%
  slice_sample(prop = 0.5) %>% mutate("(Intercept)" = 1) %>% ungroup()
num_synthetic = 2
unimpaired_betas = unimpaired_betas
unimpaired_cov = unimpaired_cov
other_betas = other_betas
other_cov = other_cov
mci_betas = mci_betas
mci_cov = mci_cov
alpha_0_dist = alpha_0_dist
prior_Sigma = prior_Sigma
prior_V_inv = prior_V_inv
prior_beta = priors_beta
nu_0 = nu_0
kappa_0 = kappa_0
contrasts_matrix = A
path_to_folder = paste0(path_to_box, "figures/simulation_study/HCAP_HRS_",
                        unique(synthetic_data_list[[1]][, "dataset_name"]),
                        "/prior_predictive_checks/")
