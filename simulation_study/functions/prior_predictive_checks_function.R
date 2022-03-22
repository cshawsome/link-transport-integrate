prior_predictive_checks <- 
  function(unimpaired_preds, other_preds, mci_preds, categorical_vars, 
           continuous_vars, id_var, dementia_var, dataset_to_copy, 
           num_synthetic, unimpaired_betas, unimpaired_cov, other_betas, 
           other_cov, mci_betas, mci_cov, alpha_0_dist, prior_Sigma, prior_V_inv, 
           prior_beta, nu_0, kappa_0, contrasts_matrix, path_to_folder){
    
    #---- create folders for results ----
    dir.create(paste0(path_to_folder, "impairment_classes/"), recursive = TRUE)
    
    for(class in c("unimpaired", "mci", "dementia", "other")){
      dir.create(paste0(path_to_folder, "continuous_vars/", 
                        tolower(class)), recursive = TRUE)
    }
    
    #---- select variables ----
    vars <- unique(c(unimpaired_preds, other_preds, mci_preds, 
                     "Unimpaired", "MCI", "Dementia", "Other"))
    
    synthetic_sample <- dataset_to_copy %>%  
      dplyr::select(all_of(id_var), all_of(vars)) %>% 
      #pre-allocate columns
      mutate("group_num" = 0, "p_unimpaired" = 0, "p_other" = 0, "p_mci" = 0)
    
    generate_data <- function(){
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
      
      #pre-allocate
      mu <- matrix(0, ncol = 4*6, nrow = 10) %>%
        set_colnames(apply(expand.grid(seq(1, 4), seq(1, 6)), 1, paste0,
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
        pi <- rdirichlet(1, alpha = as.numeric(unlist(prior_counts[, 1]))*
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
                                               nrow = length(Z)))
        
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
        mu[, paste0(i, ":", seq(1, 6))] <-
          t(contrasts_matrix %*% 
              matrix(beta_Sigma_Y, nrow = ncol(contrasts_matrix), 
                     ncol = nrow(Z), byrow = FALSE))
        
        #---- **draw data ----
        #reformat contingency table
        table <- contingency_table %>% as.data.frame() %>%
          cbind(do.call(cbind, list(
            #Black              #Hispanic           
            rep(c(1, 0, 0), 2), rep(c(0, 1, 0), 2),
            #Stroke
            c(rep(0, 3), rep(1, 3))))) %>% set_colnames(c("Count", W))
        
        for(j in 1:nrow(table)){
          if(table[j, "Count"] == 0){next}
          if(j == 1){
            index = 1
          } else{
            index = sum(table[1:(j - 1), "Count"]) + 1
          }
          #Z (continuous data)
          if(table[j, "Count"] == 1){
            subset[index:(index - 1 + table[j, "Count"]), Z[, "var"]] <-
              t(as.matrix(mvrnorm(n = table[j, "Count"],
                                  mu = mu[, paste0(i, ":", j)], Sigma = sig_Y)))
          } else{
            subset[index:(index - 1 + table[j, "Count"]), Z[, "var"]] <-
              mvrnorm(n = table[j, "Count"],
                      mu = mu[, paste0(i, ":", j)], Sigma = sig_Y)
          }
        }
        assign(paste0("Z_", i), subset[, all_of(Z[, "var"])])
      }
      
      #---- **return ----
      return(list("Group" = synthetic_sample$Group,
                  "Z_unimpaired" = Z_1 %>% mutate("color" = "#00a389"), 
                  "Z_other" = Z_2 %>% mutate("color" = "#28bed9"), 
                  "Z_mci" = Z_3 %>% mutate("color" = "#fdab00"),
                  "Z_dementia" = Z_4 %>% mutate("color" = "#ff0000")))
    }
    
    #---- multiruns ----
    #start <- Sys.time()
    runs = num_synthetic
    synthetic <- replicate(num_synthetic, generate_data(), simplify = FALSE) 
    #stop <- Sys.time() - start
    
    #---- plots ----
    #---- **dem class ----
    to_copy_dementia_plot_data <- 
      as.data.frame(table(dataset_to_copy[, dementia_var])) %>% 
      mutate("prop" = Freq/sum(Freq))
    
    dem_sub <- lapply(synthetic, "[[", "Group") %>% do.call(cbind, .) %>% 
      set_colnames(seq(1, runs)) %>% as.data.frame() %>%
      pivot_longer(everything()) %>% 
      mutate("Group_label" = case_when(value == 1 ~ "Unimpaired", 
                                       value == 2 ~ "Other", 
                                       value == 3 ~ "MCI", 
                                       TRUE ~ "Dementia"))
    
    synthetic_dementia_plot_data <- 
      dem_sub %>% dplyr::count(name, Group_label) %>%
      group_by(name) %>%
      mutate(prop = n/sum(n)) %>% 
      mutate_at("name", as.numeric) %>% 
      mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                                 Group_label == "Other" ~ "#28bed9", 
                                 Group_label == "MCI" ~ "#fdab00", 
                                 TRUE ~ "#ff0000"))
    
    for(class in unique(synthetic_dementia_plot_data$Group_label)){
      subset <- synthetic_dementia_plot_data %>% filter(Group_label == class)
      
      #original data only
      ggplot(data = subset) + 
        geom_histogram(aes(x = prop),
                       color = "white", fill = "white") +
        theme_minimal() + ggtitle(class) + xlab("Proportion") + ylab("Count") +
        geom_vline(xintercept = to_copy_dementia_plot_data[
          which(to_copy_dementia_plot_data$Var1 == class), "prop"], size = 2, 
          color = unique(subset$color))
      
      ggsave(filename = paste0(path_to_folder, "impairment_classes/", class, 
                               "_line_only.jpeg"), 
             height = 3, width = 5, units = "in")
      
      #Prior predictive check
      ggplot(data = subset) + 
        geom_histogram(aes(x = prop), 
                       color = unique(subset$color), 
                       fill = unique(subset$color)) + 
        theme_minimal() + ggtitle(class) + xlab("Proportion") + ylab("Count") +
        geom_vline(xintercept = to_copy_dementia_plot_data[
          which(to_copy_dementia_plot_data$Var1 == class), "prop"], size = 2)
      
      ggsave(filename = paste0(path_to_folder, "impairment_classes/", class, 
                               ".jpeg"), 
             height = 3, width = 5, units = "in")
    }
    
    #---- **class-specific continuous ----
    for(class in unlist(unique(dataset_to_copy[, dementia_var]))){
      true_data <- dataset_to_copy %>% 
        dplyr::select(c(all_of(Z[, "var"]), all_of(dementia_var))) %>% 
        filter(!!as.symbol(dementia_var) == class) %>% mutate("color" = "black")
      
      continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
      
      for(i in 1:length(continuous_list)){
        continuous_list[[i]] <- continuous_list[[i]] %>% 
          mutate("run" = i, "type" = "synthetic") %>% 
          rbind(., true_data %>% dplyr::select(-dementia_var) %>% 
                  mutate("run" = i, "type" = "True Data"))
      }
      
      continuous_list %<>% do.call(rbind, .) %>% as.data.frame() 
      
      for(var in Z[, "var"]){
        data <- continuous_list[, c(var, "run", "type", "color")]
        #unstandardize
        data[, var] <- data[, var]*orig_sds[var] + orig_means[var]
        
        synthetic_subset <- data %>% filter(type == "synthetic")
        # data[which(data$type == "synthetic"), var] <- 
        #   synthetic_subset[, var]*orig_sds[var] + orig_means[var]
        
        continuous_plot <- 
          ggplot(data = data, aes(color = type, fill = type)) + 
          geom_density(aes(x = data[, 1]), alpha = 0.5) + 
          theme_minimal() + xlab(Z[which(Z[, "var"] == var), "label"]) + 
          scale_color_manual(values = rev(unique(data$color))) + 
          scale_fill_manual(values = rev(unique(data$color))) + 
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
id_var = "HHIDPN" 
#dementia_var 
dataset_to_copy = dataset_to_copy 
num_synthetic = 10
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
path_to_folder = paste0(path_to_box, "figures/simulation_study/", 
                        "HCAP_normal_250_unimpaired/prior_predictive_checks/")
