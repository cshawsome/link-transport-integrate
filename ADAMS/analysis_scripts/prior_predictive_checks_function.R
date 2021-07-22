prior_predictive_checks <- 
  function(normal_preds, other_preds, mci_preds, categorical_vars, 
           continuous_vars, id_var, extra_vars = NULL, dataset_to_copy, 
           num_synthetic, orig_means, orig_sds){
    #---- select variables ----
    vars <- unique(c(normal_preds, other_preds, mci_preds, extra_vars))
    
    synthetic_sample <- dataset_to_copy %>% 
      dplyr::select(id_var, all_of(vars)) %>% 
      #pre-allocate columns
      mutate("Group" = 0, "p_normal" = 0, "p_other" = 0, "p_mci" = 0)
    
    
    
    generate_data <- function(){
      #---- latent class ----
      group = 1
      synthetic_sample[, "Group"] <- 0
      
      for(model in c("normal", "other", "mci")){
        subset_index <- which(synthetic_sample$Group == 0)
        random_draw <- sample(seq(1, 10000), size = 1)
        
        prior_betas <- as.vector(get(paste0(model, "_betas"))[, random_draw])
        prior_cov <- matrix(unlist(get(paste0(model, "_cov"))[, random_draw]), 
                            nrow = nrow(prior_betas))
        
        betas <- mvrnorm(n = 1, mu = t(prior_betas), Sigma = prior_cov)
        
        synthetic_sample[subset_index, paste0("p_", model)] <- 
          expit(as.matrix(synthetic_sample[subset_index, 
                                           get(paste0(model, "_preds"))]) %*% 
                  as.matrix(betas))
        
        synthetic_sample[subset_index, "Group"] <- 
          rbernoulli(n = length(subset_index), 
                     p = synthetic_sample[subset_index, paste0("p_", model)])*group
        
        group = group + 1
      }
      synthetic_sample[which(synthetic_sample$Group == 0), "Group"] <- 4
      
      #pre-allocate
      mu <- matrix(0, ncol = 4*6, nrow = 10) %>%
        set_colnames(apply(expand.grid(seq(1, 4), seq(1, 6)), 1, paste0,
                           collapse = ":"))
      
      for(i in 1:4){
        #---- **contingency cells ----
        subset <- synthetic_sample %>% filter(Group == i)
        prior_counts <- alpha_0_dist[, c(sample(seq(1, 10000), size = 1), 
                                         ncol(alpha_0_dist))] %>% 
          filter(group_number == i)
        
        #---- **p(contingency table cell) ----
        pi <- rdirichlet(1, alpha = as.numeric(unlist(prior_counts[, 1])))
        
        #---- **contingency table count ----
        contingency_table <- rmultinom(n = 1, size = nrow(subset), prob = pi)
        UtU <- diag(contingency_table[, 1])
        
        #---- **draw new UtU if needed ----
        while(det(t(A) %*% UtU %*% A) < 1e-9){
          random_draw <- sample(seq(1, 10000), size = 1)
          new_counts <- alpha_0_dist[, c(random_draw, ncol(alpha_0_dist))] %>% 
            filter(group_number == i)
          
          UtU <- diag(unlist(new_counts[, 1]))
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
        random_draw <- sample(seq(1, 10000), size = 1)
        Sigma_prior <- prior_Sigma[, c(random_draw, ncol(prior_Sigma))] %>% 
          filter(group_number == i)
        sig_Y <- riwish(v = (nu_0), S = matrix(unlist(Sigma_prior[, 1]), 
                                               nrow = nrow(Z)))
        
        #---- **beta_0 ----
        V_0_inv <- prior_V_inv[, c(random_draw, ncol(prior_V_inv))] %>% 
          filter(group_number == i)
        beta_0 <- priors_beta[, c(random_draw, ncol(priors_beta))] %>% 
          filter(group_number == i)
        
        #as matrices
        V_0_inv <- matrix(unlist(V_0_inv[, 1]), nrow = 4)
        beta_0 <- matrix(unlist(beta_0[, 1]), nrow = nrow(V_0_inv))
        
        #---- **draw beta | Sigma----
        beta_Sigma_Y <- matrix.normal(beta_0, solve(V_0_inv), sig_Y/kappa_0[i])
        
        #---- **compute mu ----
        mu[, paste0(i, ":", seq(1, 6))] <-
          t(A %*% matrix(beta_Sigma_Y, nrow = ncol(A), ncol = nrow(Z),
                         byrow = FALSE))
        
        #---- **draw data ----
        #reformat contingency table
        table <- contingency_table %>% as.data.frame() %>%
          cbind(do.call(cbind, list(
            #Black              #Hispanic           #Stroke
            rep(c(1, 0, 0), 2), rep(c(0, 1, 0), 2), c(rep(0, 3), rep(1, 3))))) %>%
          set_colnames(c("Count", W))
        
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
                      mu = mu[, paste0(i, ":", j)],
                      Sigma = sig_Y)
          }
        }
        assign(paste0("Z_", i), subset[, all_of(Z[, "var"])])
      }
      
      #---- **return ----
      return(list("Group" = synthetic_sample$Group,
                  "Z_normal" = Z_1 %>% mutate("color" = "#00a389"), 
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
    ADAMS_data[which(ADAMS_data$Adem_dx_cat == "Normal"), "Adem_dx_cat"] <- 
      "Unimpaired"
    ADAMS_dementia_plot_data <- as.data.frame(table(ADAMS_data$Adem_dx_cat)) %>% 
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
      
      #ADAMS data only
      ggplot(data = subset) + 
        geom_histogram(aes(x = prop),
                       color = "white", fill = "white") +
        theme_minimal() + ggtitle(class) + xlab("Proportion") + ylab("Count") +
        geom_vline(xintercept = ADAMS_dementia_plot_data[
          which(ADAMS_dementia_plot_data$Var1 == class), "prop"], size = 2, 
          color = unique(subset$color))
      
      ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                               "priors/impairment_classes/", class, 
                               "_line_only.jpeg"), 
             height = 3, width = 5, units = "in")
      
      #Prior predictive check
      ggplot(data = subset) + 
        geom_histogram(aes(x = prop), 
                       color = unique(subset$color), fill = unique(subset$color)) + 
        theme_minimal() + ggtitle(class) + xlab("Proportion") + ylab("Count") +
        geom_vline(xintercept = ADAMS_dementia_plot_data[
          which(ADAMS_dementia_plot_data$Var1 == class), "prop"], size = 2)
      
      ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                               "priors/impairment_classes/", class, ".jpeg"), 
             height = 3, width = 5, units = "in")
    }
    
    #---- **class-specific continuous ----
    ADAMS_data[which(ADAMS_data$Adem_dx_cat == "Unimpaired"), "Adem_dx_cat"] <- 
      "Normal"
    
    for(class in unique(ADAMS_train$Adem_dx_cat)){
      true_data <- ADAMS_data %>% 
        dplyr::select(c(all_of(Z[, "var"]), "Adem_dx_cat")) %>% 
        filter(Adem_dx_cat == class) %>% mutate("color" = "black")
      
      continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
      
      for(i in 1:length(continuous_list)){
        continuous_list[[i]] <- continuous_list[[i]] %>% 
          mutate("run" = i, "type" = "synthetic") %>% 
          rbind(., true_data %>% dplyr::select(-"Adem_dx_cat") %>% 
                  mutate("run" = i, "type" = "ADAMS"))
      }
      
      continuous_list %<>% do.call(rbind, .) %>% as.data.frame() 
      
      for(var in Z[, "var"]){
        data <- continuous_list[, c(var, "run", "type", "color")]
        #unstandardize
        synthetic_subset <- data %>% filter(type == "synthetic")
        data[which(data$type == "synthetic"), var] <- 
          synthetic_subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
        
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
        
        anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/",
                                    "priors/continuous_vars/", tolower(class), "/", 
                                    var, ".gif"),
                  animation = last_animation(),
                  renderer = gifski_renderer())
      }
    }
  }


