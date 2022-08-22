prior_predictive_checks <- 
  function(dataset_to_copy, calibration_sample = FALSE, calibration_prop = NA, 
           calibration_sample_name = NA, path_to_raw_prior_sample, path_to_data, 
           path_to_output_folder, continuous_check_test = FALSE, 
           continuous_check = c("Unimpaired", "MCI", "Dementia", "Other"), 
           categorical_vars, continuous_vars, variable_labels, color_palette,
           contrasts_matrix, kappa_0_mat, nu_0_mat, num_synthetic){
    
    #---- update path to output folder ----
    if(!calibration_sample){
      path_to_output_folder <- 
        paste0(path_to_output_folder, "no_calibration_sample/")
    } else{
      path_to_output_folder <- 
        paste0(path_to_output_folder, "calibration_", calibration_sample_name, 
               "/")
    }
    
    #---- create folders for results ----
    if(!dir.exists(paste0(path_to_output_folder, "impairment_classes/"))){
      dir.create(paste0(path_to_output_folder, "impairment_classes/"), 
                 recursive = TRUE) 
    }
    
    for(class in c("unimpaired", "mci", "dementia", "other")){
      if(!dir.exists(paste0(path_to_output_folder, "continuous_vars/", 
                            tolower(class)))){
        dir.create(paste0(path_to_output_folder, "continuous_vars/", 
                          tolower(class)), recursive = TRUE)
      }
      
      if(!dir.exists(paste0(path_to_output_folder, "continuous_vars/combined"))){
        dir.create(paste0(path_to_output_folder, "continuous_vars/combined"), 
                   recursive = TRUE)
      }
      
      if(!dir.exists(paste0(path_to_output_folder, "continuous_vars/error_set/", 
                            tolower(class)))){
        dir.create(paste0(path_to_output_folder, "continuous_vars/error_set/", 
                          tolower(class)), recursive = TRUE)
      }
      
      if(!dir.exists(paste0(path_to_output_folder, "continuous_vars/error_set/", 
                            "combined"))){
        dir.create(paste0(path_to_output_folder, "continuous_vars/error_set/", 
                          "combined"), recursive = TRUE)
      }
      
      if(continuous_check_test & 
         !dir.exists(paste0(path_to_output_folder, "continuous_vars/test_set/", 
                            tolower(class)))){
        
        dir.create(paste0(path_to_output_folder, "continuous_vars/test_set/", 
                          tolower(class)), recursive = TRUE)
        
      }
      
      if(continuous_check_test & 
         !dir.exists(paste0(path_to_output_folder, "continuous_vars/test_set/", 
                            "combined"))){
        
        dir.create(paste0(path_to_output_folder, "continuous_vars/test_set/", 
                          "combined"), recursive = TRUE)
        
      }
    }
    
    #---- set prior data ----
    if(!calibration_sample){
      #---- **latent classes ----
      for(group in c("unimpaired", "mci", "other")){
        assign(paste0(group, "_betas"), 
               vroom(paste0(path_to_data, "data/prior_data/latent_class_", 
                            group, "_betas.csv"), delim = ","))
        assign(paste0(group, "_cov"), 
               readRDS(paste0(path_to_data, "data/prior_data/latent_class_", 
                              group, "_cov")))
        
        assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
      }
      
      #---- **contingency cells ----
      alpha_0_dist <- 
        readRDS(paste0(path_to_data, "data/prior_data/imputation_cell_props")) 
      
      #--- **beta and sigma ----
      priors_beta <- 
        readRDS(paste0(path_to_data, "data/prior_data/priors_beta")) 
      prior_V_inv <- 
        readRDS(paste0(path_to_data, "data/prior_data/priors_V_inv"))  
      prior_Sigma <- 
        readRDS(paste0(path_to_data, "data/prior_data/priors_Sigma")) 
    } else{
      #---- read in raw prior sample ----
      prior_imputed_clean <- readRDS(path_to_raw_prior_sample) %>%
        lapply(function(x) mutate_at(x, "HHIDPN", as.numeric)) 
      
      #---- selected vars ----
      selected_vars <- 
        read_csv(paste0(path_to_data, "data/variable_selection/", 
                        "model_coefficients.csv")) 
      
      #unimpaired model predictors
      unimpaired_preds <- 
        c("(Intercept)", selected_vars %>% 
            filter(data_label != "Intercept" & Unimpaired != 0) %>% 
            dplyr::select(data_label) %>% unlist())
      
      #other model predictors
      other_preds <- c("(Intercept)", selected_vars %>% 
                         filter(data_label != "Intercept" & Other != 0) %>% 
                         dplyr::select(data_label) %>% unlist())
      
      #mci model predictors
      mci_preds <- c("(Intercept)", selected_vars %>% 
                       filter(data_label != "Intercept" & MCI != 0) %>% 
                       dplyr::select(data_label) %>% unlist())
      
      #---- cell ID key ----
      cell_ID_key <- read_csv(paste0(path_to_data, "data/cell_ID_key.csv"))
      
      #---- calibration subset ----
      calibration_var <- paste0("calibration_", calibration_prop*100)
      calibration_subset <- 
        dataset_to_copy %>% filter(!!sym(calibration_var) == 1)
    }
    
    #---- select variables ----
    vars <- unique(c(unimpaired_preds, other_preds, mci_preds))
    
    synthetic_sample <- dataset_to_copy %>% 
      dplyr::select("HHIDPN", all_of(vars)) %>% 
      #pre-allocate columns
      mutate("group_num" = 0, "p_unimpaired" = 0, "p_other" = 0, "p_mci" = 0)
    
    if(calibration_sample){
      #---- **split sample ----
      synthetic_sample <- 
        anti_join(synthetic_sample, calibration_subset, by = "HHIDPN")
    }
    
    #---- max index ----
    if(calibration_sample){
      max_index <- length(prior_imputed_clean)
    } else{
      max_index <- length(priors_beta)  
    }
    
    generate_data <- function(color_palette){
      #---- index for random draws ----
      random_draw <- sample(seq(1, max_index), size = 1)
      
      #---- latent class ----
      group_num = 1
      for(model in c("unimpaired", "other", "mci")){
        subset_index <- which(synthetic_sample$group_num == 0)
        
        if(calibration_sample){
          if(model == "mci"){
            class_name = "MCI"
          } else{
            class_name = str_to_sentence(model)
          }
          
          latent_class_model <- 
            glm(formula(paste(class_name, " ~ ", 
                              paste(get(paste0(model, "_preds"))[-1], 
                                    collapse = " + "), 
                              collapse = "")), family = "binomial", 
                #don't select (Intercept) variable
                data = rbind(prior_imputed_clean[[random_draw]][, vars[-1]], 
                             calibration_subset[, vars[-1]]))
          
          prior_betas <- coefficients(latent_class_model)
          prior_cov <- vcov(latent_class_model)
          
        } else{
          prior_betas <- get(paste0(model, "_betas"))[, random_draw]
          prior_cov <- get(paste0(model, "_cov"))[[random_draw]]
        }
        
        betas <- 
          mvnfast::rmvn(n = 1, mu = unlist(prior_betas), sigma = prior_cov)
        
        synthetic_sample[subset_index, paste0("p_", model)] <- 
          expit(as.matrix(synthetic_sample[subset_index, 
                                           get(paste0(model, "_preds"))]) %*% 
                  t(as.matrix(betas)))
        
        synthetic_sample[subset_index, "group_num"] <- 
          rbernoulli(n = length(subset_index), 
                     p = synthetic_sample[subset_index, 
                                          paste0("p_", model)])*group_num
        
        group_num = group_num + 1
      }
      
      synthetic_sample[which(synthetic_sample$group_num == 0), "group_num"] <- 4
      
      synthetic_sample[, "Group"] <- 
        case_when(synthetic_sample$group_num == 1 ~ "Unimpaired", 
                  synthetic_sample$group_num == 2 ~ "Other", 
                  synthetic_sample$group_num == 3 ~ "MCI", 
                  synthetic_sample$group_num == 4 ~ "Dementia")
      
      #pre-allocate: ncol = num impairement groups * num contingency cells
      # these are contingency-cell specific means for continuous variables by
      # impairment class
      mu <- matrix(0, ncol = 4*nrow(contrasts_matrix), 
                   nrow = length(continuous_vars)) %>%
        set_colnames(apply(
          expand.grid(c("Unimpaired", "MCI", "Dementia", "Other"), 
                      seq(1, nrow(contrasts_matrix))), 1, paste0,
          collapse = ":"))
      
      #---- hyperparameters ----
      kappa_0 <- kappa_0_mat[1, ]
      nu_0 <- nu_0_mat[, ]  
      
      for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
        
        #---- **contingency cells ----
        subset <- synthetic_sample %>% filter(Group == class)
        if(nrow(subset) == 0){
          assign(paste0("Z_", class), 
                 matrix(rep(NA, length(continuous_vars)), nrow = 1) %>% 
                   as.data.frame() %>% set_colnames(continuous_vars))
          next
        }
        
        if(calibration_sample){
          prior_counts <- 
            rbind(prior_imputed_clean[[random_draw]][, c(categorical_vars, class)], 
                  calibration_subset[, c(categorical_vars, class)]) %>% 
            filter(!!sym(class) == 1) %>% 
            unite("cell_ID", all_of(categorical_vars), sep = "") %>% 
            dplyr::select("cell_ID") %>% table() %>% as.data.frame() %>% 
            set_colnames(c("cell_ID", "Freq"))
          
          if(nrow(prior_counts) < nrow(cell_ID_key)){
            prior_counts <- 
              left_join(cell_ID_key, prior_counts, by = "cell_ID") %>% 
              dplyr::select(c("cell_ID", "Freq")) 
            
            prior_counts[which(is.na(prior_counts$Freq)), "Freq"] <- 0
          }
          
          prior_counts %<>% mutate("prop" = Freq/sum(Freq))
          prior_UtU <- diag(prior_counts$Freq)
          
          #update counts for this particular subset
          prior_counts <- prior_counts$prop*nrow(subset)
          
          
        } else{
          prior_counts <- 
            alpha_0_dist[[random_draw]][[class]][, "props"]*nrow(subset)
        }
        
        #---- **p(contingency table cell) ----
        pi <- 
          MCMCpack::rdirichlet(1, alpha = as.numeric(unlist(prior_counts)))
        
        #---- **contingency table count ----
        contingency_table <- rmultinom(n = 1, size = nrow(subset), prob = pi)
        UtU <- diag(contingency_table[, 1])
        
        # #---- **draw new UtU if needed ----
        # while(det(t(contrasts_matrix) %*% UtU %*% contrasts_matrix) < 1e-9){
        #   #---- ****new random draw index ----
        #   random_draw <- sample(seq(1, max_index), size = 1)
        #   
        #   if(calibration_sample){
        #     new_counts <- 
        #       rbind(prior_imputed_clean[[random_draw]][, c(categorical_vars, class)], 
        #             calibration_subset[, c(categorical_vars, class)]) %>% 
        #       filter(!!as.symbol(class) == 1) %>% 
        #       unite("cell_ID", all_of(categorical_vars), sep = "") %>% 
        #       dplyr::select("cell_ID") %>% table() %>% as.data.frame()
        #     
        #     if(nrow(new_counts) < nrow(cell_ID_key)){
        #       new_counts <- left_join(cell_ID_key, new_counts) %>% 
        #         dplyr::select(c("cell_ID", "Freq")) 
        #       
        #       new_counts[which(is.na(new_counts$Freq)), "Freq"] <- 0
        #     }
        #     
        #     new_counts %<>% mutate("prop" = Freq/sum(Freq))
        #     
        #     new_counts <- new_counts$prop*nrow(subset)
        #     
        #     UtU <- diag(new_counts)
        #     
        #   } else{
        #     new_counts <- 
        #       alpha_0_dist[[random_draw]][[class]][, "props"]*nrow(subset)
        #     
        #     UtU <- diag(unlist(new_counts[, 1]))
        #   }
        # }
        
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
        
        #---- **beta_0 + Sigma_0 ----
        if(calibration_sample){
          
          #---- ****prior U matrix ----
          prior_U <- matrix(0, nrow = sum(prior_UtU), ncol = nrow(prior_UtU))
          
          for(j in 1:ncol(prior_UtU)){
            if(sum(prior_UtU[, j]) == 0){next}
            if(j == 1){
              index = 1
            } else{
              index = sum(prior_UtU[, 1:(j - 1)]) + 1
            }
            prior_U[index:(index - 1 + sum(prior_UtU[, j])), j] <- 1
          }
          
          continuous_covariates <- 
            rbind(prior_imputed_clean[[random_draw]][, c(continuous_vars, class)], 
                  calibration_subset[, c(continuous_vars, class)]) %>% 
            filter(!!sym(class) == 1) %>% 
            dplyr::select(all_of(continuous_vars)) %>% as.matrix()
          
          V_0_inv <- t(A) %*% prior_UtU %*% A
          V_0 <- solve(V_0_inv)
          
          beta_0 <- V_0 %*% t(A) %*% t(prior_U) %*% continuous_covariates
          
          residual <- continuous_covariates - 
            prior_U %*% contrasts_matrix %*% beta_0
          
          Sigma_prior <- t(residual) %*% residual
          
        } else{
          V_0_inv <- 
            as.matrix(
              prior_V_inv[[random_draw]][[class]][, seq(1, ncol(A))])
          
          beta_0 <- 
            as.matrix(
              priors_beta[[random_draw]][[
                class]][, seq(1, length(continuous_vars))])
          
          Sigma_prior <-
            as.matrix(
              prior_Sigma[[random_draw]][[
                class]][, seq(1, length(continuous_vars))])
        }
        
        #Sigma_prior <- diag(1, nrow = length(continuous_vars))
        
        redraws = 0
        while(is.character(tryCatch(sig_Y <- 
                                    MCMCpack::riwish(
                                      v = as.numeric(nu_0[, class]), 
                                      S = Sigma_prior), 
                                    error = function(e) "error")) & 
              redraws <= 100){
          
          Sigma_prior = Sigma_prior + 0.001
          redraws = redraws + 1
          #print(paste0("class = ", class, "; b = ", b))            
        }
        
        #---- **draw beta | Sigma----
        beta_Sigma_Y <- 
          rmatrixnorm(M = beta_0, U = as.positive.definite(solve(V_0_inv)), 
                      V = as.positive.definite(sig_Y/unlist(kappa_0[, class])))
        
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
                    mutate("Group" = "Unimpaired") %>% 
                    left_join(color_palette, by = "Group"), 
                  "Z_other" = Z_Other %>% 
                    mutate("Group" = "Other") %>% 
                    left_join(color_palette, by = "Group"), 
                  "Z_mci" = Z_MCI %>% 
                    mutate("Group" = "MCI") %>% 
                    left_join(color_palette, by = "Group"), 
                  "Z_dementia" = Z_Dementia %>% 
                    mutate("Group" = "Dementia") %>% 
                    left_join(color_palette, by = "Group")))
    }
    
    #---- multiruns ----
    runs = num_synthetic
    synthetic <- replicate(num_synthetic, generate_data(color_palette), 
                           simplify = FALSE) 
    
    #---- plots ----
    #---- **dem class ----
    dem_sub <- lapply(synthetic, "[[", "Group") %>% do.call(cbind, .) %>% 
      set_colnames(seq(1, runs)) %>% as.data.frame() %>%
      pivot_longer(everything()) 
    
    synthetic_dementia_plot_data <- 
      dem_sub %>% dplyr::count(name, value) %>%
      group_by(name) %>%
      mutate_at("name", as.numeric) %>% 
      left_join(., color_palette, by = c("value" = "Group"))
    
    #---- ****individual ----
    for(class in unique(synthetic_dementia_plot_data$value)){
      subset <- synthetic_dementia_plot_data %>% filter(value == class)
      
      if(nrow(subset) == 0){
        next 
      } else{
        #Prior predictive check
        ggplot(data = subset) + 
          geom_histogram(aes(x = n), 
                         color = unique(subset$Color), 
                         fill = unique(subset$Color)) + 
          theme_minimal() + ggtitle(class) + xlab("Count") + ylab("Frequency") 
        
        ggsave(filename = paste0(path_to_output_folder, "impairment_classes/", 
                                 class, ".jpeg"), 
               height = 3, width = 5, units = "in")
      }
    }
    
    #---- ****combined ----
    synthetic_dementia_plot_data$value <- 
      factor(synthetic_dementia_plot_data$value, 
             levels = c("Unimpaired", "MCI", "Dementia", "Other"))
    
    colors <- synthetic_dementia_plot_data %>% arrange(value) %>% 
      group_by(value) %>% slice_head(n = 1) %>% ungroup() %>% 
      dplyr::select(Color) %>% unlist() %>% unname()
    
    ggplot(data = synthetic_dementia_plot_data) + 
      geom_histogram(aes(x = n, color = value, fill = value)) + 
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_minimal() + theme(legend.title=element_blank()) + 
      ggtitle("Impairment Classes") + xlab("Count") + ylab("Frequency") 
    
    ggsave(filename = paste0(path_to_output_folder, 
                             "impairment_classes/all_classes.jpeg"), 
           height = 3, width = 5, units = "in")
    
    #---- **continuous ----
    #---- ****class-specific ----
    for(class in continuous_check){
      
      continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
      
      for(i in 1:length(continuous_list)){
        continuous_list[[i]] <- continuous_list[[i]] %>% 
          dplyr::select(-c("Group")) %>%
          mutate("run" = i, "Type" = "Synthetic") 
      }
      
      continuous_list %<>% do.call(rbind, .) %>% as.data.frame() 
      
      for(var in continuous_vars){
        data <- continuous_list[, c(var, "run", "Type", "Color")] 
        
        if(is.na(sum(data[, var])) | continuous_check_test){
          continuous_plot <- 
            ggplot(data = data, aes(color = Type, fill = Type)) + 
            geom_density(aes(x = data[, 1]), alpha = 0.5) + 
            theme_minimal() + 
            xlab(variable_labels[variable_labels$data_label == var, 
                                 "figure_label"]) + 
            scale_color_manual(values = rev(unique(data$Color))) + 
            scale_fill_manual(values = rev(unique(data$Color)))
          
          if(continuous_check_test){
            ggsave(filename = paste0(path_to_output_folder, 
                                     "continuous_vars/test_set/", 
                                     tolower(class), "/", var, ".jpeg"), 
                   height = 3, width = 5, units = "in")
          } else{
            ggsave(filename = paste0(path_to_output_folder, 
                                     "continuous_vars/error_set/", 
                                     tolower(class), "/", var, ".jpeg"), 
                   height = 3, width = 5, units = "in")
          }
          
        } else{
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
          
          anim_save(filename = paste0(path_to_output_folder, "continuous_vars/", 
                                      tolower(class), "/", var, ".gif"),
                    animation = last_animation(),
                    renderer = gifski_renderer())
        }
      }
    }
    
    #---- ****combined ----
    true_data <- dataset_to_copy %>%
      dplyr::select(c(all_of(continuous_vars))) %>% mutate("Color" = "black")
    
    for(class in unique(synthetic_dementia_plot_data$value)){
      
      continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
      
      for(i in 1:length(continuous_list)){
        continuous_list[[i]] <- continuous_list[[i]] %>% 
          dplyr::select(-c("Group")) %>%
          mutate("run" = i, "Type" = "Synthetic", "Color" = "#f2caaa") %>% 
          rbind(., true_data %>% mutate("run" = i, "Type" = "Observed"))
      }
      
      if(exists("all_continuous_list")){
        all_continuous_list %<>% append(., continuous_list)
      } else{
        all_continuous_list <- continuous_list
      }
    }
    
    all_continuous_list %<>% do.call(rbind, .) %>% as.data.frame() 
    
    for(var in continuous_vars){
      data <- all_continuous_list[, c(var, "run", "Type", "Color")] 
      
      if(is.na(sum(data[, var])) | continuous_check_test){
        all_continuous_plot <- 
          ggplot(data = data, aes(color = Type, fill = Type)) + 
          geom_density(aes(x = data[, 1]), alpha = 0.5) + 
          theme_minimal() + 
          xlab(variable_labels[variable_labels$data_label == var, 
                               "figure_label"]) + 
          scale_color_manual(values = rev(unique(data$Color))) + 
          scale_fill_manual(values = rev(unique(data$Color)))
        
        if(continuous_check_test){
          ggsave(filename = paste0(path_to_output_folder, 
                                   "continuous_vars/test_set/combined/", 
                                   var, ".jpeg"), 
                 height = 3, width = 5, units = "in")
        } else{
          ggsave(filename = paste0(path_to_output_folder, 
                                   "continuous_vars/error_set/combined/", 
                                   var, ".jpeg"), 
                 height = 3, width = 5, units = "in")
        }
        
      } else{
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
        
        anim_save(filename = paste0(path_to_output_folder, "continuous_vars/", 
                                    "combined/", var, ".gif"),
                  animation = last_animation(),
                  renderer = gifski_renderer())
      }
    }
  }

# #---- test function ----
# unimpaired_preds = unimpaired_preds
# other_preds = other_preds
# mci_preds = mci_preds
# categorical_vars = W
# continuous_vars = Z
# variable_labels = variable_labels
# color_palette = color_palette
# dataset_to_copy = synthetic_HCAP_list[[7]]
# num_synthetic = 1000
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
# nu_0 = nu_0[, unique(synthetic_data_list[[1]]$dataset_name)]
# kappa_0_mat = kappa_0_mat
# contrasts_matrix = A
# path_to_output_folder = paste0(path_to_box, "figures/simulation_study/HCAP_HRS_",
#                         unique(synthetic_data_list[[1]][, "dataset_name"]),
#                         "/prior_predictive_checks/")
