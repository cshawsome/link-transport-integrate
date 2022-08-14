generate_synthetic <- 
  function(warm_up, run_number, starting_props, dataset_to_copy, 
           calibration_sample = FALSE, calibration_prop = NA, 
           calibration_sample_name = NA, path_to_raw_prior_sample, path_to_data, 
           path_to_analyses_folder, path_to_figures_folder, categorical_vars, 
           continuous_vars, id_var, variable_labels, cell_ID_key, color_palette, 
           contrasts_matrix, kappa_0_mat, nu_0_mat, num_synthetic, 
           data_only = FALSE){
    
    #---- check subfolders for results ----
    if(!calibration_sample){
      if(!data_only){
        if(!dir.exists(paste0(path_to_analyses_folder, "synthetic_data/", 
                              "no_calibration_sample/run_", run_number))){
          dir.create(paste0(path_to_analyses_folder, "synthetic_data/", 
                            "no_calibration_sample/run_", run_number), 
                     recursive = TRUE)
        }
        
        if(!dir.exists(paste0(path_to_figures_folder, "diagnostics/", 
                              "no_calibration_sample/run_", run_number))){
          dir.create(paste0(path_to_figures_folder, "diagnostics/", 
                            "no_calibration_sample/run_", run_number), 
                     recursive = TRUE)
        }
        
        if(!dir.exists(paste0(path_to_analyses_folder, "diagnostics_data/", 
                              "no_calibration_sample/run_", run_number))){
          dir.create(paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "no_calibration_sample/run_", run_number), 
                     recursive = TRUE)
        }
      }
    } else{
      if(!data_only){
        if(!dir.exists(paste0(path_to_analyses_folder, "synthetic_data/", 
                              "calibration_", calibration_sample_name, "/run_", 
                              run_number))){
          dir.create(paste0(path_to_analyses_folder, "synthetic_data/", 
                            "calibration_", calibration_sample_name, "/run_", 
                            run_number), recursive = TRUE)
        }
        
        if(!dir.exists(paste0(path_to_figures_folder, "diagnostics/", 
                              "calibration_", calibration_sample_name, "/run_", 
                              run_number))){
          dir.create(paste0(path_to_figures_folder, "diagnostics/", 
                            "calibration_", calibration_sample_name, "/run_", 
                            run_number), recursive = TRUE)
        }
        
        if(!dir.exists(paste0(path_to_analyses_folder, "diagnostics_data/", 
                              "calibration_", calibration_sample_name, "/run_", 
                              run_number))){
          dir.create(paste0(path_to_analyses_folder, "diagnostics_data/", 
                            "calibration_", calibration_sample_name, "/run_", 
                            run_number), recursive = TRUE)
        }
      } 
    }
    
    #---- set prior distributions ----
    if(!calibration_sample){
      #---- **latent classes ----
      for(group in c("unimpaired", "mci", "other")){
        assign(paste0(group, "_betas"), 
               vroom(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                            "latent_class_", group, "_betas.csv"), delim = ","))
        assign(paste0(group, "_cov"), 
               readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                              "latent_class_", group, "_cov")))
        
        assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
      }
      
      #---- **contingency cells ----
      alpha_0_dist <- 
        readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                       "imputation_cell_props")) 
      
      #--- **beta and sigma ----
      priors_beta <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                                    "prior_data/priors_beta")) 
      prior_V_inv <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                                    "prior_data/priors_V_inv"))  
      prior_Sigma <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                                    "prior_data/priors_Sigma")) 
    } else{
      #---- read in raw prior sample ----
      prior_imputed_clean <- readRDS(path_to_raw_prior_sample) %>%
        lapply(function(x) mutate_at(x, "HHIDPN", as.numeric)) 
      
      #---- **prep sample ----
      variable_labels_ADAMS <- variable_labels %>% 
        filter(ADAMS %in% colnames(prior_imputed_clean[[1]]))
      
      prior_imputed_clean <- 
        lapply(prior_imputed_clean, 
               function(x) rename_at(x, vars(variable_labels_ADAMS$ADAMS), ~ 
                                       variable_labels_ADAMS$data_label)) 
      
      #---- selected vars ----
      selected_vars <- 
        read_csv(paste0(path_to_data, 
                        "analyses/simulation_study/variable_selection/", 
                        "model_coefficients.csv")) 
      
      #unimpaired model predictors
      unimpaired_preds <- c("(Intercept)", selected_vars %>% 
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
      
      #---- calibration subset ----
      calibration_var <- paste0("calibration_", calibration_prop*100)
      calibration_subset <- 
        dataset_to_copy %>% filter(!!sym(calibration_var) == 1)
    }
    
    #---- select variables ----
    vars <- unique(c(unimpaired_preds, other_preds, mci_preds, 
                     "Unimpaired", "MCI", "Dementia", "Other"))
    
    #---- sampling counts ----
    warm_up = warm_up
    synthetic_sets = num_synthetic
    B = warm_up + synthetic_sets
    
    if(calibration_sample){
      #---- **split sample ----
      dataset_to_copy <- 
        anti_join(dataset_to_copy, calibration_subset, by = "HHIDPN")
    }
    
    #---- count contingency cells ----
    cross_class_label <- dataset_to_copy[, categorical_vars] %>%
      unite("cell_ID", everything(), sep = "") %>% table() %>% 
      as.data.frame() %>% set_colnames(c("cell_ID", "count"))
    
    if(nrow(cross_class_label) < nrow(cell_ID_key)){
      missing <- which(!cell_ID_key$cell_ID %in% cross_class_label$cell_ID)
      new_table <- cell_ID_key[, "cell_ID"] %>% left_join(., cross_class_label)
      new_table[is.na(new_table$count), "count"] <- 0
      cross_class_label <- new_table
    }
    
    #---- chain storage ----
    model_gamma_chain <- 
      matrix(nrow = sum(length(unimpaired_preds), length(other_preds), 
                        length(mci_preds)), ncol = B) %>% as.data.frame() %>%
      mutate("model" = c(rep("unimpaired", length(unimpaired_preds)), 
                         rep("other", length(other_preds)), 
                         rep("mci", length(mci_preds))), 
             "pred" = c(unimpaired_preds, mci_preds, other_preds))
    
    pi_chain <- matrix(nrow = nrow(cross_class_label), ncol = 4*B) %>% 
      set_colnames(gsub(" ", "", 
                        apply(expand.grid(
                          c("Unimpaired", "MCI", "Dementia", "Other"), 
                          seq(1, B)), 1, paste, collapse = ":"))) %>% 
      set_rownames(cross_class_label$cell_ID)
    
    mu_chain <-
      matrix(nrow = length(continuous_vars), 
             ncol = 4*nrow(cross_class_label)*B) %>%
      set_colnames(gsub(" ", "", 
                        apply(expand.grid(
                          c("Unimpaired", "MCI", "Dementia", "Other"),
                          seq(1:nrow(cross_class_label)), seq(1, B)), 1, paste, 
                          collapse = ":"))) %>% 
      set_rownames(unlist(variable_labels[variable_labels$data_label %in% 
                                            continuous_vars, "figure_label"]))
    
    if(!data_only){
      latent_class_chain <- matrix(nrow = 4, ncol = B) %>% 
        set_rownames(c("Unimpaired", "Other", "MCI", "Dementia"))
      
      Sigma_chain <- matrix(nrow = length(continuous_vars), ncol = 4*B) %>%
        set_colnames(gsub(" ", "", 
                          apply(expand.grid(
                            c("Unimpaired", "MCI", "Dementia", "Other"), 
                            seq(1, B)), 1, paste, collapse = ":"))) %>% 
        set_rownames(unlist(variable_labels[variable_labels$data_label %in% 
                                              continuous_vars, "figure_label"]))
    }
    
    #---- max index ----
    if(calibration_sample){
      max_index <- length(prior_imputed_clean)
    } else{
      max_index <- length(priors_beta)  
    }
    
    #---- nu_0 and kappa_0 hyperparameters ----
    if(!calibration_sample){
      kappa_0 <- 
        kappa_0_mat[
          which(kappa_0_mat$dataset_name == 
                  unlist(unique(dataset_to_copy[, "dataset_name"])) & 
                  kappa_0_mat$calibration_name == "no_calibration"), ]
      
      nu_0 <- 
        nu_0_mat[
          which(nu_0_mat$dataset_name == 
                  unlist(unique(dataset_to_copy[, "dataset_name"])) & 
                  nu_0_mat$calibration_name == "no_calibration"), ]  
    } else{
      kappa_0 <- 
        kappa_0_mat[
          which(kappa_0_mat$dataset_name == 
                  unlist(unique(dataset_to_copy[, "dataset_name"])) & 
                  kappa_0_mat$calibration_name == 
                  paste0("calibration_", 100*calibration_prop)), ]
      
      nu_0 <- 
        nu_0_mat[
          which(nu_0_mat$dataset_name == 
                  unlist(unique(dataset_to_copy[, "dataset_name"])) & 
                  nu_0_mat$calibration_name == 
                  paste0("calibration_", 100*calibration_prop)), ]  
    }
    
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
          random_draw <- sample(seq(1, max_index), size = 1)
          
          if(!calibration_sample){
            prior_betas <- get(paste0(model, "_betas"))[, random_draw]
            prior_cov <- get(paste0(model, "_cov"))[[random_draw]]
          } else{
            
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
          }
          
          model_gamma_chain[which(model_gamma_chain$model == model), b] <- 
            t(mvnfast::rmvn(n = 1, mu = unlist(prior_betas), sigma = prior_cov))
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
      
      dataset_to_copy[which(dataset_to_copy$group_num == 0), "group_num"] <- 4
      
      dataset_to_copy[, "Group"] <- 
        case_when(dataset_to_copy$group_num == 1 ~ "Unimpaired", 
                  dataset_to_copy$group_num == 2 ~ "Other", 
                  dataset_to_copy$group_num == 3 ~ "MCI", 
                  dataset_to_copy$group_num == 4 ~ "Dementia")
      
      #---- ****group: summary ----
      summary <- table(dataset_to_copy$group_num)/nrow(dataset_to_copy) 
      
      if(length(summary) < 4){
        missing <- which(!seq(1, 4) %in% names(summary))
        new_summary <- vector(length = 4)
        new_summary[missing] <- 0
        new_summary[-missing] <- summary
        summary <- new_summary
      } 
      
      if(!data_only){
        latent_class_chain[, b] <- summary 
      }
      
      for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
        subset <- dataset_to_copy %>% filter(Group == class) 
        if(nrow(subset) == 0){
          next
        } else{
          
          observed_count <- subset[, categorical_vars] %>% 
            unite("cell_ID", sep = "") %>% table() 
          
          if(length(observed_count) < nrow(cell_ID_key)){
            missing <- 
              which(!cell_ID_key$cell_ID %in% names(observed_count))
            new_counts <- vector(length = nrow(cell_ID_key))
            new_counts[missing] <- 0
            new_counts[-missing] <- observed_count
            observed_count <- new_counts
          }
          
          random_draw <- sample(seq(1, max_index), size = 1)
          
          if(!calibration_sample){
            prior_count <- 
              alpha_0_dist[[random_draw]][[class]][, "props"]*nrow(subset)  
          } else{
            prior_count <- 
              rbind(prior_imputed_clean[[random_draw]][, c(categorical_vars, class)], 
                    calibration_subset[, c(categorical_vars, class)]) %>% 
              filter(!!sym(class) == 1) %>% 
              unite("cell_ID", all_of(categorical_vars), sep = "") %>% 
              dplyr::select("cell_ID") %>% table() %>% as.data.frame() %>% 
              set_colnames(c("cell_ID", "Freq"))
            
            if(nrow(prior_count) < nrow(cell_ID_key)){
              prior_count <- left_join(cell_ID_key, prior_count) %>% 
                dplyr::select(c("cell_ID", "Freq")) 
              
              prior_count[which(is.na(prior_count$Freq)), "Freq"] <- 0
            }
            
            prior_count %<>% mutate("prop" = Freq/sum(Freq))
            
            prior_UtU <- diag(prior_count$Freq)
            
            prior_count <- prior_count$prop*nrow(subset)
          }
          
          posterior_count <- prior_count + observed_count
          
          #---- **p(contingency table cell) ----
          pi_chain[, paste0(class, ":", b)] <- 
            MCMCpack::rdirichlet(1, alpha = as.numeric(unlist(posterior_count)))
          
          #---- ****contingency table count ----
          contingency_table <- 
            rmultinom(n = 1, size = nrow(subset), 
                      prob = pi_chain[, paste0(class, ":", b)])
          
          UtU <- diag(contingency_table[, 1])
          
          #---- ****draw new UtU if needed ----
          # while(det(t(A) %*% UtU %*% A) < 1e-10){
          #   
          #   random_draw <- sample(seq(1, max_index), size = 1)
          #   
          #   posterior_first_count <- 
          #     alpha_0_dist[[random_draw]][[class]][, "props"]*nrow(subset)
          #   
          #   posterior_count <- 
          #     posterior_first_count + posterior_second_count
          #   
          #   pi_chain[, paste0(class, ":", b)] <- 
          #     MCMCpack::rdirichlet(1, 
          #                          alpha = as.numeric(unlist(posterior_count)))
          #   
          #   contingency_table <- 
          #     rmultinom(n = 1, size = nrow(subset), 
          #               prob = pi_chain[, paste0(class, ":", b)])
          #   
          #   UtU <- diag(contingency_table[, 1])
          # }
          
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
          continuous_covariates <- subset[, continuous_vars] %>% as.matrix
          
          V_inv <- t(A) %*% UtU %*% A 
          random_draw <- sample(seq(1, max_index), size = 1)
          
          if(!calibration_sample){
            V_0_inv <- 
              as.matrix(
                prior_V_inv[[random_draw]][[class]][, seq(1, ncol(V_inv))])
            
            beta_0 <- 
              as.matrix(
                priors_beta[[random_draw]][[
                  class]][, seq(1, length(continuous_vars))])
          } else{
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
            
            prior_continuous_covariates <- 
              rbind(prior_imputed_clean[[random_draw]][, c(continuous_vars, class)], 
                    calibration_subset[, c(continuous_vars, class)]) %>% 
              filter(!!sym(class) == 1) %>% 
              dplyr::select(all_of(continuous_vars)) %>% as.matrix()
            
            V_0_inv <- t(A) %*% prior_UtU %*% A
            V_0 <- solve(V_0_inv)
            
            beta_0 <- V_0 %*% t(A) %*% t(prior_U) %*% prior_continuous_covariates
          }
          
          M <- solve(V_inv + unlist(kappa_0[, class])*V_0_inv)
          m <-  t(A) %*% t(U) %*% continuous_covariates + 
            unlist(kappa_0[, class])*V_0_inv %*% beta_0
          
          Mm <- M %*% m
          
          #---- ****draw Sigma | Y ----
          ZtZ <- t(continuous_covariates) %*% continuous_covariates
          third_term <- 
            unlist(kappa_0[, class])*t(beta_0) %*% V_0_inv %*% beta_0
          fourth_term <- t(m) %*% M %*% m
          
          # random_draw <- sample(seq(1, max_index), size = 1)
          # Sigma_prior <- 
          #   as.matrix(
          #     prior_Sigma[[random_draw]][[
          #       class]][, seq(1, length(continuous_vars))])
          
          Sigma_prior <- diag(1, nrow = ncol(continuous_covariates))
          
          sigma_mat <- Sigma_prior + ZtZ + third_term - fourth_term
          redraws = 0
          
          while(is.character(tryCatch(sig_Y <- 
                                      MCMCpack::riwish(
                                        v = (as.numeric(nu_0[, class]) + 
                                             nrow(subset)), 
                                        S = sigma_mat), 
                                      error = function(e) "error")) & 
                redraws <= 100){
            
            sigma_mat = sigma_mat + 0.001
            redraws = redraws + 1
            #print(paste0("class = ", class, "; b = ", b))            
          }
          
          if(!data_only){
            Sigma_chain[, paste0(class, ":", b)] <- diag(sig_Y)
          }
          
          #---- ****draw beta | Sigma, Y ----
          beta_Sigma_Y <- 
            rmatrixnorm(M = Mm, U = as.positive.definite(M), 
                        V = as.positive.definite(sig_Y/unlist(kappa_0[, class])))
          
          #---- ****compute mu ----
          mu_chain[, paste0(class, ":", 
                            seq(1, nrow(cross_class_label)), ":", b)] <- 
            t(A %*% beta_Sigma_Y)
          
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
                     continuous_vars] <- 
                mvnfast::rmvn(n = contingency_table[j, "Count"],
                              mu = mu_chain[, paste0(class, ":", j, ":", b)], 
                              sigma = as.positive.definite(sig_Y))
            } else{
              subset[index:(index - 1 + contingency_table[j, "Count"]), 
                     continuous_vars] <- 
                mvnfast::rmvn(n = contingency_table[j, "Count"],
                              mu = mu_chain[, paste0(class, ":", j, ":", b)], 
                              sigma = as.positive.definite(sig_Y))
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
                          c(categorical_vars, continuous_vars)] <- 
            subset[, c(categorical_vars, continuous_vars)]
        }
      }
      
      #---- ****correct "White" column ----
      dataset_to_copy %<>% 
        mutate("White" = if_else(black == 0 & hispanic == 0, 1, 0))
      
      #---- ****save synthetic sample ----
      if(b > warm_up){
        if(!exists("dataset_list")){
          dataset_list <- list()
        }
        dataset_list[[b - warm_up]] <- dataset_to_copy
      }
    }
    
    if(data_only){
      return(dataset_list)
    } else{
      if(!calibration_sample){
        saveRDS(dataset_list,
                file = paste0(path_to_analyses_folder, "synthetic_data/", 
                              "no_calibration_sample/run_", run_number, 
                              "/synthetic_dataset_list"))
      } else{
        saveRDS(dataset_list,
                file = paste0(path_to_analyses_folder, "synthetic_data/", 
                              "calibration_", calibration_sample_name, "/run_", 
                              run_number, "/synthetic_dataset_list"))
      }
      
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
                   colour = figure_label)) + 
        geom_line(aes(group = figure_label)) + 
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
      
      if(!calibration_sample){
        ggsave(filename = "gamma_chain.jpeg", plot = gamma_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "no_calibration_sample/run_", run_number), 
               height = 7, width = 12, units = "in", device = "jpeg")
      } else{
        ggsave(filename = "gamma_chain.jpeg", plot = gamma_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "calibration_", calibration_sample_name, "/run_", 
                             run_number), height = 7, width = 12, units = "in", 
               device = "jpeg")
      }
      
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
      
      if(!calibration_sample){
        ggsave(filename = "latent_class_chain.jpeg", 
               plot = latent_class_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "no_calibration_sample/run_", run_number), 
               height = 7, width = 10, units = "in", device = "jpeg")
      } else{
        ggsave(filename = "latent_class_chain.jpeg", 
               plot = latent_class_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "calibration_", calibration_sample_name, "/run_", 
                             run_number), height = 7, width = 10, units = "in", 
               device = "jpeg")
      }
      
      #---- ****pi chain ----
      pi_chain_data <- pi_chain %>% as.data.frame() %>% 
        rownames_to_column("Cell") %>% 
        pivot_longer(-c("Cell"), names_to = c("Group", "Run"), names_sep = ":", 
                     values_to = "probability") %>% 
        arrange(desc(probability)) %>% 
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
      
      if(!calibration_sample){
        ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "no_calibration_sample/run_", run_number), 
               height = 7, width = 10, units = "in", device = "jpeg")
      } else{
        ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "calibration_", calibration_sample_name, "/run_", 
                             run_number), height = 7, width = 10, units = "in", 
               device = "jpeg")
      }
      
      #---- ****Sigma chain ----
      Sigma_chain_data <- Sigma_chain %>% as.data.frame() %>% 
        rownames_to_column("Z") %>% 
        pivot_longer(-c("Z"), names_to = c("Group", "Run"), names_sep = ":", 
                     values_to = "variance") %>% 
        mutate_at("Run", as.numeric) %>% mutate_if(is.character, as.factor) 
      
      Sigma_chain_plot <- 
        ggplot(data = Sigma_chain_data, 
               aes(x = Run, y = variance, colour = Z)) +       
        geom_line(aes(group = Z)) + 
        geom_vline(xintercept = warm_up, size = 1) +
        xlab("Run") + ylab("Variance") +  
        scale_color_manual(values = 
                             rev(colorRampPalette(wes_palette("Darjeeling1"))(
                               nrow(Sigma_chain)))) + 
        scale_x_continuous(breaks = seq(0, B, by = 100)) + 
        facet_grid(rows = vars(factor(Group, 
                                      levels = c("Unimpaired", "MCI", "Dementia", 
                                                 "Other"))), 
                   scales = "free") + theme_bw() + 
        theme(legend.position = "bottom")
      
      if(!calibration_sample){
        ggsave(filename = "Sigma_chain.jpeg", plot = Sigma_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "no_calibration_sample/run_", run_number), 
               height = 7, width = 12, units = "in", device = "jpeg")
      } else{
        ggsave(filename = "Sigma_chain.jpeg", plot = Sigma_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "calibration_", calibration_sample_name, "/run_", 
                             run_number), height = 7, width = 12, units = "in", 
               device = "jpeg")
      }
      
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
        facet_grid(rows = vars(factor(Z)), 
                   cols = vars(factor(cell_name)), 
                   scales = "free") + theme_bw() + 
        theme(legend.position = "bottom")
      
      if(!calibration_sample){
        ggsave(filename = "mu_chain.jpeg", plot = mu_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "no_calibration_sample/run_", run_number), 
               height = 15, width = 12, units = "in", device = "jpeg")
      } else{
        ggsave(filename = "mu_chain.jpeg", plot = mu_chain_plot, 
               path = paste0(path_to_figures_folder, "diagnostics/", 
                             "calibration_", calibration_sample_name, "/run_", 
                             run_number), height = 15, width = 12, units = "in", 
               device = "jpeg")
      }
      
      #---- save datasets ----
      if(!calibration_sample){
        write_csv(gamma_plot_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "no_calibration_sample/run_", run_number, 
                                "/gamma_plot_data.csv"))
        
        write_csv(latent_class_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "no_calibration_sample/run_", run_number, 
                                "/latent_class_data.csv"))
        
        write_csv(pi_chain_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "no_calibration_sample/run_", run_number, 
                                "/pi_chain_data.csv"))
        
        write_csv(Sigma_chain_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "no_calibration_sample/run_", run_number, 
                                "/Sigma_chain_data.csv"))
        
        write_csv(mu_chain_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "no_calibration_sample/run_", run_number, 
                                "/mu_chain_data.csv"))
      } else{
        write_csv(gamma_plot_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "calibration_", calibration_sample_name, "/run_", 
                                run_number, "/gamma_plot_data.csv"))
        
        write_csv(latent_class_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "calibration_", calibration_sample_name, "/run_", 
                                run_number, "/latent_class_data.csv"))
        
        write_csv(pi_chain_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "calibration_", calibration_sample_name, "/run_", 
                                run_number, "/pi_chain_data.csv"))
        
        write_csv(Sigma_chain_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "calibration_", calibration_sample_name, "/run_", 
                                run_number, "/Sigma_chain_data.csv"))
        
        write_csv(mu_chain_data, 
                  file = paste0(path_to_analyses_folder, "diagnostics_data/", 
                                "calibration_", calibration_sample_name, "/run_", 
                                run_number, "/mu_chain_data.csv"))
      }
    }
  }

# #---- test function ----
# warm_up = 100
# run_number = 1
# starting_props = c(0.25, 0.25, 0.25, 0.25)
# unimpaired_preds
# other_preds
# mci_preds
# categorical_vars = W
# continuous_vars = Z
# id_var = "HHIDPN"
# variable_labels
# dataset_to_copy = synthetic_HCAP_list[[1]]
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
# nu_0
# kappa_0_mat
# contrasts_matrix = A
# path_to_analyses_folder =
#   paste0(path_to_box, "analyses/simulation_study/HCAP_HRS_",
#          unique(synthetic_HCAP_list[[1]][, "dataset_name"]),
#          "/")
# path_to_figures_folder =
#   paste0(path_to_box,
#          "figures/simulation_study/HCAP_HRS_",
#          unique(synthetic_HCAP_list[[1]][, "dataset_name"]),
#          "/")
# data_only = FALSE
