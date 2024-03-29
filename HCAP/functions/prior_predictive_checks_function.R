prior_predictive_checks <- 
  function(dataset_to_copy, calibration_sample = FALSE, calibration_prop = NA, 
           calibration_sample_name = NA, path_to_data, path_to_output_folder, 
           continuous_check_test = FALSE, 
           continuous_check = c("Unimpaired", "MCI", "Dementia", "Other"), 
           categorical_vars, continuous_vars, variable_labels, color_palette,
           contrasts_matrix, weights_matrix, kappa_0_mat, nu_0_mat, 
           num_synthetic){
    
    #---- check calibration parameter ----
    if(calibration_sample){
      if(calibration_prop > 1 | calibration_prop < 0) 
        stop("Calibration prop must be between 0 and 1")
    }
    
    #---- update path to output folder ----
    if(!calibration_sample){
      path_to_output_folder <- 
        paste0(path_to_output_folder, "ADAMS_prior/")
    } else{
      path_to_output_folder <- 
        paste0(path_to_output_folder, calibration_sample_name, "/")
    }
    
    #---- create folders for results ----
    if(!dir.exists(paste0(path_to_output_folder, "cell_counts"))){
      dir.create(paste0(path_to_output_folder, "cell_counts"),
                 recursive = TRUE)
    }
    
    if(!dir.exists(paste0(path_to_output_folder, "continuous_vars/combined"))){
      dir.create(paste0(path_to_output_folder, "continuous_vars/combined"), 
                 recursive = TRUE)
    }
    
    if(!dir.exists(paste0(path_to_output_folder, "continuous_vars/error_set/",
                          "combined"))){
      dir.create(paste0(path_to_output_folder, "continuous_vars/error_set/",
                        "combined"), recursive = TRUE)
    }
    
    if(continuous_check_test &
       !dir.exists(paste0(path_to_output_folder, "continuous_vars/test_set/",
                          "combined"))){
      
      dir.create(paste0(path_to_output_folder, "continuous_vars/test_set/",
                        "combined"), recursive = TRUE)
    }
    
    #---- cell ID key ----
    cell_ID_key <- read_csv(paste0(path_to_data, "data/cell_ID_key.csv"))
    
    #---- set prior data ----
    if(!calibration_sample){
      #---- **latent classes ----
      for(group in c("unimpaired", "mci", "other", "dementia")){
        assign(paste0(group, "_betas"), 
               vroom(paste0(path_to_data, "data/prior_data/latent_class_", 
                            group, "_betas.csv"), delim = ","))
        # assign(paste0(group, "_cov"), 
        #        readRDS(paste0(path_to_data, "data/prior_data/latent_class_", 
        #                       group, "_cov")))
        
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
      # #---- read in raw prior sample ----
      # prior_imputed_clean <- readRDS(path_to_raw_prior_sample) %>%
      #   lapply(function(x) mutate_at(x, "HHIDPN", as.numeric)) 
      
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
      
      #dementia model predictors
      dementia_preds <- c("(Intercept)", selected_vars %>% 
                            filter(data_label != "Intercept" & Dementia != 0) %>% 
                            dplyr::select(data_label) %>% unlist())
      
      #---- calibration subset ----
      if(calibration_scenario == "calibration_100"){
        calibration_subset <- dataset_to_copy
      } else{
        calibration_subset <- 
          dataset_to_copy %>% filter(!!sym(calibration_scenario) == 1)
        
        # #test weighting
        # calibration_subset %<>% 
        #   mutate("Group" = case_when(Unimpaired == 1 ~ "Unimpaired", 
        #                              MCI == 1 ~ "MCI", 
        #                              Dementia == 1 ~ "Dementia", 
        #                              Other == 1 ~ "Other")) %>%
        #   unite("cell_ID", c("black", "hispanic", "stroke"), remove = FALSE, 
        #         sep = "")
        # 
        # cell_ID_key_long <- cell_ID_key %>% 
        #   pivot_longer(-c("cell_order", "cell_ID", "cell_name"), 
        #                names_to = c("dataset_name", "Group"), 
        #                values_to = "weights", 
        #                names_sep = "_IPW_") %>% filter()
        # 
        # calibration_subset %<>% 
        #   left_join(., cell_ID_key_long, 
        #             by = c("cell_ID", "Group", "dataset_name"))
        
      }
    }
    
    #---- select variables ----
    vars <- unique(c(unimpaired_preds, other_preds, mci_preds))
    
    synthetic_sample <- dataset_to_copy %>% 
      dplyr::select("HHIDPN", all_of(vars)) %>% 
      #pre-allocate columns
      mutate("p_unimpaired" = 0, "p_other" = 0, "p_mci" = 0, "p_dementia" = 0, 
             "Unimpaired" = 0, "MCI" = 0, "Dementia" = 0, "Other" = 0)
    
    #---- max index ----
    max_index <- length(priors_beta)  
    
    #---- DEBUG ----
    #for(run in 1:1000){
    generate_data <- function(color_palette){
      if(!calibration_sample){
        #---- index for random draws ----
        random_draw <- sample(seq(1, max_index), size = 1)  
      }
      
      #---- latent class ----
      for(model in c("unimpaired", "other", "mci", "dementia")){
        
        if(calibration_sample){
          if(model == "mci"){
            class_name = "MCI"
          } else{
            class_name = str_to_sentence(model)
          }
          
          # latent_class_model <- 
          #   suppressWarnings(glm(formula(
          #     paste(class_name, " ~ ", 
          #           paste(get(paste0(model, "_preds"))[-1], collapse = " + "), 
          #           collapse = "")), family = "binomial", 
          #     #don't select (Intercept) variable
          #     data = slice_sample(calibration_subset[, vars[-1]], 
          #                         n = nrow(calibration_subset), 
          #                         replace = TRUE)))
          # 
          # prior_betas <- coefficients(latent_class_model)
          
          bootstrap_sample <- 
            slice_sample(calibration_subset[, vars[-1]], 
                         n = nrow(calibration_subset), replace = TRUE) %>% 
            mutate("selected" = 1) %>% 
            rbind(synthetic_sample[, vars[-1]] %>% mutate("selected" = 0))
          
          latent_class_model <- 
            suppressWarnings(glm(formula(
              paste(class_name, " ~ ", 
                    paste(get(paste0(model, "_preds"))[-1], collapse = " + "), 
                    collapse = "")), family = "binomial", 
              #don't select (Intercept) variable
              data = bootstrap_sample %>% filter(selected == 1)))
          
          prior_betas <- coefficients(latent_class_model)
          
          while(sum(is.na(prior_betas)) > 0){
            # latent_class_model <- 
            #   suppressWarnings(glm(formula(
            #     paste(class_name, " ~ ", 
            #           paste(get(paste0(model, "_preds"))[-1], collapse = " + "), 
            #           collapse = "")), family = "binomial", 
            #     #don't select (Intercept) variable
            #     data = slice_sample(calibration_subset[, vars[-1]], 
            #                         n = nrow(calibration_subset), 
            #                         replace = TRUE)))
            # 
            # prior_betas <- coefficients(latent_class_model)
            
            bootstrap_sample <- 
              slice_sample(calibration_subset[, vars[-1]], 
                           n = nrow(calibration_subset), replace = TRUE) %>% 
              mutate("selected" = 1) %>% 
              rbind(synthetic_sample[, vars[-1]] %>% mutate("selected" = 0))
            
            latent_class_model <- 
              suppressWarnings(glm(formula(
                paste(class_name, " ~ ", 
                      paste(get(paste0(model, "_preds"))[-1], collapse = " + "), 
                      collapse = "")), family = "binomial", 
                #don't select (Intercept) variable
                data = bootstrap_sample %>% filter(selected == 1)))
            
            prior_betas <- coefficients(latent_class_model)
          }
          
          # prior_cov <- vcov(latent_class_model)
          # 
          # prior_betas <-
          #   t(mvnfast::rmvn(n = 1, mu = unlist(prior_betas), sigma = prior_cov))
          
        } else{
          prior_betas <- get(paste0(model, "_betas"))[, random_draw]
          #prior_cov <- get(paste0(model, "_cov"))[[random_draw]]
        }
        
        # betas <- 
        #   mvnfast::rmvn(n = 1, mu = unlist(prior_betas), sigma = prior_cov)
        
        synthetic_sample[, paste0("p_", model)] <- 
          expit(as.matrix(synthetic_sample[, 
                                           get(paste0(model, "_preds"))]) %*% 
                  as.matrix(prior_betas))
      }
      
      #---- **assign class ----
      #choose Unimpaired if that has the highest probability
      synthetic_sample[, "Group"] <-
        apply(synthetic_sample[, c("p_unimpaired", "p_mci", "p_dementia", 
                                   "p_other")], 1,
              function(x) str_remove(names(which.max(x)), "p_"))
      
      #randomly choose other classes based on predicted probabilities
      synthetic_sample[synthetic_sample$Group != "unimpaired", "Group"] <- 
        apply(synthetic_sample[synthetic_sample$Group != "unimpaired", 
                               c("p_mci", "p_dementia", "p_other")], 1, 
              function(x) 
                sample(c("mci", "dementia", "other"), size = 1, prob = x))
      
      #reformat group labels
      synthetic_sample %<>% 
        mutate("Group" = case_when(Group == "mci" ~ str_to_upper(Group), 
                                   TRUE ~ str_to_sentence(Group)))
      
      synthetic_sample %<>% 
        mutate("Unimpaired" = ifelse(Group == "Unimpaired", 1, 0), 
               "MCI" = ifelse(Group == "MCI", 1, 0), 
               "Dementia" = ifelse(Group == "Dementia", 1, 0), 
               "Other" = ifelse(Group == "Other", 1, 0)) 
      
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
      
      nu_0 <- nu_0_mat[1, ]  
      
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
          prior_counts <- calibration_subset[, c(categorical_vars, class)] %>% 
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
          
          if(str_detect(calibration_sample_name, "SRS_race") | 
             str_detect(calibration_sample_name, "design")){
            
            #Make column for observed sampled counts
            prior_counts$Observed <- prior_counts$Freq
            
            prior_counts %<>% 
              left_join(., 
                        cell_ID_key[, c("cell_ID", calibration_sample_name)] %>% 
                          rename("weights" = calibration_sample_name), 
                        by = "cell_ID")
            prior_counts[is.na(prior_counts)] <- 0
            
            #Full observed count
            prior_counts %<>% mutate("Freq" = Freq*weights)
            prior_UtU <- diag(prior_counts$Observed)
            
          } else{
            prior_UtU <- diag(prior_counts$Freq)
          }
          
          prior_counts %<>% mutate("prop" = Freq/sum(Freq))
          
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
            calibration_subset[, c(categorical_vars, continuous_vars, class)] %>% 
            filter(!!sym(class) == 1) %>% 
            arrange(black, hispanic, stroke) %>%
            dplyr::select(all_of(continuous_vars)) %>% as.matrix()
          
          V_0_inv <- t(A) %*% prior_UtU %*% A
          
          while(is.character(tryCatch(V_0 <- solve(V_0_inv), 
                                      error = function(e) "error"))){
            
            prior_UtU = prior_UtU + 
              diag(1, nrow = nrow(prior_UtU), ncol = ncol(prior_UtU)) 
            
            V_0_inv <- t(A) %*% prior_UtU %*% A
            V_0 <- solve(V_0_inv)
          }
          
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
          cbind(., contrasts_matrix[, -1]) %>% dplyr::rename("Count" = "V1")
        
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
        #---- **save categorical ----
        if(calibration_sample){
          assign(paste0("W_", class), rbind(
            subset[, all_of(categorical_vars)], 
            calibration_subset[calibration_subset[, class] == 1, 
                               all_of(categorical_vars)]))
        } else{
          assign(paste0("W_", class), subset[, all_of(categorical_vars)])
        }
        
        #---- **save continuous ----
        if(calibration_sample){
          assign(paste0("Z_", class), rbind(
            subset[, all_of(continuous_vars)], 
            calibration_subset[calibration_subset[, class] == 1, 
                               all_of(continuous_vars)]))
        } else{
          assign(paste0("Z_", class), subset[, all_of(continuous_vars)])
        }
      }
      #}
      # #---- END DEBUG ----
      if(calibration_sample){
        group <- c(synthetic_sample$Group, 
                   calibration_subset %>% 
                     dplyr::select(
                       c("Unimpaired", "MCI", "Dementia", "Other")) %>% 
                     pivot_longer(everything()) %>% filter(value == 1) %>% 
                     dplyr::select("name") %>% unlist() %>% unname())
      } else{
        group <- synthetic_sample$Group
      }
      
      #---- **return ----
      return(list("Group" = group,
                  "W_unimpaired" = W_Unimpaired %>% 
                    mutate("Group" = "Unimpaired"),
                  "W_other" = W_Other %>% 
                    mutate("Group" = "Other"), 
                  "W_mci" = W_MCI %>% 
                    mutate("Group" = "MCI"), 
                  "W_dementia" = W_Dementia %>% 
                    mutate("Group" = "Dementia"),
                  
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
    #---- **categorical ----
    #pre-allocate results matrix
    synthetic_counts <- 
      matrix(0, nrow = nrow(contrasts_matrix), 
             ncol = (num_synthetic + 1)) %>% as.data.frame() %>% 
      set_colnames(c(seq(1, num_synthetic), "cell"))
    
    cells <- contrasts_matrix[, -1] %>% as.data.frame() %>% 
      unite("cell_ID", sep = "") %>% left_join(cell_ID_key, by = "cell_ID")
    
    synthetic_counts[, "cell"] <- cells$cell_name
    
    #counts from synthetic datasets
    for(num in 1:num_synthetic){
      for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
        counts <- synthetic[[num]][[paste0("W_", tolower(group))]] 
        
        if(nrow(counts) == 0){
          synthetic_counts[, num] <- 0
        } else{
          counts %<>% dplyr::select(-one_of("Group")) %>%
            unite("cell_ID", sep = "") %>% table() %>% as.data.frame() %>% 
            set_colnames(c("cell_ID", "Freq")) %>% 
            left_join(cell_ID_key, ., by = "cell_ID")
          
          counts[which(is.na(counts$Freq)), "Freq"] <- 0
          
          synthetic_counts[
            which(synthetic_counts$cell %in% counts$cell_name), num] <- 
            synthetic_counts[
              which(synthetic_counts$cell %in% counts$cell_name), num] +
            counts$Freq
        }
      }
    }
    
    synthetic_count_plot_data <- synthetic_counts %>% mutate("truth" = 0) 
    
    #true counts
    counts <- dataset_to_copy %>%  
      dplyr::select(all_of(categorical_vars)) %>% 
      unite("cell_ID", sep = "") %>% table() %>% as.data.frame() %>% 
      set_colnames(c("cell_ID", "Freq")) %>% 
      left_join(cell_ID_key, ., by = "cell_ID")
    
    counts[which(is.na(counts$Freq)), "Freq"] <- 0
    
    synthetic_count_plot_data[
      which(synthetic_count_plot_data$cell %in% counts$cell_name), 
      "truth"] <- counts$Freq
    
    synthetic_count_plot_data %<>% 
      pivot_longer(-c("cell", "truth")) 
    
    ggplot(data = synthetic_count_plot_data , aes(x = value)) + 
      geom_histogram(fill = "black", color = "black") + theme_bw() + 
      xlab("Contingency Cell Count") + ylab("") + 
      facet_wrap(facets = vars(cell), ncol = 2, scales = "free") +
      geom_vline(aes(xintercept = truth), color = "#f2caaa",
                 size = 1) +
      theme(text = element_text(size = 8), strip.text = element_text(size = 8))  
    
    ggsave(filename = paste0(path_to_output_folder, "cell_counts/", 
                             "overall_count.jpeg"), 
           width = 3, height = 4, units = "in")
    
    #---- **continuous ----
    true_data <- dataset_to_copy %>%
      dplyr::select(c(all_of(continuous_vars))) %>% mutate("Color" = "black")
    
    for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
      
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
# dataset_to_copy = HCAP_analytic
# calibration_sample = FALSE
# calibration_prop = NA
# calibration_sample_name = NA
# path_to_raw_prior_sample = NA
# path_to_data = path_to_box
# path_to_output_folder =
#   paste0(path_to_box,
#          "figures/chapter_6/prior_predictive_checks/")
# continuous_check_test = TRUE
# continuous_check = c("Unimpaired", "MCI", "Dementia", "Other")
# categorical_vars = W
# continuous_vars = Z
# variable_labels = variable_labels
# color_palette = color_palette
# contrasts_matrix = A
# weights_matrix = NA
# kappa_0_mat = kappa_0_mat
# nu_0_mat = nu_0_mat
# num_synthetic = 1000