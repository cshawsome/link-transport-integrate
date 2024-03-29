generate_synthetic_continuous <- 
  function(data, sample_size, dementia_prop = 0.25, mci_prop = 0.25, 
           other_prop = 0.25, dist, parameters, selected_vars_estimates, 
           path_to_results, scenario_name = NA){
    #---- check parameters ----
    if(sum(dementia_prop, mci_prop, other_prop) > 1){
      stop("Impairment proportions sum to more than 1. Please check values.") 
    }
    
    if(!dir.exists(path_to_results)){
      dir.create(path_to_results, recursive = TRUE)
    }
    
    #--- bootstrap data ----
    synthetic_data <- sample_n(data, size = sample_size, replace = TRUE) %>% 
      mutate("HHIDPN" = seq(1, sample_size))
    
    #---- assign impairment category ----
    unimpaired_prop <- 1 - sum(dementia_prop, mci_prop, other_prop)
    
    for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
      if(group != "Unimpaired"){
        synthetic_data %<>% mutate(!!paste0(group, "_prob") := 0)
      }
      synthetic_data %<>% mutate({{group}} := 0)
    }
    
    #---- **predicted probabilities ----
    for(group in c("Dementia", "MCI", "Other")){
      X <- as.matrix(synthetic_data[, selected_vars_estimates$data_label])
      beta <- as.matrix(selected_vars_estimates[, group])
      
      synthetic_data[, paste0(group, "_prob")] <- exp(X %*% beta)
    }
    
    #---- **group assignment ----
    for(group in c("Dementia", "MCI", "Other", "Unimpaired")){
      num_people <- round(get(paste0(tolower(group), "_prop"))*sample_size)
      
      synthetic_data %<>% 
        mutate("group_membership" = 
                 rowSums(synthetic_data[, c("Unimpaired", "Other", "MCI", 
                                            "Dementia")])) 
      
      subset <- synthetic_data %>% filter(group_membership == 0)
      
      if(group == "Unimpaired"){
        indices <- subset %>% dplyr::select("HHIDPN") %>% unlist()
        
      } else{
        indices <- 
          sample_n(subset, size = num_people, 
                   weight = unlist(subset[, paste0(group, "_prob")])) %>% 
          dplyr::select("HHIDPN") %>% unlist()
          
      }
      
      synthetic_data[indices, group] <- 1
    }
    
    #---- **get rid of extra column ----
    synthetic_data %<>% dplyr::select(-one_of("group_membership"))
    
    #---- scenario name ----
    if(is.na(scenario_name)){
      scenario_name <- tolower(names(which.max(
        colSums(synthetic_data[, c("Unimpaired", "MCI", "Dementia", "Other")]))))
    }
    
    #---- reformat synthetic data ----
    synthetic_data %<>% as.matrix()  
    
    #---- generate synthetic data ----
    #---- **shell for continuous vars ----
    synthetic_data %<>% 
      cbind(., matrix(nrow = sample_size, 
                      ncol = ncol(parameters$Unimpaired$beta_center)) %>%
              set_colnames(colnames(parameters$Unimpaired$beta_center)))
    
    #---- **draw continuous variables ----
    for(i in 1:nrow(synthetic_data)){
      group <- names(
        which(
          synthetic_data[i, c("Unimpaired", "MCI", "Dementia", "Other")] == 1))
      
      #---- **sigma | Y ----
      sigma <- riwish(v = normal_parameter_list[[group]]$sigma_dof, 
                      S = normal_parameter_list[[group]]$sigma_center)
      
      #---- **beta | sigma, Y ----
      beta <- 
        rmatrixnorm(M = normal_parameter_list[[group]]$beta_center, 
                    U = as.positive.definite(normal_parameter_list[[group]]$row_cov), 
                    V = as.positive.definite(sigma))
      
      #---- **predicted Y ----
      X = synthetic_data[i, rownames(beta)]
      
      #---- ****normal ----
      if(dist == "normal"){
        synthetic_data[i, colnames(beta)] = 
          X %*% beta + rnorm(n = length(colnames(beta)), mean = 0, sd = 0.5)
        #---- ****lognormal ----
      } else if(dist == "lognormal"){
        synthetic_data[i, colnames(beta)] = 
          X %*% beta + 
          exp(rnorm(n = length(colnames(beta)), mean = 0, sd = 0.5))
        #---- ****bathtub ----
      } else if(dist == "bathtub"){
        #subset beta-- some of the variables are reverse coded so we want to 
        # flip the skew for them
        beta_no_flip <- beta[, which(!colnames(beta) %in% 
                                       c("age_Z", "adl_Z", "trailsA_Z"))]
        beta_flip <- beta[, which(colnames(beta) %in% 
                                    c("age_Z", "adl_Z", "trailsA_Z"))] 
        
        if(group %in% c("MCI", "Dementia")){
          synthetic_data[i, colnames(beta_no_flip)] = 
            X %*% beta_no_flip - 
            2*rbeta(n = length(colnames(beta_no_flip)), shape1 = 5, shape2 = 1)
          
          synthetic_data[i, colnames(beta_flip)] = 
            X %*% beta_flip + 
            2*rbeta(n = length(colnames(beta_flip)), shape1 = 5, shape2 = 1)
        } else{
          synthetic_data[i, colnames(beta_no_flip)] = 
            X %*% beta_no_flip + 
            2*rbeta(n = length(colnames(beta_no_flip)), shape1 = 5, shape2 = 1)
          
          synthetic_data[i, colnames(beta_flip)] = 
            X %*% beta_flip - 
            2*rbeta(n = length(colnames(beta_flip)), shape1 = 5, shape2 = 1)
        }
      } else{
        stop(paste0("Invalid distribution argument. Please choose from: 'normal'", 
                    ", 'lognormal', or 'bathtub'."))
      }
    }
    #---- return values ----
    write_csv(as.data.frame(synthetic_data), 
              paste0(path_to_results, dist, "_", sample_size, "_", 
                     scenario_name, ".csv"))
  }

# #---- testing ----
# data <- HRS_imputed
# sample_size <- 1000
# dementia_prop = 0.35
# mci_prop = 0.10
# other_prop = 0.20
# dist <- "normal"
# scenario_name = NA
# parameters <- normal_parameter_list
# selected_vars_estimates <- fixed_betas
# path_to_results = paste0(path_to_box, "analyses/simulation_study/",
#                          "test_superpopulations/")
# 
# generate_synthetic_continuous(data, sample_size, dementia_prop,
#                               mci_prop, other_prop, dist, 
#                               selected_vars_estimates, parameters, 
#                               path_to_results)


