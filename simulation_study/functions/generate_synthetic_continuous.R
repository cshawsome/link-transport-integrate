generate_synthetic_continuous <- 
  function(data, sample_size, unimpaired_prop = 0.25, mci_prop = 0.25, 
           dementia_prop = 0.25, dist, parameters, path_to_results){
    #---- check parameters ----
    if(sum(unimpaired_prop, mci_prop, dementia_prop) > 1){
      stop("Impairment proportions sum to more than 1. Please check values.") 
    }
    
    if(!dir.exists(path_to_results)){
      dir.create(path_to_results, recursive = TRUE)
    }
    
    #--- bootstrap data ----
    synthetic_data <- sample_n(data, size = sample_size, replace = TRUE)
    
    #---- assign impairment category ----
    other_prop <- 1 - sum(unimpaired_prop, mci_prop, dementia_prop)
    
    for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
      synthetic_data %<>% mutate({{group}} := 0)
    }
    
    #---- **Unimpaired ----
    synthetic_data[1:round(unimpaired_prop*sample_size), "Unimpaired"] <- 1
    
    #---- **MCI ----
    first_1 <- max(which(synthetic_data$Unimpaired == 1)) + 1
    last_1 <- first_1 + round(mci_prop*sample_size)
    synthetic_data[first_1:last_1, "MCI"] <- 1
    
    #---- **Dementia ----
    first_1 <- max(which(synthetic_data$MCI == 1)) + 1
    last_1 <- first_1 + round(dementia_prop*sample_size)
    synthetic_data[first_1:last_1, "Dementia"] <- 1
    
    #---- **Other ----
    first_1 <- max(which(synthetic_data$Dementia == 1)) + 1
    last_1 <- sample_size
    synthetic_data[first_1:last_1, "Other"] <- 1
    
    #---- scenario name ----
    scenario_name <- tolower(names(which.max(
      colSums(synthetic_data[, c("Unimpaired", "MCI", "Dementia", "Other")]))))
    
    #---- reformat synthetic data ----
    synthetic_data %<>% as.matrix()  
    
    #---- generate synthetic data ----
    #---- **normal distribution ----
    if(dist == "normal"){
      #---- ****shell for continuous vars ----
      synthetic_data %<>% 
        cbind(., matrix(nrow = sample_size, 
                        ncol = ncol(parameters$Unimpaired$beta_center)) %>%
                set_colnames(colnames(parameters$Unimpaired$beta_center)))
      
      #---- ****draw continuous variables ----
      for(i in 1:nrow(synthetic_data)){
        group <- names(
          which(
            synthetic_data[i, c("Unimpaired", "MCI", "Dementia", "Other")] == 1))
        
        #---- ****sigma | Y ----
        sigma <- riwish(v = normal_parameter_list[[group]]$sigma_dof, 
                        S = normal_parameter_list[[group]]$sigma_center)
        
        #---- ****beta | sigma, Y ----
        beta <- matrix.normal(M = normal_parameter_list[[group]]$beta_center, 
                              U = normal_parameter_list[[group]]$row_cov, 
                              V = sigma)
        
        #---- ****predicted Y ----
        X = synthetic_data[i, rownames(beta)]
        synthetic_data[i, colnames(beta)] = 
          X%*%beta + rnorm(n = length(colnames(beta)), mean = 0, sd = 0.5)
      }
    } else{
      stop("Invalid distribution argument. Please choose from: \"normal\".")
    }
    
    #---- return values ----
    write_csv(as.data.frame(synthetic_data), 
              paste0(path_to_results, "synthetic_", dist, "_", sample_size, 
                     "_", scenario_name, ".csv"))
  }

# #---- testing ----
# data <- HRS_analytic
# sample_size <- 1000
# unimpaired_prop = 0.35
# mci_prop = 0.10
# dementia_prop = 0.35
# dist <- "normal"
# parameters <- normal_parameter_list
# path_to_results = paste0(path_to_box, "analyses/simulation_study/", 
#                          "synthetic_data/")
# 
# generate_synthetic_continuous(data, sample_size, unimpaired_prop,
#                               mci_prop, dementia_prop, dist,
#                               parameters, path_to_results)


