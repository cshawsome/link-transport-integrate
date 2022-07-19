specify_priors <- 
  function(prior_data = "ADAMS", calibration_sample = NA, calibration_prop = NA, 
           path_to_box){
    if(is.na(calibration_sample)){
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
      
      #---- ****contingency cells ----
      alpha_0_dist <- 
        readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                       "imputation_cell_props")) 
      
      #--- ****beta and sigma ----
      priors_beta <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                                    "prior_data/priors_beta")) 
      prior_V_inv <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                                    "prior_data/priors_V_inv"))  
      prior_Sigma <- readRDS(paste0(path_to_box, "analyses/simulation_study/",
                                    "prior_data/priors_Sigma")) 
    }
    
    #---- return ----
    return("unimpaired_betas" = unimpaired_betas)
  }
