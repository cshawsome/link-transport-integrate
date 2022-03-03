generate_synthetic_continuous <- 
  function(data, sample_size, normal_prop = 0.25, mci_prop = 0.25, 
           dementia_prop = 0.25, dist, parameters){
    #---- check parameters ----
    if(sum(normal_prop, mci_prop, dementia_prop) > 1){
      stop("Impairment proportions sum to more than 1. Please check values.") 
    }
    
    #--- bootstrap data ----
    synthetic_data <- sample_n(data, size = sample_size, replace = TRUE)
    
    #---- assign impairment category ----
    
    return(synthetic_data)
  }

#Testing
data <- HRS_analytic
sample_size <- 2000
normal_prop = 0.5
mci_prop = 0.3
dementia_prop = 0.3
dist <- "normal"
parameters <- normal_parameter_list

test <- 
  generate_synthetic_continuous(data, sample_size, normal_prop = 0.3, 
                                mci_prop = 0.3, dementia_prop = 0.3, dist, 
                                parameters)


