generate_synthetic_continuous <- 
  function(data, sample_size, unimpaired_prop = 0.25, mci_prop = 0.25, 
           dementia_prop = 0.25, dist, parameters){
    #---- check parameters ----
    if(sum(unimpaired_prop, mci_prop, dementia_prop) > 1){
      stop("Impairment proportions sum to more than 1. Please check values.") 
    }
    
    #--- bootstrap data ----
    synthetic_data <- sample_n(data, size = sample_size, replace = TRUE)
    
    #---- assign impairment category ----
    other_prop <- 1 - sum(unimpaired_prop, mci_prop, dementia_prop)
    
    for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
      synthetic_data %<>% mutate({{group}} := 0)
    }
    
    synthetic_data[1:round(unimpaired_prop*sample_size), "Unimpaired"] <- 1
    synthetic_data[1:round(unimpaired_prop*sample_size), "MCI"] <- 1
    
    return(synthetic_data)
  }

#Testing
data <- HRS_analytic
sample_size <- 2000
unimpaired_prop = 0.3
mci_prop = 0.3
dementia_prop = 0.3
dist <- "normal"
parameters <- normal_parameter_list

test <- 
  generate_synthetic_continuous(data, sample_size, unimpaired_prop = 0.3, 
                                mci_prop = 0.3, dementia_prop = 0.3, dist, 
                                parameters)


