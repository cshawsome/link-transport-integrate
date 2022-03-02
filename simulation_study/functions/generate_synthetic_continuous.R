generate_synthetic_continuous <- function(data, sample_size, dist, parameters){
  #--- bootstrap data ----
  synthetic_data <- sample_n(data, size = sample_size, replace = TRUE)
  
  return(synthetic_data)
}

#Testing
data <- HRS_analytic
sample_size <- 2000
dist <- "normal"
parameters <- normal_parameter_list
