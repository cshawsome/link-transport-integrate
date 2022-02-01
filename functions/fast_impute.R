fast_impute <- 
  function(predictor_matrix, data, method, m, maxit, save = "no"){
    
    #---- where matrix ----
    where <- is.na(data)*1  
    impute_vars <- rownames(predictor_matrix)
    
    #---- imputation 0: mean imputation ----
    # #check var types: only HHIDPN and Adem_dx_label_collapsed are non-numeric
    # which(lapply(data, class) != "numeric")
    avgs <- colMeans(data[, -(which(lapply(data, class) != "numeric"))], 
                     na.rm = TRUE)
    
    for(var in impute_vars){
      data[where[, var] == 1 , var] <- avgs[var]
    }
    
    #---- add intercept ----
    data[, "intercept"] <- 1
    
    #---- pre-allocate chain storage ---- 
    if(save == "yes"){
      #colnames are imputation number:iteration number:stat number
      #stat number: 1 = mean; 2 = sd
      trace_data <- matrix(nrow = length(impute_vars), ncol = 2*m*maxit) %>% 
        set_colnames(apply(expand.grid(seq(1:m), seq(1:maxit), seq(1, 2)), 
                           1, paste, collapse = ":")) %>% 
        set_rownames(impute_vars)
    }
    
    #---- pre-allocate list of imputed datasets ----
    impute_list <- list()
    
    #---- imputation loop ----
    for(run in 1:m){
      imputed_data <- data
      for(iter in 1:maxit){
        for(var in impute_vars){
          imputed_data[where[, var] == 1, var] <- NA
          #---- **PMM ----
          if(method == "PMM"){
            #fastMice will only do PMM if the outcome variable is factored
            imputed_data[, var] <- factor(unlist(imputed_data[, var]))
            
            imputed_data[, var] <- 
              as.numeric(as.character(
                fill_NA_N(imputed_data, model = "pmm", posit_y = var, k = 10,
                          posit_x = c("intercept", names(
                            predictor_matrix[
                              var, which(predictor_matrix[var, ] == 1)])))))
          } else{
            imputed_data[, var] <- 
              as.numeric(fill_NA(imputed_data, model = "lm_bayes", 
                                 posit_y = var, 
                                 posit_x = c("intercept", names(
                                   predictor_matrix[
                                     var, 
                                     which(predictor_matrix[var, ] == 1)]))))
          }
          
          if(exists("trace_data")){
            #in the third slot: 1 = mean, 2 = sd
            trace_data[var, paste0(c(run, iter, 1), collapse = ":")] <- 
              mean(unlist(imputed_data[where[, var] == 1, var]))
            trace_data[var, paste0(c(run, iter, 2), collapse = ":")] <- 
              sd(unlist(imputed_data[where[, var] == 1, var]))
          }
        }
      }
      impute_list[[run]] <- imputed_data
    }
    
    #---- save results ----
    if(save == "yes"){
      #where matrix
      write_csv(as.data.frame(where), 
                file = here::here("MI datasets", 
                                  paste0("where_", tolower(method), "_", 
                                         tolower(mechanism), 
                                         as.numeric(sub("%","", 
                                                        mask_percent)), 
                                         ".csv")))
      
      #trace_data plots data
      write_csv(as.data.frame(trace_data), 
                file = here::here("MI datasets", 
                                  paste0("trace_data_", tolower(method), "_", 
                                         tolower(mechanism), 
                                         as.numeric(sub("%","", 
                                                        mask_percent)),
                                         ".csv")))
    }
    
    #---- return ----
    return(impute_list)
  }

#---- function testing ----
predictor_matrix <- predict
data <- ADAMS_analytic

test <- fast_impute(predictor_matrix = predict, data = ADAMS_analytic, 
                    method = "PMM", m = 2, maxit = 5, save = "no")
















