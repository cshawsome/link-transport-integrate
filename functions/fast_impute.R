fast_impute <- 
  function(predictor_matrix, data, path_for_output, method, m, maxit){
    
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
    #colnames are imputation number:iteration number:stat number
    trace_data <- matrix(nrow = length(impute_vars), ncol = 2*m*maxit) %>% 
      set_colnames(apply(expand.grid(seq(1:m), seq(1:maxit), c("mean", "sd")), 
                         1, paste0, collapse = ":")) %>% 
      set_rownames(impute_vars)
    colnames(trace_data) <- gsub(" ", "", colnames(trace_data))  
    
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
          
          trace_data[var, paste0(c(run, iter, "mean"), collapse = ":")] <- 
            mean(unlist(imputed_data[where[, var] == 1, var]))
          trace_data[var, paste0(c(run, iter, "sd"), collapse = ":")] <- 
            sd(unlist(imputed_data[where[, var] == 1, var]))
        }
      }
      impute_list[[run]] <- imputed_data
    }
    
    #---- save results ----
    #create directory for results
    dir.create(here::here(path_for_output, "MI"))
    
    #---- **where matrix ----
    write_csv(as.data.frame(where), 
              file = paste0(path_for_output, "MI/where.csv"))
    
    #---- **trace plots data ----
    trace_data %<>% as.data.frame() %>% 
      rownames_to_column(var = "impute_vars")
    
    write_csv(trace_data, file = paste0(path_for_output, "MI/trace_data.csv"))
    
    #---- **trace plots ----
    plot_data <- trace_data %>% 
      pivot_longer(-"impute_vars", names_to = c("run", "iteration", "stat"), 
                   names_sep = ":") %>% 
      mutate_at(.vars = c("iteration"), as.integer) %>%
      mutate_at(.vars = c("run"), 
                function(x) factor(x, levels = seq(1, maxit)))
    
    trace_plots <- 
      ggplot(plot_data, aes(x = iteration, y = value, color = run)) +
      geom_line(aes(group = run)) + theme_bw() + 
      facet_wrap_paginate(impute_vars~stat, ncol = 2, nrow = 6, page = 1, 
                          scales = "free", strip.position = "top")
    
    n = n_pages(trace_plots)
    
    pdf(paste0(path_for_output, "MI/trace_plots.pdf"), paper = "letter", 
        height = 10.5, width = 8)
    
    for(i in 1:n){
      print(ggplot(plot_data, aes(x = iteration, y = value, color = run)) +
              geom_line(aes(group = run)) + theme_bw() + 
              facet_wrap_paginate(impute_vars~stat, ncol = 2, nrow = 6, 
                                  page = i, scales = "free", 
                                  strip.position = "top"))
    }
    
    dev.off()
    
    #---- **imputed data ----
    saveRDS(impute_list, file = paste0(path_for_output, "MI/MI_datasets"))
  }

#---- function testing ----
start <- Sys.time()
fast_impute(predictor_matrix = predict, data = ADAMS_analytic, 
            path_for_output = paste0(path_to_box, "data/ADAMS/cleaned/"), 
            method = "PMM", m = 2, maxit = 15)
end <- Sys.time() - start

#test output
test_where <- read_csv(file = here::here("ADAMS", "MI", "where.csv"))
test_trace_data <- read_csv(file = here::here("ADAMS", "MI", "trace_data.csv"))
test_impute_data <- readRDS(file = here::here("ADAMS", "MI", "MI_datasets"))












