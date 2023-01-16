analysis_function <- 
  function(warm_up, starting_props, dataset_to_copy, orig_means,
           orig_sds, calibration_sample = FALSE, calibration_prop = NA, 
           calibration_sample_name = NA, path_to_data, categorical_vars, 
           continuous_vars, id_var, variable_labels, cell_ID_key, color_palette, 
           contrasts_matrix, kappa_0_mat, nu_0_mat, num_synthetic){
    
    #---- pre-allocated results ----
    result_names <- 
      c("mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other",
        "LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other",
        "UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other", 
        "dementia_logOR_age5", "dementia_SE_age5", 
        "mci_logOR_age5", "mci_SE_age5",
        "dementia_logOR_edu", "dementia_SE_edu",
        "mci_logOR_edu", "mci_SE_edu",
        "dementia_logOR_female", "dementia_SE_female",
        "mci_logOR_female", "mci_SE_female",
        "dementia_logOR_black", "dementia_SE_black", 
        "dementia_logOR_hispanic", "dementia_SE_hispanic", 
        "mci_logOR_black", "mci_SE_black", 
        "mci_logOR_hispanic", "mci_SE_hispanic")
    
    results <- matrix(ncol = length(result_names), nrow = 1) %>% 
      set_colnames(all_of(result_names))
    
    if(calibration_sample){
      stop("This script is not equipped to handle calibration subsets yet.")
    } else{
      calibration_sample_name <- "ADAMS_prior"
    }
    
    #---- generate synthetic data ----
    synthetic_HCAP <- 
      generate_synthetic(warm_up, run_number = NA, 
                         starting_props = starting_props, 
                         dataset_to_copy = dataset_to_copy, 
                         orig_means = orig_means, orig_sds = orig_sds,
                         calibration_sample = calibration_sample, 
                         calibration_prop = calibration_prop, 
                         calibration_sample_name = calibration_sample_name,
                         path_to_data = path_to_data, 
                         path_to_analyses_folder = NA, 
                         path_to_figures_folder = NA, 
                         categorical_vars = categorical_vars, 
                         continuous_vars = continuous_vars, id_var = id_var, 
                         variable_labels = variable_labels, 
                         cell_ID_key = cell_ID_key,
                         color_palette = color_palette, 
                         contrasts_matrix = contrasts_matrix, 
                         kappa_0_mat = kappa_0_mat, 
                         nu_0_mat = nu_0_mat, num_synthetic = num_synthetic, 
                         data_only = TRUE)
    
    #---- function to clean missing counts ----
    all_classes <- c("Unimpaired", "MCI", "Dementia", "Other")
    
    clean_counts <- function(counts, all_classes){
      missing_class <- setdiff(all_classes, names(counts))
      
      if(length(missing_class) > 0){counts[missing_class] <- 0}
      
      return(counts[all_classes])
    } 
    
    #---- synthetic impairment class counts ----
    counts <- 
      lapply(synthetic_HCAP, function(x){
        temp <- table(x[, c("Unimpaired", "MCI", "Dementia", "Other")]) %>% 
          as.data.frame() %>% filter(Freq != 0) %>% 
          pivot_longer(c("Unimpaired", "MCI", "Dementia", "Other"), 
                       names_to = "Group") %>% filter(value == 1) %>% 
          dplyr::select(-one_of("value"))
        
        vec <- temp$Freq %>% set_names(temp$Group)
        
        return(vec)}) %>% 
      lapply(., function(x) clean_counts(x, all_classes)) %>% 
      do.call(rbind, .)
    
    #---- **mean counts ----
    mean_counts <- round(colMeans(counts))
    results[, paste0("mean_", names(mean_counts))] <- mean_counts
    
    #---- **CI ----
    results[, paste0("LCI_", colnames(counts))] <- 
      apply(counts, 2, function(x) quantile(x, 0.025))
    
    results[, paste0("UCI_", colnames(counts))] <- 
      apply(counts, 2, function(x) quantile(x, 0.975))
    
    #---- age associations ----
    #---- **create 5-year age bands ----
    synthetic_HCAP %<>% lapply(., function(data) 
      data %<>% mutate("age5" = age/5))
    
    #---- **fit models ----
    dementia_results_list <- lapply(synthetic_HCAP, function(x){
      model <- glm(Dementia ~ age5, family = binomial(link = "logit"), data = x)
      return(summary(model)$coefficients["age5", c("Estimate", "Std. Error")])
    })
    
    mci_results_list <- lapply(synthetic_HCAP, function(x){
      model <- glm(MCI ~ age5, family = binomial(link = "logit"), data = x)
      return(summary(model)$coefficients["age5", c("Estimate", "Std. Error")])
    })
    
    #---- **summarize results ----
    results[, c("dementia_logOR_age5", "dementia_SE_age5")] <- 
      do.call(rbind, dementia_results_list) %>% colMeans()
    
    results[, c("mci_logOR_age5", "mci_SE_age5")] <- 
      do.call(rbind, mci_results_list) %>% colMeans()
    
    #---- education associations ----
    #---- **fit models ----
    dementia_results_list <- lapply(synthetic_HCAP, function(x){
      model <- glm(Dementia ~ edyrs, family = binomial(link = "logit"), data = x)
      return(summary(model)$coefficients["edyrs", c("Estimate", "Std. Error")])
    })
    
    mci_results_list <- lapply(synthetic_HCAP, function(x){
      model <- glm(MCI ~ edyrs, family = binomial(link = "logit"), data = x)
      return(summary(model)$coefficients["edyrs", c("Estimate", "Std. Error")])
    })
    
    #---- **summarize results ----
    results[, c("dementia_logOR_edu", "dementia_SE_edu")] <- 
      do.call(rbind, dementia_results_list) %>% colMeans()
    
    results[, c("mci_logOR_edu", "mci_SE_edu")] <- 
      do.call(rbind, mci_results_list) %>% colMeans()
    
    #---- race/ethnicity associations ----
    #---- **fit models ----
    dementia_results_list <- lapply(synthetic_HCAP, function(x){
      model <- 
        glm(Dementia ~ black + hispanic, family = binomial(link = "logit"), 
            data = x)
      return(
        c(summary(model)$coefficients["black", c("Estimate", "Std. Error")], 
          summary(model)$coefficients["hispanic", c("Estimate", "Std. Error")]) %>% 
          set_names(c("black", "black_stderror", "hispanic", "hispanic_stderror")))
    })
    
    mci_results_list <- lapply(synthetic_HCAP, function(x){
      model <- 
        glm(MCI ~ black + hispanic, family = binomial(link = "logit"), 
            data = x)
      return(
        c(summary(model)$coefficients["black", c("Estimate", "Std. Error")], 
          summary(model)$coefficients["hispanic", c("Estimate", "Std. Error")]) %>% 
          set_names(c("black", "black_stderror", "hispanic", "hispanic_stderror")))
    })
    
    #---- **summarize results ----
    results[, c("dementia_logOR_black", "dementia_SE_black", 
                "dementia_logOR_hispanic", "dementia_SE_hispanic")] <- 
      do.call(rbind, dementia_results_list) %>% colMeans()
    
    results[, c("mci_logOR_black", "mci_SE_black", 
                "mci_logOR_hispanic", "mci_SE_hispanic")] <- 
      do.call(rbind, mci_results_list) %>% colMeans()
    
    #---- sex/gender associations ----
    #---- **fit models ----
    dementia_results_list <- lapply(synthetic_HCAP, function(x){
      model <- glm(Dementia ~ female, family = binomial(link = "logit"), data = x)
      return(summary(model)$coefficients["female", c("Estimate", "Std. Error")])
    })
    
    mci_results_list <- lapply(synthetic_HCAP, function(x){
      model <- glm(MCI ~ female, family = binomial(link = "logit"), data = x)
      return(summary(model)$coefficients["female", c("Estimate", "Std. Error")])
    })
    
    #---- **summarize results ----
    results[, c("dementia_logOR_female", "dementia_SE_female")] <- 
      do.call(rbind, dementia_results_list) %>% colMeans()
    
    results[, c("mci_logOR_female", "mci_SE_female")] <- 
      do.call(rbind, mci_results_list) %>% colMeans()
    
    #---- return results ----
    return(results)
  }

# #---- test function ----
# warm_up = 100
# run_number = 1
# starting_props = c(0.25, 0.25, 0.25, 0.25)
# dataset_to_copy = HCAP_imputed[[1]]
# orig_means = dataset_to_copy %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
#   colMeans() %>% t() %>% as.data.frame()
# orig_sds = dataset_to_copy %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
#   apply(., 2, sd) %>% t() %>% as.data.frame()
# calibration_sample = FALSE
# calibration_prop = NA
# calibration_sample_name = NA
# path_to_data = paste0(path_to_box,"data/")
# categorical_vars = W
# continuous_vars = Z
# id_var = "HHIDPN"
# variable_labels = variable_labels
# cell_ID_key = cell_ID_key
# color_palette = color_palette
# contrasts_matrix = A
# kappa_0_mat = read_csv(paste0(path_to_box, "data/tuning/kappa_0_matrix_HCAP.csv"))
# nu_0_mat = read_csv(paste0(path_to_box, "data/tuning/nu_0_matrix_HCAP.csv"))
# num_synthetic = 1000
# data_only = FALSE
# 
