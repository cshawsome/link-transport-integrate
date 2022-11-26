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
        "UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other")
    
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
    
    return(results)
  }

# #---- test function ----
# warm_up = 100
# run_number = 1
# starting_props = c(0.25, 0.25, 0.25, 0.25)
# dataset_to_copy = HCAP_analytic
# orig_means = HCAP_analytic %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
#   colMeans() %>% t() %>% as.data.frame()
# orig_sds = HCAP_analytic %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
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
