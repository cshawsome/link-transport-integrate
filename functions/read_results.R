read_results <- function(path, skip = 0){
  #---- truth table results + synthetic datasets ----
  if(str_detect(path, "superpopulations")){
    dataset_name <- str_split(path, pattern = "superpopulations/")[[1]][2]
    dataset_name <- str_remove(dataset_name, pattern = ".csv")
    return(data.table::fread(path, fill = TRUE, skip = skip) %>% 
             mutate("dataset_name" = dataset_name))
  }
  
  #---- simulation study results ----
  if(str_detect(path, "simulation_study/results")){
    # dataset_name <- str_split(path, pattern = "results/")[[1]][2]
    # dataset_name <- str_remove(dataset_name, pattern = ".csv")
    
    # if(str_detect(dataset_name, "ADAMS_HCAP")){
    #   data <- data.table::fread(path, fill = TRUE, skip = skip) %>% 
    #     na.omit() %>% mutate("dataset_name" = dataset_name) %>% 
    #     separate(dataset_name, into = c("name", "calibration"), 
    #              sep = "ADAMS_", remove = FALSE) %>%
    #     separate(name, 
    #              into = c("Distribution", "sample_size", "prior_props"), 
    #              sep = "_", remove = FALSE) %>% 
    #     mutate("prior_props" = "ADAMS") %>%
    #     mutate_at("Distribution", str_to_title) %>% 
    #     mutate("color" = case_when(Distribution == "Normal" ~ "#135467", 
    #                                Distribution == "Lognormal" ~ "#f0824f", 
    #                                Distribution == "Bathtub" ~ "b51661")) %>% 
    #     dplyr::select(-one_of("name")) %>%
    #     mutate_at(vars(-one_of("time", "dataset_name", "Distribution", 
    #                            "prior_props", "calibration", "color")), 
    #               as.numeric) 
    # } else{
    
    data <- data.table::fread(path, fill = TRUE, skip = skip) %>% 
      na.omit("time")
    
    if(str_detect(unique(data$dataset_name), "ADAMS_prior")){
      data %<>% mutate("dataset_name" = paste0(dataset_name, "_NA_NA"))
    } 
    
    if(str_detect(unique(data$dataset_name), "SRS") & 
       !str_detect(unique(data$dataset_name), "race")){
      data %<>% mutate("dataset_name" = paste0(dataset_name, "_NA"))
    }
    
    data %<>% 
      separate(dataset_name, 
               into = c("HRS_text", "sample_size", "sample_text", "HCAP_prop", 
                        "calibration1", "calibration2", "calibration_sampling", 
                        "sampling_strata"), 
               sep = "_", remove = FALSE) %>% 
      dplyr::select(-one_of("HRS_text", "sample_text")) %>%
      unite("calibration", c("calibration1", "calibration2")) %>%
      mutate_at(vars(-one_of("time", "dataset_name", "calibration", 
                             "calibration_sampling", "sampling_strata")), 
                as.numeric)
    #}
    
    return(data)
  }
  
  #---- extra_analyses ----
  if(str_detect(path, "HRS_synthetic")){
    dataset_name <- str_split(path, pattern = "HRS_synthetic_")[[1]][2]
    dataset_name <- str_remove(dataset_name, pattern = ".csv")
    return(data.table::fread(path, fill = TRUE, skip = skip) %>% 
             mutate("dataset_name" = dataset_name))
  }
  
  #---- OLD ----
  # #---- simulation study results ----
  # if(str_detect(path, "HCAP_")){
  #   dataset_name <- str_split(path, pattern = "HCAP_")[[1]][2]
  #   dataset_name <- str_split(dataset_name, pattern = "/")[[1]][1]
  #   
  #   sim_name <- str_split(path, pattern = "run_")[[1]][2]
  #   sim_name <- str_split(sim_name, pattern = "/")[[1]][1]
  #   
  #   dataset_number <- str_split(path, pattern = "run_")[[1]][2]
  #   dataset_number <- str_split(dataset_number, pattern = "/")[[1]][2]
  #   dataset_number <- str_remove(dataset_number, pattern = ".csv")
  #   
  #   return(data.table::fread(path, fill = TRUE) %>% 
  #     mutate("dataset_name" = dataset_name, 
  #            "sim_name" = sim_name, 
  #            "dataset_number" = dataset_number))
  # }
}
