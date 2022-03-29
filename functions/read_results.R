read_results <- function(path){
  #---- truth table results ----
  if(str_detect(path, "HRS_synthetic")){
    dataset_name <- str_split(path, pattern = "HRS_synthetic_")[[1]][2]
    dataset_name <- str_remove(dataset_name, pattern = ".csv")
    return(data.table::fread(path, fill = TRUE) %>% 
      mutate("dataset_name" = dataset_name))
  }
  
  #---- simulation study results ----
  if(str_detect(path, "HCAP_")){
    dataset_name <- str_split(path, pattern = "HCAP_")[[1]][2]
    dataset_name <- str_split(dataset_name, pattern = "/")[[1]][1]
    
    sim_name <- str_split(path, pattern = "run_")[[1]][2]
    sim_name <- str_split(sim_name, pattern = "/")[[1]][1]
    
    dataset_number <- str_split(path, pattern = "run_")[[1]][2]
    dataset_number <- str_split(dataset_number, pattern = "/")[[1]][2]
    dataset_number <- str_remove(dataset_number, pattern = ".csv")
    
    return(data.table::fread(path, fill = TRUE) %>% 
      mutate("dataset_name" = dataset_name, 
             "sim_name" = sim_name, 
             "dataset_number" = dataset_number))
  }
}
