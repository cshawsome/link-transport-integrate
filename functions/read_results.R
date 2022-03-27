read_results <- function(path){
  #---- truth table results ----
  if(str_detect(path, "HRS_synthetic")){
    dataset_name <- str_split(path, pattern = "HRS_synthetic_")[[1]][2]
    dataset_name <- str_remove(dataset_name, pattern = ".csv")
  }
  
  #---- simulation study results ----
  if(str_detect(path, "HCAP_")){
    dataset_name <- str_split(path, pattern = "HCAP_")[[1]][2]
    dataset_name <- str_split(dataset_name, pattern = "/")[[1]][1]
  }
  
  data.table::fread(path, fill = TRUE) %>% mutate("dataset_name" = dataset_name)
}
