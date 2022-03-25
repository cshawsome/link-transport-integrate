read_results <- function(path){
  dataset_name <- str_split(path, pattern = "HRS_synthetic_")[[1]][2]
  dataset_name <- str_remove(dataset_name, pattern = ".csv")
  data.table::fread(path, fill = TRUE) %>% mutate("dataset_name" = dataset_name)
}
