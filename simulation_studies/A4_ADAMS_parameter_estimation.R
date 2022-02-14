#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned")) 

#---- read in variable selection results ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv")) %>% 
  dplyr::select("Variable") %>% unlist()

#---- cell IDs ----
#Merge Black, Hispanic, and Stroke (ever/never)
# Ex: 001 is a white participant with a history of stroke
sample <- ADAMS_imputed_clean[[1]]
cell_IDs <- 
  as.data.frame(table(sample$Black, sample$Hispanic, sample$Astroke)) %>%
  set_colnames(c("Black", "Hispanic", "Stroke", "Count")) %>% 
  filter(!(Black == 1 & Hispanic == 1)) %>% 
  unite("cell_ID", c("Black", "Hispanic", "Stroke"), sep = "") %>% 
  dplyr::select("cell_ID") %>% unlist()

#add column to all MI datasets
ADAMS_imputed_clean %<>% 
  lapply(., function(x) x %<>% 
           unite("cell_ID", c("Black", "Hispanic", "Astroke"), sep = "", 
                 remove = FALSE))

# #Sanity check
# colnames(ADAMS_imputed_clean[[1]])
# head(ADAMS_imputed_clean[[1]]$cell_ID)
# colnames(ADAMS_imputed_clean[[25]])
# head(ADAMS_imputed_clean[[25]]$cell_ID)

#---- parameter estimation ----
continuous_vars <- selected_vars[str_detect(selected_vars, "_Z")] %>% 
  str_remove_all(., "_Z")

normal_parameter_list <- list() 

#---- **Normal distribution ----
for(cell in cell_IDs){
  filtered_data <- 
    lapply(ADAMS_imputed_clean, function(x) x %>% filter(cell_ID == cell))
  
  #counts are going to differ across datasets because stroke is imputed
  min_count <- lapply(filtered_data, nrow) %>% unlist() %>% min()
  filtered_data <- 
    lapply(filtered_data, function(x) 
      sample_n(x, size = min_count, replace = FALSE) %>% 
        dplyr::select(all_of(continuous_vars)))
  
  #---- **mean matrix ----
  M = Reduce("+", filtered_data)/length(filtered_data)
  
  #---- **row covariance ----
  #independent because independent observations
  U = diag(1, nrow = min_count, ncol = min_count)
  
  #---- **column covariance ----
  #center matrices and take (matrix)^T(matrix)
  centered_data <- lapply(filtered_data, function(x) x - M)
  prod_data <- lapply(centered_data, function(x) 
    t(as.matrix(x)) %*% as.matrix(x))
  V = Reduce("+", prod_data)/(min_count*length(prod_data))
  
  #---- **save values ----
  normal_parameter_list[[cell]] <- list("M" = M, "U" = U, "V" = V)
}
