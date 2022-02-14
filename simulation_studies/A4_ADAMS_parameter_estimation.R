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
continuous_vars <- c("AAGE", "EDYRS", "ANMSETOT_norm", "ANBNTTOT", "ANIMMCR", 
                     "ANDELCOR", "ANSER7T", "ANAFTOT", "")

#---- **Normal distribution ----

