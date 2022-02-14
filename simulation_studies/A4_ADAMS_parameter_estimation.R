#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

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


