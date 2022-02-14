#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS_analytic <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv")) %>%
  #complete-case for now
  na.omit()

#---- **continuous distribution parameters ----
normal_parameter_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/", 
                 "continuous_distribution_parameters/", 
                 "normal_parameters"))

#---- cell IDs ----
#Merge Black, Hispanic, and Stroke (ever/never)
# Ex: 001 is a white participant with a history of stroke
cell_IDs <- 
  as.data.frame(table(HRS_analytic$Black, HRS_analytic$Hispanic, 
                      HRS_analytic$r13stroke)) %>%
  set_colnames(c("Black", "Hispanic", "Stroke", "Count")) %>% 
  filter(!(Black == 1 & Hispanic == 1)) %>% 
  unite("cell_ID", c("Black", "Hispanic", "Stroke"), sep = "") %>% 
  dplyr::select("cell_ID") %>% unlist()

#add column to dataset
HRS_analytic %<>% 
  unite("cell_ID", c("Black", "Hispanic", "r13stroke"), sep = "", remove = FALSE)

#Sanity check
head(HRS_analytic$cell_ID)
