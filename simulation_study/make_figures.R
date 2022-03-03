#---- Figure X: comparing ADAMS with synthetic HRS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- ****ADAMS imputed data ----
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned"))

#---- ****synthetic HRS ----
synthetic_normal_1000 <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/", "synthetic_data/", 
                  "synthetic_normal_1000.csv"))
