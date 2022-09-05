#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HRS analytic dataset ----
HRS <- read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

# #Sanity check: only imputed memory scores should have missingness
# colMeans(is.na(HRS_analytic))[which(colMeans(is.na(HRS_analytic)) > 0)]

#---- **HCAP analytic dataset ----
HCAP <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv"))

#---- **imputation matrix ----
hotdeck_vars_mat <- 
  read_csv(paste0(path_to_box, 
                  "data/superpopulations/hotdeck_impute_mat.csv")) %>% 
  column_to_rownames("var_names")

# #---- **fixed betas ----
# fixed_betas <- 
#   read_csv(paste0(path_to_box, "data/variable_selection/", 
#                   "fixed_model_coefficients.csv"))

#---- source functions ----
source(here("simulation_study", "functions", "hotdeck_function.R"))

#---- synthetic superpopulation ----
set.seed(20220905)

#About XX hours for superpop
start <- Sys.time()
superpop_size <- 1000000
superpop <- sample_n(HRS_analytic, size = superpop_size, replace = TRUE) %>% 
  mutate("HHIDPN" = seq(1, superpop_size))

#add columns for neuropsych
superpop[, rownames(hotdeck_vars_mat)] <- NA

superpop %<>% 
  hotdeck(dataset_to_impute = ., hotdeck_dataset = HCAP, 
          imputation_mat = hotdeck_vars_mat, binary_vars = NA)

end <- Sys.time() - start

#Sanity check
colMeans(is.na(superpop))[which(colMeans(is.na(superpop)) > 0)]

#---- **impairment classes ----

#---- **QC superpop ----


#---- synthetic HRS ----
#---- **source functions ----
source(here::here("functions", "read_results.R"))

#---- **read in superpop data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **data paths ----
superpop_data_paths <- 
  list.files(path = paste0(path_to_box, 
                           "analyses/simulation_study/superpopulations"), 
             full.names = TRUE, pattern = "*.csv")

superpop_data_list <- lapply(superpop_data_paths, read_results)

#---- **create one set of synthetic HRS ----
create_HRS_datasets <- function(superpop, n){
  sample_n(superpop, size = n) %>% 
    separate(dataset_name, sep = "_", into = c("dist", "size", "prior")) %>% 
    mutate_at("size", as.numeric) %>% mutate(size = n) %>% 
    unite("dataset_name", c("dist", "size", "prior"), sep = "_")
}

set.seed(20220507)

for(n in c(500, 1000, 2000, 4000, 8000)){
  if(!exists("synthetic_HRS_list")){
    synthetic_HRS_list <- 
      lapply(superpop_data_list, function(x) create_HRS_datasets(x, n))
  } else{
    synthetic_HRS_list <- 
      append(synthetic_HRS_list, lapply(superpop_data_list, function(x) 
        create_HRS_datasets(x, n)))
  }
}

#---- **save data ----
saveRDS(synthetic_HRS_list, 
        file = paste0(path_to_box, 
                      "analyses/simulation_study/synthetic_HRS_list"))

#---- synthetic HCAP ----
#---- **read in synthetic HRS data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_HRS_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HRS_list"))

#---- **create one set of synthetic HCAP ----
set.seed(20220507)

synthetic_HCAP_list <- 
  lapply(synthetic_HRS_list, 
         function(x) 
           x %>% group_by(married_partnered) %>% slice_sample(prop = 0.5) %>% 
           mutate("(Intercept)" = 1, 
                  "calibration_50" = 0) %>% ungroup())

#---- **flag calibration subsample ----
for(i in 1:length(synthetic_HCAP_list)){
  synthetic_HCAP_list[[i]][sample(seq(1, nrow(synthetic_HCAP_list[[i]])), 
                                  size = 0.5*nrow(synthetic_HCAP_list[[i]]), 
                                  replace = FALSE), "calibration_50"] <- 1
}

#---- **save data ----
saveRDS(synthetic_HCAP_list, 
        file = paste0(path_to_box, 
                      "analyses/simulation_study/synthetic_HCAP_list"))

#---- summary stats ----
#---- **read in synthetic HCAP data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

synthetic_HCAP_list <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/synthetic_HCAP_list"))

#---- **filter to normal list ----
dataset_names <- 
  unlist(lapply(synthetic_HCAP_list, function(x) unique(x$dataset_name)))

indices <- which(dataset_names %in% 
                   paste0("normal_", c(500, 1000, 2000, 4000, 8000), "_ADAMS"))

synthetic_HCAP_list <- synthetic_HCAP_list[indices]

#---- **summarize race/ethnicity x dementia ---- 
test <- synthetic_HCAP_list[[5]]

table(test$White, test$black, test$hispanic, test$Dementia) %>% 
  as.data.frame %>% filter(!Freq == 0) %>% 
  set_colnames(c("White", "Black", "Hispanic", "Dementia", "Freq"))