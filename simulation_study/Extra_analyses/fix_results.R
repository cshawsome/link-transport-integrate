#---- source functions ----
source(here::here("functions", "read_results.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

#---- **format headers ----
result_names <- 
  c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other",
    paste0(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
           "_white"),
    paste0(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
           "_black"),
    paste0(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
           "_hispanic"),
    "mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other",
    paste0(c("mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other"), 
           "_white"), 
    paste0(c("mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other"), 
           "_black"), 
    paste0(c("mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other"), 
           "_hispanic"),
    "bias_Unimpaired", "bias_MCI", "bias_Dementia", "bias_Other",
    paste0(c("bias_Unimpaired", "bias_MCI", "bias_Dementia", "bias_Other"), 
           "_white"), 
    paste0(c("bias_Unimpaired", "bias_MCI", "bias_Dementia", "bias_Other"), 
           "_black"),
    paste0(c("bias_Unimpaired", "bias_MCI", "bias_Dementia", "bias_Other"), 
           "_hispanic"), 
    "LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other",
    paste0(c("LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other"), 
           "_white"),
    paste0(c("LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other"), 
           "_black"),
    paste0(c("LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other"), 
           "_hispanic"),
    "UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other",
    paste0(c("UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other"), 
           "_white"),
    paste0(c("UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other"), 
           "_black"),
    paste0(c("UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other"), 
           "_hispanic"),
    "Unimpaired_coverage", "MCI_coverage", "Dementia_coverage", 
    "Other_coverage", 
    paste0(c("Unimpaired_coverage", "MCI_coverage", "Dementia_coverage", 
             "Other_coverage"), "_white"),
    paste0(c("Unimpaired_coverage", "MCI_coverage", "Dementia_coverage", 
             "Other_coverage"), "_black"),
    paste0(c("Unimpaired_coverage", "MCI_coverage", "Dementia_coverage", 
             "Other_coverage"), "_hispanic"),
    "black_beta", "black_se", "black_LCI", "black_UCI",
    "hispanic_beta", "hispanic_se", "hispanic_LCI", "hispanic_UCI",
    "black_coverage", "hispanic_coverage", "time", "seed", "dataset_name")

#---- fix files ----
for(path in results_paths[5]){
  read_results(path, skip = 1001) %>%
    set_colnames(c(all_of(result_names),
                   "Distribution", "sample_size", "prior_props", "color")) %>%
    write_csv(path)
}



