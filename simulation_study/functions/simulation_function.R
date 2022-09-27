simulation_function <- 
  function(warm_up, starting_props, categorical_vars, continuous_vars, id_var, 
           variable_labels, scenario, superpopulation, orig_means, orig_sds, 
           all_scenarios_list, cell_ID_key, color_palette, num_synthetic, 
           contrasts_matrix, kappa_0_mat, nu_0_mat, truth, seed, 
           path_to_raw_prior_sample, path_to_data, path_to_results){
    
    #---- pre-allocated results ----
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
        paste0("true_dem_prev_", c("white", "black", "hispanic")),
        paste0("mean_dem_prev_", c("white", "black", "hispanic")),
        paste0("LCI_dem_prev_", c("white", "black", "hispanic")),
        paste0("UCI_dem_prev_", c("white", "black", "hispanic")),
        paste0("dem_prev_coverage_", c("white", "black", "hispanic")),
        paste0("true_PR_", c("black", "hispanic")),
        paste0("mean_PR_", c("black", "hispanic")),
        paste0("LCI_PR_", c("black", "hispanic")),
        paste0("UCI_PR_", c("black", "hispanic")),
        paste0("PR_coverage_", c("black", "hispanic")),
        "time", "seed", "dataset_name")
    
    results <- matrix(ncol = length(result_names), nrow = 1) %>% 
      set_colnames(all_of(result_names))
    
    #---- scenario name ----
    all_sim_scenarios %<>% mutate("sample_text" = "sample")
    
    scenario_name <- 
      as.character(unite(
        all_sim_scenarios[scenario, c("sample_size", "sample_text", "HCAP_prop", 
                                      "calibration")], 
        col = "name", sep = "_"))
    
    #---- start time ----
    start <- Sys.time()
    
    #---- create synthetic HRS ----
    synthetic_HRS <- superpop %>% 
      slice_sample(n = unlist(all_sim_scenarios[scenario, "sample_size"]), 
                   replace = FALSE)
    
    #---- create synthetic HCAP ----
    dataset_to_copy <- synthetic_HRS %>%
      group_by(married_partnered) %>% 
      slice_sample(prop = unlist(all_sim_scenarios[scenario, "HCAP_prop"])/100, 
                   replace = FALSE) %>% 
      mutate("(Intercept)" = 1) %>% ungroup() %>% 
      mutate("dataset_name" = paste0("HRS_", scenario_name))
    
    #---- **calibration status ----
    calibration_status <- 
      ifelse(all_sim_scenarios[scenario, "calibration"] == "no_calibration", 
             FALSE, TRUE)
    
    #---- **true impairment class counts ----
    if(calibration_status){
      #---- ****flag calibration subsample ----
      calibration_sample_name <- 
        as.character(all_sim_scenarios[scenario, "calibration"])
      
      calibration_prop <- readr::parse_number(calibration_sample_name)/100
      
      dataset_to_copy[sample(seq(1, nrow(dataset_to_copy)), 
                             size = calibration_prop*nrow(dataset_to_copy), 
                             replace = FALSE), 
                      paste0("calibration_", calibration_prop*100)] <- 1
      
      #set unselected to 0
      not_selected <- 
        which(is.na(
          dataset_to_copy[, paste0("calibration_", calibration_prop*100)]))
      
      dataset_to_copy[not_selected, 
                      paste0("calibration_", calibration_prop*100)] <- 0
      
      results[, 
              c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other")] <- 
        colSums(dataset_to_copy[not_selected, 
                                c("Unimpaired", "MCI", "Dementia", "Other")])
      
      #---- ****stratified ----
      for(race in c("White", "black", "hispanic")){
        subset <- dataset_to_copy %>% 
          filter(!!sym(race) == 1 & 
                   !!sym(paste0("calibration_", calibration_prop*100)) == 0)
        
        results[, paste0(c("true_Unimpaired", "true_MCI", "true_Dementia", 
                           "true_Other"), "_", tolower(race))] <- 
          colSums(subset[, c("Unimpaired", "MCI", "Dementia", "Other")])
      }
    } else{
      results[, 
              c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other")] <- 
        colSums(dataset_to_copy[, c("Unimpaired", "MCI", "Dementia", "Other")])
      
      #---- ****stratified ----
      for(race in c("White", "black", "hispanic")){
        results[, paste0(c("true_Unimpaired", "true_MCI", "true_Dementia", 
                           "true_Other"), "_", tolower(race))] <- 
          colSums(dataset_to_copy[which(
            dataset_to_copy[, race] == 1), 
            c("Unimpaired", "MCI", "Dementia", "Other")])
      }
    }
    
    #---- generate synthetic data ----
    synthetic_HCAP <- 
      generate_synthetic(warm_up, run_number = NA, 
                         starting_props = starting_props, 
                         dataset_to_copy = dataset_to_copy, 
                         orig_means = orig_means, orig_sds = orig_sds,
                         calibration_sample = calibration_status, 
                         calibration_prop = calibration_prop, 
                         calibration_sample_name = calibration_sample_name,
                         path_to_raw_prior_sample = path_to_raw_prior_sample, 
                         path_to_data = path_to_data, 
                         path_to_analyses_folder = NA, 
                         path_to_figures_folder = NA, 
                         categorical_vars = categorical_vars, 
                         continuous_vars = continuous_vars, id_var = id_var, 
                         variable_labels = variable_labels, 
                         cell_ID_key = cell_ID_key,
                         color_palette = color_palette, 
                         contrasts_matrix = contrasts_matrix, 
                         kappa_0_mat = kappa_0_mat, 
                         nu_0_mat = nu_0_mat, num_synthetic = num_synthetic, 
                         data_only = TRUE)
    
    #---- function to clean missing counts ----
    all_classes <- c("Unimpaired", "MCI", "Dementia", "Other")
    
    clean_counts <- function(counts, all_classes){
      missing_class <- setdiff(all_classes, names(counts))
      
      if(length(missing_class) > 0){counts[missing_class] <- 0}
      
      return(counts[all_classes])
    } 
    
    #---- synthetic impairment class counts ----
    counts <- 
      lapply(synthetic_HCAP, function(x) table(x[, "Group"])) %>%
      lapply(., function(x) clean_counts(x, all_classes)) %>% 
      do.call(rbind, .)
    
    #---- **mean counts ----
    mean_counts <- round(colMeans(counts))
    results[, paste0("mean_", names(mean_counts))] <- mean_counts
    
    #---- **CI ----
    results[, paste0("LCI_", colnames(counts))] <- 
      apply(counts, 2, function(x) quantile(x, 0.025))
    
    results[, paste0("UCI_", colnames(counts))] <- 
      apply(counts, 2, function(x) quantile(x, 0.975))
    
    #---- **bias ----
    results[, paste0("bias_", colnames(counts))] <- 
      mean_counts - results[, paste0("true_", colnames(counts))]
    
    #---- **coverage ----
    for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
      results[, paste0(class, "_coverage")] <- 
        (results[, paste0("true_", class)] >= results[, paste0("LCI_", class)])*
        (results[, paste0("true_", class)] <= results[, paste0("UCI_", class)])
    }
    
    #---- stratified class counts ----
    for(race in c("White", "black", "hispanic")){
      strat_datasets <- 
        lapply(synthetic_HCAP, function(x) x %>% filter(!!as.symbol(race) == 1))
      
      strat_counts <- 
        lapply(strat_datasets, function(x) table(x[, "Group"])) %>% 
        lapply(., function(x) clean_counts(x, all_classes)) %>% 
        do.call(rbind, .)
      
      #---- **mean counts ----
      strat_mean_counts <- round(colMeans(strat_counts))
      results[, paste0("mean_", names(strat_mean_counts), "_", tolower(race))] <- 
        strat_mean_counts
      
      #---- **CI ----
      results[, paste0("LCI_", colnames(strat_counts), "_", tolower(race))] <- 
        apply(strat_counts, 2, function(x) quantile(x, 0.025))
      
      results[, paste0("UCI_", colnames(strat_counts), "_", tolower(race))] <- 
        apply(counts, 2, function(x) quantile(x, 0.975))
      
      #---- **bias ----
      results[, paste0("bias_", colnames(strat_counts), "_", tolower(race))] <- 
        strat_mean_counts - 
        results[, paste0("true_", colnames(strat_counts), "_", tolower(race))]
      
      #---- **coverage ----
      for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
        results[, paste0(class, "_coverage_", tolower(race))] <- 
          (results[, paste0("true_", class, "_", tolower(race))] >= 
             results[, paste0("LCI_", class, "_", tolower(race))])*
          (results[, paste0("true_", class, "_", tolower(race))] <= 
             results[, paste0("UCI_", class, "_", tolower(race))])
      }
    }
    
    #---- standardized dem prevalences ----
    dem_estimates <- lapply(synthetic_HCAP, standardized_dem_estimates, 
                            standard_data = synthetic_HRS) %>% do.call(rbind, .)
    
    #---- **true dem prevs ----
    results[, paste0("true_dem_prev_", c("white", "black", "hispanic"))] <- 
      colSums(truth[, paste0(c("expected_white", "expected_black", 
                               "expected_hispanic"), "_dem_count")])/
      nrow(superpop)
    
    #---- **mean dem prevs ----
    results[, paste0("mean_dem_prev_", c("white", "black", "hispanic"))] <- 
      colMeans(dem_estimates[, paste0(c("white", "black", "hispanic"), "_risk")])
    
    #---- **prev CI ----
    results[, paste0("LCI_dem_prev_", c("white", "black", "hispanic"))] <- 
      apply(dem_estimates[, paste0(c("white", "black", "hispanic"), "_risk")], 2, 
            function(x) quantile(x, 0.025))
    
    results[, paste0("UCI_dem_prev_", c("white", "black", "hispanic"))] <- 
      apply(dem_estimates[, paste0(c("white", "black", "hispanic"), "_risk")], 2, 
            function(x) quantile(x, 0.975))
    
    #---- **prev coverage ----
    for(race in c("white", "black", "hispanic")){
      results[, paste0("dem_prev_coverage_", race)] <- 
        (results[, paste0("true_dem_prev_", race)] >= 
           results[, paste0("LCI_dem_prev_", race)])*
        (results[, paste0("true_dem_prev_", race)] <= 
           results[, paste0("UCI_dem_prev_", race)])
    }
    
    #---- **true PR ----
    results[, paste0("true_PR_", c("black", "hispanic"))] <- 
      c(results[, "true_dem_prev_black"]/results[, "true_dem_prev_white"], 
        results[, "true_dem_prev_hispanic"]/results[, "true_dem_prev_white"])
    
    #---- **mean PR ----
    results[, paste0("mean_PR_", c("black", "hispanic"))] <- 
      colMeans(dem_estimates[, paste0("PR_", c("black", "hispanic"))])
    
    #---- **RR CI ----
    results[, paste0("LCI_PR_", c("black", "hispanic"))] <- 
      apply(dem_estimates[, paste0("PR_", c("black", "hispanic"))], 2, 
            function(x) quantile(x, 0.025))
    
    results[, paste0("UCI_PR_", c("black", "hispanic"))] <- 
      apply(dem_estimates[, paste0("PR_", c("black", "hispanic"))], 2, 
            function(x) quantile(x, 0.975))
    
    #---- **PR coverage ----
    for(race in c("black", "hispanic")){
      results[, paste0("PR_coverage_", race)] <- 
        (results[, paste0("true_PR_", race)] >= 
           results[, paste0("LCI_PR_", race)])*
        (results[, paste0("true_PR_", race)] <= 
           results[, paste0("UCI_PR_", race)])
    }
    
    # #---- models ----
    # #---- **add weights ----
    # #the only variable that contributed to selection was marital status
    # # we selected half of married people and half of unmarried people, thus the 
    # # weight for every observation should be 2
    # 
    # synthetic_HCAP %<>% lapply(function(x) x %<>% mutate("weight" = 2))
    # 
    # #---- **predicted dementia status ----
    # synthetic_HCAP %<>% 
    #   lapply(function(x) x %<>% 
    #            mutate("predicted_Dementia" = ifelse(Group == "Dementia", 1, 0)))
    # 
    # #---- **run models ----
    # models <- 
    #   lapply(synthetic_HCAP, 
    #          function(dataset) glm(predicted_Dementia ~ 
    #                                  age_Z + female + black + hispanic, 
    #                                family = "poisson", data = dataset)) 
    # 
    # #---- **pooled ----
    # pooled_model <- summary(mice::pool(models))
    # 
    # #---- ****truth table ----
    # truth_table <- 
    #   truth[which(truth$dataset_name == superpop_name), ] 
    # 
    # for(race_eth in c("black", "hispanic")){
    #   #---- ****estimates ----
    #   results[, paste0(race_eth, "_beta")] <- 
    #     pooled_model[which(pooled_model$term == race_eth), "estimate"]
    #   
    #   results[, paste0(race_eth, "_se")] <- 
    #     pooled_model[which(pooled_model$term == race_eth), "std.error"]
    #   
    #   results[, paste0(race_eth, "_LCI")] <- 
    #     results[, paste0(race_eth, "_beta")] - 
    #     1.96*results[, paste0(race_eth, "_se")]
    #   
    #   results[, paste0(race_eth, "_UCI")] <- 
    #     results[, paste0(race_eth, "_beta")] + 
    #     1.96*results[, paste0(race_eth, "_se")]
    #   
    #   #---- ****coverage ----
    #   results[, paste0(race_eth, "_coverage")] <- 
    #     (truth_table[which(truth_table$term == race_eth), "estimate"] >= 
    #        results[, paste0(race_eth, "_LCI")])*
    #     (truth_table[which(truth_table$term == race_eth), "estimate"] <= 
    #        results[, paste0(race_eth, "_UCI")])
    # }
    # 
    #---- end time ----
    results[, "time"] <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    
    #---- seed ----
    results[, "seed"] <- seed
    
    #---- dataset name ----
    results[, "dataset_name"] <- paste0("HRS_", scenario_name)
    
    #---- write results ----
    file_path <- 
      paste0(path_to_results, scenario_name, ".csv")
    
    if(file.exists(file_path)){
      write_csv(as.data.frame(results), file = file_path, append = TRUE)
    } else{
      write_csv(as.data.frame(results), file = file_path, col_names = TRUE)
    }
  }

#---- test function ----
library("tidyverse")
library("DirichletReg")
library("magrittr")
library("MCMCpack")
library("locfit")
library("vroom")
library("mvnfast")
library("mice")
library("LaplacesDemon")

path_to_RScripts <- here::here("simulation_study", "functions", "/")
source(here::here("functions", "read_results.R"))
source(paste0(path_to_RScripts, "generate_synthetic_function.R"))
source(paste0(path_to_RScripts, "standardized_dem_estimates.R"))

path_to_data <- paste0("/Users/crystalshaw/Library/CloudStorage/Box-Box/",
                       "Dissertation/data/")
superpop <-
  read_results(paste0(path_to_data, "superpopulations/superpop_1000000.csv"))
truth <- read_csv(paste0(path_to_data,
                         "superpopulations/agesex_standardized_prevs.csv"))
variable_labels <-
  read_csv(paste0(path_to_data, "variable_crosswalk.csv"))
cell_ID_key <- read_csv(paste0(path_to_data, "cell_ID_key.csv")) %>%
  mutate_all(as.character)
color_palette <- read_csv(paste0(path_to_data, "color_palette.csv"))
all_sim_scenarios <- read_csv(paste0(path_to_data, "sim_study_scenarios.csv"))

warm_up = 100
starting_props = rep(0.25, 4)
categorical_vars = W = c("black", "hispanic", "stroke")
continuous_vars = Z = colnames(superpop)[str_detect(colnames(superpop), "_Z")]
id_var = "HHIDPN"
scenario = scenario_num = 404 #calibration 50 SRS sample size 2000, HCAP prop 25
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
superpopulation <- superpop
orig_means = means <-
  read_csv(paste0(path_to_box, "data/superpopulations/superpop_means.csv"))
orig_sds = sds <-
  read_csv(paste0(path_to_box, "data/superpopulations/superpop_sds.csv"))

all_scenarios_list = all_sim_scenarios

num_synthetic = 1000
nu_0_mat <- read_csv(paste0(path_to_data, "tuning/nu_0_matrix.csv"))
kappa_0_mat <- read_csv(paste0(path_to_data, "tuning/kappa_0_matrix.csv"))
contrasts_matrix = A =
  read_csv(paste0(path_to_data, "contrasts_matrix.csv")) %>% as.matrix()
path_to_results <- paste0(path_to_box, "analyses/simulation_study/results/")
seed = 1

set.seed(20220512)

replicate(2,
          simulation_function(warm_up = 100, starting_props = rep(0.25, 4),
                              categorical_vars = W, continuous_vars = Z,
                              id_var = "HHIDPN",
                              variable_labels = variable_labels,
                              scenario = scenario_num,
                              superpopulation = superpop, orig_means = means,
                              orig_sds = sds,
                              all_scenarios_list = all_sim_scenarios,
                              cell_ID_key = cell_ID_key,
                              color_palette = color_palette,
                              num_synthetic = 1000, contrasts_matrix = A,
                              kappa_0_mat = kappa_0_mat, nu_0_mat = nu_0_mat,
                              truth = truth, seed = seed,
                              path_to_raw_prior_sample =
                                paste0(path_to_data,
                                       "prior_data/MI/MI_datasets_cleaned"),
                              path_to_data = path_to_data,
                              path_to_results =
                                paste0(path_to_box,
                                       "analyses/simulation_study/results/")))

