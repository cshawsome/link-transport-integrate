#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "vroom", "locfit")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **HRS analytic dataset ----
HRS <- read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

# #Sanity check: only imputed memory scores should have missingness
# colMeans(is.na(HRS))[which(colMeans(is.na(HRS)) > 0)]

#---- **HCAP analytic dataset ----
HCAP <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv"))

#---- **imputation matrix ----
hotdeck_vars_mat <- 
  read_csv(paste0(path_to_box, 
                  "data/superpopulations/hotdeck_impute_mat.csv")) %>% 
  column_to_rownames("var_names")

#---- **ADAMS variable selection results ----
selected_vars_mat <- 
  read_csv(paste0(path_to_box, 
                  "data/variable_selection/model_coefficients.csv"))

# #---- **fixed betas ----
# fixed_betas <- 
#   read_csv(paste0(path_to_box, "data/variable_selection/", 
#                   "fixed_model_coefficients.csv"))

#---- source functions ----
source(here("simulation_study", "functions", "hotdeck_function.R"))

#---- synthetic superpopulation ----
set.seed(20220905)

start <- Sys.time()
superpop_size <- 1000000
superpop <- sample_n(HRS, size = superpop_size, replace = TRUE) %>% 
  mutate("HHIDPN" = seq(1, superpop_size))

#add columns for neuropsych
superpop[, rownames(hotdeck_vars_mat)] <- NA

superpop %<>% 
  hotdeck(dataset_to_impute = ., hotdeck_dataset = HCAP, 
          imputation_mat = hotdeck_vars_mat, binary_vars = NA)

end <- Sys.time() - start

# #Sanity check
# colMeans(is.na(superpop))[colMeans(is.na(superpop)) > 0]

#---- **standardize continuous vars ----
standardize_vars <- str_remove(unique(selected_vars_mat$data_label)[
  str_detect(unique(selected_vars_mat$data_label), "_Z")], "_Z")

Z_score <- function(data, vars){
  subset <- data %>% dplyr::select(all_of(vars)) %>% 
    mutate_all(scale) %>% mutate_all(as.numeric) %>%
    set_colnames(paste0(all_of(vars), "_Z"))
  
  data %<>% cbind(., subset)
  
  return(data)
}

superpop %<>% Z_score(., standardize_vars)

#---- **store superpop means and variances ----
orig_mean_sd <- superpop %>% 
  summarize_at(standardize_vars, .funs = c("mean" = mean, "sd" = sd))

means <- orig_mean_sd %>% dplyr::select(contains("mean")) %>% 
  set_colnames(standardize_vars)

sds <- orig_mean_sd %>% dplyr::select(contains("sd")) %>% 
  set_colnames(standardize_vars)

# #Sanity check
# test_subset <- head(superpop[, paste0(standardize_vars, "_Z")])
# test_subset*
#   matrix(rep(as.numeric(sds), nrow(test_subset)), nrow = nrow(test_subset),
#          byrow = TRUE) +
#   matrix(rep(as.numeric(means), nrow(test_subset)), nrow = nrow(test_subset),
#          byrow = TRUE)
# View(head(superpop[, standardize_vars]))

#save
write_csv(means, paste0(path_to_box, "data/superpopulations/superpop_means.csv"))
write_csv(sds, paste0(path_to_box, "data/superpopulations/superpop_sds.csv"))

#---- **impairment classes ----
#---- ****predict values in superpop ----
superpop %<>% mutate("Intercept" = 1)

for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  filtered_vars_mat <- selected_vars_mat[, c("data_label", class)] %>% 
    filter(!!sym(class) != 0)
  
  preds <- filtered_vars_mat$data_label
  beta <- unlist(filtered_vars_mat[, class])
  
  superpop[, paste0("p_", class)] <- 
    as.numeric(expit(as.matrix(superpop[, preds]) %*% as.matrix(beta)))
}

# #Sanity check
# View(head(superpop[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")]))

# #---- ****draw impairment class assignments ----
# superpop %<>%
#   mutate("Unimpaired" =
#            rbernoulli(nrow(superpop), p = superpop$p_Unimpaired)*1,
#          "Other" =
#            ifelse(Unimpaired == 0, rbernoulli(n = 1, p = p_Other)*1, 0),
#          "MCI" = ifelse(Unimpaired == 0 & Other == 0,
#                         rbernoulli(n = 1, p = p_MCI)*1 , 0),
#          "Dementia" = ifelse(Unimpaired == 0 & Other == 0 & MCI == 0, 1, 0)) %>%
#   mutate("num_classes" = Unimpaired + Other + MCI + Dementia)
# 
# #Sanity check
# table(superpop$num_classes)

# #---- **assign impairment class ----
# superpop[, "dem_class"] <- 
#   apply(superpop[, c("p_Unimpaired", "p_Other", "p_MCI", "p_Dementia")], 1, 
#         function(x) str_remove(names(which.max(x)), "p_"))
# 
# superpop %<>% 
#   mutate("Unimpaired" = ifelse(dem_class == "Unimpaired", 1, 0), 
#          "MCI" = ifelse(dem_class == "MCI", 1, 0), 
#          "Dementia" = ifelse(dem_class == "Dementia", 1, 0), 
#          "Other" = ifelse(dem_class == "Other", 1, 0)) %>% 
#   mutate("num_classes" = Unimpaired + MCI + Dementia + Other)
# 
# 
# #Sanity check
# table(superpop$num_classes)

# #---- **draw dem class by using a weighted vector ----
# superpop[, "dem_class"] <- 
#   apply(superpop[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")], 1, 
#         function(x) 
#           sample(c("Unimpaired", "MCI", "Dementia", "Other"), size = 1, 
#                  prob = x))
# 
# superpop %<>% 
#   mutate("Unimpaired" = ifelse(dem_class == "Unimpaired", 1, 0), 
#          "MCI" = ifelse(dem_class == "MCI", 1, 0), 
#          "Dementia" = ifelse(dem_class == "Dementia", 1, 0), 
#          "Other" = ifelse(dem_class == "Other", 1, 0)) %>% 
#   mutate("num_classes" = Unimpaired + MCI + Dementia + Other)
# 
# # #Sanity check
# # table(superpop$num_classes)

#---- **draw impaired categories only using weighted vector ----
superpop[, "dem_class"] <-
  apply(superpop[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")], 1,
        function(x) str_remove(names(which.max(x)), "p_"))

superpop[superpop$dem_class != "Unimpaired", "dem_class"] <- 
  apply(superpop[superpop$dem_class != "Unimpaired", 
                 c("p_MCI", "p_Dementia", "p_Other")], 1, 
        function(x) 
          sample(c("MCI", "Dementia", "Other"), size = 1, prob = x))

superpop %<>% 
  mutate("Unimpaired" = ifelse(dem_class == "Unimpaired", 1, 0), 
         "MCI" = ifelse(dem_class == "MCI", 1, 0), 
         "Dementia" = ifelse(dem_class == "Dementia", 1, 0), 
         "Other" = ifelse(dem_class == "Other", 1, 0)) %>% 
  mutate("num_classes" = Unimpaired + MCI + Dementia + Other)

# #Sanity check
# table(superpop$num_classes)

#---- **QC superpop ----
#---- ****overall summaries ----
#U: 37.0%, M: 16.1%, D: 26.4%, O: 20.4%
colMeans(superpop[, c("Unimpaired", "MCI", "Dementia", "Other")])

#More dementia among women? Yes (M: 25.2%, W: 27.2%)
mean(superpop[superpop$female == 1, "Dementia"])
mean(superpop[superpop$female == 0, "Dementia"])

#More dementia among racial/ethnic minorities? Kind of? 
# (w: 25.7%, b: 33.3%, h: 22.3%)
mean(superpop[superpop$White == 1, "Dementia"])
mean(superpop[superpop$black == 1, "Dementia"])
mean(superpop[superpop$hispanic == 1, "Dementia"])

#Age: increased risk by year (PR = 1.03)
exp(coefficients(glm(Dementia ~ age, data = superpop, family = "poisson")))

#Years of Education: decreased risk by higher education (PR = 0.98)
exp(coefficients(glm(Dementia ~ edyrs, data = superpop, family = "poisson")))

#Stroke: increased risk for yes vs. no (PR = 1.70)
exp(coefficients(glm(Dementia ~ stroke, data = superpop, family = "poisson")))

#Diabetes: no increased risk for yes vs. no (PR = 1.00)
exp(coefficients(glm(Dementia ~ diabe, data = superpop, family = "poisson")))

#Diabetes: increased risk for any impairment yes vs. no (PR = 1.10)
superpop %<>% mutate("any_impairment" = Dementia + MCI)
exp(coefficients(glm(any_impairment ~ diabe, data = superpop, family = "poisson")))

#---- ****age and sex-standardized estimates by race ----
#---- ******create age strata ----
superpop %<>% 
  mutate("age_cat" = cut(age, breaks = c(70, 75, 80, 85, 90, max(superpop$age)), 
                         include.lowest = TRUE, right = FALSE))

# #Sanity check
# table(superpop$age_cat, useNA = "ifany")

#---- ******standardization tables ----
agesex_totals <- 
  superpop %>% group_by(female, age_cat) %>% count() %>% arrange(female) %>% 
  set_colnames(c("female", "age", "superpop_count")) %>% 
  dplyr::select("superpop_count", everything())

for(race in c("White", "black", "hispanic")){
  assign(paste0(tolower(race), "_dem_risk_table"), 
         superpop %>% filter(!!sym(race) == 1) %>% 
           group_by(female, age_cat) %>% 
           summarise("dem_prop" = mean(Dementia)) %>% 
           arrange(female) %>% 
           set_colnames(c("female", "age", 
                          paste0(tolower(race), "_dem_risk"))))
}

# #Sanity check
# test_data <- superpop %>% filter(White == 1)
# table(test_data$female, test_data$age_cat, test_data$Dementia)

#merge tables
agesex_standardized <- 
  left_join(agesex_totals, white_dem_risk_table, by = c("female", "age")) %>% 
  left_join(., black_dem_risk_table, by = c("female", "age")) %>% 
  left_join(., hispanic_dem_risk_table, by = c("female", "age"))

#expected counts
for(race in c("white", "black", "hispanic")){
  agesex_standardized %<>% 
    mutate(!!paste0("expected_", race, "_dem_count") := 
             !!sym(paste0(race, "_dem_risk"))*superpop_count)
}

# #Sanity check
# agesex_standardized$white_dem_risk*agesex_standardized$superpop_count
# agesex_standardized$hispanic_dem_risk*agesex_standardized$superpop_count

#---- ******estimates ----
#0.254
white_risk <- 
  sum(agesex_standardized$expected_white_dem_count)/nrow(superpop)
#0.343
black_risk <- 
  sum(agesex_standardized$expected_black_dem_count)/nrow(superpop)
#0.225
hispanic_risk <- 
  sum(agesex_standardized$expected_hispanic_dem_count)/nrow(superpop)

#PR compared to white
#1.35
PR_black <- black_risk/white_risk
#0.88
PR_hispanic <- hispanic_risk/white_risk

#---- **save superpop data ----
write_csv(superpop, 
          paste0(path_to_box, "data/superpopulations/superpop_1000000.csv"))

#---- **save impairment class props ----
write_csv(colMeans(superpop[, c("Unimpaired", "MCI", "Dementia", "Other")]) %>% 
            as.matrix() %>% as.data.frame() %>% set_colnames("prop") %>% 
            rownames_to_column("Group"), 
          paste0(path_to_box, 
                 "data/superpopulations/impairment_class_props.csv"))

#---- **save age and sex-standardized table ----
write_csv(agesex_standardized, 
          paste0(path_to_box, 
                 "data/superpopulations/agesex_standardized_prevs.csv"))

#---- synthetic HRS ----
#create one set of synthetic HRS for tuning 
create_HRS_datasets <- function(superpop, n){
  sample_n(superpop, size = n) %>% 
    mutate("dataset_name" = paste0("HRS_", n), 
           "(Intercept)" = 1)
}

set.seed(20220507)

for(n in c(2000, 4000, 8000)){
  if(!exists("synthetic_HRS_list")){
    synthetic_HRS_list <- list(create_HRS_datasets(superpop, n))
  } else{
    synthetic_HRS_list <- 
      append(synthetic_HRS_list, list(create_HRS_datasets(superpop, n)))
  }
}

#---- **save data ----
saveRDS(synthetic_HRS_list, 
        file = paste0(path_to_box, "data/HRS/synthetic_HRS_list"))

#---- synthetic HCAP ----
#---- **create one set of synthetic HCAP ----
set.seed(20220507)

sample_props <- c(0.25, 0.50)

for(sample_prop in sample_props){
  if(exists("synthetic_HCAP_list")){
    synthetic_HCAP_list <- 
      append(synthetic_HCAP_list, 
             lapply(synthetic_HRS_list, 
                    function(x) 
                      x %>% group_by(married_partnered) %>% 
                      slice_sample(prop = sample_prop) %>% 
                      mutate("calibration_25_SRS" = 0, 
                             "calibration_50_SRS" = 0, 
                             "calibration_25_design" = 0, 
                             "calibration_50_design" = 0) %>% ungroup() %>% 
                      mutate("dataset_name" = 
                               paste0(dataset_name, "_sample_", 
                                      sample_prop*100))))
  } else{
    synthetic_HCAP_list <- 
      lapply(synthetic_HRS_list, 
             function(x) 
               x %>% group_by(married_partnered) %>% 
               slice_sample(prop = sample_prop) %>% 
               mutate("calibration_25_SRS" = 0, 
                      "calibration_50_SRS" = 0, 
                      "calibration_25_design" = 0, 
                      "calibration_50_design" = 0) %>% ungroup() %>% 
               mutate("dataset_name" = 
                        paste0(dataset_name, "_sample_", sample_prop*100)))
  }
}

#---- **flag calibration subsamples ----
#---- ****SRS calibration ----
for(i in 1:length(synthetic_HCAP_list)){
  synthetic_HCAP_list[[i]][sample(seq(1, nrow(synthetic_HCAP_list[[i]])), 
                                  size = 0.25*nrow(synthetic_HCAP_list[[i]]), 
                                  replace = FALSE), "calibration_25_SRS"] <- 1
  
  synthetic_HCAP_list[[i]][sample(seq(1, nrow(synthetic_HCAP_list[[i]])), 
                                  size = 0.50*nrow(synthetic_HCAP_list[[i]]), 
                                  replace = FALSE), "calibration_50_SRS"] <- 1
}

#---- ****design calibration ----
filtered_vars_mat <- selected_vars_mat[, c("data_label", "Dementia")] %>% 
  filter(Dementia != 0)

preds <- filtered_vars_mat$data_label
beta <- unlist(filtered_vars_mat[, "Dementia"])

for(i in 1:length(synthetic_HCAP_list)){
  synthetic_HCAP_list[[i]][, "p_dementia"] <- 
    as.numeric(expit(as.matrix(synthetic_HCAP_list[[i]][, preds]) %*% 
                       as.matrix(beta)))
  
  for(calibration_prop in c(0.25, 0.50)){
    #sample 20% group within each race/ethnicity based on p(dementia)
    impaired_sample_IDs <- synthetic_HCAP_list[[i]] %>% 
      dplyr::select("HHIDPN", "White", "black", "hispanic", "p_dementia") %>% 
      unite("race_eth_code", c("White", "black", "hispanic")) %>% 
      group_by(race_eth_code) %>% arrange(desc(p_dementia)) %>%
      slice_head(n = 0.20*calibration_prop*nrow(synthetic_HCAP_list[[i]])/3) %>% 
      ungroup() %>% dplyr::select("HHIDPN") %>% unlist()
    
    #sample 80% randomly within each race/ethnicity
    random_sample_IDs <- 
      tryCatch(synthetic_HCAP_list[[i]] %>% 
                 filter(!HHIDPN %in% impaired_sample_IDs) %>% 
                 dplyr::select("HHIDPN", "White", "black", "hispanic") %>% 
                 unite("race_eth_code", c("White", "black", "hispanic")) %>% 
                 group_by(race_eth_code) %>%
                 slice_sample(n = 0.80*calibration_prop*
                                nrow(synthetic_HCAP_list[[i]])/3) %>% 
                 ungroup() %>% dplyr::select("HHIDPN") %>% unlist(), 
               error = function(e) {
                 #how to handle lack of people in race/ethnic category
                 subset <- synthetic_HCAP_list[[i]] %>% 
                   filter(!HHIDPN %in% impaired_sample_IDs) %>% 
                   dplyr::select("HHIDPN", "White", "black", "hispanic") %>% 
                   unite("race_eth_code", c("White", "black", "hispanic"))
                 
                 counts <- subset %>% group_by(race_eth_code) %>% 
                   count(race_eth_code)
                 
                 #how many do we need to sample
                 num_to_sample <- 
                   0.80*calibration_prop*nrow(synthetic_HCAP_list[[i]])/3
                 
                 counts %<>% mutate("missing" = num_to_sample - n)
                 
                 #take everyone from groups that have fewer than the number we 
                 #  need
                 for(code in counts$race_eth_code){
                   if(counts[which(counts$race_eth_code == code), "missing"] > 0){
                     selected <- subset %>% filter(race_eth_code == code) %>% 
                       dplyr::select("HHIDPN") %>% unlist() %>% unname()
                   } else{
                     selected <- subset %>% filter(race_eth_code == code) %>% 
                       slice_sample(n = num_to_sample) %>% dplyr::select("HHIDPN") %>% 
                       unlist() %>% unname()
                   }
                   
                   if(exists("ids_vec")){
                     ids_vec <- c(ids_vec, selected)
                   } else{
                     ids_vec <- selected
                   }
                 }
                 
                 #fill in the rest with random sample from the unsampled 
                 #  observations
                 missing_obs <- num_to_sample*3 - length(ids_vec)
                 
                 if(missing_obs > 0){
                   #still unselected
                   subset %<>% filter(!HHIDPN %in% ids_vec)  
                   ids_vec <- 
                     c(ids_vec, subset %>% 
                         slice_sample(n = missing_obs) %>% 
                         dplyr::select("HHIDPN") %>% unlist() %>% unname())
                 }
                 return(ids_vec)
               })
    
    synthetic_HCAP_list[[i]][
      synthetic_HCAP_list[[i]]$HHIDPN %in% 
        c(impaired_sample_IDs, random_sample_IDs), 
      paste0("calibration_", calibration_prop*100, "_design")] <- 1
    
    #---- ****calculate sampling weight ----
    for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
      possible_counts <- 
        synthetic_HCAP_list[[i]] %>% filter(!!sym(group) == 1) %>% 
        dplyr::select("black", "hispanic", "stroke") %>% 
        unite("cell_code", c("black", "hispanic", "stroke"), sep = "") %>% 
        table() %>% as.data.frame() %>% set_colnames(c("cell_ID", "Freq"))
      
      selected_counts <- 
        synthetic_HCAP_list[[i]] %>% 
        filter(!!sym(group) == 1 &
                 !!sym(paste0("calibration_", calibration_prop*100, 
                              "_design")) == 1) %>% 
        dplyr::select("black", "hispanic", "stroke") %>% 
        unite("cell_code", c("black", "hispanic", "stroke"), sep = "") %>% 
        table() %>% as.data.frame() %>% set_colnames(c("cell_ID", "Freq")) %>% 
        left_join(possible_counts, by = "cell_ID") %>% 
        set_colnames(c("cell_ID", "selected", "possible"))
      
      cell_ID_key[which(cell_ID_key$cell_ID %in% selected_counts$cell_ID), 
                  paste0(unique(synthetic_HCAP_list[[i]]$dataset_name), "_", 
                         "calibration_", calibration_prop*100, "_design_IPW_", 
                         group)] <-
        selected_counts$possible/selected_counts$selected
    }
  }
}

#---- **save data ----
saveRDS(synthetic_HCAP_list, 
        file = paste0(path_to_box, "data/HCAP/synthetic_HCAP_list"))

#---- **save updated cell ID key ----
cell_ID_key[is.na(cell_ID_key)] <- 0
write_csv(cell_ID_key, paste0(path_to_box, "data/cell_ID_key.csv")) 
