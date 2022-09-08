#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "vroom", "locfit")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

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

#About XX hours for superpop
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
colMeans(superpop[, c("Unimpaired", "MCI", "Dementia", "Other")])

#More dementia among women? Yes (M: 25.5%, W: 27.6%)
mean(superpop[superpop$female == 1, "Dementia"])
mean(superpop[superpop$female == 0, "Dementia"])

#More dementia among racial/ethnic minorities? Kind of? 
# (w: 26.0%, b: 33.6%, h: 22.6%)
mean(superpop[superpop$White == 1, "Dementia"])
mean(superpop[superpop$black == 1, "Dementia"])
mean(superpop[superpop$hispanic == 1, "Dementia"])

#Age: increased risk by year (RR = 1.03)
exp(coefficients(glm(Dementia ~ age, data = superpop, family = "poisson")))

#Years of Education: decreased risk by higher education (RR = 0.98)
exp(coefficients(glm(Dementia ~ edyrs, data = superpop, family = "poisson")))

#Stroke: increased risk for yes vs. no (RR = 1.69)
exp(coefficients(glm(Dementia ~ stroke, data = superpop, family = "poisson")))

#Diabetes: no increased risk for yes vs. no (RR = 1.00)
exp(coefficients(glm(Dementia ~ diabe, data = superpop, family = "poisson")))

#Diabetes: increased risk for any impairment yes vs. no (RR = 1.10)
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
white_risk <- 
  sum(agesex_standardized$expected_white_dem_count)/nrow(superpop)
black_risk <- 
  sum(agesex_standardized$expected_black_dem_count)/nrow(superpop)
hispanic_risk <- 
  sum(agesex_standardized$expected_hispanic_dem_count)/nrow(superpop)

#RR compared to white
RR_black <- black_risk/white_risk
RR_hispanic <- hispanic_risk/white_risk

#---- **save superpop data ----
write_csv(superpop, 
          paste0(path_to_box, "data/superpopulations/superpop_1000000.csv"))

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

for(n in c(500, 1000, 2000, 4000, 8000)){
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

synthetic_HCAP_list <- 
  lapply(synthetic_HRS_list, 
         function(x) 
           x %>% group_by(married_partnered) %>% slice_sample(prop = 0.5) %>% 
           mutate("calibration_50" = 0) %>% ungroup())

#---- **flag calibration subsample ----
for(i in 1:length(synthetic_HCAP_list)){
  synthetic_HCAP_list[[i]][sample(seq(1, nrow(synthetic_HCAP_list[[i]])), 
                                  size = 0.5*nrow(synthetic_HCAP_list[[i]]), 
                                  replace = FALSE), "calibration_50"] <- 1
}

#---- **save data ----
saveRDS(synthetic_HCAP_list, 
        file = paste0(path_to_box, "data/HCAP/synthetic_HCAP_list"))