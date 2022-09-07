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
    mutate_all(scale) %>% set_colnames(paste0(all_of(vars), "_Z"))
  
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
    expit(as.matrix(superpop[, preds]) %*% as.matrix(beta))
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

#---- **draw dem class by using a weighted vector ----
superpop[, "dem_class"] <- 
  apply(superpop[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")], 1, 
        function(x) 
          sample(c("Unimpaired", "MCI", "Dementia", "Other"), size = 1, 
                 prob = x))

superpop %<>% 
  mutate("Unimpaired" = ifelse(dem_class == "Unimpaired", 1, 0), 
         "MCI" = ifelse(dem_class == "MCI", 1, 0), 
         "Dementia" = ifelse(dem_class == "Dementia", 1, 0), 
         "Other" = ifelse(dem_class == "Other", 1, 0)) %>% 
  mutate("num_classes" = Unimpaired + MCI + Dementia + Other)

#Sanity check
table(superpop$num_classes)

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

#Sanity check
table(superpop$num_classes)

#---- **QC superpop ----
#---- ****overall summaries ----
colMeans(superpop[, c("Unimpaired", "MCI", "Dementia", "Other")])

#More dementia among women? Yes (M: 25.5%, W: 27.6%)
mean(superpop[superpop$female == 1, "Dementia"])
mean(superpop[superpop$female == 0, "Dementia"])

#More dementia among racial/ethnic minorities? Kind of? 
# (w: 26.0%, b: 33.7%, h: 22.6%)
mean(superpop[superpop$White == 1, "Dementia"])
mean(superpop[superpop$black == 1, "Dementia"])
mean(superpop[superpop$hispanic == 1, "Dementia"])

#Age: increased risk by year (RR = 1.03)
exp(coefficients(glm(Dementia ~ age, data = superpop, family = "poisson")))

#Years of Education: decreased risk by higher education (RR = 0.98)
exp(coefficients(glm(Dementia ~ edyrs, data = superpop, family = "poisson")))

#Stroke: increased risk for yes vs. no (RR = 1.69)
exp(coefficients(glm(Dementia ~ stroke, data = superpop, family = "poisson")))

#Diabetes: no increased risk for yes vs. no (RR = 1)
exp(coefficients(glm(Dementia ~ diabe, data = superpop, family = "poisson")))

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
  set_colnames(c("Female", "Age Strata", "Superpop Count")) %>% 
  dplyr::select("Superpop Count", everything())

for(race in c("White", "black", "hispanic")){
  assign(paste0(tolower(race), "_dem_risk_table"), 
         superpop %>% filter(!!sym(race) == 1) %>% 
           group_by(female, age_cat) %>% 
           summarise("dem_prop" = mean(Dementia)) %>% 
           arrange(female) %>% 
           set_colnames(c("Female", "Age Strata", paste0(race, " Dem Risk"))))
}

# #Sanity check
# test_data <- superpop %>% filter(White == 1)
# table(test_data$female, test_data$age_cat, test_data$Dementia)


#---- **save superpop data ----

#---- OLD ----
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
