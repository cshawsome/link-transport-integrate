#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "vroom", "locfit", "miceFast", "ggforce")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_at("cell_ID", as.character) %>% 
  dplyr::select(c("cell_order", "cell_ID", "cell_name"))

#---- **imputation mat ----
hotdeck_impute_mat <- read_csv(paste0(path_to_box, "data/superpopulations/", 
                                      "hotdeck_impute_mat.csv")) %>% 
  column_to_rownames("var_names")

#---- **HRS analytic dataset ----
HRS <- read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv")) 

# #Sanity check: should be no missing data
# colMeans(is.na(HRS))[which(colMeans(is.na(HRS)) > 0)]

#---- **HCAP analytic dataset ----
HCAP <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv")) %>% 
  rename("no_drinking" = "r13no_drinking")

# #Sanity check: there should be no missing data
# colMeans(is.na(HCAP))[which(colMeans(is.na(HCAP)) > 0)]

# JZ: I'm a little confused -- should I expect to see missingness in these variables? 

#---- **ADAMS variable selection results ----
selected_vars_mat <- 
  read_csv(paste0(path_to_box, 
                  "data/variable_selection/model_coefficients.csv"))

#---- source functions ----
source(here::here("simulation_study", "functions", "hotdeck_function.R"))

#---- draw synthetic superpopulation ----
set.seed(20220905)

superpop_size <- 1000000
superpop <- sample_n(HRS %>% dplyr::select(-one_of("r13proxy", "HCAP_SELECT")), 
                     size = superpop_size, replace = TRUE) %>% 
  mutate("HHIDPN" = seq(1, superpop_size))

#add columns for neuropsych
neuropsych_cols <- colnames(HCAP)[!colnames(HCAP) %in% colnames(HRS)]
neuropsych_cols <- 
  neuropsych_cols[-c(which(neuropsych_cols == "intercept"), 
                     which(str_detect(neuropsych_cols, "cat")), 
                     which(str_detect(neuropsych_cols, "pool")))]

superpop[, neuropsych_cols] <- NA

#---- hotdeck imputation ----
superpop_imputed <- hotdeck(superpop, HCAP, hotdeck_impute_mat)

# #Sanity check
# colMeans(is.na(superpop_imputed))[which(colMeans(is.na(superpop_imputed)) > 0)]

#---- **standardize non-interaction continuous vars ----
standardize_vars <- str_remove(unique(selected_vars_mat$data_label)[
  str_detect(unique(selected_vars_mat$data_label), "_Z")], "_Z")
standardize_vars <- standardize_vars[!grepl("\\*", standardize_vars)]

Z_score <- function(data, vars){
  subset <- data %>% dplyr::select(all_of(vars)) %>% 
    mutate_all(scale) %>% mutate_all(as.numeric) %>%
    set_colnames(paste0(all_of(vars), "_Z"))
  
  data %<>% cbind(., subset)
  
  return(data)
}

superpop_imputed %<>% Z_score(., standardize_vars)

# JZ: Warning message about deprecated use of 'all_of()' outside of 
# a selecting function -- not a big problem

#---- **store superpop means and variances ----
orig_mean_sd <- superpop_imputed %>% 
  summarize_at(standardize_vars, .funs = c("mean" = mean, "sd" = sd))

means <- orig_mean_sd %>% dplyr::select(contains("mean")) %>% 
  set_colnames(standardize_vars)

sds <- orig_mean_sd %>% dplyr::select(contains("sd")) %>% 
  set_colnames(standardize_vars)

# #Sanity check
# JZ: this sanity check cannot be run 
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
#---- **** add vars ----
superpop_imputed %<>% mutate("Intercept" = 1)

# interaction_vars <- 
#   selected_vars_mat$data_label[grepl("\\*", selected_vars_mat$data_label)]
# 
# for(var in interaction_vars){
#   split_var <- str_split(var, pattern = "[*]")
#   
#   superpop_imputed %<>% 
#     mutate(!!sym(var) := !!sym(split_var[[1]][1])*!!sym(split_var[[1]][2]))
# }

#---- ****predict values in superpop ----
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  filtered_vars_mat <- selected_vars_mat[, c("data_label", class)] %>% 
    filter(!!sym(class) != 0)
  
  preds <- filtered_vars_mat$data_label
  beta <- unlist(filtered_vars_mat[, class])
  
  superpop_imputed[, paste0("p_", class)] <- 
    as.numeric(expit(as.matrix(superpop_imputed[, preds]) %*% as.matrix(beta)))
}

# #Sanity check
# View(head(superpop_imputed[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")]))

# #---- **draw impaired categories using highest prob ----
# superpop_imputed[, "dem_class"] <-
#   apply(superpop_imputed[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")], 1,
#         function(x) str_remove(names(which.max(x)), "p_"))

# #---- **draw all categories using weighted vector ----
# superpop_imputed[, "dem_class"] <-
#   apply(superpop_imputed[,
#                  c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")], 1,
#         function(x)
#           sample(c("Unimpaired", "MCI", "Dementia", "Other"), size = 1, prob = x))

#---- **draw impaired categories only using weighted vector ----
superpop_imputed[, "dem_class"] <-
  apply(superpop_imputed[, c("p_Unimpaired", "p_MCI", "p_Dementia", "p_Other")], 1,
        function(x) str_remove(names(which.max(x)), "p_"))

superpop_imputed[superpop_imputed$dem_class != "Unimpaired", "dem_class"] <-
  apply(superpop_imputed[superpop_imputed$dem_class != "Unimpaired",
                         c("p_MCI", "p_Dementia", "p_Other")], 1,
        function(x)
          sample(c("MCI", "Dementia", "Other"), size = 1, prob = x))

superpop_imputed %<>% 
  mutate("Unimpaired" = ifelse(dem_class == "Unimpaired", 1, 0), 
         "MCI" = ifelse(dem_class == "MCI", 1, 0), 
         "Dementia" = ifelse(dem_class == "Dementia", 1, 0), 
         "Other" = ifelse(dem_class == "Other", 1, 0)) %>% 
  mutate("num_classes" = Unimpaired + MCI + Dementia + Other)

# #Sanity check
# table(superpop_imputed$num_classes)

#---- **QC superpop_imputed ----
#---- ****overall summaries ----
#hotdeck: U: 37.3%, M: 16.4%, D: 26.0%, O: 20.4%
colMeans(superpop_imputed[, c("Unimpaired", "MCI", "Dementia", "Other")])

#More dementia among women? 
#hotdeck: Yes (M: 24.5%, W: 26.9%)
mean(superpop_imputed[superpop_imputed$female == 1, "Dementia"])
mean(superpop_imputed[superpop_imputed$female == 0, "Dementia"])

#More dementia among racial/ethnic minorities?  
# hotdeck: Yes (w: 24.3%, b: 33.1%, h: 29.6%)
# JZ: b should be 32.1%
mean(superpop_imputed[superpop_imputed$White == 1, "Dementia"])
mean(superpop_imputed[superpop_imputed$black == 1, "Dementia"])
mean(superpop_imputed[superpop_imputed$hispanic == 1, "Dementia"])

#Age: increased risk by year (hotdeck: PR = 1.03)
exp(coefficients(glm(Dementia ~ age, data = superpop_imputed, 
                     family = "poisson")))

#Years of Education: decreased risk by higher education (hotdeck: PR = 0.95)
exp(coefficients(glm(Dementia ~ edyrs, data = superpop_imputed, 
                     family = "poisson")))

#Stroke: increased risk for yes vs. no (hotdeck: PR = 1.72)
exp(coefficients(glm(Dementia ~ stroke, data = superpop_imputed, 
                     family = "poisson")))

#Diabetes: no increased risk for yes vs. no (hotdeck: PR = 1.03)
exp(coefficients(glm(Dementia ~ diabe, data = superpop_imputed, 
                     family = "poisson")))

#Diabetes: increased risk for any impairment yes vs. no (hotdeck: PR = 1.11)
superpop_imputed %<>% mutate("any_impairment" = Dementia + MCI)
exp(coefficients(glm(any_impairment ~ diabe, data = superpop_imputed, 
                     family = "poisson")))

#---- ****age and sex-standardized estimates by race ----
#---- ******create age strata ----
superpop_imputed %<>% 
  mutate("age_cat" = 
           cut(age, breaks = c(70, 75, 80, 85, 90, max(superpop_imputed$age)), 
               include.lowest = TRUE, right = FALSE))

# #Sanity check
# table(superpop_imputed$age_cat, useNA = "ifany")

#---- ******standardization tables ----
agesex_totals <- 
  superpop_imputed %>% group_by(female, age_cat) %>% count() %>% 
  arrange(female) %>% 
  set_colnames(c("female", "age", "superpop_imputed_count")) %>% 
  dplyr::select("superpop_imputed_count", everything())

for(race in c("White", "black", "hispanic")){
  assign(paste0(tolower(race), "_dem_risk_table"), 
         superpop_imputed %>% filter(!!sym(race) == 1) %>% 
           group_by(female, age_cat) %>% 
           summarise("dem_prop" = mean(Dementia)) %>% 
           arrange(female) %>% 
           set_colnames(c("female", "age", 
                          paste0(tolower(race), "_dem_risk"))))
}

# #Sanity check
# test_data <- superpop_imputed %>% filter(White == 1)
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
             !!sym(paste0(race, "_dem_risk"))*superpop_imputed_count)
}

# #Sanity check
# agesex_standardized$white_dem_risk*agesex_standardized$superpop_imputed_count
# agesex_standardized$hispanic_dem_risk*agesex_standardized$superpop_imputed_count

#---- ******estimates ----
#hotdeck: 0.241
white_risk <- 
  sum(agesex_standardized$expected_white_dem_count)/nrow(superpop_imputed)
#hotdeck: 0.331
black_risk <- 
  sum(agesex_standardized$expected_black_dem_count)/nrow(superpop_imputed)
#hotdeck: 0.302
hispanic_risk <- 
  sum(agesex_standardized$expected_hispanic_dem_count)/nrow(superpop_imputed)

#PR compared to white
#hotdeck: 1.37
PR_black <- black_risk/white_risk
#hotdeck: 1.25
PR_hispanic <- hispanic_risk/white_risk

#---- **save superpop_imputed data ----
write_csv(superpop_imputed, 
          paste0(path_to_box, "data/superpopulations/superpop_1000000.csv"))

#---- **save impairment class props ----
write_csv(colMeans(superpop_imputed[, c("Unimpaired", "MCI", "Dementia", "Other")]) %>% 
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
create_HRS_datasets <- function(superpop_imputed, n){
  sample_n(superpop_imputed, size = n) %>% 
    mutate("dataset_name" = paste0("HRS_", n), 
           "(Intercept)" = 1)
}

set.seed(20220507)

for(n in c(2000, 4000, 8000)){
  if(!exists("synthetic_HRS_list")){
    synthetic_HRS_list <- list(create_HRS_datasets(superpop_imputed, n))
  } else{
    synthetic_HRS_list <- 
      append(synthetic_HRS_list, list(create_HRS_datasets(superpop_imputed, n)))
  }
}

#---- **save data ----
saveRDS(synthetic_HRS_list, 
        file = paste0(path_to_box, "data/HRS/synthetic_HRS_list"))

#---- synthetic HCAP ----
#---- **create one set of synthetic HCAP ----
set.seed(20220507)

sample_props <- c(0.25, 0.50)
# JZ: in the diagram (Fig 3), are we only doing the version with sample prop of 0.5? 

for(sample_prop in sample_props){
  if(exists("synthetic_HCAP_list")){
    synthetic_HCAP_list <- 
      append(synthetic_HCAP_list, 
             lapply(synthetic_HRS_list, 
                    function(x) 
                      x %>% group_by(married_partnered) %>% 
                      slice_sample(prop = sample_prop) %>% 
                      mutate("calibration_20_SRS" = 0,
                             "calibration_35_SRS" = 0,
                             "calibration_50_SRS" = 0,
                             "calibration_20_SRS_race" = 0,
                             "calibration_35_SRS_race" = 0,
                             "calibration_50_SRS_race" = 0) %>% ungroup() %>% 
                      mutate("dataset_name" = 
                               paste0(dataset_name, "_sample_", 
                                      sample_prop*100))))
  } else{
    synthetic_HCAP_list <- 
      lapply(synthetic_HRS_list, 
             function(x) 
               x %>% group_by(married_partnered) %>% 
               slice_sample(prop = sample_prop) %>% 
               mutate("calibration_20_SRS" = 0,
                      "calibration_35_SRS" = 0,
                      "calibration_50_SRS" = 0, 
                      "calibration_20_SRS_race" = 0,
                      "calibration_35_SRS_race" = 0,
                      "calibration_50_SRS_race" = 0) %>% ungroup() %>% 
               mutate("dataset_name" = 
                        paste0(dataset_name, "_sample_", sample_prop*100)))
  }
}

#---- **flag calibration subsamples ----
#---- ****SRS calibration ----
for(i in 1:length(synthetic_HCAP_list)){
  for(calibration_prop in c(0.20, 0.35, 0.50)){
    synthetic_HCAP_list[[i]][
      sample(seq(1, nrow(synthetic_HCAP_list[[i]])), 
             size = calibration_prop*nrow(synthetic_HCAP_list[[i]]), 
             replace = FALSE), 
      paste0("calibration_", calibration_prop*100,"_SRS")] <- 1 
  }
}

# #Sanity check
# lapply(synthetic_HCAP_list, function(x) 
#   colMeans(x[, paste0("calibration_", c(20, 35, 50),"_SRS")]))

#---- ****SRS race calibration ----
for(i in 1:length(synthetic_HCAP_list)){
  sample_race_props <- 
    colMeans(synthetic_HCAP_list[[i]][, c("White", "black", "hispanic")])
  
  for(calibration_prop in c(0.20, 0.35, 0.50)){
    #set Black and Hispanic proportions
    B = 0.60
    H = 0.60
    
    #calculate p(selection) for White participants
    W = (calibration_prop - B*sample_race_props["black"] - 
           H*sample_race_props["hispanic"])/sample_race_props["White"]
    
    #select IDs
    if(exists("race_IDs")){rm(race_IDs)}
    
    for(race in c("White", "black", "hispanic")){
      race_prop = case_when(race == "White" ~ W, 
                            race == "black" ~ B, 
                            race == "hispanic" ~ H)
      
      subset_IDs <- synthetic_HCAP_list[[i]] %>% filter(!!sym(race) == 1) %>% 
        slice_sample(prop = race_prop) %>% dplyr::select("HHIDPN") %>% 
        unlist() %>% unname()
      
      if(exists("race_IDs")){
        race_IDs <- c(race_IDs, subset_IDs)
      } else{
        race_IDs <- subset_IDs
      }
    }
    
    #flag sample
    synthetic_HCAP_list[[i]][
      synthetic_HCAP_list[[i]]$HHIDPN %in% race_IDs, 
      paste0("calibration_", calibration_prop*100,"_SRS_race")] <- 1
    
    
    #store weights
    cell_ID_key %<>% 
      mutate(!!sym(paste0("calibration_", calibration_prop*100,"_SRS_race")) := 
               case_when(str_detect(cell_name, "White") ~ 1/W, 
                         str_detect(cell_name, "Black") ~ 1/B, 
                         str_detect(cell_name, "Hispanic") ~ 1/H))
  }
}

# #Calculations
# colMeans(synthetic_HCAP_list[[1]][, c("White", "black", "hispanic")])
# colMeans(synthetic_HCAP_list[[6]][, c("White", "black", "hispanic")])
# 
# 0.5 = W*0.765 + B*0.145 + H*0.090
# Set B and H = 0.6 to sample 60% of Black and Hispanic participants
# Then solve for W

# #Sanity check
# lapply(synthetic_HCAP_list, function(x)
#   colMeans(x[, paste0("calibration_", c(20, 35, 50),"_SRS_race")]))

# #---- ****design calibration ----
# for(i in 1:length(synthetic_HCAP_list)){
# 
#   synthetic_HCAP_list[[i]] %<>%
#     mutate("impaired" = ifelse(Dementia == 1 | MCI == 1, 1, 0))
# 
#   synthetic_HCAP_list[[i]] %<>%
#     unite("sample_cells", c("impaired", "black", "hispanic", "stroke"), sep = "",
#           remove = FALSE)
# 
#   # #test cells
#   # table(synthetic_HCAP_list[[i]]$sample_cells) %>% as.data.frame() %>%
#   #   mutate("prop" = Freq/sum(Freq))
# 
#   for(calibration_prop in c(0.50)){
#     # #to start fresh for troubleshooting
#     # synthetic_HCAP_list[[i]]$calibration_50_design <- 0
# 
#     num_impaired <- round(0.60*calibration_prop*nrow(synthetic_HCAP_list[[i]]))
#     num_unimpaired <-
#       round(calibration_prop*nrow(synthetic_HCAP_list[[i]]) - num_impaired)
# 
#     #sample 60% within each race/ethnicity x stroke cell
#     prop_bh <- 0.60
# 
#     bh_sample_IDs <-
#       synthetic_HCAP_list[[i]] %>% filter(White == 0) %>%
#       group_by(sample_cells) %>% slice_sample(prop = prop_bh)
# 
#     #sample remaining from white participants
#     num_bh_impaired <- sum(bh_sample_IDs$impaired)
#     num_bh_unimpaired <- nrow(bh_sample_IDs) - num_bh_impaired
# 
#     num_w_impaired <- num_impaired - num_bh_impaired
#     num_w_unimpaired <- num_unimpaired - num_bh_unimpaired
# 
#     # #Sanity check
#     # num_unimpaired == num_bh_unimpaired + num_w_unimpaired
#     # num_impaired == num_bh_impaired + num_w_impaired
# 
#     w_impaired_sample_IDs <-
#       synthetic_HCAP_list[[i]] %>% filter(White == 1 & impaired == 1) %>%
#       slice_sample(n = num_w_impaired)
# 
#     w_unimpaired_sample_IDs <-
#       synthetic_HCAP_list[[i]] %>% filter(White == 1 & impaired == 0) %>%
#       slice_sample(n = num_w_unimpaired)
# 
#     synthetic_HCAP_list[[i]][
#       synthetic_HCAP_list[[i]]$HHIDPN %in%
#         c(bh_sample_IDs$HHIDPN, w_impaired_sample_IDs$HHIDPN,
#           w_unimpaired_sample_IDs$HHIDPN),
#       paste0("calibration_", calibration_prop*100, "_design")] <- 1
# 
#     #Sanity check props
#     print("entire sample:")
#     synthetic_HCAP_list[[i]] %>%
#       unite("race_cells", c(impaired, black, hispanic), sep = "") %>%
#       dplyr::select("race_cells") %>% table()/nrow(synthetic_HCAP_list[[i]])
# 
#     print("random sample:")
#     synthetic_HCAP_list[[i]] %>%
#       filter(!!sym(paste0("calibration_", calibration_prop*100, "_SRS")) == 1) %>%
#       unite("race_cells", c(impaired, black, hispanic), sep = "") %>%
#       dplyr::select("race_cells") %>%
#       table()/sum(synthetic_HCAP_list[[i]][, paste0("calibration_",
#                                                     calibration_prop*100, "_SRS")])
# 
#     print("design sample:")
#     synthetic_HCAP_list[[i]] %>%
#       filter(!!sym(paste0("calibration_", calibration_prop*100, "_design")) == 1) %>%
#       unite("race_cells", c(impaired, black, hispanic), sep = "") %>%
#       dplyr::select("race_cells") %>%
#       table()/sum(synthetic_HCAP_list[[i]][, paste0("calibration_",
#                                                     calibration_prop*100, "_design")])
# 
#     #---- ****calculate sampling weight ----
#     for(group in c("Unimpaired", "MCI", "Dementia", "Other")){
#       #calculate weight for white participants
#       possible_counts <-
#         synthetic_HCAP_list[[i]] %>% filter(!!sym(group) == 1 & White == 1) %>%
#         dplyr::select("black", "hispanic", "stroke") %>%
#         unite("cell_code", c("black", "hispanic", "stroke"), sep = "") %>%
#         table() %>% as.data.frame() %>% set_colnames(c("cell_ID", "Freq"))
# 
#       selected_counts <-
#         synthetic_HCAP_list[[i]] %>%
#         filter(!!sym(group) == 1 & White == 1 &
#                  !!sym(paste0("calibration_", calibration_prop*100,
#                               "_design")) == 1) %>%
#         dplyr::select("black", "hispanic", "stroke") %>%
#         unite("cell_code", c("black", "hispanic", "stroke"), sep = "") %>%
#         table() %>% as.data.frame() %>% set_colnames(c("cell_ID", "Freq")) %>%
#         left_join(possible_counts, by = "cell_ID") %>%
#         set_colnames(c("cell_ID", "selected", "possible"))
# 
#       cell_ID_key[which(cell_ID_key$cell_ID %in% selected_counts$cell_ID),
#                   paste0(unique(synthetic_HCAP_list[[i]]$dataset_name), "_",
#                          "calibration_", calibration_prop*100, "_design_IPW_",
#                          group)] <-
#         selected_counts$possible/selected_counts$selected
# 
#       #same weight for all black and hispanic participants
#       cell_ID_key[which(!cell_ID_key$cell_ID %in% selected_counts$cell_ID),
#                   paste0(unique(synthetic_HCAP_list[[i]]$dataset_name), "_",
#                          "calibration_", calibration_prop*100, "_design_IPW_",
#                          group)] <- 1/prop_bh
#     }
#   }
# }

#---- **save data ----
saveRDS(synthetic_HCAP_list, 
        file = paste0(path_to_box, "data/HCAP/synthetic_HCAP_list"))

#---- **save updated cell ID key ----
cell_ID_key[is.na(cell_ID_key)] <- 0
write_csv(cell_ID_key, paste0(path_to_box, "data/cell_ID_key.csv")) 
