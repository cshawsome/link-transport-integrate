#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "here", "ggbreak")
#see if I actually need these
#devtools
#meta
#install_github("thomasp85/patchwork")

#---- source functions ----
source(here::here("functions", "read_results.R"))

#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- ****results ----
results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

results <- do.call(rbind, lapply(results_paths, read_results)) %>% 
  group_by(dataset_name) %>% slice_head(n = 1000) %>% 
  mutate("HRS_sample_size" = sample_size) %>% 
  mutate("HCAP_sample_size" = HCAP_prop/100*HRS_sample_size) %>% 
  dplyr::select(-one_of("sample_size")) %>% 
  filter(calibration == "ADAMS_prior" & HCAP_prop == 50)

#---- ******check number of simulation runs ----
#there should be 1000 runs of each scenario
table(results$dataset_name, useNA = "ifany")

#---- ****superpop ----
superpop_impairment_props <- 
  read_csv(paste0(path_to_box, 
                  "data/superpopulations/impairment_class_props.csv"))

superpop_impairment_props$Group <- 
  factor(superpop_impairment_props$Group, 
         levels = c("Unimpaired", "MCI", "Dementia", "Other"))

#---- extra calcs ----
# #---- **log PR measures ----
# PR_cols <- 
#   expand_grid(c("mean", "LCI", "UCI"), "PR",
#               c("black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# results %<>% 
#   cbind(., results %>% ungroup() %>% dplyr::select(all_of(PR_cols)) %>% 
#           mutate_all(.funs = log) %>% 
#           set_colnames(str_replace(PR_cols, "_PR_", "_log_PR_")))
# 
# #---- **squared SEs ----
# #on the count scale
# for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
#   for(race in c("overall", "white", "black", "hispanic"))
#     if(race == "overall"){
#       results %<>%
#         mutate(!!sym(paste0("SE_", class)) :=
#                  (rowMeans(cbind(
#                    abs(!!sym(paste0("mean_", class)) -
#                          !!sym(paste0("LCI_", class))),
#                    abs(!!sym(paste0("mean_", class)) -
#                          !!sym(paste0("UCI_", class)))))/1.96)^2)
#     } else{
#       results %<>%
#         mutate(!!sym(paste0("SE_", class, "_", race)) :=
#                  (rowMeans(cbind(
#                    abs(!!sym(paste0("mean_", class, "_", race)) -
#                          !!sym(paste0("LCI_", class, "_", race))),
#                    abs(!!sym(paste0("mean_", class, "_", race)) -
#                          !!sym(paste0("UCI_", class, "_", race)))))/1.96)^2)
#     }
# }
# 
# for(measure in c("dem_prev", "log_PR")){
#   for(race in c("white", "black", "hispanic")){
#     if(measure == "log_PR" & race == "white"){next}
#     
#     results %<>%
#       mutate(!!sym(paste0("SE_", measure, "_", race)) :=
#                (rowMeans(cbind(
#                  abs(!!sym(paste0("mean_", measure, "_", race)) -
#                        !!sym(paste0("LCI_", measure, "_", race))),
#                  abs(!!sym(paste0("mean_", measure, "_", race)) -
#                        !!sym(paste0("UCI_", measure, "_", race)))))/1.96)^2)
#   }
# }
# 
# #---- format data ----
# results %<>% 
#   mutate("prior_sample" = 
#            case_when(calibration_sampling == "NA" ~ "ADAMS", 
#                      calibration_sampling == "SRS" & 
#                        calibration == "calibration_50" & 
#                        sampling_strata == "NA" ~ 
#                        "HCAP 50% SRS Adjudication",
#                      calibration_sampling == "SRS" & 
#                        calibration == "calibration_50" & 
#                        sampling_strata == "race" ~ 
#                        "HCAP 50% Race-stratified SRS Adjudication",
#                      calibration_sampling == "SRS" & 
#                        calibration == "calibration_35" & 
#                        sampling_strata == "NA" ~ 
#                        "HCAP 35% SRS Adjudication", 
#                      calibration_sampling == "SRS" & 
#                        calibration == "calibration_35" & 
#                        sampling_strata == "race" ~ 
#                        "HCAP 35% Race-stratified SRS Adjudication", 
#                      calibration_sampling == "SRS" & 
#                        calibration == "calibration_20" & 
#                        sampling_strata == "NA" ~ 
#                        "HCAP 20% SRS Adjudication", 
#                      calibration_sampling == "SRS" & 
#                        calibration == "calibration_20" & 
#                        sampling_strata == "race" ~ 
#                        "HCAP 20% Race-stratified SRS Adjudication"))
# 
# #---- match ADAMS data ----
# # like a cbind but making sure correct rows are matched
# #create place-holder columns
# overall_impairment_cols <- 
#   expand_grid(c("mean", "LCI", "UCI"),
#               c("Unimpaired", "MCI", "Dementia", "Other")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# by_race_impairment_cols <- 
#   expand_grid(c("mean", "LCI", "UCI"),
#               c("Unimpaired", "MCI", "Dementia", "Other"),
#               c("white", "black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# dem_prev_cols <- 
#   expand_grid(c("mean"), "dem_prev", c("white", "black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# PR_cols <- 
#   expand_grid(c("mean"), "log_PR", c("black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# impairment_SE_cols <- 
#   expand_grid(c("SE"),
#               c("Unimpaired", "MCI", "Dementia", "Other")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# impairment_SE_by_race_cols <- 
#   expand_grid(c("SE"),
#               c("Unimpaired", "MCI", "Dementia", "Other"),
#               c("white", "black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname() 
# 
# dem_prev_SE_cols <- 
#   expand_grid(c("SE"), c("dem_prev"),
#               c("white", "black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname() 
# 
# PR_SE_cols <- 
#   expand_grid(c("SE"), c("log_PR"),
#               c("black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname() 
# 
# pt_est_cols_cores <- c(overall_impairment_cols, by_race_impairment_cols, 
#                        dem_prev_cols, PR_cols)
# SE_cols_cores <- c(impairment_SE_cols, impairment_SE_by_race_cols, 
#                    dem_prev_SE_cols, PR_SE_cols)
# 
# for(prior_group in unique(results$prior_sample)[-which(
#   unique(results$prior_sample) == "ADAMS")]){
#   for(HRS_n in unique(results$HRS_sample_size)){
#     for(sample_prop in unique(results$HCAP_prop)){
#       rows <-
#         which(results$prior_sample == prior_group &
#                 results$HRS_sample_size == HRS_n &
#                 results$HCAP_prop == sample_prop)
#       
#       results[rows, c(paste0("ADAMS_", all_of(pt_est_cols_cores)), 
#                       paste0("ADAMS_", all_of(SE_cols_cores)))] <- 
#         results %>% filter(prior_sample == "ADAMS" & HRS_sample_size == HRS_n & 
#                              HCAP_prop == sample_prop) %>%
#         slice_head(n = length(rows)) %>% ungroup() %>%
#         dplyr::select(c(all_of(pt_est_cols_cores), 
#                         all_of(SE_cols_cores)))
#     }
#   }
# }
# 
# #---- **combined estimates ----
# combined_results <- results %>% filter(!prior_sample == "ADAMS")
# 
# #---- ****impairment classes ----
# for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
#   for(race in c("overall", "white", "black", "hispanic")){
#     if(race == "overall"){
#       #---- ******mean estimates ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_mean_", class)) :=
#                  rowMeans(cbind(!!sym(paste0("mean_", class)),
#                                 !!sym(paste0("ADAMS_mean_", class)))))
#       
#       #---- ******within variance ----
#       combined_results %<>%
#         mutate(!!sym(paste0("within_var_", class)) :=
#                  rowMeans(cbind(!!sym(paste0("SE_", class)),
#                                 !!sym(paste0("ADAMS_SE_", class)))))
#       
#       #---- ******between variance ----
#       combined_results %<>%
#         mutate(!!sym(paste0("between_var_", class)) :=
#                  apply(cbind(!!sym(paste0("mean_", class)),
#                              !!sym(paste0("ADAMS_mean_", class))),
#                        1, var))
#       
#       #---- ******LCIs ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_LCI_", class)) :=
#                  !!sym(paste0("combined_mean_", class)) -
#                  1.96*sqrt(!!sym(paste0("within_var_", class)) +
#                              !!sym(paste0("between_var_", class)) +
#                              !!sym(paste0("between_var_", class))/2))
#       
#       #---- ******UCIs ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_UCI_", class)) :=
#                  !!sym(paste0("combined_mean_", class)) +
#                  1.96*sqrt(!!sym(paste0("within_var_", class)) +
#                              !!sym(paste0("between_var_", class)) +
#                              !!sym(paste0("between_var_", class))/2))
#       
#       #---- ******coverage ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_", class, "_coverage")) :=
#                  (!!sym(paste0("true_", class)) >=
#                     !!sym(paste0("combined_LCI_", class)))*
#                  (!!sym(paste0("true_", class)) <=
#                     !!sym(paste0("combined_UCI_", class))))
#     } else{
#       #---- ******mean estimates ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_mean_", class, "_", race)) :=
#                  rowMeans(cbind(
#                    !!sym(paste0("mean_", class, "_", race)),
#                    !!sym(paste0("ADAMS_mean_", class, "_", race)))))
#       
#       #---- ******within variance ----
#       combined_results %<>%
#         mutate(!!sym(paste0("within_var_", class, "_", race)) :=
#                  rowMeans(cbind(!!sym(paste0("SE_", class, "_", race)),
#                                 !!sym(paste0("ADAMS_SE_", class, "_", race)))))
#       
#       #---- ******between variance ----
#       combined_results %<>%
#         mutate(!!sym(paste0("between_var_", class, "_", race)) :=
#                  apply(cbind(!!sym(paste0("mean_", class, "_", race)),
#                              !!sym(paste0("ADAMS_mean_", class, "_", race))), 1,
#                        var))
#       
#       #---- ******LCIs ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_LCI_", class, "_", race)) :=
#                  !!sym(paste0("combined_mean_", class, "_", race)) -
#                  1.96*sqrt(!!sym(paste0("within_var_", class, "_", race)) +
#                              !!sym(paste0("between_var_", class, "_", race)) +
#                              !!sym(paste0("between_var_", class, "_", race))/2))
#       
#       #---- ******UCIs ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_UCI_", class, "_", race)) :=
#                  !!sym(paste0("combined_mean_", class, "_", race)) +
#                  1.96*sqrt(!!sym(paste0("within_var_", class, "_", race)) +
#                              !!sym(paste0("between_var_", class, "_", race)) +
#                              !!sym(paste0("between_var_", class, "_", race))/2))
#       
#       #---- ******coverage ----
#       combined_results %<>%
#         mutate(!!sym(paste0("combined_", class, "_coverage", "_", race)) :=
#                  (!!sym(paste0("true_", class, "_", race)) >=
#                     !!sym(paste0("combined_LCI_", class, "_", race)))*
#                  (!!sym(paste0("true_", class, "_", race)) <=
#                     !!sym(paste0("combined_UCI_", class, "_", race))))
#     }
#   }
# }
# 
# #---- ****dementia prevalence + PR ----
# for(measure in c("dem_prev", "log_PR")){
#   for(race in c("white", "black", "hispanic")){
#     if(measure == "log_PR" & race == "white"){next}
#     
#     #---- ******mean estimates ----
#     combined_results %<>%
#       mutate(!!sym(paste0("combined_mean_", measure, "_", race)) :=
#                rowMeans(cbind(
#                  !!sym(paste0("mean_", measure, "_", race)),
#                  !!sym(paste0("ADAMS_mean_", measure, "_", race)))))
#     
#     #---- ******within variance ----
#     combined_results %<>%
#       mutate(!!sym(paste0("within_var_", measure, "_", race)) :=
#                rowMeans(cbind(!!sym(paste0("SE_", measure, "_", race)),
#                               !!sym(paste0("ADAMS_SE_", measure, "_", race)))))
#     
#     #---- ******between variance ----
#     combined_results %<>%
#       mutate(!!sym(paste0("between_var_", measure, "_", race)) :=
#                apply(cbind(!!sym(paste0("mean_", measure, "_", race)),
#                            !!sym(paste0("ADAMS_mean_", measure, "_", race))), 1,
#                      var))
#     
#     #---- ******LCIs ----
#     combined_results %<>%
#       mutate(!!sym(paste0("combined_LCI_", measure, "_", race)) :=
#                !!sym(paste0("combined_mean_", measure, "_", race)) -
#                1.96*sqrt(!!sym(paste0("within_var_", measure, "_", race)) +
#                            !!sym(paste0("between_var_", measure, "_", race)) +
#                            !!sym(paste0("between_var_", measure, "_", race))/2))
#     
#     #---- ******UCIs ----
#     combined_results %<>%
#       mutate(!!sym(paste0("combined_UCI_", measure, "_", race)) :=
#                !!sym(paste0("combined_mean_", measure, "_", race)) +
#                1.96*sqrt(!!sym(paste0("within_var_", measure, "_", race)) +
#                            !!sym(paste0("between_var_", measure, "_", race)) +
#                            !!sym(paste0("between_var_", measure, "_", race))/2))
#   }
# }
# 
# #---- ******convert back to PR ----
# cols_to_convert <- 
#   expand_grid(c("combined"), c("mean", "LCI", "UCI"), c("log_PR"),
#               c("black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname() 
# 
# combined_results %<>% 
#   cbind(., combined_results %>% ungroup() %>% 
#           dplyr::select(all_of(cols_to_convert)) %>% mutate_all(.funs = exp) %>% 
#           set_colnames(str_replace(cols_to_convert, "_log_PR_", "_PR_")))
# 
# #---- ******coverage ----
# for(measure in c("dem_prev", "PR")){
#   for(race in c("white", "black", "hispanic")){
#     if(measure == "PR" & race == "white"){next}
#     combined_results %<>%
#       mutate(!!sym(paste0("combined_", measure, "_coverage", "_", race)) :=
#                (!!sym(paste0("true_", measure, "_", race)) >=
#                   !!sym(paste0("combined_LCI_", measure, "_", race)))*
#                (!!sym(paste0("true_", measure, "_", race)) <=
#                   !!sym(paste0("combined_UCI_", measure, "_", race))))
#   }
# }
# 
# #----- ****bind results ----
# combined_results %<>% 
#   mutate("prior_sample" = paste0(prior_sample, " + ADAMS")) %>% 
#   dplyr::select(contains("combined"), contains("true"), "prior_sample", 
#                 "dataset_name", "HCAP_prop", "HRS_sample_size", 
#                 "HCAP_sample_size") 
# 
# combined_results %<>% 
#   set_colnames(str_remove(colnames(combined_results), "combined_"))
# 
# results %<>% plyr::rbind.fill(., combined_results)

#---- **proportions ----
results %<>% 
  mutate(true_Unimpaired_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "Unimpaired", "prop"]), 
         true_MCI_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "MCI", "prop"]), 
         true_Dementia_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "Dementia", "prop"]), 
         true_Other_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "Other", "prop"]))

#defining column names
overall_impairment_cols <-
  expand_grid(c("mean", "LCI", "UCI"),
              c("Unimpaired", "MCI", "Dementia", "Other")) %>%
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

by_race_impairment_cols <-
  expand_grid(c("mean", "LCI", "UCI"),
              c("Unimpaired", "MCI", "Dementia", "Other"),
              c("white", "black", "hispanic")) %>%
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

dem_prev_cols <-
  expand_grid(c("mean"), "dem_prev", c("white", "black", "hispanic")) %>%
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

PR_cols <-
  expand_grid(c("mean"), "log_PR", c("black", "hispanic")) %>%
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

pt_est_cols_cores <- c(overall_impairment_cols, by_race_impairment_cols, 
                       dem_prev_cols, PR_cols)

#selecting correct columns
calc_props <- 
  pt_est_cols_cores[-c(which(str_detect(pt_est_cols_cores, "dem_prev")), 
                       which(str_detect(pt_est_cols_cores, "log_PR")))]

results[, paste0(calc_props, "_prop")] <- 
  results[, calc_props]/results$HCAP_sample_size

#---- Figure XX: mean and 95% CI impairment class proportions ----
#---- **plot data ----
plot_data <- results %>% group_by(HRS_sample_size) %>% 
  summarise_at(paste0(
    c("mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other", 
      "LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other", 
      "UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other"), "_prop"), 
    mean) %>% 
  pivot_longer(paste0(
    c("mean_Unimpaired", "mean_MCI", "mean_Dementia", "mean_Other", 
      "LCI_Unimpaired", "LCI_MCI", "LCI_Dementia", "LCI_Other", 
      "UCI_Unimpaired", "UCI_MCI", "UCI_Dementia", "UCI_Other"), "_prop"),
    names_to = c(".value", "Group"), names_pattern = "(.*?)_(.*)") %>% 
  mutate_at("Group", function(x) str_remove(x, "_prop")) %>% 
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Group", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) 

#---- **plot ----
ggplot(data = plot_data, aes(x = mean, y = factor(HRS_sample_size))) +
  geom_vline(data = superpop_impairment_props, aes(xintercept = prop)) +
  geom_point(size = 3, shape = rep(c(19, 15, 17, 18), 3)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.3) + theme_bw() + 
  facet_grid(cols = vars(Group)) + 
  scale_x_continuous(breaks = seq(0.10, 0.40, by = 0.05)) +
  xlab("Impairment class proportion") + ylab("HRS sample size") + 
  theme(text = element_text(size = 24))  

ggsave(filename = 
         paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                "figureXX_mean_CI_impairment_class.jpeg"), 
       dpi = 300, width = 14.75, height = 3.25, units = "in")

#---- Figure XX: 95% CI coverage impairment classes ----
#---- **plot data ----
plot_data <- results %>% group_by(HRS_sample_size) %>% 
  summarise_at(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"), 
               mean) %>% 
  pivot_longer(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"),
               names_to = c("class", "coverage"), names_sep = "_") %>% 
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate(class = 
           factor(class, 
                  levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>%
  mutate("value_percent" = value*100)

#---- **plot ----
ggplot(data = plot_data, 
       aes(x = as.factor(HRS_sample_size), y = value_percent, group = class, 
           shape = class)) + 
  geom_line(size = 1, aes(linetype = class)) + geom_point(size = 3) + 
  geom_hline(yintercept = 95, lty = "dashed") +
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, seq(90, 100, by = 2))) + 
  scale_y_cut(breaks = c(89), space = 0.2, which = c(1, 2), scales = c(5, 0.5)) +
  scale_shape_manual(values = rep(c(19, 15, 17, 18), 3)) +
  theme_bw() + ylab("95% interval coverage") + xlab("HRS Sample Size") +
  labs(linetype = "Impairment Class", shape = "Impairment Class") +
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = 
         paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                "figureXX_impairment_class_coverage.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#----- Figure XX: percent bias overall + race-stratified ----
#---- **read in data ----
#---- **plot data ----
#pull the correct columns
by_race_cols <- expand_grid(c("mean", "true"), 
                            c("Unimpaired", "MCI", "Dementia", "Other"), 
                            c("white", "black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

overall_cols <- 
  paste0(expand_grid(c("mean", "true"), 
                     c("Unimpaired", "MCI", "Dementia", "Other")) %>% 
           unite("names", everything(), sep = "_") %>% unlist() %>% unname(), 
         "_overall")

cols <- c(overall_cols, by_race_cols)

plot_data <- results %>% ungroup() %>%
  dplyr::select("HRS_sample_size", "HCAP_sample_size",
                #"_overall" suffix is not in the results data
                all_of(str_remove(cols, "_overall"))) %>% 
  #rename the overall cols
  set_colnames(c("HRS_sample_size", "HCAP_sample_size", cols))

for(race in c("white", "black", "hispanic", "overall")){
  plot_data %<>% mutate(!!paste0("mean_total_", race) := 
                          !!sym(paste0("mean_Unimpaired_", race)) + 
                          !!sym(paste0("mean_MCI_", race)) + 
                          !!sym(paste0("mean_Dementia_", race)) + 
                          !!sym(paste0("mean_Other_", race)))
  
  plot_data %<>% mutate(!!paste0("true_total_", race) := 
                          !!sym(paste0("true_Unimpaired_", race)) + 
                          !!sym(paste0("true_MCI_", race)) + 
                          !!sym(paste0("true_Dementia_", race)) + 
                          !!sym(paste0("true_Other_", race)))
}

for(race in c("white", "black", "hispanic", "overall")){
  for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
    for(measure in c("mean", "true")){
      plot_data %<>% 
        mutate(!!paste(measure, class, race, sep = "_") := 
                 !!sym(paste(measure, class, race, sep = "_"))/
                 !!sym(paste(measure, "total", race, sep = "_")))
    }
  }  
}

plot_data %<>%
  pivot_longer(c(all_of(cols)),
               names_to = c(".value", "class", "race"), 
               names_sep = "_") %>%
  mutate("error" = mean - true) %>%
  mutate("percent_error" = error/true*100) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("race", function(x) 
    factor(x, levels = c("Overall", "White", "Black", "Hispanic"))) %>% 
  group_by(HRS_sample_size, class, race) %>% 
  summarize_at(c("error", "percent_error", "squared_error"), mean) %>% 
  rename(c("bias" = "error", "percent_bias" = "percent_error")) %>% 
  mutate("RMSE" = sqrt(squared_error)) 

#---- **plot ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = percent_bias, group = class, 
           shape = class)) + 
  geom_line(size = 1, aes(linetype = class)) + geom_point(size = 3) + 
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent bias") + 
  facet_grid(cols = vars(race)) + 
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_shape_manual(values = rep(c(19, 15, 17, 18), 3)) +
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) +
  labs(linetype = "Impairment Class", shape = "Impairment Class") +
  theme(text = element_text(size = 24), legend.position = "bottom")

ggsave(filename = paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                         "figureXX_impairment_class_percent_bias.jpeg"), 
       dpi = 300, width = 13.5, height = 4.5, units = "in")

#---- Figure XX: RMSE overall + race-stratified ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = RMSE, group = class, shape = class)) + 
  geom_line(size = 1, aes(linetype = class)) + geom_point(size = 3) + 
  theme_bw() + ylab("RMSE") + 
  facet_grid(cols = vars(race)) + 
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_shape_manual(values = rep(c(19, 15, 17, 18), 3)) +
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  labs(linetype = "Impairment Class", shape = "Impairment Class") +
  theme(text = element_text(size = 24), legend.position = "bottom")

ggsave(filename = paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                         "figureXX_impairment_class_RMSE.jpeg"), 
       dpi = 300, width = 13.5, height = 4.5, units = "in")

#---- Figure XX: dementia prevalence ----
#---- **plot data ----
truth <- 
  data.frame("Race" = c("White", "Black", "Hispanic"), 
             "prev" = c(unique(na.omit(results$true_dem_prev_white)), 
                        unique(na.omit(results$true_dem_prev_black)), 
                        unique(na.omit(results$true_dem_prev_hispanic)))) %>% 
  mutate_at("prev", function(x) round(x, 2)) %>% 
  mutate_at("Race", 
            function(x) factor(x, levels = c("White", "Black", "Hispanic")))

plot_data <- results %>% 
  group_by(HRS_sample_size) %>% 
  summarise_at(
    c("mean_dem_prev_white", "mean_dem_prev_black", "mean_dem_prev_hispanic", 
      "LCI_dem_prev_white", "LCI_dem_prev_black", "LCI_dem_prev_hispanic", 
      "UCI_dem_prev_white", "UCI_dem_prev_black", "UCI_dem_prev_hispanic"), 
    mean) %>% 
  pivot_longer(
    c("mean_dem_prev_white", "mean_dem_prev_black", "mean_dem_prev_hispanic", 
      "LCI_dem_prev_white", "LCI_dem_prev_black", "LCI_dem_prev_hispanic", 
      "UCI_dem_prev_white", "UCI_dem_prev_black", "UCI_dem_prev_hispanic"),
    names_to = c(".value", "Race"), names_pattern = "(.*?)_(.*)") %>% 
  mutate_at("Race", function(x) str_remove(x, "dem_prev_")) %>% 
  mutate_at("Race", str_to_sentence) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) 

#---- **plot ----
ggplot(data = plot_data, aes(x = mean, y = HRS_sample_size)) + 
  geom_vline(aes(xintercept = prev), data = truth) +
  geom_point(size = 3) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1) +
  facet_grid(cols = vars(Race)) + theme_bw() + 
  xlab("Prevalence of dementia") + ylab("HRS sample size") + 
  scale_x_continuous(limits = c(0, 0.45), breaks = seq(0.0, 0.45, by = 0.1)) + 
  theme(text = element_text(size = 24))  

ggsave(filename = paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                         "figureXX_mean_CI_dem_prev.jpeg"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

#---- Figure XX: 95% CI coverage dementia prevalence ----
#---- **plot data ----
plot_data <- results %>% 
  group_by(HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(paste0("dem_prev_coverage_", c("white", "black", "hispanic")), 
               mean) %>% 
  pivot_longer(paste0("dem_prev_coverage_", c("white", "black", "hispanic")),
               names_to = c(".value", "Race"), 
               names_pattern = "(.*?)_(.*)") %>%
  mutate_at("Race", function(x) str_remove(x, "prev_coverage_")) %>% 
  mutate_at("Race", str_to_sentence) %>%
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("coverage_percent" = dem*100)

#---- **plot ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = coverage_percent, group = Race)) + 
  geom_line() + geom_point(size = 3) + facet_grid(cols = vars(Race)) +
  geom_hline(yintercept = 95, lty = "dashed") +
  scale_y_continuous(limits = c(0, 102), breaks = c(0, seq(80, 100, by = 5))) + 
  scale_y_cut(breaks = c(79), space = 0.2, which = c(1, 2), scales = c(5, 0.5)) +
  theme_bw() + ylab("95% interval coverage") + xlab("HRS Sample Size") +
  theme(text = element_text(size = 24))      

ggsave(filename = paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                         "figureXX_dem_prev_coverage.jpeg"), 
       dpi = 300, width = 13.5, height = 4, units = "in") 

#---- Appendix Figure XX: bias dementia prevalence ----
#---- **plot data ----
cols_by_race <- expand_grid(c("mean", "true"), 
                            c("dem_prev"), 
                            c("white", "black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% ungroup() %>%
  dplyr::select("HRS_sample_size", "HCAP_sample_size",  all_of(cols_by_race)) %>%
  pivot_longer(all_of(cols_by_race),
               names_to = c(".value", "measure1", "measure2", "race"), 
               names_sep = "_") %>% 
  mutate("error" = mean - true) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>%
  mutate_at("race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>%
  group_by(HRS_sample_size, race) %>% 
  summarize_at(c("error", "squared_error"), mean) %>% 
  rename(c("bias" = "error")) %>% 
  mutate("RMSE" = sqrt(squared_error)) 

#---- **plot ----
ggplot(data = plot_data, aes(x = HRS_sample_size, y = bias, group = race)) + 
  geom_line() + geom_point(size = 3) + 
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("Bias") + xlab("HRS Sample Size") +
  facet_grid(cols = vars(race)) +
  theme(text = element_text(size = 24))      

ggsave(filename = paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                         "appendix_figureXX_dem_prev_bias.jpeg"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

#---- Figure 4.21: RMSE dementia prevalence no HCAP adjudication ----
ggplot(data = plot_data, aes(x = HRS_sample_size, y = RMSE, group = race)) + 
  geom_line() + geom_point(size = 3) + 
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(cols = vars(race)) + ylim(0, 0.08) +
  theme(text = element_text(size = 24))      

ggsave(filename = paste0(path_to_box, "papers/paper1_model_methods/figures/", 
                         "appendix_figureXX_dem_prev_rmse.jpeg"), 
       dpi = 300, width = 13.5, height = 4, units = "in") 

#---- Figure 5.18: RMSE dementia prevalence HCAP adjudication ----
ggplot(data = plot_data %>% filter(!str_detect(prior_sample, "\\+")), 
       aes(x = HRS_sample_size, y = RMSE, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black",
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(race)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.18_dem_prev_RMSE_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 7, units = "in")

#---- Figure 5.32: RMSE dementia prevalence HCAP adjudication + ADAMS ----
ggplot(data = plot_data %>% filter(str_detect(prior_sample, "\\+")), 
       aes(x = HRS_sample_size, y = RMSE, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(rep(1, 6))) +
  scale_color_manual(values = c(#green, pink, blue
    "#61bbb6", "#f35f5f","#288fb4",   
    "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(race)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 2, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 2, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.32_dem_prev_RMSE_HCAP_adjudication_plus_ADAMS.jpeg"), 
       dpi = 300, width = 23, height = 7, units = "in")

#---- Figure 4.22 + 5.19 + 5.33: PR dementia ----
#---- **plot data ----
truth <- 
  data.frame("Comparison" = c("Black vs. White", "Hispanic vs. White"), 
             "PR" = c(na.omit(unique(results$true_dem_prev_black)/
                                unique(results$true_dem_prev_white)), 
                      na.omit(unique(results$true_dem_prev_hispanic)/
                                unique(results$true_dem_prev_white))), 
             "PD" = c(na.omit(unique(results$true_dem_prev_black)-
                                unique(results$true_dem_prev_white)), 
                      na.omit(unique(results$true_dem_prev_hispanic)-
                                unique(results$true_dem_prev_white)))) %>% 
  mutate_at(c("PR", "PD"), function(x) round(x, 2)) %>% 
  mutate_at("Comparison", 
            function(x) 
              factor(x, levels = c("Black vs. White", "Hispanic vs. White")))

plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(
    c("mean_PR_black", "mean_PR_hispanic", "LCI_PR_black", "LCI_PR_hispanic", 
      "UCI_PR_black", "UCI_PR_hispanic", 
      "mean_PD_black", "mean_PD_hispanic", "LCI_PD_black", "LCI_PD_hispanic", 
      "UCI_PD_black", "UCI_PD_hispanic"), mean) %>% 
  pivot_longer(
    c("mean_PR_black", "mean_PR_hispanic", "LCI_PR_black", "LCI_PR_hispanic", 
      "UCI_PR_black", "UCI_PR_hispanic", 
      "mean_PD_black", "mean_PD_hispanic", "LCI_PD_black", "LCI_PD_hispanic", 
      "UCI_PD_black", "UCI_PD_hispanic"),
    names_to = c(".value", "measure", "Race"), names_sep = "_") %>%
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>%
  mutate_at("HCAP_prop", 
            function(x) factor(x, levels = c("HCAP Proportion\n50% of HRS", 
                                             "HCAP Proportion\n25% of HRS")))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 50% SRS Adjudication", "HCAP 35% SRS Adjudication", 
                    "HCAP 20% SRS Adjudication",
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 20% Race-stratified SRS Adjudication", 
                    "HCAP 50% SRS Adjudication + ADAMS",
                    "HCAP 35% SRS Adjudication + ADAMS",
                    "HCAP 20% SRS Adjudication + ADAMS",
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS",
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS",
                    "HCAP 20% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.22a: PR no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PR"), 
       aes(x = mean, y = HRS_sample_size)) + 
  geom_vline(aes(xintercept = PR), data = truth, color = "#ff0000", size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 3, shape = 1) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1) +
  facet_grid(cols = vars(Comparison), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("Prevalence Ratio (PR)") + ylab("HRS sample size") + 
  theme(text = element_text(size = 24))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.22a_mean_CI_PR_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 4.22b: PD no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PD"), 
       aes(x = mean, y = HRS_sample_size)) + 
  geom_vline(aes(xintercept = PD), data = truth, color = "#ff0000", size = 2) +
  geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
  geom_point(size = 3, shape = 1) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1) +
  facet_grid(cols = vars(Comparison), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("Prevalence Difference (PD)") + ylab("HRS sample size") + 
  theme(text = element_text(size = 24))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.22b_mean_CI_PD_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.19a: HCAP adjudication ----
ggplot(data = plot_data %>% filter(!str_detect(prior_sample, "\\+") & 
                                     measure == "PR"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PR), data = truth, size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black",
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + xlab("Prevalence Ratio (PR)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", 
                              nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), 
         nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.19a_mean_CI_PR_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in") 

#---- **plot 5.19b: HCAP adjudication ----
ggplot(data = plot_data %>% filter(!str_detect(prior_sample, "\\+") & 
                                     measure == "PD"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PD), data = truth, size = 2) +
  geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black",
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + xlab("Prevalence Difference (PD)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", 
                              nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), 
         nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.19b_mean_CI_PD_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in") 

#---- **defense plot 5.19a v1: HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PR"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PR), data = truth, size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 4, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values = c("black")) +
  theme_bw() + xlab("Prevalence Ratio (PR)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 32), legend.position = "none") 

ggsave(filename = paste0(path_to_box, "presentations/Defense/figures/", 
                         "figure5.19a_v1.jpeg"), 
       dpi = 300, width = 18, height = 7.25, units = "in") 

#---- **defense plot 5.19a v2: HCAP adjudication ----
ggplot(data = plot_data %>% 
         filter(prior_sample %in% 
                  c("ADAMS", "HCAP 20% SRS Adjudication", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication") & 
                  measure == "PR"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PR), data = truth, size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 4, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 3))) +
  scale_color_manual(values = c("black", 
                                "#61bbb6", "#f35f5f","#288fb4")) +
  theme_bw() + xlab("Prevalence Ratio (PR)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 32), legend.position = "none") 

ggsave(filename = paste0(path_to_box, "presentations/Defense/figures/", 
                         "figure5.19a_v2.jpeg"), 
       dpi = 300, width = 18, height = 7.25, units = "in") 

#---- **defense plot 5.19a v3: HCAP adjudication ----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & 
                  measure == "PR"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PR), data = truth, size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 4, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                "#61bbb6", "#f35f5f","#288fb4", 
                                "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + xlab("Prevalence Ratio (PR)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 32), legend.position = "none") 

ggsave(filename = paste0(path_to_box, "presentations/Defense/figures/", 
                         "figure5.19a_v3.jpeg"), 
       dpi = 300, width = 18, height = 7.25, units = "in") 

#---- **defense plot 5.19b v1: HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PD"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PD), data = truth, size = 2) +
  geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
  geom_point(size = 4, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1)) +
  scale_color_manual(values = c("black")) +
  theme_bw() + xlab("Prevalence Difference (PD)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 32), legend.position = "none")

ggsave(filename = paste0(path_to_box, "presentations/Defense/figures/", 
                         "figure5.19b_v1.jpeg"), 
       dpi = 300, width = 18, height = 7.25, units = "in")

#---- **defense plot 5.19b v2: HCAP adjudication ----
ggplot(data = plot_data %>% 
         filter(prior_sample %in% 
                  c("ADAMS", "HCAP 20% SRS Adjudication", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication") & 
                  measure == "PD"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PD), data = truth, size = 2) +
  geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
  geom_point(size = 4, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 3))) +
  scale_color_manual(values = c("black", 
                                "#61bbb6", "#f35f5f","#288fb4")) +
  theme_bw() + xlab("Prevalence Difference (PD)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 32), legend.position = "none") 

ggsave(filename = paste0(path_to_box, "presentations/Defense/figures/", 
                         "figure5.19b_v2.jpeg"), 
       dpi = 300, width = 18, height = 7.25, units = "in") 

#---- **defense plot 5.19b v3: HCAP adjudication ----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & 
                  measure == "PD"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PD), data = truth, size = 2) +
  geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
  geom_point(size = 4, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                "#61bbb6", "#f35f5f","#288fb4", 
                                "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + xlab("Prevalence Difference (PD)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 32), legend.position = "none") 

ggsave(filename = paste0(path_to_box, "presentations/Defense/figures/", 
                         "figure5.19b_v3.jpeg"), 
       dpi = 300, width = 18, height = 7.25, units = "in")

#---- **plot 5.33a: PR HCAP adjudication + ADAMS ----
ggplot(data = plot_data %>% 
         filter(str_detect(prior_sample, "\\+") & measure == "PR"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PR), data = truth, size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(rep(1, 6))) +
  scale_color_manual(values = c(#green, pink, blue
    "#61bbb6", "#f35f5f","#288fb4",   
    "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + xlab("Prevalence Ratio (PR)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 2, byrow = TRUE), 
         color = guide_legend(title = "Prior", nrow = 2, byrow = TRUE))

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.33a_mean_CI_PR_HCAP_adjudication_plus_ADAMS.jpeg"), 
       dpi = 300, width = 23, height = 8, units = "in")

#---- **plot 5.33b: PD HCAP adjudication + ADAMS ----
ggplot(data = plot_data %>% 
         filter(str_detect(prior_sample, "\\+") & measure == "PD"), 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PD), data = truth, size = 2) +
  geom_vline(aes(xintercept = 0), color = "black", lty = "dashed") +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(rep(1, 6))) +
  scale_color_manual(values = c(#green, pink, blue
    "#61bbb6", "#f35f5f","#288fb4",   
    "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + xlab("Prevalence Difference (PD)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 2, byrow = TRUE), 
         color = guide_legend(title = "Prior", nrow = 2, byrow = TRUE))

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.33b_mean_CI_PD_HCAP_adjudication_plus_ADAMS.jpeg"), 
       dpi = 300, width = 23, height = 8, units = "in")

#---- Figure 4.23 + 5.20 + 5.34: 95% CI coverage PR ----
#---- **plot data ----
plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(c(paste0("PR_coverage_", c("black", "hispanic")), 
                 paste0("PD_coverage_", c("black", "hispanic"))), mean) %>% 
  pivot_longer(c(paste0("PR_coverage_", c("black", "hispanic")), 
                 paste0("PD_coverage_", c("black", "hispanic"))),
               names_to = c("measure", ".value", "Race"), 
               names_sep = "_") %>%
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>%
  mutate_at("HCAP_prop", 
            function(x) factor(x, levels = c("HCAP Proportion\n50% of HRS", 
                                             "HCAP Proportion\n25% of HRS")))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 50% SRS Adjudication", "HCAP 35% SRS Adjudication", 
                    "HCAP 20% SRS Adjudication",
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 20% Race-stratified SRS Adjudication", 
                    "HCAP 50% SRS Adjudication + ADAMS",
                    "HCAP 35% SRS Adjudication + ADAMS",
                    "HCAP 20% SRS Adjudication + ADAMS",
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS",
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS",
                    "HCAP 20% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.23a: PR no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PR"), 
       aes(x = HRS_sample_size, y = coverage, group = Comparison, color = Comparison)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.23a_PR_coverage_no_HCAP_adjudiation.jpeg"), 
       dpi = 300, width = 13.25, height = 6.25, units = "in")

#---- **plot 4.23b: PD no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PD"), 
       aes(x = HRS_sample_size, y = coverage, group = Comparison, color = Comparison)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.23b_PD_coverage_no_HCAP_adjudiation.jpeg"), 
       dpi = 300, width = 13.25, height = 6.25, units = "in")

#---- **plot 5.20a: PR HCAP adjudication ----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & measure == "PR"), 
       aes(x = HRS_sample_size, y = coverage, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", 
                              nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), 
         nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.20a_PR_coverage_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in") 

#---- **plot 5.20b: PD HCAP adjudication ----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & measure == "PD"), 
       aes(x = HRS_sample_size, y = coverage, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", 
                              nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), 
         nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.20b_PD_coverage_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in") 

#---- **plot 5.34: HCAP adjudication + ADAMS ----
ggplot(data = plot_data %>% filter(str_detect(prior_sample, "\\+")), 
       aes(x = HRS_sample_size, y = PR, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(rep(1, 6))) +
  scale_color_manual(values = c(#green, pink, blue
    "#61bbb6", "#f35f5f","#288fb4",   
    "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 2, byrow = TRUE), 
         color = guide_legend(title = "Prior",nrow = 2, byrow = TRUE))

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.34_PR_coverage_HCAP_adjudication_plus_ADAMS.jpeg"), 
       dpi = 300, width = 23, height = 8, units = "in") 

#---- Figure 4.24 + 5.21 + 5.35: PR bias ----
#---- **plot data ----
cols_by_race <- expand_grid(c("mean", "true"), 
                            c("PR", "PD"), 
                            c("black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% ungroup() %>%
  dplyr::select("prior_sample", "HCAP_prop", "HRS_sample_size", 
                "HCAP_sample_size",  all_of(cols_by_race)) %>%
  pivot_longer(all_of(cols_by_race),
               names_to = c(".value", "measure", "race"), 
               names_sep = "_") %>% 
  mutate("error" = mean - true) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>% 
  mutate("Comparison" = case_when(race == "Black" ~ "Black vs. White", 
                                  race == "Hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>%
  group_by(prior_sample, HCAP_prop, HRS_sample_size, Comparison, measure) %>% 
  summarize_at(c("error", "squared_error"), mean) %>% 
  rename(c("bias" = "error")) %>% 
  mutate("RMSE" = sqrt(squared_error)) %>% 
  mutate_at("HCAP_prop", function(x) 
    factor(x, levels = c("HCAP Proportion\n50% of HRS", 
                         "HCAP Proportion\n25% of HRS")))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 50% SRS Adjudication", "HCAP 35% SRS Adjudication", 
                    "HCAP 20% SRS Adjudication",
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 20% Race-stratified SRS Adjudication", 
                    "HCAP 50% SRS Adjudication + ADAMS",
                    "HCAP 35% SRS Adjudication + ADAMS",
                    "HCAP 20% SRS Adjudication + ADAMS",
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS",
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS",
                    "HCAP 20% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.24a: PR bias ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PR"), 
       aes(x = HRS_sample_size, y = bias, group = Comparison)) + 
  geom_line(aes(color = Comparison), size = 1.5) + 
  geom_point(aes(color = Comparison), size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("Bias") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.24a_PR_bias_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in") 

#---- **plot 4.24b: PD bias ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PD"), 
       aes(x = HRS_sample_size, y = bias, group = Comparison)) + 
  geom_line(aes(color = Comparison), size = 1.5) + 
  geom_point(aes(color = Comparison), size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("Bias") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.24b_PD_bias_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in") 

#---- **plot 5.21a: PR bias ----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & measure == "PR"), 
       aes(x = HRS_sample_size, y = bias, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("Bias") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.21a_bias_PR_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in")

#---- **plot 5.21b: PD bias ----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & measure == "PD"), 
       aes(x = HRS_sample_size, y = bias, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("Bias") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.21b_bias_PD_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in")

#---- **plot 5.35 ----
ggplot(data = plot_data %>% filter(str_detect(prior_sample, "\\+")), 
       aes(x = HRS_sample_size, y = bias, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(rep(1, 6))) +
  scale_color_manual(values = c(#green, pink, blue
    "#61bbb6", "#f35f5f","#288fb4",   
    "#449187", "#cc435f", "#1d556f")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("Bias") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 2, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 2, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.35_bias_PR_HCAP_adjudication_plus_ADAMS.jpeg"), 
       dpi = 300, width = 23, height = 8, units = "in")

#---- Figure 4.25a: PR RMSE no HCAP adjudication----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PR"), 
       aes(x = HRS_sample_size, y = RMSE, group = Comparison)) + 
  geom_line(aes(color = Comparison), size = 1.5) + 
  geom_point(aes(color = Comparison), size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.25a_PR_RMSE_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in") 

#---- Figure 4.25b: PD RMSE no HCAP adjudication----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS" & measure == "PD"), 
       aes(x = HRS_sample_size, y = RMSE, group = Comparison)) + 
  geom_line(aes(color = Comparison), size = 1.5) + 
  geom_point(aes(color = Comparison), size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0, lty = "dashed") +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.25b_PD_RMSE_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in") 

#---- Figure 5.22a: PR RMSE HCAP adjudication----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & measure == "PR"), 
       aes(x = HRS_sample_size, y = RMSE, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.22a_RMSE_PR_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in")

#---- Figure 5.22b: PD RMSE HCAP adjudication----
ggplot(data = plot_data %>% 
         filter(!str_detect(prior_sample, "\\+") & measure == "PD"), 
       aes(x = HRS_sample_size, y = RMSE, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 6))) +
  scale_color_manual(values = c("black", 
                                #green, pink, blue
                                "#61bbb6", "#f35f5f","#288fb4",   
                                "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 3, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 3, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.22b_RMSE_PD_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 20, height = 8, units = "in")

#---- Figure 5.36: PR RMSE HCAP adjudication----
ggplot(data = plot_data %>% filter(str_detect(prior_sample, "\\+")), 
       aes(x = HRS_sample_size, y = RMSE, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(rep(1, 6))) +
  scale_color_manual(values = c(#green, pink, blue
    "#61bbb6", "#f35f5f","#288fb4",   
    "#449187", "#cc435f", "#1d556f")) +
  theme_bw() + ylab("RMSE") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", nrow = 2, byrow = TRUE), 
         color = guide_legend(title = "Prior"), nrow = 2, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.36_RMSE_PR_HCAP_adjudication_plus_ADAMS.jpeg"), 
       dpi = 300, width = 23, height = 8, units = "in")

#---- Sanity check ----
test <- results %>% 
  filter(prior_sample == "HCAP 20% Race-stratified SRS Adjudication" & 
           HRS_sample_size == 8000 & HCAP_sample_size == 4000)

View(test %>% 
       dplyr::select(c("true_PR_hispanic", "mean_PR_hispanic", "LCI_PR_hispanic", 
                       "UCI_PR_hispanic", "PR_coverage_hispanic")))

