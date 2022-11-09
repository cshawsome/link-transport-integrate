#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "devtools", "here", "meta")
install_github("thomasp85/patchwork")

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
  dplyr::select(-one_of("sample_size"))

#---- ****superpop ----
superpop_impairment_props <- 
  read_csv(paste0(path_to_box, 
                  "data/superpopulations/impairment_class_props.csv"))
superpop_impairment_props$Group <- 
  factor(superpop_impairment_props$Group, 
         levels = c("Unimpaired", "MCI", "Dementia", "Other"))

#---- check number of simulation runs ----
#Number of missing runs noted in simulation study log
table(results$dataset_name, useNA = "ifany")

#---- extra calcs ----
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

overall_cols <- expand_grid(c("mean", "LCI", "UCI"),
                            c("Unimpaired", "MCI", "Dementia", "Other")) %>%
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

by_race_cols <- expand_grid(c("mean", "LCI", "UCI"),
                            c("Unimpaired", "MCI", "Dementia", "Other"),
                            c("white", "black", "hispanic")) %>%
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

calc_props <- c(overall_cols, by_race_cols)

results[, paste0(calc_props, "_prop")] <- 
  results[, calc_props]/results$HCAP_sample_size

#---- **SEs ----
#on the prop scale
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  for(race in c("overall", "white", "black", "hispanic"))
    if(race == "overall"){
      results %<>%
        mutate(!!sym(paste0("SE_", class)) :=
                 rowMeans(cbind(
                   abs(!!sym(paste0("mean_", class, "_prop")) -
                         !!sym(paste0("LCI_", class, "_prop"))),
                   abs(!!sym(paste0("mean_", class, "_prop")) -
                         !!sym(paste0("UCI_", class, "_prop")))))/1.96)
    } else{
      results %<>%
        mutate(!!sym(paste0("SE_", class, "_", race)) :=
                 rowMeans(cbind(
                   abs(!!sym(paste0("mean_", class, "_", race, "_prop")) -
                         !!sym(paste0("LCI_", class, "_", race, "_prop"))),
                   abs(!!sym(paste0("mean_", class, "_", race, "_prop")) -
                         !!sym(paste0("UCI_", class, "_", race, "_prop")))))/1.96)
    }
}

#---- format data ----
results %<>% 
  mutate("prior_sample" = 
           case_when(calibration_sampling == "NA" ~ "ADAMS", 
                     calibration_sampling == "SRS" & 
                       calibration == "calibration_50" & 
                       sampling_strata == "NA" ~ 
                       "HCAP 50% SRS Adjudication",
                     calibration_sampling == "SRS" & 
                       calibration == "calibration_50" & 
                       sampling_strata == "race" ~ 
                       "HCAP 50% Race-stratified SRS Adjudication",
                     calibration_sampling == "SRS" & 
                       calibration == "calibration_35" & 
                       sampling_strata == "NA" ~ 
                       "HCAP 35% SRS Adjudication", 
                     calibration_sampling == "SRS" & 
                       calibration == "calibration_35" & 
                       sampling_strata == "race" ~ 
                       "HCAP 35% Race-stratified SRS Adjudication"))

#fill in matching ADAMS data (like a cbind but making sure correct rows are
# matched)
#create place-holder columns
new_cols_cores <- c(overall_cols, by_race_cols)

for(prior_group in unique(results$prior_sample)[-which(
  unique(results$prior_sample) == "ADAMS")]){
  for(HRS_n in unique(results$HRS_sample_size)){
    for(sample_prop in unique(results$HCAP_prop)){
      rows <-
        which(results$prior_sample == prior_group &
                results$HRS_sample_size == HRS_n &
                results$HCAP_prop == sample_prop)
      
      results[rows, c(paste0("ADAMS_", all_of(new_cols_cores)), 
                      paste0("ADAMS_", all_of(new_cols_cores), "_prop"))] <- 
        results %>% filter(prior_sample == "ADAMS" & HRS_sample_size == HRS_n & 
                             HCAP_prop == sample_prop) %>%
        slice_head(n = length(rows)) %>% ungroup() %>%
        dplyr::select(c(all_of(new_cols_cores)), 
                      paste0(all_of(new_cols_cores), "_prop"))
    }
  }
}

#---- **meta analysis ----
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  for(i in 1:nrow(results %>% filter(prior_sample != "ADAMS"))){
    m.gen <- metagen(TE = log_transform, seTE = SE, studlab = prior_sample,
                     data = results[i, ] %>% mutate("log_transform" = log()), 
                     method.tau = "PM", sm = "PLN", random = TRUE, hakn = TRUE)
  }
}




# for(prior_group in unique(results$prior_sample)[-which(
#   unique(results$prior_sample) == "ADAMS")]){
#   for(HRS_n in unique(results$HRS_sample_size)){
#     for(sample_prop in unique(results$HCAP_prop)){
#       for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
#         subset <- results %>% ungroup() %>%
#           filter(prior_sample %in% c("ADAMS", prior_group) & 
#                    HRS_sample_size == HRS_n & HCAP_prop == sample_prop) %>% 
#           dplyr::select(c("prior_sample", "HCAP_prop", "HRS_sample_size", 
#                           paste0("mean_", class, "_prop"), 
#                           paste0("SE_", class))) %>% 
#           set_colnames(c(
#             "prior_sample", "HCAP_prop", "HRS_sample_size", "mean", "SE")) %>%
#           mutate("log_transform" = log(mean)) %>% 
#           group_by(prior_sample, HCAP_prop, HRS_sample_size) %>% 
#           summarize_at(c("mean", "SE", "log_transform"), mean)
#         
#         m.gen <- metagen(TE = log_transform, seTE = SE, studlab = prior_sample,
#                          data = subset, method.tau = "PM", sm = "PLN",
#                          random = TRUE, hakn = TRUE)
#         
#         summary(m.gen)
#       }
#     }
#   }
# }

#---- Figure 4.XX + 5.XX: mean and 95% CI impairment class proportions ----
#---- **read in data ----
#---- ****color palette ----
color_palette <- read_csv(here::here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HCAP_sample_size, HRS_sample_size) %>% 
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
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS"))

#---- **plot 4.XX: no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = mean, y = factor(HRS_sample_size))) +
  geom_vline(data = superpop_impairment_props, aes(xintercept = prop), 
             size = 2, color = rep(color_palette$Color, 2)) +
  geom_point(size = 4, shape = 1) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.3, size = 1) + theme_bw() + 
  facet_grid(cols = vars(Group), rows = vars(HCAP_prop), scales = "free_y") + 
  scale_x_continuous(breaks = seq(0.10, 0.40, by = 0.05)) +
  xlab("Impairment class proportion") + ylab("HRS sample size") + 
  theme(text = element_text(size = 24))  

ggsave(filename = 
         paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                "figure4.XX_mean_CI_impairement_class_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 14.75, height = 6.5, units = "in")

#---- **plot 5.XX: HCAP calibration ----
#add data for 100% adjudication
HCAP_all_adjudicated <- results %>% ungroup() %>%
  dplyr::select(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other", 
                  "HCAP_sample_size", "HRS_sample_size", "HCAP_prop")) %>% 
  group_by(HCAP_prop, HCAP_sample_size, HRS_sample_size) %>% na.omit() %>%
  summarise_at(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
               mean) %>% 
  pivot_longer(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
               names_to = "Group", values_to = "mean") %>% 
  mutate_at("Group", function(x) str_remove(x, "true_")) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>% 
  mutate(mean = mean/HCAP_sample_size) %>%
  mutate("calibration" = "calibration_100", 
         "prior_sample" = "HCAP 100% Adjudication",
         "LCI" = NA, "UCI" = NA) %>% mutate_at("HRS_sample_size", as.factor)

plot_data %<>% rbind(., HCAP_all_adjudicated) 

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication"))
                    # "HCAP 35% SRS Adjudication + ADAMS", 
                    # "HCAP 50% SRS Adjudication + ADAMS", 
                    # "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    # "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))
plot_data$Group <- 
  factor(plot_data$Group, levels = c("Unimpaired", "MCI", "Dementia", "Other"))

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, color = prior_sample, 
           shape = prior_sample)) +
  geom_vline(data = superpop_impairment_props, aes(xintercept = prop), 
             size = 1.5) +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  scale_shape_manual(values = c(15, 1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("#fbb040", "black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) + theme_bw() + 
  facet_grid(cols = vars(Group), rows = vars(HCAP_prop), scales = "free_y") + 
  scale_x_continuous(breaks = seq(0.00, 0.70, by = 0.10)) +
  xlab("Impairment class proportion") + ylab("HRS sample size") + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Prior", title.vjust = 0.5, nrow = 3, 
                              byrow = TRUE), 
         color = guide_legend(title = "Prior", title.vjust = 0.5, nrow = 3, 
                              byrow = TRUE))  

ggsave(filename = 
         paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                "figure5.XX_mean_CI_impairement_class_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 17, height = 11, units = "in")

#---- Figure 4.XX + 5.XX: 95% CI coverage impairment classes ----
#---- **read in data ----
#---- ****color palette ----
color_palette <- read_csv(here::here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HCAP_sample_size, HRS_sample_size) %>% 
  summarise_at(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"), 
               mean) %>% 
  pivot_longer(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"),
               names_to = c("class", "coverage"), names_sep = "_") %>% 
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>%
  mutate(class = factor(class, 
                        levels = c("Unimpaired", "MCI", "Dementia", "Other")))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XX: no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = as.factor(HRS_sample_size), y = value, group = class, 
           color = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), scales = "free_x") + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = 
         paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                "figure4.XX_impairement_class_coverage_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.XX: HCAP calibration ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = value, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) + 
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(class)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = 
         paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                "figure5.XX_impairement_class_coverage_HCAP_adjudication.jpeg"), 
       dpi = 300, height = 7.5, width = 17, units = "in")

#----- Figure 4.XXa-c + 5.XXa-c: bias, percent bias, RMSE impairment class props ----
#---- **read in data ----
#---- ****color palette ----
color_palette <- read_csv(here::here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% ungroup() %>%
  dplyr::select("prior_sample", "HCAP_prop", "HCAP_sample_size", 
                "HRS_sample_size", 
                paste0("mean_", c("Unimpaired", "MCI", "Dementia", "Other")), 
                paste0("true_", c("Unimpaired", "MCI", "Dementia", "Other"))) %>%
  pivot_longer(c(paste0("mean_", c("Unimpaired", "MCI", "Dementia", "Other")), 
                 paste0("true_", c("Unimpaired", "MCI", "Dementia", "Other"))),
               names_to = c(".value", "class"), 
               names_sep = "_") %>% 
  mutate("mean" = mean/HCAP_sample_size, "true" = true/HCAP_sample_size) %>%
  mutate("error" = mean - true) %>% 
  mutate("squared_error" = error^2) %>%
  mutate("percent_error" = error/true*100) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>% 
  group_by(prior_sample, HRS_sample_size, HCAP_prop, class) %>% 
  summarize_at(c("error", "percent_error", "squared_error"), mean) %>% 
  rename(c("bias" = "error", "percent_bias" = "percent_error")) %>% 
  mutate("RMSE" = sqrt(squared_error))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XXa: bias + no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = bias, group = class, color = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + ylab("Bias") + 
  facet_grid(rows = vars(HCAP_prop), scales = "free_x") + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24)) + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")     

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XXa_impairement_class_bias_prop_no_HCAP_", 
                         "adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.XXa: bias + HCAP calibration ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = bias, color = prior_sample, 
           group = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Bias") + 
  facet_grid(rows = vars(HCAP_prop), cols = vars(class), scales = "free_x") + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24), legend.position = "bottom")  + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XXa_impairement_class_bias_prop_HCAP_", 
                         "adjudication.jpeg"), 
       dpi = 300, width = 17, height = 7.50, units = "in")

#---- **plot 4.XXb: percent bias + no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = percent_bias, group = class, 
           color = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent bias") + 
  facet_grid(rows = vars(HCAP_prop), scales = "free_x") + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24)) + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XXb_impairement_class_percent_bias_prop_no_", 
                         "HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.XXb: percent bias + HCAP calibration ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = percent_bias, color = prior_sample, 
           group = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent bias") + 
  facet_grid(rows = vars(HCAP_prop), cols = vars(class), scales = "free_x") + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24), legend.position = "bottom")  + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XXb_impairement_class_percent_bias_prop_HCAP_", 
                         "adjudication.jpeg"), 
       dpi = 300, width = 17, height = 7.50, units = "in")

#---- **plot 4.XXc: RMSE + no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = RMSE, group = class, color = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + theme_bw() + 
  scale_color_manual(values = group_colors) +
  ylab("RMSE") + 
  facet_grid(rows = vars(HCAP_prop), scales = "free_x") + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24)) + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")   

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XXc_impairement_class_RMSE_prop_no_HCAP_", 
                         "adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.XXc: RMSE + HCAP calibration ----
ggplot(data = plot_data, aes(x = HRS_sample_size, y = RMSE, color = prior_sample, 
                             group = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("RMSE") + 
  facet_grid(rows = vars(HCAP_prop), cols = vars(class), scales = "free_x") + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24), legend.position = "bottom")  + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XXc_RMSE_bias_prop_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 17, height = 7.50, units = "in")

#---- Figure 4.XXa-c + 5.XXa-c: race-stratified bias, percent bias, RMSE impairment class props ----
#---- **read in data ----
#---- ****color palette ----
color_palette <- read_csv(here::here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
cols_by_race <- expand_grid(c("mean", "true"), 
                            c("Unimpaired", "MCI", "Dementia", "Other"), 
                            c("white", "black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% ungroup() %>%
  dplyr::select("prior_sample", "HCAP_prop", "HRS_sample_size", 
                "HCAP_sample_size",  all_of(cols_by_race)) %>%
  pivot_longer(all_of(cols_by_race),
               names_to = c(".value", "class", "race"), 
               names_sep = "_") %>% 
  mutate("mean" = mean/HCAP_sample_size, "true" = true/HCAP_sample_size) %>%
  mutate("error" = mean - true) %>%
  mutate("percent_error" = error/true*100) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) %>% 
  group_by(prior_sample, HCAP_prop, HRS_sample_size, class, race) %>% 
  summarize_at(c("error", "percent_error", "squared_error"), mean) %>% 
  rename(c("bias" = "error", "percent_bias" = "percent_error")) %>% 
  mutate("RMSE" = sqrt(squared_error))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XXa: bias by race + no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = bias, color = class, group = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Bias") + 
  facet_grid(cols = vars(race), rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24), legend.position = "bottom")    

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XXa_impairement_class_bias_prop_by_race_no_", 
                         "HCAP_calibration.jpeg"), 
       dpi = 300, width = 13.5, height = 7.25, units = "in")

#---- **plot 5.XXa: by race + HCAP calibration ----
for(race_subset in unique(plot_data$race)){
  ggplot(data = plot_data %>% filter(race == race_subset), 
         aes(x = HRS_sample_size, y = bias, color = prior_sample, 
             group = prior_sample, shape = prior_sample)) + 
    geom_line(size = 1.5) + 
    geom_point(size = 3) + 
    scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
    scale_color_manual(values = c("black", 
                                  "#288fb4", "#1d556f", 
                                  "#f35f5f", "#cc435f", 
                                  "#288fb4", "#1d556f", 
                                  "#f35f5f", "#cc435f")) +
    geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
    ylab("Bias") + 
    facet_grid(rows = vars(HCAP_prop), cols = vars(class), scales = "free_x") + 
    scale_x_discrete(name = "HRS Sample Size", 
                     breaks = unique(plot_data$HRS_sample_size)) + 
    theme(text = element_text(size = 24), legend.position = "bottom")  + 
    guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                                nrow = 5, byrow = FALSE), 
           color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
           nrow = 5, byrow = FALSE) + ggtitle(str_to_sentence(race_subset))
  
  ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                           "figure5.XXa_impairement_class_bias_prop_", 
                           race_subset, "_HCAP_adjudication.jpeg"), 
         dpi = 300, width = 17, height = 7.50, units = "in")
  
}

#---- **plot 4.XXb: percent bias by race + no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = percent_bias, color = class, group = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent bias") + 
  facet_grid(cols = vars(race), rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24), legend.position = "bottom")    

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XXb_impairement_class_percent_bias_prop_by_",
                         "race_no_HCAP_calibration.jpeg"), 
       dpi = 300, width = 13.5, height = 7.25, units = "in")

#---- **plot 5.XXb: percent by race + HCAP calibration ----
for(race_subset in unique(plot_data$race)){
  ggplot(data = plot_data %>% filter(race == race_subset), 
         aes(x = HRS_sample_size, y = percent_bias, color = prior_sample, 
             group = prior_sample, shape = prior_sample)) + 
    geom_line(size = 1.5) + 
    geom_point(size = 3) + 
    scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
    scale_color_manual(values = c("black", 
                                  "#288fb4", "#1d556f", 
                                  "#f35f5f", "#cc435f", 
                                  "#288fb4", "#1d556f", 
                                  "#f35f5f", "#cc435f")) +
    geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
    ylab("Percent bias") + 
    facet_grid(rows = vars(HCAP_prop), cols = vars(class), scales = "free_x") + 
    scale_x_discrete(name = "HRS Sample Size", 
                     breaks = unique(plot_data$HRS_sample_size)) + 
    theme(text = element_text(size = 24), legend.position = "bottom")  + 
    guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                                nrow = 5, byrow = FALSE), 
           color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
           nrow = 5, byrow = FALSE) + ggtitle(str_to_sentence(race_subset))
  
  ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                           "figure5.XXb_impairement_class_percent_bias_prop_", 
                           race_subset, "_HCAP_adjudication.jpeg"), 
         dpi = 300, width = 17, height = 7.50, units = "in")
}

#---- **plot 4.XXc: RMSE by race + no HCAP calibration ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = RMSE, color = class, group = class)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("RMSE") + 
  facet_grid(cols = vars(race), rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  theme(text = element_text(size = 24), legend.position = "bottom")    

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XXc_impairement_class_rmse_by_race_no_HCAP_", 
                         "calibration.jpeg"), 
       dpi = 300, width = 13.5, height = 7.25, units = "in")

#---- **plot 5.XXC: RMSE by race + HCAP calibration ----
for(race_subset in unique(plot_data$race)){
  ggplot(data = plot_data %>% filter(race == race_subset), 
         aes(x = HRS_sample_size, y = RMSE, color = prior_sample, 
             group = prior_sample, shape = prior_sample)) + 
    geom_line(size = 1.5) + 
    geom_point(size = 3) + 
    scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
    scale_color_manual(values = c("black", 
                                  "#288fb4", "#1d556f", 
                                  "#f35f5f", "#cc435f", 
                                  "#288fb4", "#1d556f", 
                                  "#f35f5f", "#cc435f")) +
    geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
    ylab("RMSE") + 
    facet_grid(rows = vars(HCAP_prop), cols = vars(class), scales = "free_x") + 
    scale_x_discrete(name = "HRS Sample Size", 
                     breaks = unique(plot_data$HRS_sample_size)) + 
    theme(text = element_text(size = 24), legend.position = "bottom")  + 
    guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                                nrow = 5, byrow = FALSE), 
           color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
           nrow = 5, byrow = FALSE) + ggtitle(str_to_sentence(race_subset))
  
  ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                           "figure5.XXC_impairement_class_RMSE_", race_subset, 
                           "_HCAP_adjudication.jpeg"), 
         dpi = 300, width = 17, height = 7.50, units = "in")
}

#---- Figure 4.XX + 5.XX: dementia prevalence ----
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
  group_by(prior_sample, HCAP_prop, HRS_sample_size, HCAP_sample_size) %>% 
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
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) 

plot_data[plot_data$LCI < 0, "LCI"] <- 0

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XX: no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = mean, y = HRS_sample_size)) + 
  geom_vline(aes(xintercept = prev), data = truth, color = "#ff0000", size = 2) +
  geom_point(size = 3, shape = 1) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1) +
  facet_grid(cols = vars(Race), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("Prevalence of dementia") + ylab("HRS sample size") + 
  scale_x_continuous(breaks = seq(0.0, 0.5, by = 0.10)) + 
  theme(text = element_text(size = 24))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XX_mean_CI_dem_prev_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.XX: HCAP adjudication ----
ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = prior_sample, 
           color = prior_sample)) + 
  geom_vline(aes(xintercept = prev), data = truth, size = 2) +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) +
  facet_grid(cols = vars(Race), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("Prevalence of dementia") + ylab("HRS sample size") + 
  scale_x_continuous(breaks = seq(0.0, 0.9, by = 0.10)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = TRUE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = TRUE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XX_mean_CI_dem_prev_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 17, height = 11, units = "in")

#---- Figure 4.XX + 5.XX: 95% CI coverage dementia prevalence ----
#---- **plot data ----
plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HRS_sample_size, HCAP_sample_size) %>% 
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
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) 

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XX: no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = dem, group = Race)) + 
  geom_line(aes(color = Race), size = 1.5) + 
  geom_point(aes(color = Race), size = 3, shape = 1) + 
  scale_color_manual(values = c("#006d9e", "#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Race/Ethnicity")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XX_dem_prev_coverage_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in") 

#---- **plot 5.XX: HCAP adjudication ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = dem, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Race)) + 
  guides(color = guide_legend(title = "Race/Ethnicity")) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XX_dem_prev_coverage_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 17, height = 8, units = "in") 

#---- Figure 4.XX + 5.XX: PR dementia ----
#---- **plot data ----
truth <- 
  data.frame("Comparison" = c("Black vs. White", "Hispanic vs. White"), 
             "PR" = c(na.omit(unique(results$true_dem_prev_black)/
                                unique(results$true_dem_prev_white)), 
                      na.omit(unique(results$true_dem_prev_hispanic)/
                                unique(results$true_dem_prev_white)))) %>% 
  mutate_at("PR", function(x) round(x, 2)) %>% 
  mutate_at("Comparison", 
            function(x) 
              factor(x, levels = c("Black vs. White", "Hispanic vs. White")))

plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(
    c("mean_PR_black", "mean_PR_hispanic", "LCI_PR_black", "LCI_PR_hispanic", 
      "UCI_PR_black", "UCI_PR_hispanic"), mean) %>% 
  pivot_longer(
    c("mean_PR_black", "mean_PR_hispanic", "LCI_PR_black", "LCI_PR_hispanic", 
      "UCI_PR_black", "UCI_PR_hispanic"),
    names_to = c(".value", "Race"), names_pattern = "(.*?)_(.*)") %>% 
  mutate_at("Race", function(x) str_remove(x, "PR_")) %>% 
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Race", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS")) 

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XX: no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = mean, y = HRS_sample_size)) + 
  geom_vline(aes(xintercept = PR), data = truth, color = "#ff0000", size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 3, shape = 1) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1) +
  facet_grid(cols = vars(Comparison), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("Prevalence Ratio (PR)") + ylab("HRS sample size") + 
  theme(text = element_text(size = 24))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XX_mean_CI_PR_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 13.5, height = 6.25, units = "in")

#---- **plot 5.XX: HCAP adjudication ----
ggplot(data = plot_data, 
       aes(y = HRS_sample_size, x = mean, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_vline(aes(xintercept = PR), data = truth, size = 2) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 3, position = position_dodge(-0.95)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.95)) +
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  theme_bw() + xlab("Prevalence Ratio (PR)") + ylab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XX_mean_CI_PR_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 17, height = 8, units = "in") 

#---- Figure 4.XX + 5.XX: 95% CI coverage PR ----
#---- **plot data ----
plot_data <- results %>% 
  group_by(prior_sample, HCAP_prop, HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(paste0("PR_coverage_", c("black", "hispanic")), mean) %>% 
  pivot_longer(paste0("PR_coverage_", c("black", "hispanic")),
               names_to = c(".value", "Race"), 
               names_pattern = "(.*?)_(.*)") %>%
  mutate_at("Race", function(x) str_remove(x, "coverage_")) %>% 
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "HCAP Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "HCAP Proportion\n50% of HRS"))

plot_data$prior_sample <- 
  factor(plot_data$prior_sample, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 35% SRS Adjudication", "HCAP 50% SRS Adjudication", 
                    "HCAP 35% Race-stratified SRS Adjudication", 
                    "HCAP 50% Race-stratified SRS Adjudication", 
                    "HCAP 35% SRS Adjudication + ADAMS", 
                    "HCAP 50% SRS Adjudication + ADAMS", 
                    "HCAP 35% Race-stratified SRS Adjudication + ADAMS", 
                    "HCAP 50% Race-stratified SRS Adjudication + ADAMS"))

#---- **plot 4.XX: no HCAP adjudication ----
ggplot(data = plot_data %>% filter(prior_sample == "ADAMS"), 
       aes(x = HRS_sample_size, y = PR, group = Comparison, color = Comparison)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3, shape = 1) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "figure4.XX_PR_coverage_no_HCAP_adjudiation.jpeg"), 
       dpi = 300, width = 13.25, height = 6.25, units = "in")

#---- **plot 5.XX: HCAP adjudication ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = PR, group = prior_sample, 
           color = prior_sample, shape = prior_sample)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) + 
  scale_shape_manual(values = c(1, rep(19, 4), rep(1, 4))) +
  scale_color_manual(values = c("black", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f", 
                                "#288fb4", "#1d556f", 
                                "#f35f5f", "#cc435f")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(rows = vars(HCAP_prop), cols = vars(Comparison)) + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Adjudicated\nPrior\nSample", 
                              nrow = 5, byrow = FALSE), 
         color = guide_legend(title = "Adjudicated\nPrior\nSample"), 
         nrow = 5, byrow = FALSE)

ggsave(filename = paste0(path_to_box, "figures/chapter_5/simulation_study/", 
                         "figure5.XX_PR_coverage_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 17, height = 8, units = "in") 

#---- OLD ----
# #---- **combined scenario ----
# impairment_by_race <- expand_grid(c("mean", "LCI", "UCI"), 
#                                   c("Unimpaired", "MCI", "Dementia", "Other"), 
#                                   c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# dem_prev_by_race <- 
#   expand_grid(c("mean_dem_prev", "LCI_dem_prev", "UCI_dem_prev"), 
#               c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# PR <- 
#   expand_grid(c("mean_PR", "LCI_PR", "UCI_PR"), 
#               c("black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# cols_by_race <- c(impairment_by_race, dem_prev_by_race, PR)
# 
# ADAMS_results <- results %>% filter(prior_sample == "ADAMS") %>% 
#   dplyr::select("dataset_name", "prior_sample", "HRS_sample_size", 
#                 "HCAP_sample_size", "HCAP_prop", 
#                 paste0("mean_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#                 paste0("LCI_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#                 paste0("UCI_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#                 all_of(cols_by_race)) %>%
#   set_colnames(c(
#     "dataset_name", "prior_sample", "HRS_sample_size", "HCAP_sample_size", 
#     "HCAP_prop", 
#     paste0("ADAMS_mean_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#     paste0("ADAMS_LCI_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#     paste0("ADAMS_UCI_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#     paste0("ADAMS_", all_of(cols_by_race)))) %>% 
#   mutate_at("dataset_name", function(x) 
#     str_split(x, pattern = "_ADAMS_")[[1]][1])
# 
# #---- ****calculate ADAMS SEs ----
# #impairment classes
# for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
#   for(race in c("overall", "white", "black", "hispanic")){
#     if(race == "overall"){
#       ADAMS_results %<>%
#         mutate(!!sym(paste0("ADAMS_SE_", class)) :=
#                  rowMeans(cbind(abs(!!sym(paste0("ADAMS_mean_", class)) -
#                                       !!sym(paste0("ADAMS_LCI_", class))),
#                                 abs(!!sym(paste0("ADAMS_mean_", class)) -
#                                       !!sym(paste0("ADAMS_UCI_", class)))))/1.96)
#     } else{
#       ADAMS_results %<>%
#         mutate(!!sym(paste0("ADAMS_SE_", class, "_", race)) :=
#                  rowMeans(cbind(
#                    abs(!!sym(paste0("ADAMS_mean_", class, "_", race)) -
#                          !!sym(paste0("ADAMS_LCI_", class, "_", race))),
#                    abs(!!sym(paste0("ADAMS_mean_", class, "_", race)) -
#                          !!sym(paste0("ADAMS_UCI_", class, "_", race)))))/1.96)
#     }
#   }
# }
# 
# #dem prev + PR
# for(measure in c("dem_prev", "PR")){
#   for(race in c("white", "black", "hispanic")){
#     if(measure == "PR" & race == "white"){next}
#     
#     ADAMS_results %<>%
#       mutate(!!sym(paste0("ADAMS_SE_", measure, "_", race)) :=
#                rowMeans(cbind(
#                  abs(!!sym(paste0("ADAMS_mean_", measure, "_", race)) -
#                        !!sym(paste0("ADAMS_LCI_", measure, "_", race))),
#                  abs(!!sym(paste0("ADAMS_mean_", measure, "_", race)) -
#                        !!sym(paste0("ADAMS_UCI_", measure, "_", race)))))/1.96)
#     
#   }
# }
# 
# combined_results <- results %>% filter(!prior_sample %in% c("ADAMS")) %>%
#   mutate_at("prior_sample", function(x) paste0(x, " + ADAMS")) %>% 
#   mutate_at("dataset_name", function(x)
#     str_split(x, pattern = "_calibration_")[[1]][1]) 
# 
# #---- ****calculate other scenario SEs ----
# #impairment classes
# for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
#   for(race in c("overall", "white", "black", "hispanic")){
#     if(race == "overall"){
#       combined_results %<>%
#         mutate(!!sym(paste0("SE_", class)) :=
#                  rowMeans(cbind(abs(!!sym(paste0("mean_", class)) -
#                                       !!sym(paste0("LCI_", class))),
#                                 abs(!!sym(paste0("mean_", class)) -
#                                       !!sym(paste0("UCI_", class)))))/1.96)
#     } else{
#       combined_results %<>%
#         mutate(!!sym(paste0("SE_", class, "_", race)) :=
#                  rowMeans(cbind(
#                    abs(!!sym(paste0("mean_", class, "_", race)) -
#                          !!sym(paste0("LCI_", class, "_", race))),
#                    abs(!!sym(paste0("mean_", class, "_", race)) -
#                          !!sym(paste0("UCI_", class, "_", race)))))/1.96)
#     }
#   }
# }
# 
# #dem prev + PR
# for(measure in c("dem_prev", "PR")){
#   for(race in c("white", "black", "hispanic")){
#     if(measure == "PR" & race == "white"){next}
#     
#     combined_results %<>%
#       mutate(!!sym(paste0("SE_", measure, "_", race)) :=
#                rowMeans(cbind(
#                  abs(!!sym(paste0("mean_", measure, "_", race)) -
#                        !!sym(paste0("LCI_", measure, "_", race))),
#                  abs(!!sym(paste0("mean_", measure, "_", race)) -
#                        !!sym(paste0("UCI_", measure, "_", race)))))/1.96)
#     
#   }
# }
# 
# 
# race_dem_prev_placeholders <- 
#   expand_grid(c("ADAMS_mean_dem_prev", "ADAMS_SE_dem_prev"), 
#               c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# PR_placeholders <- 
#   expand_grid(c("ADAMS_mean_PR", "ADAMS_SE_PR"), 
#               c("black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# race_placeholders <- c(race_impairment_placeholders, race_dem_prev_placeholders, 
#                        PR_placeholders)
# 
# combined_results[, c(paste0("ADAMS_mean_", 
#                             c("Unimpaired", "MCI", "Dementia", "Other")), 
#                      paste0("ADAMS_SE_", 
#                             c("Unimpaired", "MCI", "Dementia", "Other")), 
#                      all_of(race_placeholders))] <- NA
# 
# 
# #impairment classes
# for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
#   for(race in c("overall", "white", "black", "hispanic")){
#     if(race == "overall"){
#       #---- ****mean estimates ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_mean_", class)) := 
#                  rowMeans(cbind(!!sym(paste0("mean_", class)), 
#                                 !!sym(paste0("ADAMS_mean_", class)))))
#       
#       #---- ****within variance ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("within_var_", class)) := 
#                  rowMeans(cbind(!!sym(paste0("SE_", class)), 
#                                 !!sym(paste0("ADAMS_SE_", class)))))
#       
#       #---- ****between variance ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("between_var_", class)) := 
#                  apply(cbind(!!sym(paste0("mean_", class)), 
#                              !!sym(paste0("ADAMS_mean_", class))), 1, var))
#       
#       #---- ****LCIs ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_LCI_", class)) := 
#                  !!sym(paste0("combined_mean_", class)) - 
#                  1.96*sqrt(!!sym(paste0("within_var_", class)) + 
#                              !!sym(paste0("between_var_", class)) + 
#                              !!sym(paste0("between_var_", class))/2))
#       
#       #---- ****UCIs ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_UCI_", class)) := 
#                  !!sym(paste0("combined_mean_", class)) + 
#                  1.96*sqrt(!!sym(paste0("within_var_", class)) + 
#                              !!sym(paste0("between_var_", class)) + 
#                              !!sym(paste0("between_var_", class))/2))
#       
#       #---- ****coverage ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_", class, "_coverage")) := 
#                  (!!sym(paste0("true_", class)) >= 
#                     !!sym(paste0("combined_LCI_", class)))* 
#                  (!!sym(paste0("true_", class)) <= 
#                     !!sym(paste0("combined_UCI_", class))))
#     } else{
#       #---- ****mean estimates ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_mean_", class, "_", race)) := 
#                  rowMeans(cbind(
#                    !!sym(paste0("mean_", class, "_", race)), 
#                    !!sym(paste0("ADAMS_mean_", class, "_", race)))))
#       
#       #---- ****within variance ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("within_var_", class, "_", race)) := 
#                  rowMeans(cbind(!!sym(paste0("SE_", class, "_", race)), 
#                                 !!sym(paste0("ADAMS_SE_", class, "_", race)))))
#       
#       #---- ****between variance ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("between_var_", class, "_", race)) := 
#                  apply(cbind(!!sym(paste0("mean_", class, "_", race)), 
#                              !!sym(paste0("ADAMS_mean_", class, "_", race))), 1, 
#                        var))
#       
#       #---- ****LCIs ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_LCI_", class, "_", race)) := 
#                  !!sym(paste0("combined_mean_", class, "_", race)) - 
#                  1.96*sqrt(!!sym(paste0("within_var_", class, "_", race)) + 
#                              !!sym(paste0("between_var_", class, "_", race)) + 
#                              !!sym(paste0("between_var_", class, "_", race))/2))
#       
#       #---- ****UCIs ----
#       combined_results %<>% 
#         mutate(!!sym(paste0("combined_UCI_", class, "_", race)) := 
#                  !!sym(paste0("combined_mean_", class, "_", race)) + 
#                  1.96*sqrt(!!sym(paste0("within_var_", class, "_", race)) + 
#                              !!sym(paste0("between_var_", class, "_", race)) + 
#                              !!sym(paste0("between_var_", class, "_", race))/2))
#       
#       #---- ****coverage ----
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
# #dem prev + PR
# for(measure in c("dem_prev", "PR")){
#   for(race in c("white", "black", "hispanic")){
#     if(measure == "PR" & race == "white"){next}
#     #---- ****mean estimates ----
#     combined_results %<>% 
#       mutate(!!sym(paste0("combined_mean_", measure, "_", race)) := 
#                rowMeans(cbind(
#                  !!sym(paste0("mean_", measure, "_", race)), 
#                  !!sym(paste0("ADAMS_mean_", measure, "_", race)))))
#     
#     #---- ****within variance ----
#     combined_results %<>% 
#       mutate(!!sym(paste0("within_var_", measure, "_", race)) := 
#                rowMeans(cbind(!!sym(paste0("SE_", measure, "_", race)), 
#                               !!sym(paste0("ADAMS_SE_", measure, "_", race)))))
#     
#     #---- ****between variance ----
#     combined_results %<>% 
#       mutate(!!sym(paste0("between_var_", measure, "_", race)) := 
#                apply(cbind(!!sym(paste0("mean_", measure, "_", race)), 
#                            !!sym(paste0("ADAMS_mean_", measure, "_", race))), 1, 
#                      var))
#     
#     #---- ****LCIs ----
#     combined_results %<>% 
#       mutate(!!sym(paste0("combined_LCI_", measure, "_", race)) := 
#                !!sym(paste0("combined_mean_", measure, "_", race)) - 
#                1.96*sqrt(!!sym(paste0("within_var_", measure, "_", race)) + 
#                            !!sym(paste0("between_var_", measure, "_", race)) + 
#                            !!sym(paste0("between_var_", measure, "_", race))/2))
#     
#     #---- ****UCIs ----
#     combined_results %<>% 
#       mutate(!!sym(paste0("combined_UCI_", measure, "_", race)) := 
#                !!sym(paste0("combined_mean_", measure, "_", race)) + 
#                1.96*sqrt(!!sym(paste0("within_var_", measure, "_", race)) + 
#                            !!sym(paste0("between_var_", measure, "_", race)) + 
#                            !!sym(paste0("between_var_", measure, "_", race))/2))
#     
#     #---- ****coverage ----
#     combined_results %<>% 
#       mutate(!!sym(paste0("combined_", measure, "_coverage_", race)) := 
#                (!!sym(paste0("true_", measure, "_", race)) >= 
#                   !!sym(paste0("combined_LCI_", measure, "_", race)))* 
#                (!!sym(paste0("true_", measure, "_", race)) <= 
#                   !!sym(paste0("combined_UCI_", measure, "_", race))))
#   }
# }
# 
# #---- ****select relevant columns ----
# impairment_class_measures_by_race <- 
#   expand_grid(
#     c("true", "combined_mean", "combined_LCI", "combined_UCI"), 
#     c("Unimpaired", "MCI", "Dementia", "Other"), 
#     c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# dem_prev_measures_by_race <- 
#   expand_grid(
#     c("combined_mean_dem_prev", "combined_LCI_dem_prev", 
#       "combined_UCI_dem_prev"), 
#     c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# PR_measures <- 
#   expand_grid(
#     c("combined_mean_PR", "combined_LCI_PR", "combined_UCI_PR"), 
#     c("black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# impairment_coverage_by_race <- 
#   expand_grid(c("combined"), 
#               c("Unimpaired", "MCI", "Dementia", "Other"), 
#               c("coverage"), 
#               c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# dem_prev_coverage_by_race <- expand_grid(c("combined_dem_prev_coverage"), 
#                                          c("white", "black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# PR_coverage <- expand_grid(c("combined_PR_coverage"), 
#                            c("black", "hispanic")) %>% 
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# cols_by_race <- c(impairment_class_measures_by_race, dem_prev_measures_by_race,
#                   PR_measures,
#                   impairment_coverage_by_race, dem_prev_coverage_by_race, 
#                   PR_coverage)
# 
# combined_results %<>% 
#   dplyr::select("dataset_name", "HRS_sample_size", "HCAP_sample_size", 
#                 "HCAP_prop", "prior_sample", 
#                 paste0("true_", c("Unimpaired", "MCI", "Dementia", "Other")),
#                 paste0("combined_mean_", 
#                        c("Unimpaired", "MCI", "Dementia", "Other")),
#                 paste0("combined_LCI_", 
#                        c("Unimpaired", "MCI", "Dementia", "Other")), 
#                 paste0("combined_UCI_", 
#                        c("Unimpaired", "MCI", "Dementia", "Other")), 
#                 paste0("combined_", c("Unimpaired", "MCI", "Dementia", "Other"), 
#                        "_coverage"), all_of(cols_by_race)) %>% 
#   set_colnames(c("dataset_name", "HRS_sample_size", "HCAP_sample_size", 
#                  "HCAP_prop", "prior_sample", 
#                  paste0("true_", c("Unimpaired", "MCI", "Dementia", "Other")),
#                  paste0("mean_", c("Unimpaired", "MCI", "Dementia", "Other")),
#                  paste0("LCI_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#                  paste0("UCI_", c("Unimpaired", "MCI", "Dementia", "Other")), 
#                  paste0(c("Unimpaired", "MCI", "Dementia", "Other"), 
#                         "_coverage"), 
#                  str_remove(all_of(cols_by_race), "combined_")))
# 
# #---- ****bind with original results ----
# results %<>% plyr::rbind.fill(., combined_results)


#---- Figure X: comparing real HRS with synthetic HRS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))

#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

#---- ****ADAMS imputed data ----
ADAMS_imputed_clean <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets_cleaned"))

#stack data and rename variables
ADAMS_imputed_stacked <- do.call(rbind, ADAMS_imputed_clean) %>% 
  rename_at(vars(variable_labels$ADAMS), ~ variable_labels$data_label)

#---- **ADAMS plots ----
plot_data <- ADAMS_imputed_stacked[1:nrow(ADAMS_imputed_clean[[1]]), ] %>% 
  dplyr::select(c(all_of(continuous_vars), 
                  "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  pivot_longer(cols = c("Unimpaired", "MCI", "Dementia", "Other"), 
               names_to = "Group") %>% filter(value == 1) %>% 
  rename_at(vars(plot_labels$data_label), ~ plot_labels$figure_label) %>% 
  dplyr::select(-c("value")) %>% 
  pivot_longer(-c("Group"), names_to = "Variable") %>% 
  left_join(color_palette) %>% 
  mutate("order" = case_when(Group == "Unimpaired" ~ 1, 
                             Group == "MCI" ~ 2, 
                             Group == "Dementia" ~ 3, 
                             Group == "Other" ~ 4))

plot_data$Color <- reorder(plot_data$Color, plot_data$order)
plot_data$Group <- reorder(plot_data$Group, plot_data$order)

ggplot(data = plot_data, aes(x = value, color = Color, fill = Color)) + 
  geom_density(alpha = 0.5) + theme_minimal() + 
  theme(text = element_text(size = 10)) + 
  scale_color_identity(guide = "legend", labels = levels(plot_data$Group)) + 
  scale_fill_identity(guide = "legend", labels = levels(plot_data$Group)) +
  theme(legend.position = "bottom") + xlab("") + 
  facet_wrap(vars(Variable), scales = "free", ncol = 4) +
  guides(fill = guide_legend(title = "Group")) +
  guides(color = guide_legend(title = "Group"))

ggsave(filename = "ADAMS_mix_Z.jpeg", plot = last_plot(), 
       path = paste0(path_to_box, "figures/simulation_study/"), 
       width = 8, height = 10, units = "in", device = "jpeg")

