#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "devtools", "here")
install_github("thomasp85/patchwork")

#---- source functions ----
source(here::here("functions", "read_results.R"))

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

#---- ****synthetic HRS ----
synthetic_normal_1000 <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/synthetic_data/", 
                  "ADAMS_props/HRS/HRS_synthetic_normal_1000_dementia.csv"))

synthetic_lognormal_1000 <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/synthetic_data/", 
                  "ADAMS_props/HRS/HRS_synthetic_lognormal_1000_dementia.csv"))

synthetic_bathtub_1000 <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/synthetic_data/", 
                  "ADAMS_props/HRS/HRS_synthetic_bathtub_1000_dementia.csv"))

#---- **define plot variables ----
continuous_vars <- colnames(synthetic_normal_1000)[
  str_detect(colnames(synthetic_normal_1000), "_Z")] 

plot_labels <- variable_labels %>% filter(data_label %in% continuous_vars) 

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

#---- **synthetic HRS plots ----
#---- ****normal ----
plot_data <- synthetic_normal_1000 %>% 
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

ggsave(filename = "synthetic_normal_mix_Z.jpeg", plot = last_plot(), 
       path = paste0(path_to_box, "figures/simulation_study/"), 
       width = 8, height = 10, units = "in", device = "jpeg")

#---- ****lognormal ----
plot_data <- synthetic_lognormal_1000 %>% 
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

ggsave(filename = "synthetic_lognormal_mix_Z.jpeg", plot = last_plot(), 
       path = paste0(path_to_box, "figures/simulation_study/"), 
       width = 8, height = 10, units = "in", device = "jpeg")

#---- ****bathtub ----
plot_data <- synthetic_bathtub_1000 %>% 
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

ggsave(filename = "synthetic_bathtub_mix_Z.jpeg", plot = last_plot(), 
       path = paste0(path_to_box, "figures/simulation_study/"), 
       width = 8, height = 10, units = "in", device = "jpeg")

#---- check number of simulation runs ----
#Missing runs:
# n = 2000, sample = 25, calibration_50_SRS: 556
# n = 2000, sample = 50, calibration_50_SRS: 27
# n = 4000, sample = 25, calibration 50 SRS: 19

#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

results <- do.call(rbind, lapply(results_paths, read_results)) %>% 
  group_by(dataset_name) %>% slice_head(n = 1000)

table(results$dataset_name, useNA = "ifany")

#---- Figure X: mean and 95% CI impairment class counts ----
#---- **truth ----
superpop_impairment_props <- 
  read_csv(paste0(path_to_box, 
                  "data/superpopulations/impairment_class_props.csv"))
superpop_impairment_props$Group <- 
  factor(superpop_impairment_props$Group, 
         levels = c("Unimpaired", "MCI", "Dementia", "Other"))

#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
#---- ****extra calcs ----
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

results %<>% mutate(sample_size = HCAP_prop/100*sample_size)

results %<>% 
  mutate(mean_Unimpaired_prop = mean_Unimpaired/sample_size, 
         mean_MCI_prop = mean_MCI/sample_size, 
         mean_Dementia_prop = mean_Dementia/sample_size, 
         mean_Other_prop = mean_Other/sample_size) %>% 
  mutate(LCI_Unimpaired_prop = LCI_Unimpaired/sample_size, 
         LCI_MCI_prop = LCI_MCI/sample_size, 
         LCI_Dementia_prop = LCI_Dementia/sample_size, 
         LCI_Other_prop = LCI_Other/sample_size) %>%
  mutate(UCI_Unimpaired_prop = UCI_Unimpaired/sample_size, 
         UCI_MCI_prop = UCI_MCI/sample_size, 
         UCI_Dementia_prop = UCI_Dementia/sample_size, 
         UCI_Other_prop = UCI_Other/sample_size)

plot_data <- results %>% 
  group_by(calibration, calibration_sampling, HCAP_prop, sample_size) %>% 
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
  mutate_at("sample_size", as.factor) %>% 
  mutate_at("Group", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "Sample Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "Sample Proportion\n50% of HRS")) %>% 
  mutate("calibration_sampling" = 
           case_when(calibration_sampling == "NA" ~ "ADAMS", 
                     calibration_sampling == "SRS" & 
                       calibration == "calibration_50" ~ 
                       "HCAP 50% SRS Adjudication"))

#---- **plot v1: no HCAP calibration ----
ggplot(data = plot_data %>% filter(calibration == "no_calibration"), 
       aes(x = mean, y = sample_size)) +
  geom_vline(data = superpop_impairment_props, aes(xintercept = prop), 
             size = 2, color = rep(color_palette$Color, 2)) +
  geom_point(size = 3) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.3, size = 1) + theme_bw() + 
  facet_grid(cols = vars(Group), rows = vars(HCAP_prop), scales = "free_y") + 
  scale_x_continuous(breaks = seq(0.10, 0.40, by = 0.05)) +
  xlab("Impairment class proportion") + ylab("HCAP sample size") + 
  theme(text = element_text(size = 24))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "mean_CI_impairement_class_no_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 14.75, height = 6.5, units = "in")

#---- **plot v2: HCAP calibration ----
#add data for 100% adjudication
HCAP_all_adjudicated <- results %>% ungroup() %>%
  dplyr::select(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other", 
                  "sample_size", "HCAP_prop")) %>% 
  group_by(HCAP_prop, sample_size) %>% 
  summarise_at(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
               mean) %>% 
  pivot_longer(c("true_Unimpaired", "true_MCI", "true_Dementia", "true_Other"), 
               names_to = "Group", values_to = "mean") %>% 
  mutate_at("Group", function(x) str_remove(x, "true_")) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "Sample Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "Sample Proportion\n50% of HRS")) %>% 
  mutate(mean = mean/sample_size) %>%
  mutate("calibration" = "calibration_100", 
         "calibration_sampling" = "HCAP 100% Adjudication",
         "LCI" = NA, "UCI" = NA) %>% mutate_at("sample_size", as.factor)

plot_data %<>% rbind(., HCAP_all_adjudicated)
plot_data$calibration_sampling <- 
  factor(plot_data$calibration_sampling, 
         levels = c("HCAP 100% Adjudication", "ADAMS", 
                    "HCAP 50% SRS Adjudication"))
plot_data$Group <- 
  factor(plot_data$Group, levels = c("Unimpaired", "MCI", "Dementia", "Other"))

ggplot(data = plot_data, 
       aes(x = mean, y = sample_size, shape = calibration_sampling)) +
  geom_vline(data = superpop_impairment_props, aes(xintercept = prop), 
             size = 2, color = rep(color_palette$Color, 2)) +
  geom_point(size = 4, position = position_dodge(-0.6)) + 
  scale_shape_manual(values = c(18, 1, 19)) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.4, size = 1, 
                position = position_dodge(-0.6)) + theme_bw() + 
  facet_grid(cols = vars(Group), rows = vars(HCAP_prop), scales = "free_y") + 
  scale_x_continuous(breaks = seq(0.00, 0.70, by = 0.10)) +
  xlab("Impairment class proportion") + ylab("HCAP sample size") + 
  theme(text = element_text(size = 24), legend.position = "bottom") + 
  guides(shape = guide_legend(title = "Adjudicated Sample for Prior"))

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "mean_CI_impairement_class_HCAP_adjudication.jpeg"), 
       dpi = 300, width = 14.75, height = 6.5, units = "in")

#---- Figure X: 95% CI coverage impairment classes ----
#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% 
  group_by(calibration, calibration_sampling, HCAP_prop, sample_size) %>%
  summarise_at(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"), 
               mean) %>% 
  pivot_longer(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"),
               names_to = c("class", "coverage"), names_sep = "_") %>% 
  mutate_at("sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ 
                                   "Sample Proportion\n25% of HRS", 
                                 HCAP_prop == 50 ~ 
                                   "Sample Proportion\n50% of HRS")) %>% 
  mutate("calibration_sampling" = 
           case_when(calibration_sampling == "NA" ~ "ADAMS", 
                     calibration_sampling == "SRS" & 
                       calibration == "calibration_50" ~ 
                       "HCAP 50% SRS Adjudication"))

#---- **plot v1: no HCAP calibration ----
ggplot(data = plot_data %>% filter(calibration == "no_calibration"), 
       aes(x = sample_size, y = value, group = class)) + 
  geom_line(aes(color = class), size = 1.5) + 
  geom_point(aes(color = class), size = 3) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HCAP Sample Size") +
  facet_grid(cols = vars(HCAP_prop), scales = "free_x") + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "impairement_class_coverage.jpeg"))

#---- **plot v2: HCAP calibration ----
ggplot(data = plot_data, aes(x = sample_size, y = value, group = class)) + 
  geom_line(aes(color = class), size = 1.5) + 
  geom_point(aes(color = class), size = 3) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HCAP Sample Size") +
  facet_grid(cols = vars(HCAP_prop), rows = vars(calibration), 
             scales = "free_x") + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 24), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "impairement_class_coverage.jpeg"))

#----- Figure X: bias impairment class counts ----
#---- **% bias ----
#---- ****overall ----
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  results %<>% 
    mutate(!!paste0("percent_increase_", class) := 
             !!sym(paste0("bias_", class))/!!sym(paste0("true_", class))*100)
}

#---- ****by race ----
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  for(race in c("white", "black", "hispanic")){
    results %<>% 
      mutate(!!paste0("percent_increase_", class, "_", race) := 
               !!sym(paste0("bias_", class, "_", race))/
               !!sym(paste0("true_", class, "_", race))*100)
  }
}

#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data: counts ----
plot_data <- results %>% 
  dplyr::select("calibration", "HCAP_prop", "sample_size", 
                paste0("percent_increase_", 
                       c("Unimpaired", "MCI", "Dementia", "Other"))) %>%
  pivot_longer(paste0("percent_increase_", 
                      c("Unimpaired", "MCI", "Dementia", "Other")),
               names_to = c("text", "class"), 
               names_sep = "_increase_") %>% 
  mutate_at("class", function(x) str_remove(x, pattern = "_prop")) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ "25% of HRS", 
                                 HCAP_prop == 50 ~ "50% of HRS"))

#---- **plot: count ----
ggplot(data = plot_data, aes(x = sample_size, y = value, color = class)) + 
  geom_boxplot(size = 1) +
  scale_color_manual(values = group_colors) + 
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent Bias") + 
  facet_grid(cols = vars(HCAP_prop), scales = "free_x") + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HCAP Sample Size", 
                   breaks = unique(plot_data$sample_size)) + 
  theme(text = element_text(size = 18), legend.position = "bottom")       

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "impairement_class_percent_bias_counts.jpeg"))

#---- **plot data: counts x race ----
cols_by_race <- expand_grid("percent_increase", 
                            c("Unimpaired", "MCI", "Dementia", "Other"), 
                            c("white", "black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% 
  dplyr::select("calibration", "HCAP_prop", "sample_size",  
                all_of(cols_by_race)) %>%
  mutate_at("sample_size", as.factor) %>%
  pivot_longer(all_of(cols_by_race), names_to = c("text", "class_race"), 
               names_sep = "_increase_") %>% 
  separate(col = "class_race", into = c("class", "race"), sep = "_") %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ "25% of HRS", 
                                 HCAP_prop == 50 ~ "50% of HRS"))

ggplot(data = plot_data, aes(x = sample_size, y = value, color = class)) + 
  geom_boxplot(size = 0.60) +
  scale_color_manual(values = group_colors) + 
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent Bias") + 
  facet_grid(rows = vars(race), cols = vars(HCAP_prop), scales = "free") + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HCAP Sample Size", 
                   breaks = unique(plot_data$sample_size)) + 
  theme(text = element_text(size = 18), legend.position = "bottom")    

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "impairement_class_percent_bias_count_by_race.jpeg"))

#---- Figure X: dementia prevalence ----
#---- **plot data ----
truth <- data.frame("Race" = c("White", "Black", "Hispanic"), 
                    "prev" = c(unique(results$true_dem_prev_white), 
                               unique(results$true_dem_prev_black), 
                               unique(results$true_dem_prev_hispanic))) %>% 
  mutate_at("prev", function(x) round(x, 2)) %>% 
  mutate_at("Race", 
            function(x) factor(x, levels = c("White", "Black", "Hispanic")))

plot_data <- results %>% 
  group_by(calibration, HCAP_prop, sample_size) %>% 
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
  mutate("sample_size" = 1/(HCAP_prop/100)*sample_size) %>%
  mutate_at("sample_size", as.factor) %>% 
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ "25% of HRS", 
                                 HCAP_prop == 50 ~ "50% of HRS"))

#---- **plot ----
ggplot(data = plot_data, aes(x = mean, y = sample_size)) + 
  geom_vline(aes(xintercept = prev), data = truth, color = "#ff0000", size = 1) +
  geom_point(size = 2) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2) +
  facet_grid(cols = vars(Race), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("Dementia prevalence") + ylab("HRS sample size") + 
  scale_x_continuous(breaks = seq(0.0, 0.5, by = 0.10)) + 
  theme(text = element_text(size = 18))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "mean_CI_dem_prev.jpeg"))

#---- Figure X: 95% CI coverage dementia prevalence ----
#---- **plot data ----
plot_data <- results %>% 
  group_by(calibration, HCAP_prop, sample_size) %>% 
  summarise_at(paste0("dem_prev_coverage_", c("white", "black", "hispanic")), 
               mean) %>% 
  pivot_longer(paste0("dem_prev_coverage_", c("white", "black", "hispanic")),
               names_to = c(".value", "Race"), 
               names_pattern = "(.*?)_(.*)") %>%
  mutate_at("Race", function(x) str_remove(x, "prev_coverage_")) %>% 
  mutate_at("Race", str_to_sentence) %>%
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>%
  mutate("sample_size" = 1/(HCAP_prop/100)*sample_size) %>%
  mutate_at("sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ "25% of HRS", 
                                 HCAP_prop == 50 ~ "50% of HRS"))

#---- **plot ----
ggplot(data = plot_data, aes(x = sample_size, y = dem, group = Race)) + 
  geom_line(aes(color = Race), size = 1) + 
  geom_point(aes(color = Race), size = 2) + 
  scale_color_manual(values = c("#006d9e", "#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(cols = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Race/Ethnicity")) + 
  theme(text = element_text(size = 18), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "dem_prev_coverage.jpeg")) 

#---- Figure X: RR dementia ----
#---- **plot data ----
truth <- data.frame("Comparison" = c("Black vs. White", "Hispanic vs. White"), 
                    "PR" = c(unique(results$true_dem_prev_black)/
                               unique(results$true_dem_prev_white), 
                             unique(results$true_dem_prev_hispanic)/
                               unique(results$true_dem_prev_white))) %>% 
  mutate_at("PR", function(x) round(x, 2)) %>% 
  mutate_at("Comparison", 
            function(x) 
              factor(x, levels = c("Black vs. White", "Hispanic vs. White")))

plot_data <- results %>% 
  group_by(calibration, HCAP_prop, sample_size) %>% 
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
  mutate("sample_size" = 1/(HCAP_prop/100)*sample_size) %>%
  mutate_at("sample_size", as.factor) %>% 
  mutate_at("Race", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ "25% of HRS", 
                                 HCAP_prop == 50 ~ "50% of HRS"))

#---- **plot ----
ggplot(data = plot_data, aes(x = mean, y = sample_size)) + 
  geom_vline(aes(xintercept = PR), data = truth, color = "#ff0000", size = 1) +
  geom_vline(aes(xintercept = 1), color = "black", lty = "dashed") +
  geom_point(size = 2) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2) +
  facet_grid(cols = vars(Comparison), rows = vars(HCAP_prop)) + theme_bw() + 
  xlab("PR") + ylab("HRS sample size") + 
  theme(text = element_text(size = 18))  

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "mean_CI_PR.jpeg"))

#---- Figure X: 95% CI coverage RR ----
#---- **plot data ----
plot_data <- results %>% 
  group_by(calibration, HCAP_prop, sample_size) %>% 
  summarise_at(paste0("PR_coverage_", c("black", "hispanic")), mean) %>% 
  pivot_longer(paste0("PR_coverage_", c("black", "hispanic")),
               names_to = c(".value", "Race"), 
               names_pattern = "(.*?)_(.*)") %>%
  mutate_at("Race", function(x) str_remove(x, "coverage_")) %>% 
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate("sample_size" = 1/(HCAP_prop/100)*sample_size) %>%
  mutate_at("sample_size", as.factor) %>% 
  mutate("HCAP_prop" = case_when(HCAP_prop == 25 ~ "25% of HRS", 
                                 HCAP_prop == 50 ~ "50% of HRS"))

#---- **plot ----
ggplot(data = plot_data, aes(x = sample_size, y = PR, group = Comparison)) + 
  geom_line(aes(color = Comparison), size = 1) + 
  geom_point(aes(color = Comparison), size = 2) + 
  scale_color_manual(values = c("#e09a3b", "#9ddfdf")) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HRS Sample Size") +
  facet_grid(cols = vars(HCAP_prop)) + 
  guides(color = guide_legend(title = "Comparison")) + 
  theme(text = element_text(size = 18), legend.position = "bottom")      

ggsave(filename = paste0(path_to_box, "figures/chapter_4/simulation_study/", 
                         "PR_coverage.jpeg"))
