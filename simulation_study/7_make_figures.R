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
# n = 2000: 1
# n = 4000: 64
# n = 8000: 274
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

results <- do.call(rbind, lapply(results_paths, read_results)) 

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
#---- ****extra calcs: delete later ----
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

results %<>% mutate(sample_size = 0.5*sample_size)

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

#---- ****end extra calcs: delete later ----
plot_data <- results %>% 
  group_by(calibration, sample_size) %>% 
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
            factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other")))

#---- **plot ----
ggplot(data = plot_data, aes(x = mean, y = sample_size)) +
  geom_vline(data = superpop_impairment_props, aes(xintercept = prop), 
             size = 1, color = color_palette$Color) +
  geom_point(size = 2) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2) + theme_bw() + 
  facet_grid(cols = vars(Group)) + 
  scale_x_continuous(breaks = seq(0.10, 0.40, by = 0.05)) +
  xlab("Proportion") + ylab("HCAP sample size") + 
  theme(text = element_text(size = 18))  

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "mean_CI_impairement_class.jpeg"))

#---- Figure X: 95% CI coverage impairment classes ----
#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% 
  group_by(calibration, sample_size) %>% 
  summarise_at(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"), 
               mean) %>% 
  pivot_longer(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"),
               names_to = c("class", "coverage"), names_sep = "_") %>% 
  mutate_at("sample_size", as.factor)

#---- **plot ----
ggplot(data = plot_data, aes(x = sample_size, y = value, group = class)) + 
  geom_line(aes(color = class), size = 1) + 
  geom_point(aes(color = class), size = 2) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + xlab("HCAP Sample Size") +
  #facet_grid(rows = vars(calibration)) + 
  guides(color = guide_legend(title = "Group")) + 
  theme(text = element_text(size = 18))       

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
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
  dplyr::select("calibration", "sample_size", 
                paste0("percent_increase_", 
                       c("Unimpaired", "MCI", "Dementia", "Other"))) %>%
  mutate_at("sample_size", as.numeric) %>% 
  mutate("sample_size" = 0.5*sample_size) %>% 
  mutate_at("sample_size", as.factor) %>%
  pivot_longer(paste0("percent_increase_", 
                      c("Unimpaired", "MCI", "Dementia", "Other")),
               names_to = c("text", "class"), 
               names_sep = "_increase_") %>% 
  mutate_at("class", function(x) str_remove(x, pattern = "_prop")) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other")))

#---- **plot: count ----
ggplot(data = plot_data, aes(x = sample_size, y = value, color = class)) + 
  geom_boxplot(size = 1) +
  scale_color_manual(values = group_colors) + 
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent Bias") + 
  #facet_grid(rows = vars(calibration)) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HCAP Sample Size", 
                   breaks = unique(plot_data$sample_size)) + 
  theme(text = element_text(size = 18))       

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "impairement_class_percent_bias_counts.jpeg"))

#---- **plot data: counts x race ----
cols_by_race <- expand_grid("percent_increase", 
                            c("Unimpaired", "MCI", "Dementia", "Other"), 
                            c("white", "black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% 
  dplyr::select("calibration", "sample_size",  
                all_of(cols_by_race)) %>%
  mutate_at("sample_size", as.numeric) %>% 
  mutate("sample_size" = 0.5*sample_size) %>% 
  mutate_at("sample_size", as.factor) %>%
  pivot_longer(all_of(cols_by_race), names_to = c("text", "class_race"), 
               names_sep = "_increase_") %>% 
  separate(col = "class_race", into = c("class", "race"), sep = "_") %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic")))

ggplot(data = plot_data, aes(x = sample_size, y = value, color = class)) + 
  geom_boxplot(size = 0.60) +
  scale_color_manual(values = group_colors) + 
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent Bias") + 
  facet_grid(cols = vars(race), scales = "free_y") + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "HCAP Sample Size", 
                   breaks = unique(plot_data$sample_size)) + 
  theme(text = element_text(size = 18))    

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "impairement_class_percent_bias_count_by_race.jpeg"))

#---- Figure X: HRS model results ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

results <- do.call(rbind, lapply(results_paths, read_results)) %>% 
  #restrict results for now
  filter(prior_props == "ADAMS")

#---- **truth ----
truth <- 
  read_csv(paste0(path_to_box, "data/truth.csv")) %>% 
  filter(term %in% c("black", "hispanic")) %>% 
  dplyr::select("term", "estimate", "dataset_name") %>% 
  separate(dataset_name, 
           into = c("Distribution", "sample_size", "prior_props"), 
           sep = "_") %>% 
  mutate_at(c("term", "Distribution"), str_to_sentence) %>% 
  rename_with(c("term", "estimate"), .fn = ~ c("race_eth", "beta")) %>% 
  mutate_at("race_eth", function(x) paste0(x, " vs. White")) %>%
  #restrict results for now
  filter(prior_props == "ADAMS")

truth$Distribution <- 
  factor(truth$Distribution, levels = c("Normal", "Lognormal", "Bathtub"))

#---- **plot data ----
results_summary <- results %>%
  group_by(calibration, Distribution, sample_size) %>%
  summarize_at(.vars = c("black_beta", "hispanic_beta", "black_se", 
                         "hispanic_se", "black_LCI", "hispanic_LCI", 
                         "black_UCI", "hispanic_UCI", "black_coverage", 
                         "hispanic_coverage"), 
               ~mean(., na.rm = TRUE)) %>% 
  pivot_longer(cols = !c("calibration", "Distribution", "sample_size"), 
               names_to = c("race_eth", ".value"),
               names_sep = "_") %>% 
  mutate_at("race_eth", str_to_sentence) %>% 
  mutate_at("race_eth", function(x) paste0(x, " vs. White")) %>%
  mutate_at("sample_size", as.numeric) %>% 
  mutate("sample_size" = 0.5*sample_size) 

results_summary$sample_size <- 
  factor(results_summary$sample_size, 
         levels = c("250", "500", "1000", "2000", "4000"))

results_summary$Distribution <- 
  factor(truth$Distribution, levels = c("Normal", "Lognormal", "Bathtub"))

#---- **color palette ----
navy <- "#135467"

#---- **plot ----
ggplot(results_summary, 
       #%>% filter(!sample_size %in% c("500", "1000")), 
       aes(x = beta, y = sample_size)) +
  geom_point(size = 4, position = position_dodge(-0.8), color = navy) + 
  # geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .3,
  #               position = position_dodge(-0.8), color = navy) +
  theme_bw() + xlab("Beta (log RR)") + ylab("\"HCAP\" Sample Size") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  geom_vline(xintercept = 0, color = "dark gray", linetype = "dashed", size = 1) + 
  geom_vline(data = truth %>% filter(Distribution == "Normal"), 
             aes(xintercept = beta), size = 1) +
  facet_grid(cols = vars(race_eth), rows = vars(calibration)) + 
  theme(text = element_text(size = 18))   


ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "HRS_model_results.jpeg"), 
       height = 6, width = 15, units = "in")

