#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "devtools", "here")
install_github("thomasp85/patchwork")

#---- source functions ----
source(here::here("functions", "read_results.R"))

#---- Figure X: comparing ADAMS with synthetic HRS ----
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

#---- Figure X: 95% CI impairment classes ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

results <- do.call(rbind, lapply(results_paths, read_results))

#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% 
  group_by(distribution, sample_size, prior_props) %>% 
  summarise_at(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"), 
               mean) %>% 
  mutate_at("sample_size", as.numeric) %>% 
  mutate("sample_size" = 0.5*sample_size) %>%
  pivot_longer(paste0(c("Unimpaired", "MCI", "Dementia", "Other"), "_coverage"),
               names_to = c("class", "coverage"), names_sep = "_") 

#---- **plot ----
ggplot(data = plot_data, aes(x = sample_size, y = value, group = class)) + 
  geom_line(aes(color = class)) + geom_point(aes(color = class)) + 
  scale_color_manual(values = group_colors) +
  geom_hline(yintercept = 0.95, lty = "dashed") +
  theme_bw() + ylab("95% CI Coverage") + 
  facet_grid(cols = vars(distribution)) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_continuous(name = "Sample Size", 
                     breaks = unique(plot_data$sample_size))

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "impairement_class_coverage.jpeg"))

#----- Figure X: bias impairment class counts ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*.csv")

results <- do.call(rbind, lapply(results_paths, read_results))

#---- ****% bias ----
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  results %<>% 
    mutate(!!paste0("percent_increase_", class) := 
             !!sym(paste0("bias_", class))/!!sym(paste0("true_", class))*100)
}

#---- **color palette ----
color_palette <- read_csv(here("color_palette.csv"))

group_colors <- color_palette$Color
names(group_colors) <- color_palette$Group

#---- **plot data ----
plot_data <- results %>% 
  dplyr::select("sample_size", "distribution", 
                paste0("percent_increase_", 
                       c("Unimpaired", "MCI", "Dementia", "Other"))) %>%
  mutate_at("sample_size", as.numeric) %>% 
  mutate("sample_size" = 0.5*sample_size) %>% 
  mutate_at("sample_size", as.factor) %>%
  pivot_longer(paste0("percent_increase_", 
                      c("Unimpaired", "MCI", "Dementia", "Other")),
               names_to = c("text", "class"), 
               names_sep = "_increase_") %>% 
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other")))

#---- **plot ----
ggplot(data = plot_data, aes(x = sample_size, y = value, color = class)) + 
  geom_boxplot() +
  scale_color_manual(values = group_colors) + 
  geom_hline(yintercept = 0, lty = "dashed") + theme_bw() + 
  ylab("Percent Bias") + facet_grid(cols = vars(distribution)) + 
  guides(color = guide_legend(title = "Group")) + 
  scale_x_discrete(name = "Sample Size", 
                   breaks = unique(plot_data$sample_size))

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "impairement_class_percent_bias.jpeg"))

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
  read_csv(paste0(path_to_box, "analyses/simulation_study/truth.csv")) %>% 
  filter(term %in% c("black", "hispanic")) %>% 
  dplyr::select("term", "estimate", "dataset_name") %>% 
  separate(dataset_name, 
                    into = c("Distribution", "sample_size", "prior_props"), 
                    sep = "_") %>% 
  mutate_at(c("term", "Distribution"), str_to_sentence) %>% 
  rename_with(c("term", "estimate"), .fn = ~ c("race_eth", "beta")) %>% 
  #restrict results for now
  filter(prior_props == "ADAMS")

truth$Distribution <- 
  factor(truth$Distribution, levels = c("Normal", "Lognormal", "Bathtub"))

#---- **plot data ----
results_summary <- results %>%
  group_by(Distribution, sample_size) %>%
  summarize_at(.vars = c("black_beta", "hispanic_beta", "black_se", 
                         "hispanic_se", "black_LCI", "hispanic_LCI", 
                         "black_UCI", "hispanic_UCI", "black_coverage", 
                         "hispanic_coverage"), 
               ~mean(., na.rm = TRUE)) %>% 
  pivot_longer(cols = !c("Distribution", "sample_size"), 
               names_to = c("race_eth", ".value"),
               names_sep = "_") %>% 
  mutate_at("race_eth", str_to_sentence) 

results_summary$Distribution <- 
  factor(results_summary$Distribution, levels = c("Normal", "Lognormal", "Bathtub"))
  
#---- color palette ----
navy <- "#135467"

#---- **plot ----
ggplot(results_summary, 
       aes(x = beta, y = sample_size)) +
  geom_point(size = 2.0, position = position_dodge(-0.8), color = navy) + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = .3,
                position = position_dodge(-0.8), color = navy) +
  theme_bw() + xlab("Beta") + ylab("Sample Size") +
  theme(legend.position = "bottom", legend.direction = "horizontal") + 
  geom_vline(xintercept = 0, color = "dark gray", linetype = "dashed") + 
  geom_vline(data = truth, aes(xintercept = beta)) +
  facet_grid(rows = vars(race_eth), cols = vars(Distribution))  

ggsave(filename = paste0(path_to_box, "figures/simulation_study/", 
                         "HRS_model_results.jpeg"))






