#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "devtools", "here")
install_github("thomasp85/patchwork")

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
  read_csv(paste0(path_to_box, "analyses/simulation_study/", "synthetic_data/", 
                  "synthetic_normal_1000.csv"))

#---- **define plot variables ----
continuous_vars <- colnames(synthetic_normal_1000)[
  str_detect(colnames(synthetic_normal_1000), "_Z")] 

plot_labels <- variable_labels %>% filter(data_label %in% continuous_vars) 

#---- **ADAMS plots ----
plot_data <- ADAMS_imputed_stacked %>% 
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


