#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **prior imputed clean ----
prior_imputed_clean <- 
  readRDS(paste0(path_to_box, "analyses/simulation_study/prior_data/MI/", 
                 "MI_datasets_cleaned")) %>%
  lapply(function(x) mutate_at(x, "HHIDPN", as.numeric))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(ADAMS %in% colnames(prior_imputed_clean[[1]]))

#---- relabel columns ----
prior_imputed_clean <- 
  lapply(prior_imputed_clean, 
         function(x) rename_at(x, vars(variable_labels$ADAMS), ~ 
                                 variable_labels$data_label)) 

#---- categorical vars ----
#notation from Schafer 1997
W <- c("black", "hispanic", "stroke")

#---- imputation props ----
get_props <- function(data, W){
  counts <- data %>% unite("cell_ID", all_of(W), sep = "", remove = FALSE) %>% 
    dplyr::select("cell_ID", "Unimpaired", "MCI", "Dementia", "Other") %>% 
    pivot_longer(-c("cell_ID")) %>% filter(value == 1) %>% group_by(name) %>% 
    table() %>% as.data.frame() %>% dplyr::select(-"value") 
  
  counts %<>% 
    left_join(., counts %>% group_by(name) %>% summarize_at("Freq", sum), 
              by = "name") %>% mutate("props" = Freq.x/Freq.y) %>% 
    dplyr::select(-contains("Freq"))
}

imputation_props <- lapply(prior_imputed_clean, get_props, W) %>% 
  reduce(left_join, by = c("cell_ID", "name")) %>% 
  set_colnames(c("cell_ID", "group", seq(1, length(prior_imputed_clean)))) %>% 
  write_csv(paste0(path_to_box, "analyses/simulation_study/prior_data/", 
                   "imputation_cell_props.csv"))

#---- plots ----
#---- **read in data ----
#---- ****cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv"))

#---- ****color palette ----
color_palette <- read_csv(here("color_palette.csv"))

#---- **join data ----
imputation_props_plot_data <- left_join(imputation_props, cell_ID_key) %>% 
  left_join(., color_palette, by = c("group" = "Group")) %>% 
  pivot_longer(-c("group", "cell_ID", "cell_name", "Color", "cell_order")) %>% 
  mutate_at("group", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("cell_name", function(x) 
    factor(x, levels = c("White | No Stroke", "White | Stroke", 
                         "Black | No Stroke", "Black | Stroke", 
                         "Hispanic | No Stroke", "Hispanic | Stroke")))

#---- **make plot ----
colors <- color_palette$Color
names(colors) <- color_palette$Group

ggplot(data = imputation_props_plot_data , aes(x = value)) + 
  geom_histogram(aes(fill = group, color = group)) + theme_bw() + 
  xlab("Proportion") + 
  facet_grid(rows = vars(cell_name), cols = vars(group), scales = "free") + 
  scale_color_manual(values = colors) + scale_fill_manual(values = colors) + 
  theme(legend.position = "none")

ggsave(filename = "prior_cell_props.jpeg", plot = last_plot(), 
       path = paste0(path_to_box, "figures/simulation_study/"), 
       width = 8, height = 8, units = "in", device = "jpeg")
