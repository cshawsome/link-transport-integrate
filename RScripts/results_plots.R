#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer", 
       "devtools")
install_github("thomasp85/patchwork")

#---- read in data ----
synthetic_ADAMS <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "analyses/results/ADAMSA/ADAMSA_synthetic.csv"), 
           col_types = cols(HHIDPN = col_character()))

ADAMS_columns <- c(colnames(synthetic_ADAMS)[
  which(!colnames(synthetic_ADAMS) %in% 
          c("(Intercept)", "Group", "Black", "Hispanic", "p_Unimpaired", 
            "p_Other", "p_MCI"))], "Adem_dx_cat")

ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character())) %>% 
  #only take those part of modeling (complete cases only for now)
  filter(HHIDPN %in% synthetic_ADAMS$HHIDPN) %>% 
  dplyr::select(all_of(ADAMS_columns))

#---- merge datasets ----
#format column names
synthetic_ADAMS %<>% 
  set_colnames(c("HHIDPN", paste0("synthetic:", colnames(synthetic_ADAMS))[-1]))

ADAMS_subset %<>% 
  set_colnames(c("HHIDPN", paste0("ADAMSA:", colnames(ADAMS_subset))[-1]))

merged_data <- left_join(synthetic_ADAMS, ADAMS_subset)

#---- data cleaning: dem group ----
# #data check
# table(merged_data$`ADAMSA:Adem_dx_cat`, useNA = "ifany")
# table(merged_data$`synthetic:Group`, useNA = "ifany")

merged_data %<>% 
  mutate("ADAMSA:group_class" = 
           case_when(`ADAMSA:Adem_dx_cat` %in% 
                       c("Dementia", "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia",
                     `ADAMSA:Adem_dx_cat` == "Normal" ~ "Unimpaired",
                     TRUE ~ `ADAMSA:Adem_dx_cat`), 
         "synthetic:group_class" = 
           case_when(`synthetic:Group` == 1 ~ "Unimpaired", 
                     `synthetic:Group` == 2 ~ "Other", 
                     `synthetic:Group` == 3 ~ "MCI", 
                     `synthetic:Group` == 4 ~ "Dementia"))

# #Sanity check
# table(merged_data$`ADAMSA:group_class`)
# table(merged_data$`synthetic:group_class`)

#---- plots: categorical vars ----
#---- **race/ethnicity ----
race_ethnicity_data <- 
  merged_data %>% dplyr::select(contains("ETHNIC_label")) %>% 
  pivot_longer(everything(), names_to = c("Data", "Var"), 
               names_pattern = "(.*):(.*)") %>% 
  count(Data, value) %>%
  group_by(Data) %>%
  mutate(n = n/sum(n))

race_ethnicity_plot <- 
  ggplot(data = race_ethnicity_data) + 
  geom_bar(mapping = aes(x = factor(value), y = n, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Data") + theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **stroke ----
stroke_data <- 
  merged_data %>% dplyr::select(contains("Astroke")) %>% 
  pivot_longer(everything(), names_to = c("Data", "Var"), 
               names_pattern = "(.*):(.*)") %>% 
  count(Data, value) %>%
  group_by(Data) %>%
  mutate(n = n/sum(n))

stroke_plot <- 
  ggplot(data = stroke_data) + 
  geom_bar(mapping = aes(x = factor(value), y = n, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Stroke") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Data") + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **stroke x race/ethnicity ----
stroke_race_ethnicity_data <- 
  merged_data %>% dplyr::select(contains(c("Astroke", "ETHNIC_label"))) %>% 
  mutate_if(is.numeric, as.character) %>%
  pivot_longer(everything(), names_to = c("Data", ".value"), 
               names_pattern = "(.*):(.*)") %>% 
  count(Data, Astroke, ETHNIC_label) %>%
  group_by(Data, ETHNIC_label) %>%
  mutate(n = n/sum(n))

stroke_race_ethnicity_plot <- 
  ggplot(data = stroke_race_ethnicity_data) + 
  geom_bar(mapping = aes(x = factor(Astroke), y = n, 
                         fill = factor(ETHNIC_label)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Stroke") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Race/Ethnicity") + 
  facet_grid(cols = vars(Data)) + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **patchwork plot ----
race_ethnicity_plot + stroke_plot + stroke_race_ethnicity_plot

ggsave(filename = "race_ethnicity_stroke.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/"), width = 14, height = 5, units = "in", 
       device = "jpeg")

#---- plots: continuous vars ----
continuous_vars <- c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", 
                     "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

for(var in continuous_vars){
  plot_data <- merged_data %>% dplyr::select(contains(var)) %>% 
    pivot_longer(everything(), names_to = c("Data", "Var"), 
                 names_pattern = "(.*):(.*)")
  
  plot <- ggplot(data = plot_data, aes(x = value, color = Data)) + 
    geom_density() + theme_minimal() + xlab(var) + 
    scale_color_manual(values = rev(wes_palette("Darjeeling1")))
  
  if(var != continuous_vars[length(continuous_vars)]){
    plot <- plot + theme(legend.position = "none") 
  }
  
  assign(paste0(var, "_plot"), plot)
}

#---- **patchwork plot ----
continuous_var_plot_names <- paste0(continuous_vars, "_plot")

((((get(continuous_var_plot_names[1]) + get(continuous_var_plot_names[2]) + 
  get(continuous_var_plot_names[3])) /
  (get(continuous_var_plot_names[4]) + get(continuous_var_plot_names[5]) + 
     get(continuous_var_plot_names[6])))) / 
  (get(continuous_var_plot_names[7]) + get(continuous_var_plot_names[8]) + 
        get(continuous_var_plot_names[9]))) / 
  get(continuous_var_plot_names[10])

ggsave(filename = "continuous_vars.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/"), width = 12, height = 12, units = "in", 
       device = "jpeg")

#---- plots: dementia classes ----
dementia_class_plot_data <- 
  merged_data %>% dplyr::select(contains("group_class")) %>% 
  pivot_longer(everything(), names_to = c("Data", "Var"), 
               names_pattern = "(.*):(.*)") %>% 
  count(Data, value) %>%
  group_by(Data) %>%
  mutate(prop = n/sum(n))

#releveling factors
dementia_class_plot_data$value <- 
  fct_relevel(dementia_class_plot_data$value, 
              c("Unimpaired", "MCI", "Dementia", "Other"))

dementia_class_plot <- 
  ggplot(data = dementia_class_plot_data) + 
  geom_bar(mapping = aes(x = factor(value), y = prop, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Impairment Class") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Data") + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "group_class.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/results/ADAMSA", 
       width = 5, height = 5, units = "in", device = "jpeg")

#---- plots: dementia classes x race ----
dementia_class_by_race_plot_data <- 
  merged_data %>% dplyr::select(contains(c("group_class", "ETHNIC_label"))) %>% 
  pivot_longer(everything(), names_to = c("Data", ".value"), 
               names_pattern = "(.*):(.*)") %>% 
  count(Data, group_class, ETHNIC_label) %>%
  group_by(Data, ETHNIC_label) %>%
  mutate(prop = n/sum(n))

#releveling factors
dementia_class_by_race_plot_data$group_class <- 
  fct_relevel(dementia_class_by_race_plot_data$group_class, 
              c("Unimpaired", "MCI", "Dementia", "Other"))

dementia_class_by_race_plot_n <- 
  ggplot(data = dementia_class_by_race_plot_data) + 
  geom_bar(mapping = aes(x = factor(group_class), y = n, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Impairment Class") + ylab("Count") + 
  labs(fill = "Data") + facet_grid(cols = vars(ETHNIC_label)) + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "group_class_by_race_n.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/results/ADAMSA", 
       width = 10, height = 7, units = "in", device = "jpeg")

dementia_class_by_race_plot_p <- 
  ggplot(data = dementia_class_by_race_plot_data) + 
  geom_bar(mapping = aes(x = factor(group_class), y = prop, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Impairment Class") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Data") + facet_grid(cols = vars(ETHNIC_label)) + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "group_class_by_race_p.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/results/ADAMSA", 
       width = 10, height = 7, units = "in", device = "jpeg")

#---- sensitivity/specificity ----
#---- **overall ----
#---- **race-stratified ----



