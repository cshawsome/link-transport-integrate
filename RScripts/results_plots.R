#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer", 
       "devtools", "scales", "plyr", "rje")
install_github("thomasp85/patchwork")

#---- read in data ----
synthetic_ADAMS <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "analyses/results/ADAMSA/double_transform/", 
                  "ADAMSA_synthetic.csv"), 
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

#---- double_transform ----
cont_vars <- c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
               "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

ADAMS_transform <- ADAMS_subset %>% 
  dplyr::select(c("HHIDPN", all_of(cont_vars))) 

ADAMS_transform[, cont_vars] <- 
  apply(ADAMS_transform[, cont_vars], 2, function(x) (x + 1)/(max(x) + 2))

# #Sanity check
# apply(ADAMS_transform[, cont_vars], 2, max)
# apply(ADAMS_transform[, cont_vars], 2, min)

#Transform for Real
ADAMS_transform[, cont_vars] <- 
  apply(ADAMS_transform[, cont_vars], 2, function(x) logit(x))

# #Sanity check
# apply(ADAMS_transform[, cont_vars], 2, max)
# apply(ADAMS_transform[, cont_vars], 2, min)

synthetic_transformed <- synthetic_ADAMS %>% 
  dplyr::select(c("HHIDPN", all_of(cont_vars)))

#Transform to [0, 1]
synthetic_transformed[, cont_vars] <- 
  apply(synthetic_transformed[, cont_vars], 2, function(x) expit(x))

max <- ADAMS_subset %>% dplyr::select(all_of(cont_vars)) %>% apply(., 2, max)

for(var in names(max)){
  synthetic_transformed[, var] <- 
    synthetic_transformed[, var]*(max[var] + 2) - 1
}

#---- merge datasets ----
#format column names
synthetic_ADAMS %<>% 
  set_colnames(c("HHIDPN", 
                 paste0("syntheticTransformed:", 
                        colnames(synthetic_ADAMS))[-1]))

synthetic_transformed %<>% 
  set_colnames(c("HHIDPN", 
                 paste0("syntheticX:", colnames(synthetic_transformed))[-1]))

ADAMS_subset %<>% 
  set_colnames(c("HHIDPN", paste0("ADAMSAX:", colnames(ADAMS_subset))[-1]))

ADAMS_transform %<>% 
  set_colnames(c("HHIDPN", paste0("ADAMSATransformed:", 
                                  colnames(ADAMS_transform))[-1]))

merged_data <- 
  join_all(list(ADAMS_subset, ADAMS_transform, synthetic_ADAMS, 
                synthetic_transformed), by = "HHIDPN", type = "left")

#---- data cleaning: dem group ----
# #data check
# table(merged_data$`ADAMSA:Adem_dx_cat`, useNA = "ifany")
# table(merged_data$`synthetic:Group`, useNA = "ifany")

merged_data %<>% 
  mutate("ADAMSA:group_class" = 
           case_when(`ADAMSAX:Adem_dx_cat` %in% 
                       c("Dementia", "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia",
                     `ADAMSAX:Adem_dx_cat` == "Normal" ~ "Unimpaired",
                     TRUE ~ `ADAMSAX:Adem_dx_cat`), 
         "synthetic:group_class" = 
           case_when(`syntheticTransformed:Group` == 1 ~ "Unimpaired", 
                     `syntheticTransformed:Group` == 2 ~ "Other", 
                     `syntheticTransformed:Group` == 3 ~ "MCI", 
                     `syntheticTransformed:Group` == 4 ~ "Dementia"))

# #Sanity check
# table(merged_data$`ADAMSA:group_class`)
# table(merged_data$`synthetic:group_class`)

#---- plots: categorical vars ----
#---- **race/ethnicity ----
race_ethnicity_data <- 
  merged_data %>% dplyr::select(contains("ETHNIC_label")) %>% 
  pivot_longer(everything(), names_to = c("Data", "Var"), 
               names_pattern = "(.*):(.*)") %>% 
  dplyr::count(Data, value) %>%
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
  dplyr::count(Data, value) %>%
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
  dplyr::count(Data, Astroke, ETHNIC_label) %>%
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
                     "ADAMSA/double_transform"), width = 14, height = 5, 
       units = "in", device = "jpeg")

#---- plots: continuous vars ----
for(var in cont_vars){
  plot_data <- merged_data %>% dplyr::select(contains(var)) %>% 
    pivot_longer(everything(), names_to = c("Data", "Var"), 
                 names_pattern = "(.*):(.*)")
  
  plotX <- ggplot(data = plot_data %>% 
                    filter(Data %in% c("ADAMSAX", "syntheticX")), 
                  aes(x = value, color = Data)) + 
    geom_density() + theme_minimal() + xlab(var) + 
    scale_color_manual(values = rev(wes_palette("Darjeeling1")))
  
  plotZ <- ggplot(data = plot_data %>% 
                    filter(Data %in% 
                             c("ADAMSATransformed", "syntheticTransformed")), 
                  aes(x = value, color = Data)) + 
    geom_density() + theme_minimal() + xlab(var) + 
    scale_color_manual(values = rev(wes_palette("Darjeeling1")))
  
  if(var != cont_vars[length(cont_vars)]){
    plotX <- plotX + theme(legend.position = "none") 
    plotZ <- plotZ + theme(legend.position = "none") 
  }
  
  assign(paste0(var, "_plotX"), plotX)
  assign(paste0(var, "_plotZ"), plotZ)
}

#---- **patchwork plot ----
continuous_var_plot_names <- paste0(cont_vars, "_plotX")

((((get(continuous_var_plot_names[1]) + get(continuous_var_plot_names[2]) + 
      get(continuous_var_plot_names[3])) /
     (get(continuous_var_plot_names[4]) + get(continuous_var_plot_names[5]) + 
        get(continuous_var_plot_names[6])))) / 
    (get(continuous_var_plot_names[7]) + get(continuous_var_plot_names[8]) + 
       get(continuous_var_plot_names[9]))) / 
  get(continuous_var_plot_names[10])

ggsave(filename = "continuous_varsX.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), width = 12, height = 12, 
       units = "in", device = "jpeg")

continuous_var_plot_names <- paste0(cont_vars, "_plotZ")

((((get(continuous_var_plot_names[1]) + get(continuous_var_plot_names[2]) + 
      get(continuous_var_plot_names[3])) /
     (get(continuous_var_plot_names[4]) + get(continuous_var_plot_names[5]) + 
        get(continuous_var_plot_names[6])))) / 
    (get(continuous_var_plot_names[7]) + get(continuous_var_plot_names[8]) + 
       get(continuous_var_plot_names[9]))) / 
  get(continuous_var_plot_names[10])

ggsave(filename = "continuous_varsTransform.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), width = 12, height = 12, 
       units = "in", device = "jpeg")

#---- plots: dementia classes ----
#---- **overall ----
dementia_class_plot_data <- 
  merged_data %>% dplyr::select(contains("group_class")) %>% 
  pivot_longer(everything(), names_to = c("Data", "Var"), 
               names_pattern = "(.*):(.*)") %>% 
  dplyr::count(Data, value) %>%
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
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), 
       width = 5, height = 5, units = "in", device = "jpeg")

#---- **race-stratified ----
dementia_class_by_race_plot_data <- 
  merged_data %>% dplyr::select(contains(c("group_class", "ETHNIC_label"))) %>%
  set_colnames(c("ADAMSA:group_class", "synthetic:group_class", 
                 "ADAMSA:ETHNIC_label", "synthetic:ETHNIC_label")) %>%
  pivot_longer(everything(), names_to = c("Data", ".value"), 
               names_pattern = "(.*):(.*)") %>% 
  dplyr::count(Data, group_class, ETHNIC_label) %>%
  dplyr::group_by(Data, ETHNIC_label) %>%
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
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), 
       width = 10, height = 7, units = "in", device = "jpeg")

dementia_class_by_race_plot_p <- 
  ggplot(data = dementia_class_by_race_plot_data) + 
  geom_bar(mapping = aes(x = factor(group_class), y = prop, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Impairment Class") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Data") + facet_grid(cols = vars(ETHNIC_label)) + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "group_class_by_race_p.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), 
       width = 10, height = 7, units = "in", device = "jpeg")

#---- plot: %change by category ----
#---- **overall ----
synthetic_counts <- dementia_class_plot_data %>% filter(Data == "synthetic")  
ADAMS_counts <- dementia_class_plot_data %>% filter(Data == "ADAMSA")

overall_change_data <- 
  tibble("change" = ((synthetic_counts$n - ADAMS_counts$n)/ADAMS_counts$n)*100, 
         "group" = synthetic_counts$value) %>% arrange(desc(change))

#releveling factors by %change
overall_change_data$group <- 
  fct_relevel(overall_change_data$group, 
              paste0(unique(overall_change_data$group)))

overall_change <- 
  ggplot(data = overall_change_data, aes(x = group, y = change)) + 
  geom_bar(stat = "identity", width = 0.7, 
           color = rev(wes_palette("Darjeeling1"))[1],
           fill = rev(wes_palette("Darjeeling1"))[1],
           position = position_dodge(width = 0.4)) + ylab("% change") + 
  theme_minimal()

ggsave(filename = "group_percent_change_overall.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), 
       width = 7, height = 5, units = "in", device = "jpeg")

#---- **race-stratified ----
synthetic_counts_by_race <- dementia_class_by_race_plot_data %>% 
  filter(Data == "synthetic")  
ADAMS_counts_by_race <- dementia_class_by_race_plot_data %>% 
  filter(Data == "ADAMSA")

overall_change_by_race_data <- 
  tibble("change" = ((synthetic_counts_by_race$n - ADAMS_counts_by_race$n)/
           ADAMS_counts_by_race$n)*100, 
         "group" = synthetic_counts_by_race$group_class, 
         "race_ethnicity" = synthetic_counts_by_race$ETHNIC_label) %>% 
  unite("group_by_race", c(group, race_ethnicity), sep = "\n", 
        remove = FALSE) %>%
  arrange(desc(change))

#releveling factors by %change
overall_change_by_race_data$group_by_race <- 
  fct_relevel(overall_change_by_race_data$group_by_race, 
              paste0(unique(overall_change_by_race_data$group_by_race)))

overall_change_by_race <- 
  ggplot(data = overall_change_by_race_data, 
         aes(x = group_by_race, y = change, 
             color = race_ethnicity, fill = race_ethnicity)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) + ylab("% change") + 
  xlab("Group:Race") + 
  scale_color_manual(values = rev(wes_palette("Darjeeling1"))) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1"))) + 
  theme_minimal() + theme(text = element_text(size = 7)) 

ggsave(filename = "group_percent_change_by_race.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/results/", 
                     "ADAMSA/double_transform"), 
       width = 7, height = 5, units = "in", device = "jpeg")


