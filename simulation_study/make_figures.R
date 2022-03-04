#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "devtools")
install_github("thomasp85/patchwork")

#---- Figure X: comparing ADAMS with synthetic HRS ----
#---- **read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))

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
continuous_vars <- colnames(ADAMS_imputed_stacked)[
  str_detect(colnames(ADAMS_imputed_stacked), "_Z")] 

#---- **ADAMS plots ----
plot_data <- ADAMS_imputed_stacked %>% 
  dplyr::select(c("age_Z", "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  pivot_longer(cols = -c("age_Z"), names_to = "Group") %>% filter(value == 1)

ggplot(data = plot_data, aes(x = age_Z, color = Group, fill = Group)) + 
  geom_density(alpha = 0.5) + theme_minimal() + 
  xlab(unlist(variable_labels[which(variable_labels$data_label == "age_Z"), 
                       "figure_label"])) + 
  theme(text = element_text(size = 8)) + 
  scale_color_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)], 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)], 
                    guide = guide_legend(reverse = TRUE)) + 
  theme(legend.position = "bottom") 

#---- **synthetic HRS plots ----
plot_data <- synthetic_normal_1000 %>% 
  dplyr::select(c("age_Z", "Unimpaired", "MCI", "Dementia", "Other")) %>% 
  pivot_longer(cols = -c("age_Z"), names_to = "Group") %>% filter(value == 1)

ggplot(data = plot_data, aes(x = age_Z, color = Group, fill = Group)) + 
  geom_density(alpha = 0.5) + theme_minimal() + 
  xlab(unlist(variable_labels[which(variable_labels$data_label == "age_Z"), 
                              "figure_label"])) + 
  theme(text = element_text(size = 8)) + 
  scale_color_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)], 
                     guide = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)], 
                    guide = guide_legend(reverse = TRUE)) + 
  theme(legend.position = "bottom") 


