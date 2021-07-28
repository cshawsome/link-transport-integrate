#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/ADAMS/cleaned/ADAMS_train.csv"))

#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")

#---- cross-class race/ethnicity x stroke ----
overall_data_counts <- as.data.frame(table(ADAMS_subset$ETHNIC_label, 
                                           ADAMS_subset$Astroke)) %>% 
  mutate("prop" = Freq/sum(Freq))

for(group in unique(ADAMS_subset$Adem_dx_cat)){
  subset <- ADAMS_subset %>% filter(Adem_dx_cat == group)
  assign(paste0(group, "_data_counts"), 
         as.data.frame(table(subset$ETHNIC_label, 
                             subset$Astroke)) %>% 
           mutate("prop" = Freq/sum(Freq)) %>% 
           unite("cell", c("Var1", "Var2"), sep = ":")) 
}

#---- bootstrap props ----
B = 10000
bootstrap_props <- matrix(0, nrow = 6*4, ncol = (B + 2)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, B), "group", "cell"))
bootstrap_props[, "group"] <- rep(unique(ADAMS_subset$Adem_dx_cat), each = 6)

cells <- 
  as.data.frame(table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke)) %>% 
  unite("cell", c("Var1", "Var2"), sep = ":")

bootstrap_props[, "cell"] <- rep(cells$cell, 4)

for(b in 1:B){
  sample <- sample_n(ADAMS_subset, size = nrow(ADAMS_subset), replace = TRUE)
  for(group in unique(ADAMS_subset$Adem_dx_cat)){
    sub_sample <- sample %>% filter(Adem_dx_cat == group) 
    
    counts <- 
      as.data.frame(table(sub_sample$ETHNIC_label, sub_sample$Astroke)) %>% 
      unite("cell", c("Var1", "Var2"), sep = ":")
    
    bootstrap_props[
      which(bootstrap_props$group == group & 
              bootstrap_props$cell %in% counts$cell), b] <- 
      counts$Freq/nrow(sub_sample)
  }
}

bootstrap_props_plot_data <- bootstrap_props %>% mutate("truth" = 0) 

for(group in unique(ADAMS_subset$Adem_dx_cat)){
  true_counts <- get(paste0(group, "_data_counts"))
  bootstrap_props_plot_data[
    which(bootstrap_props_plot_data$group == group & 
            bootstrap_props_plot_data$cell %in% true_counts$cell), "truth"] <- 
    true_counts$prop
}

bootstrap_props_plot_data %<>% 
  mutate("cat" = rep(c("Black + No Stroke", "Hispanic + No Stroke", 
                       "White + No Stroke", "Black + Stroke", 
                       "Hispanic + Stroke", "White + Stroke"), 4)) %>% 
  pivot_longer(-c("group", "cell", "truth", "cat")) %>% 
  mutate("color" = case_when(group == "Unimpaired" ~ "#00a389", 
                             group == "Other" ~ "#28bed9", 
                             group == "MCI" ~ "#fdab00", 
                             group == "Dementia" ~ "#ff0000"))

#---- **plots ----
#---- ****create directories ----
for(dem_group in unique(ADAMS_subset$Adem_dx_cat)){
 dir.create(paste0("/Users/CrystalShaw/Box/Dissertation/figures/ADAMS_test/", 
                   "prior_predictive_checks/cell_props/group_specific/", 
                   tolower(dem_group)), recursive = TRUE) 
}

#---- ****create plot ----
for(dem_group in unique(bootstrap_props_plot_data$group)){
  for(category in unique(bootstrap_props_plot_data$cat)){
    subset <- bootstrap_props_plot_data %>% 
      filter(group == dem_group & cat == category)
    ggplot(data = subset , aes(x = value)) + 
      geom_histogram(fill = "black", color = "black") + theme_minimal() + 
      xlab("Proportion") + ggtitle(category) + 
      geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                 size = 2)
    
    ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                             "ADAMS_test/prior_predictive_checks/cell_props/", 
                             "group_specific/", tolower(dem_group), "/", 
                             tolower(dem_group), "_", category, "_prop.jpeg"), 
           width = 5, height = 3, units = "in")
  } 
}

bootstrap_props %>% as.data.frame() %>% 
  mutate("group_number" = case_when(group == "Unimpaired" ~ 1, 
                                    group == "Other" ~ 2, 
                                    group == "MCI" ~ 3, 
                                    group == "Dementia" ~ 4)) %>%
  write_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/ADAMS_test/", 
                   "bootstrap_cell_props.csv"))






