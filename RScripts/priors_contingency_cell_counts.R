#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS/ADAMS_train.csv"))

#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- cbind(c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
             "proxy_cog", "ANDELCOR", "Aiadla", "Abmi"), 
           c("Age", "Total MMSE", "Serial 7s", "Immediate Word Recall", 
             "Wordlist Recall (Yes)", "Story Recall I", "Proxy Cognition (Avg)", 
             "Delayed Word Recall", "IADLs", "BMI")) %>% 
  set_colnames(c("var", "label"))

#---- cross-class race/ethnicity x stroke ----
overall_data_counts <- as.data.frame(table(ADAMS_subset$ETHNIC_label, 
                                           ADAMS_subset$Astroke)) %>% 
  mutate("percent" = round((Freq/sum(Freq))*100, 1))

for(group in unique(ADAMS_subset$Adem_dx_cat)){
  subset <- ADAMS_subset %>% filter(Adem_dx_cat == group)
  assign(paste0(group, "_data_counts"), 
         as.data.frame(table(subset$ETHNIC_label, 
                             subset$Astroke)) %>% 
           mutate("percent" = round((Freq/sum(Freq))*100, 1)) %>% 
           unite("cell", c("Var1", "Var2"), sep = ":")) 
}

#---- bootstrap counts ----
B = 10000
bootstrap_counts <- matrix(0, nrow = 6*4, ncol = (B + 2)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, B), "group", "cell"))
bootstrap_counts[, "group"] <- rep(unique(ADAMS_subset$Adem_dx_cat), each = 6)

cells <- 
  as.data.frame(table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke)) %>% 
  unite("cell", c("Var1", "Var2"), sep = ":")

bootstrap_counts[, "cell"] <- rep(cells$cell, 4)

for(b in 1:B){
  sample <- sample_n(ADAMS_subset, size = nrow(ADAMS_subset), replace = TRUE)
  for(group in unique(ADAMS_subset$Adem_dx_cat)){
    sub_sample <- sample %>% filter(Adem_dx_cat == group) 
    
    counts <- 
      as.data.frame(table(sub_sample$ETHNIC_label, sub_sample$Astroke)) %>% 
      unite("cell", c("Var1", "Var2"), sep = ":")
    
    bootstrap_counts[
      which(bootstrap_counts$group == group & 
              bootstrap_counts$cell %in% counts$cell), b] <- counts$Freq
  }
}

# bootstrap_percents <- 
#   round(bootstrap_counts/colSums(bootstrap_counts)*100, 1) %>% 
#   as.data.frame() %>%
#   mutate("truth" = data_counts$percent, 
#          "cat" = c("Black + No Stroke", "Hispanic + No Stroke", 
#                    "White + No Stroke", "Black + Stroke", "Hispanic + Stroke", 
#                    "White + Stroke")) %>% pivot_longer(-c("truth", "cat"))

bootstrap_count_plot_data <- bootstrap_counts %>% mutate("truth" = 0) 

for(group in unique(ADAMS_subset$Adem_dx_cat)){
  true_counts <- get(paste0(group, "_data_counts"))
  bootstrap_count_plot_data[
    which(bootstrap_count_plot_data$group == group & 
            bootstrap_count_plot_data$cell %in% counts$cell), "truth"] <- 
    true_counts$Freq
}

bootstrap_count_plot_data %<>% 
  mutate("cat" = rep(c("Black + No Stroke", "Hispanic + No Stroke", 
                       "White + No Stroke", "Black + Stroke", 
                       "Hispanic + Stroke", "White + Stroke"), 4)) %>% 
  pivot_longer(-c("group", "cell", "truth", "cat")) %>% 
  mutate("color" = case_when(group == "Normal" ~ "#00a389", 
                             group == "Other" ~ "#28bed9", 
                             group == "MCI" ~ "#fdab00", 
                             group == "Dementia" ~ "#ff0000"))

#---- **plots ----
# #percent
# for(category in unique(bootstrap_percents$cat)){
#   subset <- bootstrap_percents %>% filter(cat == category)
#   ggplot(data = subset , aes(x = value)) + 
#     geom_histogram(fill = "black", color = "black") + theme_minimal() + 
#     xlab("Percent") + ggtitle(category) + 
#     geom_vline(xintercept = subset$truth, color = "#f2caaa", size = 2)
#   
#   ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                            "priors/cell_counts/overall/", category, 
#                            "_percent.jpeg"), 
#          width = 5, height = 3, units = "in")
# }

#count
for(dem_group in unique(bootstrap_count_plot_data$group)){
  for(category in unique(bootstrap_count_plot_data$cat)){
    subset <- bootstrap_count_plot_data %>% 
      filter(group == dem_group & cat == category)
    ggplot(data = subset , aes(x = value)) + 
      geom_histogram(fill = "black", color = "black") + theme_minimal() + 
      xlab("Count") + ggtitle(category) + 
      geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                 size = 2)
    
    ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                             "priors/cell_counts/group_specific/", 
                             tolower(dem_group), "/", tolower(dem_group), "_", 
                             category, "_count.jpeg"), 
           width = 5, height = 3, units = "in")
  } 
}

bootstrap_counts %>% as.data.frame() %>% 
  mutate("group_number" = case_when(group == "Normal" ~ 1, 
                                    group == "Other" ~ 2, 
                                    group == "MCI" ~ 3, 
                                    group == "Dementia" ~ 4)) %>%
  write_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                   "bootstrap_cell_counts.csv"))






