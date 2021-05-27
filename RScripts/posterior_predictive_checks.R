#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer")

#---- read in data ----
#---- **ADAMS ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"), 
                        col_types = cols(HHIDPN = col_character())) 

#---- **synthetic ----
for(run in 1:10){
  if(run == 1){
    synthetic_sample <- 
      read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                      "results/ADAMSA/standard_normal/ADAMSA_synthetic_", run, 
                      ".csv")) %>% mutate("sample" = run)
  } else{
    synthetic_sample %<>% 
      rbind(., 
            read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/ADAMSA_synthetic_", 
                            run, ".csv")) %>% mutate("sample" = run))
  }
}

num_samples <- max(synthetic_sample$sample)

#---- categorical checks ----
#---- **race/ethnicity x stroke ----
for(group in unique(ADAMS_train$Adem_dx_cat)){
  subset <- ADAMS_train %>% filter(Adem_dx_cat == group)
  assign(paste0(group, "_data_counts"), 
         as.data.frame(table(subset$ETHNIC_label, 
                             subset$Astroke)) %>% 
           mutate("percent" = round((Freq/sum(Freq))*100, 1)) %>% 
           unite("cell", c("Var1", "Var2"), sep = ":")) 
}

synthetic_counts <- matrix(0, nrow = 6*4, ncol = (num_samples + 2)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, num_samples), "group", "cell"))
synthetic_counts[, "group"] <- rep(unique(ADAMS_train$Adem_dx_cat), each = 6)

cells <- 
  as.data.frame(table(ADAMS_train$ETHNIC_label, ADAMS_train$Astroke)) %>% 
  unite("cell", c("Var1", "Var2"), sep = ":")

synthetic_counts[, "cell"] <- rep(cells$cell, 4)

#counts from synthetic datasets
for(num in 1:num_samples){
  subsample <- synthetic_sample %>% filter(sample == num)
  
  for(group in unique(ADAMS_train$Adem_dx_cat)){
    subset <- subsample %>% filter(Adem_dx_cat == group) 
    
    counts <- 
      as.data.frame(table(subset$ETHNIC_label, subset$Astroke)) %>% 
      unite("cell", c("Var1", "Var2"), sep = ":")
    
    synthetic_counts[
      which(synthetic_counts$group == group & 
              synthetic_counts$cell %in% counts$cell), num] <- counts$Freq
  }
}

synthetic_count_plot_data <- synthetic_counts %>% mutate("truth" = 0) 

for(group in unique(ADAMS_train$Adem_dx_cat)){
  true_counts <- get(paste0(group, "_data_counts"))
  synthetic_count_plot_data[
    which(synthetic_count_plot_data$group == group & 
            synthetic_count_plot_data$cell %in% true_counts$cell), "truth"] <- 
    true_counts$Freq
}

synthetic_count_plot_data %<>% 
  mutate("cat" = rep(c("Black + No Stroke", "Hispanic + No Stroke", 
                       "White + No Stroke", "Black + Stroke", 
                       "Hispanic + Stroke", "White + Stroke"), 4)) %>% 
  pivot_longer(-c("group", "cell", "truth", "cat")) %>% 
  mutate("color" = case_when(group == "Normal" ~ "#00a389", 
                             group == "Other" ~ "#28bed9", 
                             group == "MCI" ~ "#fdab00", 
                             group == "Dementia" ~ "#ff0000"))

#---- ****plot ----
for(dem_group in unique(synthetic_count_plot_data$group)){
  for(category in unique(synthetic_count_plot_data$cat)){
    subset <- synthetic_count_plot_data %>% 
      filter(group == dem_group & cat == category)
    ggplot(data = subset , aes(x = value)) + 
      geom_histogram(fill = "black", color = "black") + theme_minimal() + 
      xlab("Count") + ggtitle(category) + 
      geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                 size = 2)
    
    ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                             "posteriors/cell_counts/group_specific/", 
                             tolower(dem_group), "/", tolower(dem_group), "_", 
                             category, "_count.jpeg"), 
           width = 5, height = 3, units = "in")
  } 
}

#---- continuous checks ----
#---- **density plots ----
#---- **median plots ----
#---- **skew plots ----

#---- impairment classification ----

