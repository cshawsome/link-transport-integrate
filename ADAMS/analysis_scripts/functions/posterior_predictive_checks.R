
#---- read in data ----
#---- **ADAMS train ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"), 
                        col_types = cols(HHIDPN = col_character())) 

#Continuous vars (notation from Schafer 1997)
Z <- cbind(c("AAGE", "ANMSETOT_norm", "ANSER7T", "ANIMMCR", "ANRECYES", 
             "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Abmi"), 
           c("Age", "Normed Total MMSE", "Serial 7s", "Immediate Word Recall", 
             "Wordlist Recall (Yes)", "Story Recall I", "Proxy Cognition (Avg)", 
             "Delayed Word Recall", "IADLs", "BMI")) %>% 
  set_colnames(c("var", "label"))

#---- **ADAMS unscaled ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character())) %>% 
  filter(HHIDPN %in% ADAMS_train$HHIDPN)

ADAMS_means <- colMeans(ADAMS_subset %>% dplyr::select(all_of(Z[, "var"])))
ADAMS_sds <- apply(ADAMS_subset %>% dplyr::select(all_of(Z[, "var"])), 2, sd)

#---- **synthetic ----
num_samples = 1000
num_chains = 5
for(sample in 1:num_samples){
  for(chain in 1:num_chains){
    if(sample == 1 & chain == 1){
      synthetic_sample <- 
        read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                        "results/ADAMSA/standard_normal/run_", chain, 
                        "/ADAMSA_synthetic_", 
                        sample, ".csv")) %>% 
        mutate("sample" = sample, "chain" = chain)
    } else{
      synthetic_sample %<>% 
        rbind(., 
              read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                              "results/ADAMSA/standard_normal/run_", chain, 
                              "/ADAMSA_synthetic_", sample, ".csv")) %>% 
                mutate("sample" = sample, "chain" = chain))
    }
  }
  
}

#---- chain convergence ----
for(chain in 1:num_chains){
  if(chain == 1){
    group_membership <- 
      read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                      "results/ADAMSA/standard_normal/diagnostics_data/", 
                      "run_", chain,  "/ADAMSA_latent_class_data.csv")) %>% 
      mutate("chain" = chain)
  } else{
    group_membership %<>% 
      rbind(., 
            read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/diagnostics_data/", 
                            "run_", chain,  
                            "/ADAMSA_latent_class_data.csv")) %>% 
              mutate("chain" = chain))
  }
}

#---- **plot ----
chain_convergence <- ggplot(data = group_membership, 
                            aes(x = run, y = prob, colour = as.factor(chain))) +       
  geom_line(aes(group = as.factor(chain))) + xlab("Run") + ylab("Proportion") +
  facet_wrap(facets = vars(factor(Group, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other"))), nrow = 2, ncol = 2) + 
  theme_bw() 

ggsave(filename = "group_membership_chains.jpeg", plot = chain_convergence, 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                     "diagnostics/standard_normal"), 
       width = 10, height = 5, units = "in", device = "jpeg")

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
                             "posteriors/run_1/cell_counts/group_specific/", 
                             tolower(dem_group), "/", tolower(dem_group), "_", 
                             category, "_count.jpeg"), 
           width = 5, height = 3, units = "in")
  } 
}

#---- continuous checks ----
#---- **density ----
#---- **median ----
#---- ****overall ----
synthetic_continuous <- 
  matrix(0, nrow = nrow(Z), ncol = (num_samples + 1)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, num_samples), "var"))
synthetic_continuous[, "var"] <- Z[, "var"]

#medians from synthetic datasets
for(var in Z[, "var"]){
  subset <- synthetic_sample %>% dplyr::select(!!as.symbol(var), "sample") 
  subset[, var] <- subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
  synthetic_continuous[which(synthetic_continuous$var == var), 
                       seq(1, num_samples)] <- 
    t(subset %>% group_by(sample) %>% summarize_all(median) %>% 
        dplyr::select(var))
}

synthetic_continuous %<>% 
  mutate("truth" = 0, 
         "label" = Z[, "label"])

for(var in Z[, "var"]){
  synthetic_continuous[which(synthetic_continuous$var == var), "truth"] <- 
    median(unlist(ADAMS_subset[, var]))
}

synthetic_continuous %<>% pivot_longer(seq(1:num_samples))

#---- ****plot ----
for(var_name in Z[, "label"]){
  subset <- synthetic_continuous %>% filter(label == var_name)
  ggplot(data = subset , aes(x = value)) + 
    geom_histogram(fill = "black", color = "black") + theme_minimal() + 
    xlab("Median") + ggtitle(var_name) + 
    geom_vline(xintercept = subset$truth, color = "#f2caaa" , size = 2)
  
  ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                           "posteriors/run_1/continuous_vars/median/overall/", 
                           var_name, ".jpeg"), 
         width = 5, height = 3, units = "in")
} 

#---- ****by class ----
synthetic_continuous <- 
  matrix(0, nrow = nrow(Z)*4, ncol = (num_samples + 2)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, num_samples), "group", "var"))
synthetic_continuous[, "group"] <- rep(seq(1, 4), each = nrow(Z))
synthetic_continuous[, "var"] <- rep(Z[, "var"], 4)

#medians from synthetic datasets
for(group in 1:4){
  subsample <- synthetic_sample %>% filter(Group == group)
  
  for(var in Z[, "var"]){
    subset <- subsample %>% dplyr::select(!!as.symbol(var), "sample") 
    subset[, var] <- subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
    synthetic_continuous[which(synthetic_continuous$group == group & 
                                 synthetic_continuous$var == var), 
                         seq(1, num_samples)] <- 
      t(subset %>% group_by(sample) %>% summarize_all(median) %>% 
          dplyr::select(var))
  }
}

synthetic_continuous %<>% 
  mutate("truth" = 0, 
         "label" = rep(Z[, "label"], 4), 
         "group_label" = case_when(group == 1 ~ "Normal", 
                                   group == 2 ~ "Other", 
                                   group == 3 ~ "MCI", 
                                   group == 4 ~ "Dementia")) %>% 
  mutate("color" = case_when(group_label == "Normal" ~ "#00a389", 
                             group_label == "Other" ~ "#28bed9", 
                             group_label == "MCI" ~ "#fdab00", 
                             group_label == "Dementia" ~ "#ff0000"))

for(group in unique(ADAMS_subset$Adem_dx_cat)){
  subsample <- ADAMS_subset %>% filter(Adem_dx_cat == group)
  
  for(var in Z[, "var"]){
    synthetic_continuous[which(synthetic_continuous$group_label == group & 
                                 synthetic_continuous$var == var), 
                         "truth"] <- median(unlist(subsample[, var]))
  }
}

synthetic_continuous %<>% pivot_longer(seq(1:num_samples))

#---- ****plot ----
for(dem_group in unique(synthetic_continuous$group_label)){
  for(var_name in Z[, "label"]){
    subset <- synthetic_continuous %>% 
      filter(group_label == dem_group & label == var_name)
    ggplot(data = subset , aes(x = value)) + 
      geom_histogram(fill = "black", color = "black") + theme_minimal() + 
      xlab("Median") + ggtitle(var_name) + 
      geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                 size = 2)
    
    ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                             "posteriors/run_1/continuous_vars/median/", 
                             tolower(dem_group), "/", tolower(dem_group), "_", 
                             var_name, ".jpeg"), 
           width = 5, height = 3, units = "in")
  } 
}

#---- **skew ----
#---- ****overall ----
synthetic_continuous <- 
  matrix(0, nrow = nrow(Z), ncol = (num_samples + 1)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, num_samples), "var"))
synthetic_continuous[, "var"] <- Z[, "var"]

#skew from synthetic datasets
for(var in Z[, "var"]){
  subset <- synthetic_sample %>% dplyr::select(!!as.symbol(var), "sample") 
  subset[, var] <- subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
  synthetic_continuous[which(synthetic_continuous$var == var), 
                       seq(1, num_samples)] <- 
    t(subset %>% group_by(sample) %>% summarize_all(skewness) %>% 
        dplyr::select(var))
}

synthetic_continuous %<>% 
  mutate("truth" = 0, 
         "label" = Z[, "label"])

for(var in Z[, "var"]){
  synthetic_continuous[which(synthetic_continuous$var == var), "truth"] <- 
    skewness(unlist(ADAMS_subset[, var]))
}

synthetic_continuous %<>% pivot_longer(seq(1:num_samples))

#---- ****plot ----
for(var_name in Z[, "label"]){
  subset <- synthetic_continuous %>% filter(label == var_name)
  ggplot(data = subset , aes(x = value)) + 
    geom_histogram(fill = "black", color = "black") + theme_minimal() + 
    xlab("Skew") + ggtitle(var_name) + 
    geom_vline(xintercept = subset$truth, color = "#f2caaa" , size = 2)
  
  ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                           "posteriors/run_1/continuous_vars/skew/overall/", 
                           var_name, ".jpeg"), 
         width = 5, height = 3, units = "in")
} 

#---- ****by class ----
synthetic_continuous <- 
  matrix(0, nrow = nrow(Z)*4, ncol = (num_samples + 2)) %>% 
  as.data.frame() %>% set_colnames(c(seq(1, num_samples), "group", "var"))
synthetic_continuous[, "group"] <- rep(seq(1, 4), each = nrow(Z))
synthetic_continuous[, "var"] <- rep(Z[, "var"], 4)

#skewness from synthetic datasets
for(group in 1:4){
  subsample <- synthetic_sample %>% filter(Group == group)
  
  for(var in Z[, "var"]){
    subset <- subsample %>% dplyr::select(!!as.symbol(var), "sample") 
    subset[, var] <- subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
    synthetic_continuous[which(synthetic_continuous$group == group & 
                                 synthetic_continuous$var == var), 
                         seq(1, num_samples)] <- 
      t(subset %>% group_by(sample) %>% summarize_all(skewness) %>% 
          dplyr::select(var))
  }
}

synthetic_continuous %<>% 
  mutate("truth" = 0, 
         "label" = rep(Z[, "label"], 4), 
         "group_label" = case_when(group == 1 ~ "Normal", 
                                   group == 2 ~ "Other", 
                                   group == 3 ~ "MCI", 
                                   group == 4 ~ "Dementia")) %>% 
  mutate("color" = case_when(group_label == "Normal" ~ "#00a389", 
                             group_label == "Other" ~ "#28bed9", 
                             group_label == "MCI" ~ "#fdab00", 
                             group_label == "Dementia" ~ "#ff0000"))

for(group in unique(ADAMS_subset$Adem_dx_cat)){
  subsample <- ADAMS_subset %>% filter(Adem_dx_cat == group)
  
  for(var in Z[, "var"]){
    synthetic_continuous[which(synthetic_continuous$group_label == group & 
                                 synthetic_continuous$var == var), 
                         "truth"] <- skewness(unlist(subsample[, var]))
  }
}

synthetic_continuous %<>% pivot_longer(seq(1:num_samples))

#---- ****plot ----
for(dem_group in unique(synthetic_continuous$group_label)){
  for(var_name in Z[, "label"]){
    subset <- synthetic_continuous %>% 
      filter(group_label == dem_group & label == var_name)
    ggplot(data = subset , aes(x = value)) + 
      geom_histogram(fill = "black", color = "black") + theme_minimal() + 
      xlab("Skew") + ggtitle(var_name) + 
      geom_vline(xintercept = subset$truth, color = unique(subset$color), 
                 size = 2)
    
    ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                             "posteriors/run_1/continuous_vars/skew/", 
                             tolower(dem_group), "/", tolower(dem_group), "_", 
                             var_name, ".jpeg"), 
           width = 5, height = 3, units = "in")
  } 
}

#---- impairment classification ----
#truth
ADAMS_train[which(ADAMS_train$Adem_dx_cat == "Normal"), "Adem_dx_cat"] <- 
  "Unimpaired"
ADAMS_dementia_plot_data <- as.data.frame(table(ADAMS_train$Adem_dx_cat)) %>% 
  mutate("prop" = Freq/sum(Freq))

#synthetic
dem_sub <- synthetic_sample[, c("Group", "sample", "chain")] %>% 
  mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                   Group == 2 ~ "Other", 
                                   Group == 3 ~ "MCI", 
                                   TRUE ~ "Dementia"))

synthetic_dementia_plot_data <- 
  dem_sub %>% dplyr::count(chain, sample, Group_label) %>%
  group_by(chain, sample) %>%
  mutate(prop = n/sum(n)) %>% 
  mutate_at("sample", as.numeric) %>% 
  mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                             Group_label == "Other" ~ "#28bed9", 
                             Group_label == "MCI" ~ "#fdab00", 
                             TRUE ~ "#ff0000"))

#---- ****combined plot ----
combined_plot_data <- synthetic_dementia_plot_data %>% ungroup() %>% 
  group_by(chain, Group_label) %>% 
  summarize_at("n", list("mean" = mean, "lower" = ~ quantile(.x, probs = 0.025), 
                         "upper" = ~ quantile(.x, probs = 0.975))) %>% 
  mutate("truth" = ADAMS_dementia_plot_data$Freq) %>% 
  mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                             Group_label == "Other" ~ "#28bed9", 
                             Group_label == "MCI" ~ "#fdab00", 
                             Group_label == "Dementia" ~ "#ff0000")) %>% 
  mutate("chain" = paste0("Chain ", chain))


combined_plot_data$Group_label <- 
  factor(combined_plot_data$Group_label, 
         levels = c("Unimpaired", "Other", "MCI", "Dementia"))

ggplot(data = combined_plot_data, 
       aes(x = Group_label, y = mean, color = Group_label)) + 
  geom_point(aes(size = 1)) + theme_bw() + 
  geom_point(aes(x = Group_label, y = truth, color = Group_label, size = 1), 
             shape = 18) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.10) + 
  xlab("") + ylab("Mean Count") + theme(legend.position = "none") + 
  scale_color_manual(values = rev(c(combined_plot_data$color))) + 
  facet_wrap(facets = as.factor(combined_plot_data$chain), nrow = 2, ncol = 3) +
  ggtitle(paste0("95% Credible intervals from ", num_samples, 
                 " synthetic datasets"))

ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                         "posteriors/impairment_classes_combined_all_runs.jpeg"), 
       height = 5, width = 10, units = "in")

# #----****plot: prop ----
# percentiles <- tibble("class" = c("Unimpaired", "Other", "MCI", "Dementia"), 
# "upper" = 0, "lower" = 0, 
# "upper_count" = 0, "lower_count" = 0)
# 
# for(class in unique(synthetic_dementia_plot_data$Group_label)){
#   #props
#   percentiles[which(percentiles$class == class), "upper"] <- 
#     synthetic_dementia_plot_data %>% filter(Group_label == class) %>% 
#     ungroup() %>% summarise_at("prop", ~ quantile(.x, probs = 0.975))
#   
#   percentiles[which(percentiles$class == class), "lower"] <- 
#     synthetic_dementia_plot_data %>% filter(Group_label == class) %>% 
#     ungroup() %>% summarise_at("prop", ~ quantile(.x, probs = 0.025))
#   
#   #counts
#   percentiles[which(percentiles$class == class), "upper_count"] <- 
#     synthetic_dementia_plot_data %>% filter(Group_label == class) %>% 
#     ungroup() %>% summarise_at("n", ~ quantile(.x, probs = 0.975))
#   
#   percentiles[which(percentiles$class == class), "lower_count"] <- 
#     synthetic_dementia_plot_data %>% filter(Group_label == class) %>% 
#     ungroup() %>% summarise_at("n", ~ quantile(.x, probs = 0.025))
# }
# 
# for(class in unique(synthetic_dementia_plot_data$Group_label)){
#   subset <- synthetic_dementia_plot_data %>% filter(Group_label == class)
#   ggplot(data = subset) + 
#     geom_histogram(aes(x = prop), 
#                    color = unique(subset$color), fill = unique(subset$color)) + 
#     theme_minimal() + ggtitle(class) + xlab("Proportion") + ylab("Frequency") +
#     geom_vline(xintercept = ADAMS_dementia_plot_data[
#       which(ADAMS_dementia_plot_data$Var1 == class), "prop"], size = 2) + 
#     geom_vline(xintercept = 
#                  unlist(percentiles[which(percentiles$class == class), 
#                                     "upper"]), size = 1, linetype = "dashed") + 
#     geom_vline(xintercept = 
#                  unlist(percentiles[which(percentiles$class == class), 
#                                     "lower"]), size = 1, 
#                linetype = "dashed")
#   
#   
#   ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                            "posteriors/run_1/impairment_classes/prop/", class, 
#                            ".jpeg"), 
#          height = 3, width = 5, units = "in")
# }
# 
# #----****plot: count ----
# for(class in unique(synthetic_dementia_plot_data$Group_label)){
#   subset <- synthetic_dementia_plot_data %>% filter(Group_label == class)
#   ggplot(data = subset) + 
#     geom_histogram(aes(x = n), 
#                    color = unique(subset$color), fill = unique(subset$color)) + 
#     theme_minimal() + ggtitle(class) + xlab("Count") + ylab("Frequency") +
#     geom_vline(xintercept = ADAMS_dementia_plot_data[
#       which(ADAMS_dementia_plot_data$Var1 == class), "Freq"], size = 2) + 
#     geom_vline(xintercept = 
#                  unlist(percentiles[which(percentiles$class == class), 
#                                     "upper_count"]), size = 1, 
#                linetype = "dashed") + 
#     geom_vline(xintercept = unlist(percentiles[which(percentiles$class == class), 
#                                                "lower_count"]), size = 1, 
#                linetype = "dashed")
#   
#   
#   ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                            "posteriors/impairment_classes/count/", class, 
#                            ".jpeg"), 
#          height = 3, width = 5, units = "in")
# }



