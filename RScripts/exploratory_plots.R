#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "wesanderson", "devtools", 
       "gmodels")
install_github("thomasp85/patchwork")

#---- read in data ----
#based on analysis in priors_latent_classes.R
#Categorical vars (notation from Schafer 1997)
W <- cbind(c("ETHNIC_label", "Astroke"), 
           c("Race/Ethnicity", "Stroke")) %>% 
  set_colnames(c("var", "label"))

#Continuous vars (notation from Schafer 1997)
Z <- cbind(c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
             "proxy_cog", "ANDELCOR", "Aiadla", "Abmi"), 
           c("Age", "Total MMSE", "Serial 7s", "Immediate Word Recall", 
             "Wordlist Recall (Yes)", "Story Recall I", "Proxy Cognition (Avg)", 
             "Delayed Word Recall", "IADLs", "BMI")) %>% 
  set_colnames(c("var", "label"))

group <- c("Adem_dx_cat")

#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character())) %>% 
  dplyr::select(c("HHIDPN", all_of(group), 
                  all_of(W[, "var"]), all_of(Z[, "var"]))) %>% 
  na.omit() %>%
  mutate("Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0),
         #Add intercept
         "(Intercept)" = 1) %>% 
  mutate("group_class" = 
           case_when(Adem_dx_cat %in% 
                       c("Dementia", "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia",
                     Adem_dx_cat == "Normal" ~ "Unimpaired",
                     TRUE ~ Adem_dx_cat))

analytical_sample <- ADAMS_subset  %>% 
  #don't standardize this
  mutate_at("Astroke", as.character) %>%
  #Z-score continuous
  mutate_if(is.numeric, scale) %>%
  #transform to correct type
  mutate_at("Astroke", as.numeric) 

#---- all-way contingency table ----
cross_class_label <- table(analytical_sample$ETHNIC_label, 
                           analytical_sample$Astroke) %>% as.data.frame()

# #How many are missing from this table?-- only 144! 
# sum(cross_class_label$Freq)

# #---- summary stats ----
# #Looking at the same sample as the impaired vs. unimpaired model
# normal_model_data <- ADAMS_subset %>% 
#   dplyr::select(c("AAGE", "ETHNIC_label", "ANMSETOT", "ANSER7T", 
#                   "ANIMMCR", "ANRECYES", "ANWM1TOT", 
#                   "proxy_cog", "Adem_dx_cat")) %>% na.omit() %>% 
#   mutate("Aunimpaired" = ifelse(Adem_dx_cat == "Normal", 1, 0))
# 
# #How many of each race/ethnicity in the sample
# CrossTable(normal_model_data$ETHNIC_label, useNA = "ifany")
# 
# #How many of each race/ethnicity classified as cognitively normal
# CrossTable(normal_model_data$ETHNIC_label, normal_model_data$Aunimpaired, 
#            useNA = "ifany", prop.chisq = FALSE)

#---- plots ----
# #---- **race x age ----
# race_by_age_bar <- 
#   ggplot(data = normal_model_data) + 
#   geom_bar(mapping = aes(x = factor(AAGE), y = ..count../sum(..count..), 
#                          fill = factor(ETHNIC_label)), 
#            position = "dodge") +
#   theme_minimal() + xlab("Age") + ylab("Proportion") +
#   guides(fill = guide_legend(title = "Race/Ethnicity")) +
#   scale_fill_manual(values = rev(wes_palette("Darjeeling1")))
# 
# race_by_age_dens <-
#   ggplot(data = normal_model_data, aes(x = AAGE, fill = ETHNIC_label)) + 
#   geom_density(color = NA, alpha = 0.4, position = 'identity') +
#   scale_fill_manual(values = rev(wes_palette("Darjeeling1"))) + 
#   theme_minimal() + xlab("Age") + ylab("Density") + 
#   guides(fill = guide_legend(title = "Race/Ethnicity")) 
# 
# #---- **race x MMSE ----
# race_by_MMSE_bar <-
#   ggplot(data = normal_model_data) +
#   geom_bar(mapping = aes(x = factor(ANMSETOT), y = ..count../sum(..count..),
#                          fill = factor(ETHNIC_label)),
#            position = "dodge") +
#   theme_minimal() + xlab("MMSE") + ylab("Proportion") +
#   guides(fill = guide_legend(title = "Race/Ethnicity")) +
#   scale_fill_manual(values = rev(wes_palette("Darjeeling1")))
# 
# race_by_MMSE_dens <-
#   ggplot(data = normal_model_data, aes(x = ANMSETOT, fill = ETHNIC_label)) + 
#   geom_density(color = NA, alpha = 0.4, position = 'identity') +
#   scale_fill_manual(values = rev(wes_palette("Darjeeling1"))) + 
#   theme_minimal() + xlab("MMSE") + ylab("Density") + 
#   guides(fill = guide_legend(title = "Race/Ethnicity")) 
# 
# #---- **race x proxy cognition ----
# race_by_proxy_cog_dens <-
#   ggplot(data = normal_model_data, aes(x = proxy_cog, fill = ETHNIC_label)) + 
#   geom_density(color = NA, alpha = 0.4, position = 'identity') +
#   scale_fill_manual(values = rev(wes_palette("Darjeeling1"))) + 
#   theme_minimal() + xlab("Proxy Cognition") + ylab("Density") + 
#   guides(fill = guide_legend(title = "Race/Ethnicity")) 
# 
# #---- **patchwork plot ----
# (race_by_age_bar + race_by_age_dens)/(race_by_MMSE_bar + race_by_MMSE_dens)/
#   race_by_proxy_cog_dens
# 
# ggsave(filename = "unimpaired_two_way_by_race.jpeg", plot = last_plot(), 
#        path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                      "prelim_analyses/latent_class_unimpaired/"), 
#        width = 10, height = 8, units = "in", device = "jpeg")

#---- mixture plots ----
merged_data <- left_join(ADAMS_subset, analytical_sample, by = "HHIDPN")

for(var in Z[, "var"]){
  label <- Z[which(Z[, "var"] == var), "label"]
  plot_data <- merged_data %>% 
    dplyr::select(contains(c(var, "group_class"))) %>% 
    pivot_longer(everything(), names_to = c(".value", "type"), 
                 names_pattern = "(.*).(.)") %>% 
    set_colnames(c("type", "value", "group")) 
  plot_data[, "value"] <- unlist(plot_data[, "value"])
  plot_data[which(plot_data$type == "y"), "type"] <- "z"
  
  plotX <- ggplot(data = plot_data %>% 
                    filter(type == "x"), 
                  aes(x = value, color = group, fill = group)) + 
    geom_density(alpha = 0.5) + theme_minimal() + xlab(label) + 
    scale_color_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)]) + 
    scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)])
  
  plotZ <- ggplot(data = plot_data %>% 
                    filter(type == "z"), 
                  aes(x = value, color = group, fill = group)) + 
    geom_density(alpha = 0.5) + theme_minimal() + xlab(label) + 
    scale_color_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)]) + 
    scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)])
  
  if(var != Z[nrow(Z), "var"]){
    plotZ <- plotZ + theme(legend.position = "none") 
    plotX <- plotX + theme(legend.position = "none") 
  }
  
  assign(paste0(var, "_plotX"), plotX)
  assign(paste0(var, "_plotZ"), plotZ)
}

#---- **patchwork plot ----
continuous_var_plot_names <- paste0(Z, "_plotX")

((((get(continuous_var_plot_names[1]) + get(continuous_var_plot_names[2]) + 
      get(continuous_var_plot_names[3])) /
     (get(continuous_var_plot_names[4]) + get(continuous_var_plot_names[5]) + 
        get(continuous_var_plot_names[6])))) / 
    (get(continuous_var_plot_names[7]) + get(continuous_var_plot_names[8]) + 
       get(continuous_var_plot_names[9]))) / 
  get(continuous_var_plot_names[10])

ggsave(filename = "ADAMS_mix_X.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                     "prelim_analyses/ADAMSA"), width = 12, height = 12, 
       units = "in", device = "jpeg")

continuous_var_plot_names <- paste0(Z, "_plotZ")

((((get(continuous_var_plot_names[1]) + get(continuous_var_plot_names[2]) + 
      get(continuous_var_plot_names[3])) /
     (get(continuous_var_plot_names[4]) + get(continuous_var_plot_names[5]) + 
        get(continuous_var_plot_names[6])))) / 
    (get(continuous_var_plot_names[7]) + get(continuous_var_plot_names[8]) + 
       get(continuous_var_plot_names[9]))) / 
  get(continuous_var_plot_names[10])

ggsave(filename = "ADAMS_mix_Z.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                     "prelim_analyses/ADAMSA"), width = 12, height = 12, 
       units = "in", device = "jpeg")

#---- continuous plots ----
for(var in Z[, "var"]){
  label <- Z[which(Z[, "var"] == var), "label"]
  plot_data <- merged_data %>% 
    dplyr::select(contains(c(var, "group_class"))) %>% 
    pivot_longer(everything(), names_to = c(".value", "type"), 
                 names_pattern = "(.*).(.)") %>% 
    set_colnames(c("type", "value", "group")) 
  plot_data[, "value"] <- unlist(plot_data[, "value"])
  plot_data[which(plot_data$type == "y"), "type"] <- "z"
  
  plotX <- ggplot(data = plot_data %>% 
                    filter(type == "x"), 
                  aes(x = value, color = group, fill = group)) + 
    geom_density() + theme_minimal() + xlab(label) + 
    scale_color_manual(values = rep("black", 4)) + 
    scale_fill_manual(values = rep("black", 4)) + 
    theme(legend.position = "none")
  
  assign(paste0(var, "_plotX"), plotX)
}

continuous_var_plot_names <- paste0(Z[, "var"], "_plotX")

((((get(continuous_var_plot_names[1]) + get(continuous_var_plot_names[2]) + 
      get(continuous_var_plot_names[3])) /
     (get(continuous_var_plot_names[4]) + get(continuous_var_plot_names[5]) + 
        get(continuous_var_plot_names[6])))) / 
    (get(continuous_var_plot_names[7]) + get(continuous_var_plot_names[8]) + 
       get(continuous_var_plot_names[9]))) / 
  get(continuous_var_plot_names[10])

ggsave(filename = "ADAMS_X.jpeg", plot = last_plot(), 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                     "prelim_analyses/ADAMSA"), width = 12, height = 12, 
       units = "in", device = "jpeg")

#---- OLD CODE ----
#---- plots ----
#Categorical Variables: Race/Ethnicity, IADLs, Ever/never stroke
#Continuous Variables: Age, MMSE, Delayed Word Recall, Immediate Word Recall, 
# Word list recog (Yes), proxy cognition, 
#Create labeled data (if not already created)
analytical_sample %<>% 
  mutate("GENDER_label" = ifelse(GENDER == 1, "Male", "Female"), 
         "ETHNIC_label" = case_when(ETHNIC == 1 ~ "White", 
                                    ETHNIC == 2 ~ "Black", 
                                    TRUE ~ "Hispanic"), 
         "IADLA_label" = case_when(IADLA == 1 ~ "None", 
                                   IADLA == 2 ~ "One", 
                                   IADLA == 3 ~ "Two", 
                                   IADLA == 4 ~ "Three")) %>% 
  mutate_at(c("GENDER_label", "ETHNIC_label", "IADLA_label"), as.factor) 
analytical_sample$IADLA_label <- fct_relevel(analytical_sample$IADLA_label, 
                                             c("None", "One", "Two", "Three"))


#---- *marginal ----
#---- **summary stats ----
sex_gender_plot <- table(analytical_sample$GENDER_label) %>% 
  as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) 

race_eth_plot <- table(analytical_sample$ETHNIC_label) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) 

IADL_plot <- table(analytical_sample$IADLA_label) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) 

Age_plot <- table(exp(analytical_sample$log_AAGE)) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate_at("Var1", as.factor)

Edyrs_plot <- table(analytical_sample$EDYRS) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate_at("Var1", as.factor)

MMSE_plot <- table(analytical_sample$ANMSETOT) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate_at("Var1", as.factor)

#---- **categorical plots ----
sex_gender <- 
  ggplot(data = sex_gender_plot) + 
  geom_bar(mapping = aes(x = factor(Var1), y = Prop, fill = factor(Var1)), 
           stat = "identity") + 
  theme_minimal() + xlab("Sex/Gender") + ylab("Proportion") + 
  theme(legend.position = "none") + ylim(c(0, 1)) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

race_eth <- 
  ggplot(data = race_eth_plot) + 
  geom_bar(mapping = aes(x = factor(Var1), y = Prop, fill = factor(Var1)), 
           stat = "identity") + ylim(c(0, 1)) +
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

IADLs <- 
  ggplot(data = IADL_plot) + 
  geom_bar(mapping = aes(x = factor(Var1), y = Prop, fill = factor(Var1)), 
           stat = "identity") + ylim(c(0, 1)) +
  theme_minimal() + xlab("Difficulty with IADLs") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **continuous plots ----
Age <- ggplot(data = Age_plot) + 
  geom_bar(mapping = aes(x = Var1, y = Prop), 
           color = rev(wes_palette("Darjeeling1"))[1],
           fill = rev(wes_palette("Darjeeling1"))[1],
           stat = "identity") + 
  theme_minimal() + xlab("Age") + ylab("Proportion") + 
  theme(legend.position = "none")

Edyrs <- ggplot(data = Edyrs_plot) + 
  geom_bar(mapping = aes(x = Var1, y = Prop), 
           color = rev(wes_palette("Darjeeling1"))[2],
           fill = rev(wes_palette("Darjeeling1"))[2],
           stat = "identity") + 
  theme_minimal() + xlab("Years of Education") + ylab("Proportion") + 
  theme(legend.position = "none")

MMSE <- ggplot(data = MMSE_plot) + 
  geom_bar(mapping = aes(x = Var1, y = Prop), 
           color = rev(wes_palette("Darjeeling1"))[4],
           fill = rev(wes_palette("Darjeeling1"))[4],
           stat = "identity") + 
  theme_minimal() + xlab("MMSE") + ylab("Proportion") + 
  theme(legend.position = "none")

#---- **patchwork plot ----
(((sex_gender + race_eth + IADLs)/Age)/Edyrs)/MMSE

ggsave(filename = "marginal_dists.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
       width = 8, height = 8, units = "in", device = "jpeg")

#---- *2-way categorical ----
#Sex/Gender by Race/Ethnicity
sex_by_race <- 
  ggplot(data = analytical_sample) + 
  geom_bar(mapping = aes(x = factor(GENDER_label), y = ..count../sum(..count..), 
                         fill = factor(ETHNIC_label)), 
           position = "dodge") + ylim(0, 1) +
  theme_minimal() + xlab("Sex/Gender") + ylab("Count") +
  guides(fill = guide_legend(title = "Race/Ethnicity")) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

sex_by_IADL <- 
  ggplot(data = analytical_sample) + 
  geom_bar(mapping = aes(x = factor(GENDER_label), y = ..count../sum(..count..), 
                         fill = factor(IADLA_label)), 
           position = "dodge") + ylim(0, 1) +
  theme_minimal() + xlab("Sex/Gender") + ylab("Count") +
  guides(fill = guide_legend(title = "IADL")) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

race_by_IADL <- ggplot(data = analytical_sample) + 
  geom_bar(mapping = aes(x = factor(ETHNIC_label), y = ..count../sum(..count..), 
                         fill = factor(IADLA_label)), 
           position = "dodge") + ylim(0, 1) +
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Count") +
  guides(fill = guide_legend(title = "IADL")) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **patchwork plot ----
sex_by_race + sex_by_IADL + race_by_IADL

ggsave(filename = "2way_cat_dists.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
       width = 12, height = 3, units = "in", device = "jpeg")

#---- *3-way categorical ----
IADL_by_sex_and_race <- 
  ggplot(data = analytical_sample, aes(x = GENDER_label, y = ETHNIC_label, 
                                       color = IADLA_label)) + 
  geom_point(position = "jitter", alpha = 0.75) + 
  theme_minimal() + ylab("Race/Ethnicity") + xlab("Sex/Gender") + 
  guides(color = guide_legend(title = "IADL")) + 
  scale_color_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "3way_cat_dists.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
       width = 8, height = 8, units = "in", device = "jpeg")

#---- *cont | categorical ----
for(i in 1:nrow(cross_class_label)){
  if(cross_class_label[i, "Freq"] != 0){
    gender <- cross_class_label[i, "Var1"]
    ethnic <- cross_class_label[i, "Var2"]
    iadl <- cross_class_label[i, "Var3"]
    count <- cross_class_label[i, "Freq"]
    
    data_subset <- analytical_sample %>% 
      filter(GENDER_label == gender & ETHNIC_label == ethnic & 
               IADLA_label == iadl)
    
    #summary stats
    age_plot <- table(exp(data_subset$log_AAGE)) %>% as.data.frame() %>% 
      mutate("Prop" = Freq/nrow(data_subset)) %>% 
      mutate_at("Var1", as.character) %>% 
      mutate_at("Var1", as.numeric)
    missing_ages <- which(!seq(min(exp(analytical_sample$log_AAGE)), 
                               max(exp(analytical_sample$log_AAGE))) %in% 
                            age_plot$Var1) + 69
    age_plot %<>% rbind(as.matrix(cbind(missing_ages, NA, NA)) %>% 
                          set_colnames(c("Var1", "Freq", "Prop"))) %>% 
      mutate_at("Var1", as.factor)
    
    edyrs_plot <- table(data_subset$EDYRS) %>% as.data.frame() %>% 
      mutate("Prop" = Freq/nrow(data_subset)) %>% 
      mutate_at("Var1", as.character) %>% 
      mutate_at("Var1", as.numeric)
    missing_edyrs <- which(!seq(min(analytical_sample$EDYRS), 
                                max(analytical_sample$EDYRS)) %in% 
                             edyrs_plot$Var1) - 1
    edyrs_plot %<>% rbind(as.matrix(cbind(missing_edyrs, NA, NA)) %>% 
                            set_colnames(c("Var1", "Freq", "Prop"))) %>% 
      mutate_at("Var1", as.factor)
    
    mmse_plot <- table(data_subset$ANMSETOT) %>% as.data.frame() %>% 
      mutate("Prop" = Freq/nrow(data_subset)) %>% 
      mutate_at("Var1", as.character) %>% 
      mutate_at("Var1", as.numeric)
    missing_mmse <- which(!seq(min(analytical_sample$ANMSETOT), 
                               max(analytical_sample$ANMSETOT)) %in% 
                            mmse_plot$Var1) - 1
    mmse_plot %<>% rbind(as.matrix(cbind(missing_mmse, NA, NA)) %>% 
                           set_colnames(c("Var1", "Freq", "Prop"))) %>% 
      mutate_at("Var1", as.factor)
    
    #plots
    age <- ggplot(data = age_plot) + 
      geom_bar(mapping = aes(x = Var1, y = Prop), 
               color = rev(wes_palette("Darjeeling1"))[1],
               fill = rev(wes_palette("Darjeeling1"))[1],
               stat = "identity") + 
      theme_minimal() + xlab("Age") + ylab("Proportion") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position = "none") + 
      ggtitle(paste(paste(gender, ethnic, iadl, sep = " | "), "| n =", count)) 
    
    edyrs <- ggplot(data = edyrs_plot) + 
      geom_bar(mapping = aes(x = Var1, y = Prop), 
               color = rev(wes_palette("Darjeeling1"))[2],
               fill = rev(wes_palette("Darjeeling1"))[2],
               stat = "identity") + 
      theme_minimal() + xlab("Years of Education") + ylab("Proportion") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position = "none") + 
      ggtitle(paste(paste(gender, ethnic, iadl, sep = " | "), "| n =", count)) 
    
    mmse <- ggplot(data = mmse_plot) + 
      geom_bar(mapping = aes(x = Var1, y = Prop), 
               color = rev(wes_palette("Darjeeling1"))[4],
               fill = rev(wes_palette("Darjeeling1"))[4],
               stat = "identity") + 
      theme_minimal() + xlab("MMSE") + ylab("Proportion") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position = "none") + 
      ggtitle(paste(paste(gender, ethnic, iadl, sep = " | "), "| n =", count)) 
    
    #---- **patchwork plot ----
    age + edyrs + mmse
    
    ggsave(filename = paste0("cont_given_cat", i, ".jpeg"), plot = last_plot(), 
           path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
           width = 12, height = 3, units = "in", device = "jpeg")
  }
}


