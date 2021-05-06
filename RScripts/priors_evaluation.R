#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "devtools", "gifski")
install_github("thomasp85/gganimate")
library(gganimate)

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character(), 
                                          #Do not standardize these
                                          Astroke = col_character())) %>% 
  #Z-score continuous
  mutate_if(is.numeric, scale) %>%
  #transform to correct type
  mutate_at("Astroke", as.numeric) %>% 
  mutate("group_class" = 
           case_when(Adem_dx_cat %in% 
                       c("Dementia", "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia",
                     Adem_dx_cat == "Normal" ~ "Unimpaired",
                     TRUE ~ Adem_dx_cat)) 

#---- priors ----
#---- **latent classes ----
Unimpaired_prior <- readRDS(here::here("priors", "normal_model_25.rds"))
Other_prior <- readRDS(here::here("priors", "other_model_25.rds"))
MCI_prior <- readRDS(here::here("priors", "MCI_model_25.rds"))

Unimpaired_preds <- names(coefficients(Unimpaired_prior))
Unimpaired_preds[which(Unimpaired_preds == "ETHNIC_labelBlack")] <- "Black"
Unimpaired_preds[which(Unimpaired_preds == "ETHNIC_labelHispanic")] <- 
  "Hispanic"

Other_preds <- names(coefficients(Other_prior))
Other_preds[which(Other_preds == "ETHNIC_labelBlack")] <- "Black"
Other_preds[which(Other_preds == "ETHNIC_labelHispanic")] <- "Hispanic"

MCI_preds <- names(coefficients(MCI_prior))
MCI_preds[which(MCI_preds == "ETHNIC_labelBlack")] <- "Black"
MCI_preds[which(MCI_preds == "ETHNIC_labelHispanic")] <- "Hispanic"

#---- **contingency cells ----
alpha_0 <- read_csv(here::here("priors", "contingency_cell_counts.csv")) %>%
  set_colnames(c("Var1", "Var2", "Freq", "4_prior_count", "3_prior_count",
                 "1_prior_count", "2_prior_count")) %>%
  mutate_at(paste0(seq(1, 4), "_prior_count"), function(x) 0.30*x)

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- unique(c(Unimpaired_preds, Other_preds, MCI_preds, "ETHNIC_label"))

synthetic_sample <- ADAMS_subset %>% 
  mutate("Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0),
         #Add intercept
         "(Intercept)" = 1) %>% 
  dplyr::select("HHIDPN", all_of(vars), "group_class") %>% 
  #use complete data for now
  na.omit() %>% 
  #pre-allocate columns
  mutate("Group" = 0, "p_Unimpaired" = 0, "p_Other" = 0, "p_MCI" = 0)

ADAMS_subset <- synthetic_sample

#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", 
       "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

generate_data <- function(){
  #---- **latent class ----
  group = 1
  synthetic_sample[, "Group"] <- 0
  
  for(model in c("Unimpaired", "Other", "MCI")){
    subset_index <- which(synthetic_sample$Group == 0)
    prior_model <- get(paste0(model, "_prior"))
    betas <- mvrnorm(n = 1, mu = coefficients(prior_model), 
                     Sigma = vcov(prior_model))
    
    synthetic_sample[subset_index, paste0("p_", model)] <- 
      expit(as.matrix(synthetic_sample[subset_index, 
                                       get(paste0(model, "_preds"))]) %*% 
              as.matrix(betas))
    
    synthetic_sample[subset_index, "Group"] <- 
      rbernoulli(n = length(subset_index), 
                 p = synthetic_sample[subset_index, paste0("p_", model)])*group
    
    group = group + 1
  }
  synthetic_sample[which(synthetic_sample$Group == 0), "Group"] <- 4
  
  #pre-allocate
  contingency_table <- matrix(rep(0, 6), ncol = 1)
    
  for(i in 1:4){
    #---- **contingency cells ----
    subset <- synthetic_sample %>% filter(Group == i) 
    prior_counts <- alpha_0[, paste0(i, "_prior_count")] 
    
    #---- ****p(contingency table cell) ----
    pi <- rdirichlet(1, alpha = as.numeric(unlist(prior_counts)))
    
    #---- ****contingency table count ----
    contingency_table <- contingency_table + 
      rmultinom(n = 1, size = nrow(subset), prob = pi)
  }
  
  #---- **return ----
  return(list("Group" = synthetic_sample$Group, 
              "contingency_table" = contingency_table))
}

#---- multiruns ----
runs = 20

synthetic <- replicate(runs, generate_data(), simplify = FALSE) 

#---- plots ----
#---- **dem class ----
dem_sub <- lapply(synthetic, "[[", "Group") %>% do.call(cbind, .) %>% 
  set_colnames(seq(1, runs)) %>% as.data.frame() %>%
  pivot_longer(everything()) %>% 
  mutate("Group_label" = case_when(value == 1 ~ "Unimpaired", 
                                   value == 2 ~ "Other", 
                                   value == 3 ~ "MCI", 
                                   TRUE ~ "Dementia"))

synthetic_dementia_plot_data <- 
  dem_sub %>% dplyr::count(name, Group_label) %>%
  group_by(name) %>%
  mutate(prop = n/sum(n)) %>% 
  mutate_at("name", as.numeric)

synthetic_dementia_class_plot <- 
  ggplot(data = synthetic_dementia_plot_data) + 
  geom_bar(mapping = aes(x = factor(Group_label, 
                                    levels = c("Unimpaired", "MCI", 
                                               "Dementia", "Other")), y = prop, 
                         fill = factor(Group_label, 
                                       levels = c("Unimpaired", "MCI", 
                                                  "Dementia", "Other"))), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(2, 3, 1, 5)]) + 
  #gganimate
  transition_states(name, transition_length = 1, state_length = 1) +
  labs(title = "Synthetic {round(frame_time)}", 
       x = "Impairment Class", y = "Proportion") + transition_time(name) + 
  ease_aes('linear')

animate(synthetic_dementia_class_plot, 
        duration = max(synthetic_dementia_plot_data$name), fps = 1, 
        height = 4, width = 4, units = "in", res = 150, 
        renderer = gifski_renderer())

anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                            "priors/synthetic_dem_class.gif"), 
          animation = last_animation(), 
          renderer = gifski_renderer())

#---- ***ADAMS ----
ADAMS_dementia_plot_data <- as.data.frame(table(ADAMS_subset$group_class)) %>% 
  mutate("prop" = Freq/sum(Freq))

ADAMS_dementia_class_plot <- 
  ggplot(data = ADAMS_dementia_plot_data) + 
  geom_bar(mapping = 
             aes(x = factor(Var1, 
                            levels = c("Unimpaired", "MCI", 
                                       "Dementia", "Other")), y = prop, 
                 fill = factor(Var1, levels = c("Unimpaired", "MCI", 
                                                "Dementia", "Other"))), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Impairment Class") + ylab("Proportion") + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(2, 3, 1, 5)]) + 
  ggtitle("ADAMS")

ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
                         "ADAMS_dem_class.png"), device = "jpeg", 
       width = 4, height = 4, units = "in")

#---- **categorical: race ----
cats <- table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke) %>% 
  as.data.frame()

cat_sub <- lapply(synthetic, "[[", "contingency_table") %>% 
  do.call(cbind, .) %>% set_colnames(seq(1, runs)) %>% as.data.frame() %>% 
  mutate("Race_Eth" = cats$Var1, "Stroke" = cats$Var2) %>%
  pivot_longer(seq(1, runs)) 

synthetic_race_plot_data <- 
  cat_sub %>% group_by(name, Race_Eth) %>% summarize_at("value", sum) %>%
  mutate(prop = value/sum(value)) %>% mutate_at("name", as.integer)

synthetic_race_plot <- 
  ggplot(data = synthetic_race_plot_data) + 
  geom_bar(mapping = aes(x = factor(Race_Eth), y = prop, 
                         fill = factor(Race_Eth)), 
                         stat = "identity", position = "dodge") + 
  theme_minimal() + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling2")) + 
  #gganimate
  transition_states(name, transition_length = 1, state_length = 1) +
  labs(title = "Synthetic {round(frame_time)}", 
       x = "Race/Ethnicity", y = "Proportion") + transition_time(name) + 
  ease_aes('linear')

animate(synthetic_race_plot, 
        duration = max(synthetic_race_plot_data$name), fps = 1, 
        height = 4, width = 4, units = "in", res = 150, 
        renderer = gifski_renderer())

anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                            "priors/synthetic_race.gif"), 
          animation = last_animation(), 
          renderer = gifski_renderer())

#---- ***ADAMS ----
ADAMS_race_plot_data <- as.data.frame(table(ADAMS_subset$ETHNIC_label)) %>% 
  mutate("prop" = Freq/sum(Freq))

ADAMS_race_plot <- 
  ggplot(data = ADAMS_race_plot_data) + 
  geom_bar(mapping = aes(x = factor(Var1), y = prop, fill = factor(Var1)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling2")) + 
  ggtitle("ADAMS")

ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
                         "ADAMS_race.png"), device = "jpeg", 
       width = 4, height = 4, units = "in")

#---- **categorical: stroke ----
synthetic_stroke_plot_data <- 
  cat_sub %>% group_by(name, Stroke) %>% summarize_at("value", sum) %>%
  mutate(prop = value/sum(value)) %>% mutate_at("name", as.integer)

synthetic_stroke_plot <- 
  ggplot(data = synthetic_stroke_plot_data) + 
  geom_bar(mapping = aes(x = factor(Stroke), y = prop, 
                         fill = factor(Stroke)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling2")) + 
  #gganimate
  transition_states(name, transition_length = 1, state_length = 1) +
  labs(title = "Synthetic {round(frame_time)}", 
       x = "Stroke", y = "Proportion") + transition_time(name) + 
  ease_aes('linear')

animate(synthetic_stroke_plot, 
        duration = max(synthetic_stroke_plot_data$name), fps = 1, 
        height = 4, width = 4, units = "in", res = 150, 
        renderer = gifski_renderer())

anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                            "priors/synthetic_stroke.gif"), 
          animation = last_animation(), 
          renderer = gifski_renderer())

#---- ***ADAMS ----
ADAMS_stroke_plot_data <- as.data.frame(table(ADAMS_subset$Astroke)) %>% 
  mutate("prop" = Freq/sum(Freq))

ADAMS_stroke_plot <- 
  ggplot(data = ADAMS_stroke_plot_data) + 
  geom_bar(mapping = aes(x = factor(Var1), y = prop, fill = factor(Var1)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Stroke") + ylab("Proportion") + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling2")) + 
  ggtitle("ADAMS")

ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
                         "ADAMS_stroke.png"), device = "jpeg", 
       width = 4, height = 4, units = "in")

#---- **categorical: race x stroke ----
synthetic_stroke_plot_data <- 
  cat_sub %>% group_by(name, Stroke) %>% summarize_at("value", sum) %>%
  mutate(prop = value/sum(value)) %>% mutate_at("name", as.integer)

synthetic_stroke_plot <- 
  ggplot(data = synthetic_stroke_plot_data) + 
  geom_bar(mapping = aes(x = factor(Stroke), y = prop, 
                         fill = factor(Stroke)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling2")) + 
  #gganimate
  transition_states(name, transition_length = 1, state_length = 1) +
  labs(title = "Synthetic {round(frame_time)}", 
       x = "Stroke", y = "Proportion") + transition_time(name) + 
  ease_aes('linear')

animate(synthetic_stroke_plot, 
        duration = max(synthetic_stroke_plot_data$name), fps = 1, 
        height = 4, width = 4, units = "in", res = 150, 
        renderer = gifski_renderer())

anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                            "priors/synthetic_stroke.gif"), 
          animation = last_animation(), 
          renderer = gifski_renderer())

#---- ***ADAMS ----
ADAMS_stroke_plot_data <- as.data.frame(table(ADAMS_subset$Astroke)) %>% 
  mutate("prop" = Freq/sum(Freq))

ADAMS_stroke_plot <- 
  ggplot(data = ADAMS_stroke_plot_data) + 
  geom_bar(mapping = aes(x = factor(Var1), y = prop, fill = factor(Var1)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Stroke") + ylab("Proportion") + 
  ylim(c(0, 1)) + theme(legend.position = "none")  + 
  scale_fill_manual(values = wes_palette("Darjeeling2")) + 
  ggtitle("ADAMS")

ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
                         "ADAMS_stroke.png"), device = "jpeg", 
       width = 4, height = 4, units = "in")

