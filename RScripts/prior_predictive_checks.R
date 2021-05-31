#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer", "devtools", "gifski", 
       "transformr")
install_github("thomasp85/gganimate")
library(gganimate)

#---- read in data ----
#---- **ADAMS ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"))

ADAMS_data <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                              "data/cleaned/ADAMS/ADAMS_subset_mixed.csv"))

#---- priors ----
#---- **latent classes ----
for(group in c("normal", "mci", "other")){
  assign(paste0(group, "_betas"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "latent_class_", group, "_cov.csv")))
  
  assign(paste0(group, "_preds"), get(paste0(group, "_betas"))$preds)
}

#---- **contingency cells ----
alpha_0_dist <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                  "bootstrap_cell_counts.csv")) 

#--- **beta and sigma ----
priors_beta <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                               "priors/priors_beta.csv")) 
prior_V_inv <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                               "priors/priors_V_inv.csv")) 
prior_Sigma <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                               "priors/priors_Sigma.csv")) 

#---- **hyperparameters ----
#DOF for inverse wishart
nu_0 <- 65
#scaling for inverse wishart as variance of Beta
kappa_0 <- 0.85

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- unique(c(normal_preds, other_preds, mci_preds, "ETHNIC_label"))

#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- cbind(c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
             "proxy_cog", "ANDELCOR", "Aiadla", "Abmi"), 
           c("Age", "Total MMSE", "Serial 7s", "Immediate Word Recall", 
             "Wordlist Recall (Yes)", "Story Recall I", "Proxy Cognition (Avg)", 
             "Delayed Word Recall", "IADLs", "BMI")) %>% 
  set_colnames(c("var", "label"))

synthetic_sample <- ADAMS_train %>% 
  dplyr::select("HHIDPN", all_of(vars), "Adem_dx_cat") %>% 
  #pre-allocate columns
  mutate("Group" = 0, "p_normal" = 0, "p_other" = 0, "p_mci" = 0)

ADAMS_data %<>% 
  dplyr::select("HHIDPN", all_of(W), all_of(Z[, "var"]), "Adem_dx_cat") %>% 
  filter(HHIDPN %in% ADAMS_train$HHIDPN)

ADAMS_means <- colMeans(ADAMS_data %>% dplyr::select(all_of(Z[, "var"])))
ADAMS_sds <- apply(ADAMS_data %>% dplyr::select(all_of(Z[, "var"])), 2, sd)

#---- contrasts matrix ----
A = do.call(cbind, list(
  #intercept
  rep(1, 6),
  #race/ethnicity main effect: Black
  rep(c(1, 0, 0), 2),
  #race/ethnicity main effect: Hispanic
  rep(c(0, 1, 0), 2),
  #stroke main effect
  rep(c(0, 1), each = 3)))

generate_data <- function(){
  #---- latent class ----
  group = 1
  synthetic_sample[, "Group"] <- 0
  
  for(model in c("normal", "other", "mci")){
    subset_index <- which(synthetic_sample$Group == 0)
    random_draw <- sample(seq(1, 10000), size = 1)
    
    prior_betas <- as.vector(get(paste0(model, "_betas"))[, random_draw])
    prior_cov <- matrix(unlist(get(paste0(model, "_cov"))[, random_draw]), 
                        nrow = nrow(prior_betas))
    
    betas <- mvrnorm(n = 1, mu = t(prior_betas), Sigma = prior_cov)
    
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
  mu <- matrix(0, ncol = 4*6, nrow = 10) %>%
    set_colnames(apply(expand.grid(seq(1, 4), seq(1, 6)), 1, paste0,
                       collapse = ":"))
  
  for(i in 1:4){
    #---- **contingency cells ----
    subset <- synthetic_sample %>% filter(Group == i)
    prior_counts <- alpha_0_dist[, c(sample(seq(1, 10000), size = 1), 
                                     ncol(alpha_0_dist))] %>% 
      filter(group_number == i)
    
    #---- **p(contingency table cell) ----
    pi <- rdirichlet(1, alpha = as.numeric(unlist(prior_counts[, 1])))
    
    #---- **contingency table count ----
    contingency_table <- rmultinom(n = 1, size = nrow(subset), prob = pi)
    UtU <- diag(contingency_table[, 1])
    
    #---- **draw new UtU if needed ----
    while(det(t(A) %*% UtU %*% A) < 1e-9){
      random_draw <- sample(seq(1, 10000), size = 1)
      new_counts <- alpha_0_dist[, c(random_draw, ncol(alpha_0_dist))] %>% 
        filter(group_number == i)
      
      UtU <- diag(unlist(new_counts[, 1]))
    }
    
    #---- **make U matrix ----
    U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
    
    for(j in 1:nrow(contingency_table)){
      if(contingency_table[j, 1] == 0){next}
      if(j == 1){
        index = 1
      } else{
        index = sum(contingency_table[1:(j - 1), 1]) + 1
      }
      U[index:(index - 1 + contingency_table[j, 1]), j] <- 1
    }
    
    #---- **draw Sigma_0----
    random_draw <- sample(seq(1, 10000), size = 1)
    Sigma_prior <- prior_Sigma[, c(random_draw, ncol(prior_Sigma))] %>% 
      filter(group_number == i)
    sig_Y <- riwish(v = (nu_0), S = matrix(unlist(Sigma_prior[, 1]), 
                                           nrow = nrow(Z)))
    
    #---- **beta_0 ----
    V_0_inv <- prior_V_inv[, c(random_draw, ncol(prior_V_inv))] %>% 
      filter(group_number == i)
    beta_0 <- priors_beta[, c(random_draw, ncol(priors_beta))] %>% 
      filter(group_number == i)
    
    #as matrices
    V_0_inv <- matrix(unlist(V_0_inv[, 1]), nrow = 4)
    beta_0 <- matrix(unlist(beta_0[, 1]), nrow = nrow(V_0_inv))
    
    #---- **draw beta | Sigma----
    beta_Sigma_Y <- matrix.normal(beta_0, solve(V_0_inv), sig_Y/kappa_0)
    
    #---- **compute mu ----
    mu[, paste0(i, ":", seq(1, 6))] <-
      t(A %*% matrix(beta_Sigma_Y, nrow = ncol(A), ncol = nrow(Z),
                     byrow = FALSE))
    
    #---- **draw data ----
    #reformat contingency table
    table <- contingency_table %>% as.data.frame() %>%
      cbind(do.call(cbind, list(
        #Black              #Hispanic           #Stroke
        rep(c(1, 0, 0), 2), rep(c(0, 1, 0), 2), c(rep(0, 3), rep(1, 3))))) %>%
      set_colnames(c("Count", W))
    
    for(j in 1:nrow(table)){
      if(table[j, "Count"] == 0){next}
      if(j == 1){
        index = 1
      } else{
        index = sum(table[1:(j - 1), "Count"]) + 1
      }
      #Z (continuous data)
      if(table[j, "Count"] == 1){
        subset[index:(index - 1 + table[j, "Count"]), Z[, "var"]] <-
          t(as.matrix(mvrnorm(n = table[j, "Count"],
                              mu = mu[, paste0(i, ":", j)], Sigma = sig_Y)))
      } else{
        subset[index:(index - 1 + table[j, "Count"]), Z[, "var"]] <-
          mvrnorm(n = table[j, "Count"],
                  mu = mu[, paste0(i, ":", j)],
                  Sigma = sig_Y)
      }
    }
    assign(paste0("Z_", i), subset[, all_of(Z[, "var"])])
  }
  
  #---- **return ----
  return(list("Group" = synthetic_sample$Group,
              "Z_normal" = Z_1 %>% mutate("color" = "#00a389"), 
              "Z_other" = Z_2 %>% mutate("color" = "#28bed9"), 
              "Z_mci" = Z_3 %>% mutate("color" = "#fdab00"),
              "Z_dementia" = Z_4 %>% mutate("color" = "#ff0000")))
}

#---- multiruns ----
start <- Sys.time()
runs = 1000
synthetic <- replicate(runs, generate_data(), simplify = FALSE) 
stop <- Sys.time() - start

#---- plots ----
#---- **dem class ----
ADAMS_data[which(ADAMS_data$Adem_dx_cat == "Normal"), "Adem_dx_cat"] <- 
  "Unimpaired"
ADAMS_dementia_plot_data <- as.data.frame(table(ADAMS_data$Adem_dx_cat)) %>% 
  mutate("prop" = Freq/sum(Freq))

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
  mutate_at("name", as.numeric) %>% 
  mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                             Group_label == "Other" ~ "#28bed9", 
                             Group_label == "MCI" ~ "#fdab00", 
                             TRUE ~ "#ff0000"))

for(class in unique(synthetic_dementia_plot_data$Group_label)){
  subset <- synthetic_dementia_plot_data %>% filter(Group_label == class)
  ggplot(data = subset) + 
    geom_histogram(aes(x = prop), 
                   color = unique(subset$color), fill = unique(subset$color)) + 
    theme_minimal() + ggtitle(class) + xlab("Proportion") + ylab("Count") +
    geom_vline(xintercept = ADAMS_dementia_plot_data[
      which(ADAMS_dementia_plot_data$Var1 == class), "prop"], size = 2)
  
  ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                           "priors/impairment_classes/", class, ".jpeg"), 
         height = 3, width = 5, units = "in")
}

#---- **class-specific continuous ----
ADAMS_data[which(ADAMS_data$Adem_dx_cat == "Unimpaired"), "Adem_dx_cat"] <- 
  "Normal"

for(class in unique(ADAMS_train$Adem_dx_cat)){
  true_data <- ADAMS_data %>% 
    dplyr::select(c(all_of(Z[, "var"]), "Adem_dx_cat")) %>% 
    filter(Adem_dx_cat == class) %>% mutate("color" = "black")
  
  continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
  
  for(i in 1:length(continuous_list)){
    continuous_list[[i]] <- continuous_list[[i]] %>% 
      mutate("run" = i, "type" = "synthetic") %>% 
      rbind(., true_data %>% dplyr::select(-"Adem_dx_cat") %>% 
              mutate("run" = i, "type" = "ADAMS"))
  }
  
  continuous_list %<>% do.call(rbind, .) %>% as.data.frame() 
  
  for(var in Z[, "var"]){
    data <- continuous_list[, c(var, "run", "type", "color")]
    #unstandardize
    synthetic_subset <- data %>% filter(type == "synthetic")
    data[which(data$type == "synthetic"), var] <- 
      synthetic_subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
    
    continuous_plot <- 
      ggplot(data = data, aes(color = type, fill = type)) + 
      geom_density(aes(x = data[, 1]), alpha = 0.5) + 
      theme_minimal() + xlab(Z[which(Z[, "var"] == var), "label"]) + 
      scale_color_manual(values = rev(unique(data$color))) + 
      scale_fill_manual(values = rev(unique(data$color))) + 
      transition_states(data$run, transition_length = 1, state_length = 1) +
      labs(title = "Synthetic {round(frame_time)}") + transition_time(run) +
      ease_aes('linear')
    
    animate(continuous_plot, fps = 2, height = 4, width = 5, units = "in", 
            res = 150, renderer = gifski_renderer())
    
    anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/",
                                "priors/continuous_vars/", tolower(class), "/", 
                                var, ".gif"),
              animation = last_animation(),
              renderer = gifski_renderer())
  }
}

#---- OLD ----

# synthetic_dementia_class_plot <- geom_bar(mapping = aes(x = factor(Group_label, 
#                                     levels = c("Unimpaired", "MCI", 
#                                                "Dementia", "Other")), y = prop, 
#                          fill = factor(Group_label, 
#                                        levels = c("Unimpaired", "MCI", 
#                                                   "Dementia", "Other"))), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + 
#   ylim(c(0, 1)) + theme(legend.position = "none")  + 
#   scale_fill_manual(values = wes_palette("Darjeeling1")[c(2, 3, 1, 5)]) + 
#   #gganimate
#   transition_states(name, transition_length = 1, state_length = 1) +
#   labs(title = "Synthetic {round(frame_time)}", 
#        x = "Impairment Class", y = "Proportion") + transition_time(name) + 
#   ease_aes('linear')
# 
# animate(synthetic_dementia_class_plot, 
#         duration = max(synthetic_dementia_plot_data$name), fps = 1, 
#         height = 4, width = 4, units = "in", res = 150, 
#         renderer = gifski_renderer())
# 
# anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                             "priors/synthetic_dem_class.gif"), 
#           animation = last_animation(), 
#           renderer = gifski_renderer())

# #---- ***ADAMS ----
# ADAMS_dementia_class_plot <- 
#   ggplot(data = ADAMS_dementia_plot_data) + 
#   geom_bar(mapping = 
#              aes(x = factor(Var1, 
#                             levels = c("Unimpaired", "MCI", 
#                                        "Dementia", "Other")), y = prop, 
#                  fill = factor(Var1, levels = c("Unimpaired", "MCI", 
#                                                 "Dementia", "Other"))), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + xlab("Impairment Class") + ylab("Proportion") + 
#   ylim(c(0, 1)) + theme(legend.position = "none")  + 
#   scale_fill_manual(values = wes_palette("Darjeeling1")[c(2, 3, 1, 5)]) + 
#   ggtitle("ADAMS")
# 
# ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
#                          "ADAMS_dem_class.png"), device = "jpeg", 
#        width = 4, height = 4, units = "in")

# #---- **categorical: race ----
# cats <- table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke) %>% 
#   as.data.frame()
# 
# cat_sub <- lapply(synthetic, "[[", "contingency_table") %>% 
#   do.call(cbind, .) %>% set_colnames(seq(1, runs)) %>% as.data.frame() %>% 
#   mutate("Race_Eth" = cats$Var1, "Stroke" = cats$Var2) %>%
#   pivot_longer(seq(1, runs)) 
# 
# synthetic_race_plot_data <- 
#   cat_sub %>% group_by(name, Race_Eth) %>% summarize_at("value", sum) %>%
#   mutate(prop = value/sum(value)) %>% mutate_at("name", as.integer)
# 
# synthetic_race_plot <- 
#   ggplot(data = synthetic_race_plot_data) + 
#   geom_bar(mapping = aes(x = factor(Race_Eth), y = prop, 
#                          fill = factor(Race_Eth)), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + 
#   ylim(c(0, 1)) + theme(legend.position = "none")  + 
#   scale_fill_manual(values = wes_palette("Darjeeling2")) + 
#   #gganimate
#   transition_states(name, transition_length = 1, state_length = 1) +
#   labs(title = "Synthetic {round(frame_time)}", 
#        x = "Race/Ethnicity", y = "Proportion") + transition_time(name) + 
#   ease_aes('linear')
# 
# animate(synthetic_race_plot, 
#         duration = max(synthetic_race_plot_data$name), fps = 10, 
#         height = 4, width = 4, units = "in", res = 150, 
#         renderer = gifski_renderer())
# 
# anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                             "priors/synthetic_race.gif"), 
#           animation = last_animation(), 
#           renderer = gifski_renderer())
# 
# #---- ***ADAMS ----
# ADAMS_race_plot_data <- as.data.frame(table(ADAMS_subset$ETHNIC_label)) %>% 
#   mutate("prop" = Freq/sum(Freq))
# 
# ADAMS_race_plot <- 
#   ggplot(data = ADAMS_race_plot_data) + 
#   geom_bar(mapping = aes(x = factor(Var1), y = prop, fill = factor(Var1)), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
#   ylim(c(0, 1)) + theme(legend.position = "none")  + 
#   scale_fill_manual(values = wes_palette("Darjeeling2")) + 
#   ggtitle("ADAMS")
# 
# ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
#                          "ADAMS_race.png"), device = "jpeg", 
#        width = 4, height = 4, units = "in")
# 
# #---- **categorical: stroke ----
# synthetic_stroke_plot_data <- 
#   cat_sub %>% group_by(name, Stroke) %>% summarize_at("value", sum) %>%
#   mutate(prop = value/sum(value)) %>% mutate_at("name", as.integer)
# 
# synthetic_stroke_plot <- 
#   ggplot(data = synthetic_stroke_plot_data) + 
#   geom_bar(mapping = aes(x = factor(Stroke), y = prop, 
#                          fill = factor(Stroke)), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + 
#   ylim(c(0, 1)) + theme(legend.position = "none")  + 
#   scale_fill_manual(values = wes_palette("Darjeeling2")) + 
#   #gganimate
#   transition_states(name, transition_length = 1, state_length = 1) +
#   labs(title = "Synthetic {round(frame_time)}", 
#        x = "Stroke", y = "Proportion") + transition_time(name) + 
#   ease_aes('linear')
# 
# animate(synthetic_stroke_plot, 
#         duration = max(synthetic_stroke_plot_data$name), fps = 1, 
#         height = 4, width = 4, units = "in", res = 150, 
#         renderer = gifski_renderer())
# 
# anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                             "priors/synthetic_stroke.gif"), 
#           animation = last_animation(), 
#           renderer = gifski_renderer())
# 
# #---- ***ADAMS ----
# ADAMS_stroke_plot_data <- as.data.frame(table(ADAMS_subset$Astroke)) %>% 
#   mutate("prop" = Freq/sum(Freq))
# 
# ADAMS_stroke_plot <- 
#   ggplot(data = ADAMS_stroke_plot_data) + 
#   geom_bar(mapping = aes(x = factor(Var1), y = prop, fill = factor(Var1)), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + xlab("Stroke") + ylab("Proportion") + 
#   ylim(c(0, 1)) + theme(legend.position = "none")  + 
#   scale_fill_manual(values = wes_palette("Darjeeling2")) + 
#   ggtitle("ADAMS")
# 
# ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
#                          "ADAMS_stroke.png"), device = "jpeg", 
#        width = 4, height = 4, units = "in")
# 
# #---- **categorical: race x stroke ----
# synthetic_race_stroke_plot_data <- 
#   cat_sub %>% group_by(name, Race_Eth, Stroke) %>% 
#   summarize_at("value", sum) %>%
#   group_by(name, Stroke) %>%
#   mutate(prop = value/sum(value)) %>% mutate_at("name", as.integer)
# 
# synthetic_race_stroke_plot <- 
#   ggplot(data = synthetic_race_stroke_plot_data) + 
#   geom_bar(mapping = aes(x = factor(Stroke), y = prop, 
#                          fill = factor(Race_Eth)), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + ylim(c(0, 1)) + theme(legend.position = "none") +
#   scale_fill_manual(values = wes_palette("Darjeeling2")) + 
#   #gganimate
#   transition_states(name, transition_length = 1, state_length = 1) +
#   labs(title = "Synthetic {round(frame_time)}", 
#        x = "Stroke", y = "Proportion") + transition_time(name) + 
#   ease_aes('linear')
# 
# animate(synthetic_race_stroke_plot, 
#         duration = max(synthetic_race_stroke_plot_data$name), fps = 1, 
#         height = 4, width = 4, units = "in", res = 150, 
#         renderer = gifski_renderer())
# 
# anim_save(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
#                             "priors/synthetic_race_stroke.gif"), 
#           animation = last_animation(), 
#           renderer = gifski_renderer())
# 
# #---- ***ADAMS ----
# ADAMS_race_stroke_plot_data <- 
#   as.data.frame(table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke)) %>% 
#   group_by(Var2) %>%
#   mutate("prop" = Freq/sum(Freq))
# 
# ADAMS_race_stroke_plot <- 
#   ggplot(data = ADAMS_race_stroke_plot_data) + 
#   geom_bar(mapping = aes(x = factor(Var2), y = prop, fill = factor(Var1)), 
#            stat = "identity", position = "dodge") + 
#   theme_minimal() + xlab("Stroke") + ylab("Proportion") + 
#   ylim(c(0, 1)) + 
#   scale_fill_manual(values = wes_palette("Darjeeling2"), 
#                     name = "Race/Ethnicity") + ggtitle("ADAMS")
# 
# ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/priors/", 
#                          "ADAMS_race_stroke.png"), device = "jpeg", 
#        width = 4, height = 4, units = "in")

