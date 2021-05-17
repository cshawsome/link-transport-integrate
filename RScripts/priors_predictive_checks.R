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
                                "data/cleaned/ADAMS/ADAMS_train.csv"))

#---- priors ----
#---- **latent classes ----
for(group in c("normal", "mci", "other")){
  assign(paste0(group, "_betas"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "latent_class_", group, "_betas.csv")))
  assign(paste0(group, "_cov"), 
         read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                         "latent_class_", group, "_cov.csv")))
}

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
alpha_0_dist <- read_csv(here::here("priors", "bootstrap_cell_counts.csv")) 

#--- **beta and sigma ----
beta_prior <- readRDS(here::here("priors", "beta.rds"))
V_inv_prior <- readRDS(here::here("priors", "V_inv.rds"))
Sigma_prior <- readRDS(here::here("priors", "Sigma.rds"))

#---- **hyperparameters ----
#DOF for inverse wishart
nu_0 <- 200
#scaling for inverse wishart as variance of Beta
kappa_0 <- 1

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- unique(c(Unimpaired_preds, Other_preds, MCI_preds, "ETHNIC_label"))

#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- cbind(c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
             "proxy_cog", "ANDELCOR", "Aiadla", "Abmi"), 
           c("Age", "Total MMSE", "Serial 7s", "Immediate Word Recall", 
             "Wordlist Recall (Yes)", "Story Recall I", "Proxy Cognition (Avg)", 
             "Delayed Word Recall", "IADLs", "BMI")) %>% 
  set_colnames(c("var", "label"))

synthetic_sample <- ADAMS_subset %>% 
  dplyr::select("HHIDPN", all_of(vars), "group_class") %>% 
  #use complete data for now
  na.omit() %>% 
  #Z-score continuous
  mutate_at(all_of(Z[, "var"]), scale) %>%
  #transform to correct type
  mutate_at("Astroke", as.numeric) %>% 
  #pre-allocate columns
  mutate("Group" = 0, "p_Unimpaired" = 0, "p_Other" = 0, "p_MCI" = 0)

ADAMS_subset %<>% dplyr::select("HHIDPN", all_of(vars), "group_class") %>% 
  #use complete data for now
  na.omit()

ADAMS_means <- colMeans(ADAMS_subset %>% dplyr::select(all_of(Z[, "var"])))
ADAMS_sds <- apply(ADAMS_subset %>% dplyr::select(all_of(Z[, "var"])), 2, sd)

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
  #---- true class props ----
  true_props <- as.data.frame(table(ADAMS_subset$group_class)) %>% 
    mutate("prop" = Freq/nrow(ADAMS_subset)) %>% 
    set_rownames(c("4", "3", "2", "1"))
  
  #---- latent class ----
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
  contingency_table <- matrix(0, ncol = 4, nrow = 6) %>% set_colnames(seq(1, 4))
  mu <- matrix(0, ncol = 4*6, nrow = 10) %>% 
    set_colnames(apply(expand.grid(seq(1, 4), seq(1, 6)), 1, paste0, 
                       collapse = ":"))
  
  for(i in 1:4){
    #---- **contingency cells ----
    subset <- synthetic_sample %>% filter(Group == i) 
    prior_counts <- 
      alpha_0_dist[, sample(seq(1, ncol(alpha_0_dist)), size = 1)]*
      true_props[paste0(i), "prop"] 
    
    #---- **p(contingency table cell) ----
    pi <- rdirichlet(1, alpha = as.numeric(unlist(prior_counts)))
    
    #---- **contingency table count ----
    contingency_table[, i] <- contingency_table[, i] + 
      rmultinom(n = 1, size = nrow(subset), prob = pi)
    
    #---- **make U matrix ----
    U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
    
    for(j in 1:nrow(contingency_table)){
      if(contingency_table[j, i] == 0){next}
      if(j == 1){
        index = 1
      } else{
        index = sum(contingency_table[1:(j - 1), i]) + 1
      }
      U[index:(index - 1 + contingency_table[j, i]), j] <- 1
    }
    
    UtU <- diag(contingency_table[, i])
    
    # if(i %in% c(2, 3)){
    #   assign(paste0("UtU_", i), UtU)
    # }
    # 
    # #---- **pool UtU if needed ----
    # if(det(t(A) %*% UtU %*% A) == 0){
    #   if(exists(paste0("UtU_", (i-1)))){
    #     UtU <- UtU + get(paste0("UtU_", (i-1)))
    #   } else{
    #     UtU <- UtU + get(paste0("UtU_", (i+1))) 
    #   }
    # }
    #
    
    #---- **draw Sigma_0----
    sig_Y <- riwish(v = (nu_0), S = Sigma_prior)
    
    #---- **beta_0 ----
    V_0_inv <- matrix(V_inv_prior[, i], nrow = 4, ncol = 4)
    beta_0 <- matrix(beta_prior[, i], nrow = nrow(V_0_inv), ncol = ncol(sig_Y))
    
    #---- **draw beta | Sigma----
    beta_Sigma_Y <- matrix.normal(beta_0, solve(V_0_inv), sig_Y/kappa_0)
    
    #---- **compute mu ----
    mu[, paste0(i, ":", seq(1, 6))] <- 
      t(A %*% matrix(beta_Sigma_Y, nrow = ncol(A), ncol = nrow(Z), 
                     byrow = FALSE))
    
    #---- **draw data ----
    #reformat contingency table
    table <- contingency_table[, i] %>% as.data.frame() %>% 
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
        subset[index:(index - 1 + table[j, "Count"]), colnames(sig_Y)] <- 
          t(as.matrix(mvrnorm(n = table[j, "Count"],
                              mu = mu[, paste0(i, ":", j)], Sigma = sig_Y)))
      } else{
        subset[index:(index - 1 + table[j, "Count"]), colnames(sig_Y)] <- 
          mvrnorm(n = table[j, "Count"],
                  mu = mu[, paste0(i, ":", j)], 
                  Sigma = sig_Y)
      }
    }
    assign(paste0("Z_", i), subset[, all_of(Z[, "var"])])
  }
  
  #---- **return ----
  return(list("Group" = synthetic_sample$Group, 
              "contingency_table" = contingency_table, 
              "Z_unimpaired" = Z_1, "Z_other" = Z_2, "Z_MCI" = Z_3, 
              "Z_dementia" = Z_4))
}

#---- multiruns ----
start <- Sys.time()
runs = 10000
synthetic <- replicate(runs, generate_data(), simplify = FALSE) 
stop <- Sys.time() - start

#---- plots ----
#---- **dem class ----
ADAMS_dementia_plot_data <- as.data.frame(table(ADAMS_subset$group_class)) %>% 
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

#---- **class-specific categorical ----
ADAMS_categorical_plot_data <- 
  as.data.frame(table(ADAMS_subset$ETHNIC_label, ADAMS_subset$Astroke, 
                      ADAMS_subset$group_class)) %>% 
  mutate("Var2_label" = case_when(Var2 == 0 ~ "No Stroke", 
                                  TRUE ~ "Stroke")) %>% 
  unite("cat", c("Var1", "Var2_label"), sep = " + ") %>% 
  set_colnames(c("cat", "Var2", "Group_label", "Truth"))

categorical_sub <- lapply(synthetic, "[[", "contingency_table") %>% 
  do.call(cbind, .) %>% as.data.frame() %>% 
  set_colnames(apply(expand_grid(seq(1, runs), seq(1, 4)), 1, 
                     paste0, collapse = ":")) %>%
  mutate("cat" = c("Black + No Stroke", "Hispanic + No Stroke", 
                   "White + No Stroke", "Black + Stroke", "Hispanic + Stroke", 
                   "White + Stroke")) %>% 
  pivot_longer(-c("cat"), names_to = c("Run", "Group"), names_sep = ":", 
               values_to = "count") %>% 
  mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                   Group == 2 ~ "Other", 
                                   Group == 3 ~ "MCI", 
                                   TRUE ~ "Dementia")) %>% 
  mutate("color" = case_when(Group_label == "Unimpaired" ~ "#00a389", 
                             Group_label == "Other" ~ "#28bed9", 
                             Group_label == "MCI" ~ "#fdab00", 
                             TRUE ~ "#ff0000")) %>% 
  left_join(., ADAMS_categorical_plot_data, by = c("cat", "Group_label"))

for(group in unique(categorical_sub$Group_label)){
  subset <- categorical_sub %>% filter(Group_label == group)
  for(category in unique(categorical_sub$cat)){
    smaller_subset <- subset %>% filter(cat == category)
    ggplot(data = smaller_subset , aes(x = count)) + 
      geom_histogram(fill = "black", color = "black") + theme_minimal() + 
      xlab("Count") + ggtitle(category) + 
      geom_vline(xintercept = smaller_subset$Truth, 
                 color = smaller_subset$color, size = 2)
    
    ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                             "priors/cell_counts/group_specific/", 
                             tolower(group), "/", tolower(group), "_",  
                             category, ".jpeg"), 
           width = 5, height = 3, units = "in")
  }
}

#---- **class-specific continuous ----
for(class in unique(ADAMS_subset$group_class)){
  ADAMS_data <- ADAMS_subset %>% 
    dplyr::select(c(all_of(Z[, "var"]), "group_class")) %>% 
    filter(group_class == class)
  
  continuous_list <- lapply(synthetic, "[[", paste0("Z_", tolower(class))) 
  
  for(i in 1:length(continuous_list)){
    continuous_list[[i]] <- continuous_list[[i]] %>% 
      mutate("run" = i, "type" = "synthetic") %>% 
      rbind(., ADAMS_data %>% dplyr::select(-"group_class") %>% 
              mutate("run" = 0, "type" = "ADAMS"))
  }
  
  continuous_list %<>% do.call(rbind, .) %>% as.data.frame()
  
  for(var in Z[, "var"]){
    data <- continuous_list[, c(var, "run", "type")]
    #unstandardize
    synthetic_subset <- data %>% filter(type == "synthetic")
    data[which(data$type == "synthetic"), var] <- 
      synthetic_subset[, var]*ADAMS_sds[var] + ADAMS_means[var]
    
    continuous_plot <- ggplot() + 
      geom_density(color = as.factor(data$type)) + 
      geom_density(data = ADAMS_data, aes(x = unlist(ADAMS_data[, var]))) +
      theme_minimal() + xlab(Z[which(Z[, "var"] == var), "label"]) + 
      transition_states(data$run, transition_length = 1, state_length = 1) +
      labs(title = "Synthetic {round(frame_time)}",
           x = var, y = "density") + transition_time(name) +
      ease_aes('linear')
  }
  
  
  
  
  ggplot(data = plot_data %>% 
           filter(type == "x"), 
         aes(x = value, color = group, fill = group)) + 
    geom_density(alpha = 0.5) + theme_minimal() + xlab(label) + 
    scale_color_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)]) + 
    scale_fill_manual(values = wes_palette("Darjeeling1")[c(1, 3, 5, 2)])
  geom_density(mapping = aes(x = factor(Group_label,
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

#---- ***ADAMS ----


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

