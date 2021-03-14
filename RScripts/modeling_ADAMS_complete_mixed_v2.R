#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "locfit", 
       "wesanderson", "RColorBrewer")

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character()))

#---- **models for priors ----
Unimpaired_prior <- readRDS(here::here("priors", "normal_model_25.rds"))
Other_prior <- readRDS(here::here("priors", "other_model_25.rds"))
MCI_prior <- readRDS(here::here("priors", "MCI_model_25.rds"))

Unimpaired_preds <- names(coefficients(Unimpaired_prior))
Unimpaired_preds[which(Unimpaired_preds == "ETHNIC_labelBlack")] <- "Black"
Unimpaired_preds[which(Unimpaired_preds == "ETHNIC_labelHispanic")] <- 
  "Hispanic"
Other_preds <- names(coefficients(Other_prior))
MCI_preds <- names(coefficients(MCI_prior))

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- unique(c(Unimpaired_preds, Other_preds, MCI_preds, "ETHNIC_label"))

synthetic_sample <- ADAMS_subset %>% 
  mutate("Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0),
         #Add intercept
         "(Intercept)" = 1) %>% 
  dplyr::select("HHIDPN", all_of(vars)) %>% 
  #use complete data for now
  na.omit() %>% 
  #pre-allocate columns
  mutate("Group" = 0, "p_Unimpaired" = 0, "p_Other" = 0, "p_MCI" = 0)

# #Sanity check
# table(analytical_sample$ETHNIC_label, analytical_sample$Black, useNA = "ifany")
# table(analytical_sample$ETHNIC_label, analytical_sample$Hispanic,
#       useNA = "ifany")
# 
# colSums(is.na(analytical_sample))

#---- all-way contingency table ----
#We have one small cell-- Hispanics who have had a stroke
cross_class_label <- table(synthetic_sample$ETHNIC_label, 
                           synthetic_sample$Astroke) %>% as.data.frame() %>% 
  mutate("Stroke" = ifelse(Var2 == 0, "No Stroke", "Stroke")) %>% 
  unite("Cell Label", c("Var1", "Stroke"), sep = " | ")

#---- Bayes Stuff ----
#---- **parameters ----
#number of runs
B = 2

#---- **chain storage ----
Unimpaired_beta_chain <- 
  matrix(nrow = length(coefficients(Unimpaired_prior)), ncol = B) %>% 
  set_rownames(paste0("Unimpaired:", Unimpaired_preds))
Other_beta_chain <- 
  matrix(nrow = length(coefficients(Other_prior)), ncol = B) %>% 
  set_rownames(paste0("Other:", Other_preds))
MCI_beta_chain <- 
  matrix(nrow = length(coefficients(MCI_prior)), ncol = B) %>% 
  set_rownames(paste0("MCI:", MCI_preds))

latent_class_chain <- matrix(nrow = 4, ncol = B) %>% 
  set_rownames(c("Unimpaired", "Other", "MCI", "Dementia"))

pi_chain <- matrix(nrow = nrow(cross_class_label), ncol = 4*B) %>% 
  set_colnames(apply(expand.grid(seq(1, B), seq(1:4)), 1, paste, 
                     collapse = ":")) %>% 
  set_rownames(cross_class_label$`Cell Label`)

#---- **priors ----
#uninformative
alpha_0 <- rep(1, nrow(cross_class_label))

#---- **sampling ----
for(b in 1:B){
  #---- ****latent class betas ----
  Unimpaired_beta_chain[, b] <- 
    mvrnorm(n = 1, mu = coefficients(Unimpaired_prior), 
            Sigma = vcov(Unimpaired_prior))
  Other_beta_chain[, b] <- 
    mvrnorm(n = 1, mu = coefficients(Other_prior), Sigma = vcov(Other_prior))
  MCI_beta_chain[, b] <- 
    mvrnorm(n = 1, mu = coefficients(MCI_prior), Sigma = vcov(MCI_prior))
  
  #---- ****group membership ----
  group = 1
  synthetic_sample[, "Group"] <- 0
  for(model in c("Unimpaired", "Other", "MCI")){
    subset_index <- which(synthetic_sample$Group == 0)
    
    synthetic_sample[subset_index, paste0("p_", model)] <- 
      expit(as.matrix(synthetic_sample[subset_index, 
                                       get(paste0(model, "_preds"))]) %*% 
              as.matrix(get(paste0(model, "_beta_chain"))[, b]))
    
    synthetic_sample[subset_index, "Group"] <- 
      rbernoulli(n = length(subset_index), 
                 p = synthetic_sample[subset_index, paste0("p_", model)])*group
    
    group = group + 1
  }
  
  synthetic_sample[which(synthetic_sample$Group == 0), "Group"] <- 4
  
  #---- ****group: summary ----
  latent_class_chain[, b] <- 
    table(synthetic_sample$Group)/sum(table(synthetic_sample$Group))
  
  #---- ****p(contingency table cell) ----
  for(i in 1:4){
    subset <- synthetic_sample %>% filter(Group == i) 
    posterior_counts <- alpha_0 + 
      table(subset$ETHNIC_label, subset$Astroke) %>% as.data.frame() %>% 
      dplyr::select("Freq") %>% unlist()
    
    pi_chain[, paste0(b, ":", i)] <- rdirichlet(1, alpha = posterior_counts)
  }
  
}

#---- **plots ----
extended_pallette14 <- colorRampPalette(wes_palette("Darjeeling1"))(14)
extended_pallette6 <- colorRampPalette(wes_palette("Darjeeling1"))(6)

#---- ****beta chains ----
beta_plot_data <- 
  do.call(cbind, list(t(Unimpaired_beta_chain), t(Other_beta_chain), 
                      t(MCI_beta_chain))) %>% as.data.frame() %>%
  mutate("run" = seq(1:B)) %>% 
  pivot_longer(-c("run"), names_to = c("Group", "Predictor"), 
               names_sep = ":", values_to = "beta")

beta_chain_plot <- 
  ggplot(data = beta_plot_data, 
         aes(x = as.factor(run), y = beta, colour = Predictor)) +       
  geom_line(aes(group = Predictor)) + 
  facet_grid(rows = vars(factor(Group, 
                                levels = c("Unimpaired", "MCI", "Other")))) + 
  theme_bw() + xlab("Run") + 
  scale_color_manual(values = rev(extended_pallette14))

ggsave(filename = "beta_chain.jpeg", plot = beta_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 6, units = "in", device = "jpeg")

#---- ****latent class chain ----
latent_class_data <- t(latent_class_chain) %>% as.data.frame() %>%
  mutate("run" = seq(1:B)) %>% 
  pivot_longer(-c("run"), names_to = c("Group"), values_to = "prob") %>% 
  arrange(desc(prob)) %>%
  mutate_at("Group", as.factor)
latent_class_data$Group <- fct_relevel(latent_class_data$Group, 
                                       paste0(unique(latent_class_data$Group)))

latent_class_chain_plot <- 
  ggplot(data = latent_class_data, 
         aes(x = as.factor(run), y = prob, colour = Group)) +       
  geom_line(aes(group = Group)) + 
  theme_minimal() + xlab("Run") + ylab("Proportion of Sample") +  
  scale_color_manual(values = rev(wes_palette("Darjeeling1")))  

ggsave(filename = "latent_class_chain.jpeg", plot = latent_class_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 3, units = "in", device = "jpeg")

#---- ****pi chain ----
pi_chain_data <- pi_chain %>% as.data.frame() %>% rownames_to_column("Cell") %>% 
  pivot_longer(-c("Cell"), names_to = c("Run", "Group"), names_sep = ":", 
               values_to = "probability") %>% 
  mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                   Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                   Group == 4 ~ "Dementia")) %>% 
  mutate_if(is.character, as.factor) 

pi_chain_plot <- ggplot(data = pi_chain_data, 
                        aes(x = as.factor(Run), y = probability, 
                            colour = Cell)) +       
  geom_line(aes(group = Cell)) + 
  theme_minimal() + xlab("Run") + ylab("Probability of cell membership") +  
  scale_color_manual(values = rev(extended_pallette6)) + 
  facet_grid(rows = vars(factor(Group_label, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other")))) + theme_bw() 

ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 5, units = "in", device = "jpeg")



