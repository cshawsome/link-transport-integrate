#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "wesanderson", "RColorBrewer")

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

#Categorical vars (notation from Schafer 1997)
W <- c("Black", "Hispanic", "Astroke")
#Continuous vars (notation from Schafer 1997)
Z <- c("AAGE", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", "ANWM1TOT", 
       "proxy_cog", "ANDELCOR", "Aiadla", "Abmi")

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
  unite("Cell Label", c("Var1", "Stroke"), sep = " | ", remove = FALSE)

#---- Bayes Stuff ----
#---- **parameters ----
#number of runs
B = 10

#categorical vars contrasts matrix
A = do.call(cbind, list(
  #intercept
  rep(1, nrow(cross_class_label)),
  #race/ethnicity first main effect contrast
  c(rep(1, 2), rep(0, 2), rep(-1, 2)),
  #race/ethnicity second main effect contrast
  c(rep(0, 2), rep(1, 2), rep(-1, 2)),
  #stroke main effect contrast
  rep(c(1, -1), 3)))

Sigma_multiplier <- c(1, 1.3, 1.3, 1.7)

#---- **chain storage ----
Unimpaired_gamma_chain <- 
  matrix(nrow = length(coefficients(Unimpaired_prior)), ncol = B) %>% 
  set_rownames(paste0("Unimpaired:", Unimpaired_preds))
Other_gamma_chain <- 
  matrix(nrow = length(coefficients(Other_prior)), ncol = B) %>% 
  set_rownames(paste0("Other:", Other_preds))
MCI_gamma_chain <- 
  matrix(nrow = length(coefficients(MCI_prior)), ncol = B) %>% 
  set_rownames(paste0("MCI:", MCI_preds))

latent_class_chain <- matrix(nrow = 4, ncol = B) %>% 
  set_rownames(c("Unimpaired", "Other", "MCI", "Dementia"))

pi_chain <- matrix(nrow = nrow(cross_class_label), ncol = 4*B) %>% 
  set_colnames(apply(expand.grid(seq(1, 4), seq(1:B)), 1, paste, 
                     collapse = ":")) %>% 
  set_rownames(cross_class_label$`Cell Label`)

Sigma_chain <- matrix(nrow = length(Z), ncol = 4*B) %>%
  set_colnames(apply(expand.grid(seq(1, 4), seq(1:B)), 1, paste,
                     collapse = ":")) %>%
  set_rownames(Z)

mu_chain <-
  matrix(nrow = length(Z), ncol = 4*nrow(cross_class_label)*B) %>%
  set_colnames(apply(
    expand.grid(seq(1, 4), seq(1:nrow(cross_class_label)), seq(1:B)), 1, paste,
    collapse = ":")) %>% set_rownames(Z)

#---- **priors ----
#uninformative
alpha_0 <- rep(1, nrow(cross_class_label))

#---- **sampling ----
for(b in 1:B){
  #---- ****latent class gammas ----
  Unimpaired_gamma_chain[, b] <- 
    mvrnorm(n = 1, mu = coefficients(Unimpaired_prior), 
            Sigma = vcov(Unimpaired_prior))
  Other_gamma_chain[, b] <- 
    mvrnorm(n = 1, mu = coefficients(Other_prior), Sigma = vcov(Other_prior))
  MCI_gamma_chain[, b] <- 
    mvrnorm(n = 1, mu = coefficients(MCI_prior), Sigma = vcov(MCI_prior))
  
  #---- ****group membership ----
  group = 1
  synthetic_sample[, "Group"] <- 0
  for(model in c("Unimpaired", "Other", "MCI")){
    subset_index <- which(synthetic_sample$Group == 0)
    
    synthetic_sample[subset_index, paste0("p_", model)] <- 
      expit(as.matrix(synthetic_sample[subset_index, 
                                       get(paste0(model, "_preds"))]) %*% 
              as.matrix(get(paste0(model, "_gamma_chain"))[, b]))
    
    synthetic_sample[subset_index, "Group"] <- 
      rbernoulli(n = length(subset_index), 
                 p = synthetic_sample[subset_index, paste0("p_", model)])*group
    
    group = group + 1
  }
  
  synthetic_sample[which(synthetic_sample$Group == 0), "Group"] <- 4
  
  #---- ****group: summary ----
  latent_class_chain[, b] <- 
    table(synthetic_sample$Group)/sum(table(synthetic_sample$Group))
  
  for(i in 1:4){
    subset <- synthetic_sample %>% filter(Group == i) 
    posterior_counts <- alpha_0 + 
      table(subset$ETHNIC_label, subset$Astroke) %>% as.data.frame() %>% 
      dplyr::select("Freq") %>% unlist()
    
    #---- ****p(contingency table cell) ----
    pi_chain[, paste0(i, ":", b)] <- rdirichlet(1, alpha = posterior_counts)
    
    #---- ****contingency table count ----
    contingency_table <- rmultinom(n = 1, size = nrow(subset), 
                                   prob = pi_chain[, paste0(i, ":", b)])
    
    #---- ****make U matrix ----
    U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
    
    for(j in 1:nrow(contingency_table)){
      if(j == 1){
        index = 1
      } else{
        index = sum(contingency_table[1:(j - 1), ]) + 1
      }
      U[index:(index - 1 + contingency_table[j, ]), j] <- 1
    }
    
    UtU <- diag(contingency_table[, 1])
    
    #---- ****beta hat ----
    continuous_covariates <- subset %>% 
      dplyr::select(all_of(Z)) %>% as.matrix
    
    V <- solve(t(A) %*% UtU %*% A)
    
    beta_hat <-  V %*% t(A) %*% t(U) %*% continuous_covariates
    
    #---- ****epsilon hat ----
    eps_hat <- continuous_covariates - U %*% A %*% beta_hat
    
    #---- ****draw Sigma | Y ----
    sig_Y <- riwish(v = nrow(subset) - nrow(beta_hat), 
                    S = solve(t(eps_hat) %*% eps_hat))*Sigma_multiplier[i]
    
    Sigma_chain[, paste0(i, ":", b)] <- diag(sig_Y)
    
    #---- ****draw beta | Sigma, Y ----
    Sigma_kron_V <- kronecker(sig_Y, V, FUN = "*", make.dimnames = TRUE)
    beta_Sigma_Y <- 
      mvrnorm(n = 1, mu = as.vector(beta_hat), Sigma = Sigma_kron_V)
    
    #---- ****compute mu ----
    mu_chain[, paste0(i, ":", seq(1, nrow(cross_class_label)), ":", b)] <- 
      t(A %*% matrix(beta_Sigma_Y, nrow = ncol(A), ncol = length(Z), 
                     byrow = FALSE))
    
    #---- ****draw data ----
    #reformat contingency table
    contingency_table %<>% cbind(do.call(cbind, list(
      #Black              #Hispanic           #Stroke
      rep(c(1, 0, 0), 2), rep(c(0, 1, 0), 2), c(rep(0, 3), rep(1, 3))))) %>% 
      set_colnames(c("Count", W))
    
    for(j in 1:nrow(contingency_table)){
      if(contingency_table[j, "Count"] == 0){next}
      if(j == 1){
        index = 1
      } else{
        index = sum(contingency_table[1:(j - 1), "Count"]) + 1
      }
      #Z (continuous data)
      subset[index:(index - 1 + contingency_table[j, "Count"]), 
             colnames(sig_Y)] <- 
        ifelse(index == (index - 1 + contingency_table[j, "Count"]), 
               t(as.matrix(mvrnorm(n = contingency_table[j, "Count"],
                                   mu = mu_chain[, paste0(i, ":", j, ":", b)], 
                                   Sigma = sig_Y))), 
               mvrnorm(n = contingency_table[j, "Count"],
                       mu = mu_chain[, paste0(i, ":", j, ":", b)], 
                       Sigma = sig_Y))
      
      #W (categorical data)
      subset[index:(index - 1 + contingency_table[j, "Count"]), 
             colnames(contingency_table)[-1]] <- 
        matrix(rep(contingency_table[j, colnames(contingency_table)[-1]], 
                   contingency_table[j, "Count"]), 
               ncol = 3, byrow = TRUE)
    }
    
    #---- ****replace synthetic data ----
    synthetic_sample[which(synthetic_sample$HHIDPN %in% subset$HHIDPN), 
                     c(W, Z)] <- subset[, c(W, Z)]
  }
}

#---- **post-processing ----
#Update race/ethnicity label to match new synthetic data
synthetic_sample %<>% 
  mutate("ETHNIC_label" = case_when(Black == 1 ~ "Black", 
                                    Hispanic == 1 ~ "Hispanic", 
                                    TRUE ~ "White"))

#---- **plots ----
extended_pallette14 <- colorRampPalette(wes_palette("Darjeeling1"))(14)
extended_pallette6 <- colorRampPalette(wes_palette("Darjeeling1"))(6)

#---- ****gamma chains ----
gamma_plot_data <- 
  do.call(cbind, list(t(Unimpaired_gamma_chain), t(Other_gamma_chain), 
                      t(MCI_gamma_chain))) %>% as.data.frame() %>%
  mutate("run" = seq(1:B)) %>% 
  pivot_longer(-c("run"), names_to = c("Group", "Predictor"), 
               names_sep = ":", values_to = "gamma")

gamma_chain_plot <- 
  ggplot(data = gamma_plot_data, 
         aes(x = reorder(run, sort(as.numeric(run))), y = gamma, 
             colour = Predictor)) +       
  geom_line(aes(group = Predictor)) + 
  facet_grid(rows = vars(factor(Group, 
                                levels = c("Unimpaired", "MCI", "Other")))) + 
  theme_bw() + xlab("Run") + 
  scale_color_manual(values = rev(extended_pallette14))

ggsave(filename = "gamma_chain.jpeg", plot = gamma_chain_plot, 
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
  ggplot(data = latent_class_data, aes(x = run, y = prob, colour = Group)) +       
  geom_line(aes(group = Group)) + theme_minimal() + xlab("Run") + 
  ylab("Proportion of Sample") +  
  scale_color_manual(values = rev(wes_palette("Darjeeling1"))) + 
  scale_x_continuous(breaks = seq(1, B))

ggsave(filename = "latent_class_chain.jpeg", plot = latent_class_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 3, units = "in", device = "jpeg")

#---- ****pi chain ----
pi_chain_data <- pi_chain %>% as.data.frame() %>% rownames_to_column("Cell") %>% 
  pivot_longer(-c("Cell"), names_to = c("Group", "Run"), names_sep = ":", 
               values_to = "probability") %>% 
  mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                   Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                   Group == 4 ~ "Dementia")) %>% 
  mutate_at("Run", as.numeric) %>%
  mutate_if(is.character, as.factor) 

pi_chain_plot <- ggplot(data = pi_chain_data, 
                        aes(x = Run, y = probability, colour = Cell)) +       
  geom_line(aes(group = Cell)) + 
  theme_minimal() + xlab("Run") + ylab("Probability of cell membership") +  
  scale_color_manual(values = rev(extended_pallette6)) + 
  scale_x_continuous(breaks = seq(1, B)) + 
  facet_grid(rows = vars(factor(Group_label, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other")))) + theme_bw() 

ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 5, units = "in", device = "jpeg")

#---- ****Sigma chain ----
Sigma_chain_data <- Sigma_chain %>% as.data.frame() %>% 
  rownames_to_column("Z") %>% 
  pivot_longer(-c("Z"), names_to = c("Group", "Run"), names_sep = ":", 
               values_to = "variance") %>% 
  mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                   Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                   Group == 4 ~ "Dementia")) %>%
  mutate_at("Run", as.numeric) %>%
  mutate_if(is.character, as.factor) 

Sigma_chain_plot <- ggplot(data = Sigma_chain_data, 
                           aes(x = Run, y = variance, colour = Z)) +       
  geom_line(aes(group = Z)) + 
  theme_minimal() + xlab("Run") + ylab("Variance") +  
  scale_color_manual(values = rev(extended_pallette14)) + 
  scale_x_continuous(breaks = seq(1, B)) + 
  facet_grid(rows = vars(factor(Group_label, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other")))) + theme_bw() 

ggsave(filename = "Sigma_chain.jpeg", plot = Sigma_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 5, units = "in", device = "jpeg")

#---- ****mu chain ----
mu_chain_data <- mu_chain %>% as.data.frame() %>% 
  rownames_to_column("Z") %>% 
  pivot_longer(-c("Z"), names_to = c("Group", "Cell", "Run"), names_sep = ":", 
               values_to = "mu") %>% 
  mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
                                   Group == 2 ~ "Other", Group == 3 ~ "MCI", 
                                   Group == 4 ~ "Dementia")) %>% 
  mutate_at("Run", as.numeric) %>%
  mutate_if(is.character, as.factor) 

mu_chain_plot <- ggplot(data = mu_chain_data, 
                        aes(x = Run, y = mu, colour = Z)) +       
  geom_line(aes(group = Z)) + 
  theme_minimal() + xlab("Run") + ylab("mu") +  
  scale_color_manual(values = rev(extended_pallette14)) +
  scale_x_continuous(breaks = seq(1, B)) + 
  facet_grid(rows = vars(factor(Group_label, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other")))) + theme_bw() 

ggsave(filename = "mu_chain.jpeg", plot = Sigma_chain_plot, 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
       width = 14, height = 5, units = "in", device = "jpeg")



