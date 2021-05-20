#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here", "MASS", "MCMCpack", 
       "locfit", "MBSP", "wesanderson", "RColorBrewer")

#---- read in data ----
#---- **ADAMS ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"), 
                        col_types = cols(HHIDPN = col_character())) 

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
kappa_0 <- 1

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
  mutate("Group" = 0, "p_Unimpaired" = 0, "p_Other" = 0, "p_MCI" = 0)

#---- all-way contingency table ----
#We have one small cell-- Hispanics who have had a stroke
cross_class_label <- table(synthetic_sample$ETHNIC_label, 
                           synthetic_sample$Astroke) %>% as.data.frame() %>% 
  mutate("Stroke" = ifelse(Var2 == 0, "No Stroke", "Stroke")) %>% 
  unite("Cell Label", c("Var1", "Stroke"), sep = " | ", remove = FALSE)

#---- arrange data ----
synthetic_sample <- arrange(synthetic_sample, 
                            Astroke, desc(Black), desc(Hispanic))

#---- Bayes Stuff ----
#---- **simulation runs ----
B = 10

#categorical vars contrasts matrix
A = do.call(cbind, list(
  #intercept
  rep(1, 6),
  #race/ethnicity main effect: Black
  rep(c(1, 0, 0), 2),
  #race/ethnicity main effect: Hispanic
  rep(c(0, 1, 0), 2),
  #stroke main effect
  rep(c(0, 1), each = 3)))

#---- **chain storage ----
model_gamma_chain <- 
  matrix(nrow = sum(length(normal_preds), length(other_preds), 
                    length(mci_preds)), ncol = B) %>% as.data.frame() %>%
  mutate("model" = c(rep("normal", length(normal_preds)), 
                     rep("other", length(other_preds)), 
                     rep("mci", length(mci_preds))), 
         "pred" = c(normal_preds, mci_preds, other_preds))
  
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

#---- START TIME ----
start <- Sys.time()
#---- **sampling ----
for(b in 1:B){
  #---- ****latent class gammas ----
  for(model in c("normal", "other", "mci")){
    random_draw <- sample(seq(1, 10000), size = 1)
    
    prior_betas <- as.vector(get(paste0(model, "_betas"))[, random_draw])
    prior_cov <- matrix(unlist(get(paste0(model, "_cov"))[, random_draw]), 
                        nrow = nrow(prior_betas))
    
    model_gamma_chain[which(model_gamma_chain$model == model), b] <- 
      mvrnorm(n = 1, mu = unlist(prior_betas), Sigma = prior_cov)
  }
  
  #---- ****group membership ----
  group = 1
  synthetic_sample[, "Group"] <- 0
  for(model in c("normal", "other", "mci")){
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
    posterior_counts <- alpha_0[, paste0(i, "_prior_count")] + 
      table(subset$ETHNIC_label, subset$Astroke) %>% as.data.frame() %>% 
      dplyr::select("Freq") %>% unlist()
    
    #---- ****p(contingency table cell) ----
    pi_chain[, paste0(i, ":", b)] <- 
      rdirichlet(1, alpha = as.numeric(unlist(posterior_counts)))
    
    #---- ****contingency table count ----
    contingency_table <- rmultinom(n = 1, size = nrow(subset), 
                                   prob = pi_chain[, paste0(i, ":", b)])
    
    #---- ****make U matrix ----
    U <- matrix(0, nrow = nrow(subset), ncol = nrow(contingency_table))
    
    for(j in 1:nrow(contingency_table)){
      if(contingency_table[j, ] == 0){next}
      if(j == 1){
        index = 1
      } else{
        index = sum(contingency_table[1:(j - 1), ]) + 1
      }
      U[index:(index - 1 + contingency_table[j, ]), j] <- 1
    }
    
    UtU <- diag(contingency_table[, 1])
    
    if(i %in% c(2, 3)){
      assign(paste0("UtU_", i), UtU)
    }
    
    #---- ****pool UtU if needed ----
    if(det(t(A) %*% UtU %*% A) == 0){
      if(exists(paste0("UtU_", (i-1)))){
        UtU <- UtU + get(paste0("UtU_", (i-1)))
      } else{
        UtU <- UtU + get(paste0("UtU_", (i+1))) 
      }
    }
    
    #---- ****Mm ----
    continuous_covariates <- subset %>% 
      dplyr::select(all_of(Z)) %>% as.matrix
    
    V_inv <- t(A) %*% UtU %*% A 
    V_0_inv <- matrix(V_inv_prior[, i], nrow = nrow(V_inv), ncol = ncol(V_inv))
    beta_0 <- matrix(beta_prior[, i], nrow = nrow(V_inv), 
                     ncol = ncol(continuous_covariates))
    
    M <- solve(V_inv + kappa_0[i]*V_0_inv)
    m <-  t(A) %*% t(U) %*% continuous_covariates - 
      kappa_0[i]*V_0_inv %*% beta_0
    
    Mm <- M %*% m
    
    #---- ****draw Sigma | Y ----
    ZtZ <- t(continuous_covariates) %*% continuous_covariates
    third_term <- kappa_0[i]*t(beta_0) %*% V_0_inv %*% beta_0
    
    sig_Y <- riwish(v = (nu_0 + nrow(subset)), 
                    S = Sigma_prior + ZtZ + third_term)
    
    Sigma_chain[, paste0(i, ":", b)] <- diag(sig_Y)
    
    #---- ****draw beta | Sigma, Y ----
    beta_Sigma_Y <- matrix.normal(Mm, M, sig_Y)
    
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
      if(contingency_table[j, "Count"] == 1){
        subset[index:(index - 1 + contingency_table[j, "Count"]), 
               colnames(sig_Y)] <- 
          t(as.matrix(mvrnorm(n = contingency_table[j, "Count"],
                              mu = mu_chain[, paste0(i, ":", j, ":", b)], 
                              Sigma = diag(sqrt(var_scale[, i])) %*% sig_Y %*% 
                                diag(sqrt(var_scale[, i])))))
      } else{
        subset[index:(index - 1 + contingency_table[j, "Count"]), 
               colnames(sig_Y)] <- 
          mvrnorm(n = contingency_table[j, "Count"],
                  mu = mu_chain[, paste0(i, ":", j, ":", b)], 
                  Sigma = diag(sqrt(var_scale[, i])) %*% sig_Y %*% 
                    diag(sqrt(var_scale[, i])))
      }
      
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
  #---- **post-processing ----
  #---- ****race/ethnicity ----
  synthetic_sample %<>% 
    mutate("ETHNIC_label" = case_when(Black == 1 ~ "Black", 
                                      Hispanic == 1 ~ "Hispanic", 
                                      TRUE ~ "White"))
}

#---- END TIME ----
stop <- Sys.time() - start

#---- **plots ----
extended_pallette14 <- colorRampPalette(wes_palette("Darjeeling1"))(14)
extended_pallette6 <- colorRampPalette(wes_palette("Darjeeling1"))(6)

#---- ****gamma chains ----
gamma_plot_data <- 
  do.call(cbind, list(t(Unimpaired_gamma_chain), t(Other_gamma_chain), 
                      t(MCI_gamma_chain))) %>% as.data.frame() %>%
  mutate("run" = seq(1:B)) %>% 
  pivot_longer(-c("run"), names_to = c("Group", "Predictor"), 
               names_sep = ":", values_to = "gamma") %>% 
  filter(Predictor != "(Intercept)")

gamma_chain_plot <- 
  ggplot(data = gamma_plot_data, 
         aes(x = reorder(run, sort(as.numeric(run))), y = gamma, 
             colour = Predictor)) +       
  geom_line(aes(group = Predictor)) + 
  facet_grid(rows = vars(factor(Group, 
                                levels = c("Unimpaired", "MCI", "Other")))) + 
  theme_bw() + xlab("Run") + 
  scale_x_discrete(breaks = seq(0, B, by = 100)) + 
  scale_color_manual(values = rev(extended_pallette14))

ggsave(filename = "gamma_chain.jpeg", plot = gamma_chain_plot, 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
                     "standard_normal"), 
       width = 7, height = 6, units = "in", device = "jpeg")

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
  scale_x_continuous(breaks = seq(0, B, by = 100)) 

ggsave(filename = "latent_class_chain.jpeg", plot = latent_class_chain_plot, 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
                     "standard_normal"), 
       width = 7, height = 3, units = "in", device = "jpeg")

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
  scale_color_manual(values = extended_pallette6) + 
  scale_x_continuous(breaks = seq(0, B, by = 100)) +
  facet_grid(rows = vars(factor(Group_label, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other")))) + theme_bw() 

ggsave(filename = "pi_chain.jpeg", plot = pi_chain_plot, 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
                     "standard_normal"), 
       width = 7, height = 5, units = "in", device = "jpeg")

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
  scale_x_continuous(breaks = seq(0, B, by = 100)) + 
  facet_grid(rows = vars(factor(Group_label, 
                                levels = c("Unimpaired", "MCI", "Dementia", 
                                           "Other")))) + theme_bw() 

ggsave(filename = "Sigma_chain.jpeg", plot = Sigma_chain_plot, 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
                     "standard_normal"), 
       width = 7, height = 5, units = "in", device = "jpeg")

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
                        aes(x = Run, y = mu, colour = Group_label)) +       
  geom_line(aes(group = Group_label), alpha = 0.75) + 
  theme_minimal() + xlab("Run") + ylab("mu") +  
  scale_color_manual(values = c(wes_palette("Darjeeling1")[1], 
                                wes_palette("Darjeeling1")[3], 
                                wes_palette("Darjeeling1")[5], 
                                wes_palette("Darjeeling1")[2])) +
  scale_x_continuous(breaks = seq(0, B, by = 100)) + 
  facet_grid(rows = vars(factor(Z))) + theme_bw() 

ggsave(filename = "mu_chain.jpeg", plot = mu_chain_plot, 
       path = paste0("/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
                     "standard_normal"), 
       width = 7, height = 14, units = "in", device = "jpeg")

# #---- ****varY chain ----
# varY_chain_data <- varY_chain %>% as.data.frame() %>% 
#   rownames_to_column("Z") %>% 
#   pivot_longer(-c("Z"), names_to = c("Group", "Run"), names_sep = ":", 
#                values_to = "variance") %>% 
#   mutate("Group_label" = case_when(Group == 1 ~ "Unimpaired", 
#                                    Group == 2 ~ "Other", Group == 3 ~ "MCI", 
#                                    Group == 4 ~ "Dementia")) %>%
#   mutate_at("Run", as.numeric) %>%
#   mutate_if(is.character, as.factor) 
# 
# varY_chain_plot <- ggplot(data = varY_chain_data, 
#                           aes(x = Run, y = variance, colour = Z)) +       
#   geom_line(aes(group = Z)) + 
#   theme_minimal() + xlab("Run") + ylab("Variance") +  
#   scale_color_manual(values = rev(extended_pallette14)) + 
#   scale_x_continuous(breaks = seq(0, B, by = 100)) + 
#   facet_grid(rows = vars(factor(Group_label, 
#                                 levels = c("Unimpaired", "MCI", "Dementia", 
#                                            "Other")))) + theme_bw() 
# 
# ggsave(filename = "varY_chain.jpeg", plot = varY_chain_plot, 
#        path = "/Users/CrystalShaw/Box/Dissertation/figures/diagnostics/", 
#        width = 7, height = 5, units = "in", device = "jpeg")

#---- save datasets ----
write_csv(synthetic_sample, 
          file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/results/", 
                        "ADAMSA/standard_normal/ADAMSA_synthetic.csv"))

write_csv(gamma_plot_data, 
          file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/results/", 
                        "ADAMSA/standard_normal/diagnostics_data/", 
                        "ADAMSA_gamma_plot_data.csv"))

write_csv(latent_class_data, 
          file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/results/", 
                        "ADAMSA/standard_normal/diagnostics_data/", 
                        "ADAMSA_latent_class_data.csv"))

write_csv(pi_chain_data, 
          file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/results/", 
                        "ADAMSA/standard_normal/diagnostics_data/", 
                        "ADAMSA_pi_chain_data.csv"))

write_csv(Sigma_chain_data, 
          file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/results/", 
                        "ADAMSA/standard_normal/diagnostics_data/", 
                        "ADAMSA_Sigma_chain_data.csv"))

write_csv(mu_chain_data, 
          file = paste0("/Users/CrystalShaw/Box/Dissertation/analyses/results/", 
                        "ADAMSA/standard_normal/diagnostics_data/", 
                        "ADAMSA_mu_chain_data.csv"))






