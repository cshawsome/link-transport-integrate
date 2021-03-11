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

analytical_sample <- ADAMS_subset %>% 
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
cross_class_label <- table(analytical_sample$ETHNIC_label, 
                           analytical_sample$Astroke) %>% as.data.frame()

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
  analytical_sample[, "Group"] <- 0
  for(model in c("Unimpaired", "Other", "MCI")){
    subset_index <- which(analytical_sample$Group == 0)
    
    analytical_sample[subset_index, paste0("p_", model)] <- 
      expit(as.matrix(analytical_sample[subset_index, 
                                        get(paste0(model, "_preds"))]) %*% 
              as.matrix(get(paste0(model, "_beta_chain"))[, b]))
    
    analytical_sample[subset_index, "Group"] <- 
      rbernoulli(n = length(subset_index), 
                 p = analytical_sample[subset_index, paste0("p_", model)])*group
    
    group = group + 1
  }
  
  analytical_sample[which(analytical_sample$Group == 0), "Group"] <- 4
  
  #---- ****group: summary ----
  latent_class_chain[, b] <- 
    table(analytical_sample$Group)/sum(table(analytical_sample$Group))
}

#---- **plots ----
extended_pallette14 <- colorRampPalette(wes_palette("Darjeeling1"))(14)

#---- **** beta chains ----
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




