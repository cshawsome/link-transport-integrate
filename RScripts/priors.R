#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr")

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset.csv"))

#---- **plot: ADAMS MMSE x dem dx ----
plot_data <- ADAMS_subset %>% 
  dplyr::select(contains(c("NMSETOT", "dem_dx_cat"))) %>%
  pivot_longer(everything(), 
               names_to = c("wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  mutate("dem_dx_new_cat" = 
           case_when(dem_dx_cat %in% 
                       c("Dementia", "Probable Dementia", 
                         "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia", 
                     TRUE ~ dem_dx_cat))

ggplot(data = plot_data %>% filter(dem_dx_new_cat != "Other"), 
       aes(x = NMSETOT, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx")) + 
  facet_grid(cols = vars(wave))

ggsave(filename = "ADAMS_MMSE_by_dem.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "figures/"), width = 15, height = 5, units = "in")

#Including "Other" category
ggplot(data = plot_data, 
       aes(x = NMSETOT, color = factor(dem_dx_new_cat), 
           fill = factor(dem_dx_new_cat))) + 
  geom_histogram(alpha = 0.4, position = "identity") + theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx")) + 
  facet_grid(cols = vars(wave))

ggsave(filename = "ADAMS_MMSE_by_dem+other.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "figures/"), width = 15, height = 5, units = "in")

#---- **table: ADAMS dem dx x wave ----
table_data <- ADAMS_subset %>% dplyr::select(contains("dem_dx_cat")) %>% 
  mutate_all(function(x) case_when(x %in% 
               c("Dementia", "Probable Dementia", 
                 "Probable/Possible AD", 
                 "Probable/Possible Vascular Dementia") ~ "Dementia", 
             TRUE ~ x)) 

for(wave in c("A", "B", "C", "D")){
  dem_dx_var <- paste0(wave, "dem_dx_cat")
  print(paste0("Wave ", wave))
  print(table(table_data[, dem_dx_var], useNA = "ifany"))
  print(sum(table(table_data[, dem_dx_var])))
  print(table(table_data[, dem_dx_var])/
          sum(!is.na(table_data[, dem_dx_var])))
}

# #Sanity check
# View(table_data)

#---- prior: density overlay with ADAMS A ----
#---- **sample sizes ----
dem_dx_cat_counts <- table(ADAMS_subset$Adem_dx_cat, useNA = "ifany")
ADAMSA_dem <- sum(dem_dx_cat_counts[c("Dementia", "Probable Dementia", 
                                      "Probable/Possible AD", 
                                      "Probable/Possible Vascular Dementia")])
ADAMSA_MCI <- dem_dx_cat_counts["MCI"]
ADAMSA_normal <- dem_dx_cat_counts["Normal"]
#---- **generate fake data ----
set.seed(20201223)
fake_normal <- rnorm(n = ADAMSA_dem, mean = 28.5, sd = 3) %>% 
  as.data.frame() %>% 
  set_colnames("NMSETOT") %>%
  mutate("dem_dx_new_cat" = "Normal") 
  
fake_MCI <- rnorm(n = ADAMSA_MCI, mean = 24, sd = 4) %>% as.data.frame() %>% 
  set_colnames("NMSETOT") %>%
  mutate("dem_dx_new_cat" = "MCI") 

fake_dem <- rnorm(n = ADAMSA_normal, mean = 15, sd = 5) %>% 
  as.data.frame() %>% set_colnames("NMSETOT") %>%
  mutate("dem_dx_new_cat" = "Dementia") 

fake_data <- do.call("rbind", list(fake_normal, fake_MCI, fake_dem)) %>% 
  mutate_at("NMSETOT", function(x) case_when(x > 30 ~ 30, 
                                             x < 0 ~ 0, 
                                             TRUE ~ x))

#---- **plot ----
plot_data <- ADAMS_subset %>% 
  dplyr::select(contains(c("NMSETOT", "dem_dx_cat"))) %>%
  pivot_longer(everything(), 
               names_to = c("wave", ".value"),
               names_pattern = "(.)(.*)") %>% 
  mutate("dem_dx_new_cat" = 
           case_when(dem_dx_cat %in% 
                       c("Dementia", "Probable Dementia", 
                         "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia", 
                     TRUE ~ dem_dx_cat)) %>% 
  filter(wave == "A")

ggplot() + 
  geom_density(data = fake_data, 
               aes(x = NMSETOT, color = factor(dem_dx_new_cat))) + 
  geom_histogram(data = plot_data %>% filter(dem_dx_new_cat != "Other"), 
                 aes(x = NMSETOT, y = ..density.., 
                     color = factor(dem_dx_new_cat), 
                     fill = factor(dem_dx_new_cat)), alpha = 0.4, 
                 position = "identity")+
  theme_minimal() + 
  xlab("MMSE") + theme(text = element_text(size = 14)) + 
  guides(fill = guide_legend(title = "Dementia Dx"), 
         color = guide_legend(title = "Dementia Dx"))

ggsave(filename = "ADAMSA_MMSE_prior_overlay.jpeg", plot = last_plot(), 
       device = "jpeg",
       path = paste0("/Users/CrystalShaw/Box/Dissertation/",
                     "figures/"), width = 7, height = 5, units = "in")

#---- prior: Dirichlet mix ----
sim_alpha <- rgamma(n = 1000, shape = 0.25, rate = 0.25)
hist(sim_alpha)
summary(sim_alpha)

sim_beta_med <- rbeta(n = 1000, shape1 = 1, shape2 = median(sim_alpha))
hist(sim_beta_med)
summary(sim_beta_med)

sim_beta_min <- rbeta(n = 1000, shape1 = 1, shape2 = min(sim_alpha))
hist(sim_beta_min)
summary(sim_beta_min)

sim_beta_max <- rbeta(n = 1000, shape1 = 1, shape2 = max(sim_alpha))
hist(sim_beta_max)
summary(sim_beta_max)

