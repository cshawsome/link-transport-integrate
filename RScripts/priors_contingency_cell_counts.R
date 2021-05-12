#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "haven", "labelled")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/cleaned/ADAMS_subset_mixed.csv")) 

#---- **RAND ----
# Year | HRS | RAND
# 2000 | G | 5
rand_waves <- 5
rand_variables <- c("hhidpn", "raracem", "rahispan", 
                    paste0("r", rand_waves, "agey_e"), 
                    paste0("r", rand_waves, "stroke"))

RAND <- read_dta(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                        "RAND_longitudinal/STATA/randhrs1992_2016v2.dta"), 
                 col_select = all_of(rand_variables)) %>% 
  mutate_at("hhidpn", as.character) %>% 
  filter(r5agey_e >= 70)

colnames(RAND)[1] <- "HHIDPN" #For merging

#Remove labeled data format
val_labels(RAND) <- NULL

#---- data cleaning: dem dx ----
ADAMS_subset %<>% 
  mutate("Adem_dx_cat_collapse" = 
           case_when(Adem_dx_cat %in% 
                       c("Dementia", "Probable/Possible AD", "Probable Dementia", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia", 
                     TRUE ~ Adem_dx_cat))

# #Sanity check
# table(ADAMS_subset$Adem_dx_cat, ADAMS_subset$Adem_dx_cat_collapse, 
#       useNA = "ifany")

dem_class_props <- table(ADAMS_subset$Adem_dx_cat_collapse)/nrow(ADAMS_subset) 

#---- data cleaning: race-ethnicity ----
RAND %<>% 
  mutate("hispanic" = ifelse(rahispan == 0 | is.na(rahispan), 0, 1)) %>% 
  mutate("Ethnicity" = case_when(hispanic == 1 ~ "Hispanic",
                                 (raracem == 1 & hispanic == 0) ~ "White", 
                                 (raracem == 2 & hispanic == 0) ~ "Black", 
                                 (raracem == 3 & hispanic == 0) ~ "Other", 
                                 TRUE ~ "Unknown")) %>% 
  filter(!Ethnicity %in% c("Unknown", "Other") & !is.na(r5stroke))

#---- cross-class race/ethnicity x stroke ----
data_counts <- as.data.frame(table(ADAMS_subset$ETHNIC_label, 
                                   ADAMS_subset$r5stroke)) %>% 
  mutate("percent" = round((Freq/sum(Freq))*100, 1))

#---- bootstrap counts ----
B = 10000
bootstrap_counts <- matrix(nrow = 6, ncol = B)
for(b in 1:B){
  sample <- sample_n(ADAMS_subset, size = nrow(ADAMS_subset), replace = TRUE) 
  #sub_sample <- sample_frac(sample, size = 0.3, replace = FALSE)
  bootstrap_counts[, b] <- 
    as.data.frame(table(sample$ETHNIC_label, sample$r5stroke))$Freq
}

bootstrap_percents <- 
  round(bootstrap_counts/colSums(bootstrap_counts)*100, 1) %>% 
  as.data.frame() %>%
  mutate("truth" = data_counts$percent, 
         "cat" = c("Black + No Stroke", "Hispanic + No Stroke", 
                   "White + No Stroke", "Black + Stroke", "Hispanic + Stroke", 
                   "White + Stroke")) %>% pivot_longer(-c("truth", "cat"))

#---- **plots ----
for(category in unique(bootstrap_percents$cat)){
  subset <- bootstrap_percents %>% filter(cat == category)
  ggplot(data = subset , aes(x = value)) + 
    geom_histogram(fill = "black", color = "black") + theme_minimal() + 
    xlab("Percent") + ggtitle(category) + 
    geom_vline(xintercept = subset$truth, color = "#f2caaa", size = 2)
  
  ggsave(filename = paste0("/Users/CrystalShaw/Box/Dissertation/figures/", 
                           "priors/cell_counts/RAND/", category, ".jpeg"), 
         width = 5, height = 3, units = "in")
}

bootstrap_counts %>% as.data.frame() %>%
  write_csv(file = here::here("priors", "bootstrap_cell_counts.csv"))






