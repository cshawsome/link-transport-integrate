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
                                 TRUE ~ "Unknown")) 

#---- cross-class race/ethnicity x stroke ----
cross_class_data <- RAND %>% filter(!Ethnicity %in% c("Unknown", "Other"))
prior_counts <- 
  as.data.frame(table(cross_class_data$Ethnicity, cross_class_data$r5stroke))

#---- stratify priors by latent class ----
for(class in names(dem_class_props)){
  prior_counts %<>% 
    mutate("{{class}}_prior_count" := 
             Freq*dem_class_props[{{class}}])
}

#---- save output ----
prior_counts %>% write_csv(file = here("priors", "contingency_cell_counts.csv"))






