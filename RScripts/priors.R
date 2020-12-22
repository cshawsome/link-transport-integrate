#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

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

#---- **table: ADAMS dem dx x wave ---
table_data <- ADAMS_subset %>% dplyr::select(contains("dem_dx_cat")) %>% 
  mutate_all(function(x) case_when(x %in% 
               c("Dementia", "Probable Dementia", 
                 "Probable/Possible AD", 
                 "Probable/Possible Vascular Dementia") ~ "Dementia", 
             TRUE ~ x)) 

for(wave in c("A", "B", "C", "D")){
  dem_dx_var <- paste0(wave, "dem_dx_cat")
  print(paste0("Wave ", wave))
  print(table(table_data[, dem_dx_var], useNA = "ifany")/
          sum(!is.na(table_data[, dem_dx_var])))
}

# #Sanity check
# View(table_data)

