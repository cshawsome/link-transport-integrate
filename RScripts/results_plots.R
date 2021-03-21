#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer")

#---- read in data ----
synthetic_ADAMS <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "analyses/results/ADAMSA/ADAMSA_synthetic.csv"), 
           col_types = cols(HHIDPN = col_character()))

ADAMS_columns <- c(colnames(synthetic_ADAMS)[
  which(!colnames(synthetic_ADAMS) %in% 
          c("(Intercept)", "Group", "Black", "Hispanic", "p_Unimpaired", 
            "p_Other", "p_MCI"))], "Adem_dx_cat")

ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character())) %>% 
  #only take those part of modeling (complete cases only for now)
  filter(HHIDPN %in% synthetic_ADAMS$HHIDPN) %>% 
  dplyr::select(all_of(ADAMS_columns))

#---- merge datasets ----
#format column names
synthetic_ADAMS %<>% 
  set_colnames(c("HHIDPN", paste0("synthetic:", colnames(synthetic_ADAMS))[-1]))

ADAMS_subset %<>% 
  set_colnames(c("HHIDPN", paste0("ADAMSA:", colnames(ADAMS_subset))[-1]))

merged_data <- left_join(synthetic_ADAMS, ADAMS_subset)

#---- data cleaning: dem group ----
#data check
table(merged_data$`ADAMSA:Adem_dx_cat`, useNA = "ifany")
table(merged_data$`synthetic:Group`, useNA = "ifany")

merged_data %<>% 
  mutate("ADAMSA:group_class" = 
           case_when(`ADAMSA:Adem_dx_cat` %in% 
                       c("Dementia", "Probable/Possible AD", 
                         "Probable/Possible Vascular Dementia") ~ "Dementia",
                     `ADAMSA:Adem_dx_cat` == "Normal" ~ "Unimpaired",
                     TRUE ~ `ADAMSA:Adem_dx_cat`), 
         "synthetic:group_class" = 
           case_when(`synthetic:Group` == 1 ~ "Unimpaired", 
                     `synthetic:Group` == 2 ~ "Other", 
                     `synthetic:Group` == 3 ~ "MCI", 
                     `synthetic:Group` == 4 ~ "Dementia"))

# #Sanity check
# table(merged_data$`ADAMSA:group_class`)
# table(merged_data$`synthetic:group_class`)

#---- plots: dementia classes ----
dementia_class_plot_data <- 
  merged_data %>% dplyr::select(contains("group_class")) %>% 
  pivot_longer(everything(), names_to = c("Data", "Var"), 
               names_pattern = "(.*):(.*)") %>% 
  count(Data, value) %>%
  group_by(Data) %>%
  mutate(n = n/sum(n))

#releveling factors
dementia_class_plot_data$value <- 
  fct_relevel(dementia_class_plot_data$value, 
              c("Unimpaired", "MCI", "Dementia", "Other"))

dementia_class_plot <- 
  ggplot(data = dementia_class_plot_data) + 
  geom_bar(mapping = aes(x = factor(value), y = n, fill = factor(Data)), 
           stat = "identity", position = "dodge") + 
  theme_minimal() + xlab("Impairment Class") + ylab("Proportion") + 
  ylim(c(0, 1)) + labs(fill = "Data") + 
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "dementia_class_overall.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/results/ADAMSA/", 
       width = 5, height = 5, units = "in", device = "jpeg")


