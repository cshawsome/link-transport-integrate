#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "wesanderson", "devtools")
install_github("thomasp85/patchwork")

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character()))

#---- **RAND ----
RAND_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/RAND_subset.csv"), 
                        col_types = cols(HHIDPN = col_character()))

#---- join data ----
all_data <- left_join(ADAMS_subset, RAND_subset, by = "HHIDPN")

#---- select variables ----
vars <- c("GENDER", "ETHNIC", "log_AAGE", "EDYRS", 
          "ANMSETOT", paste0("r", seq(5, 7), "iadla_cat"))

analytical_sample <- all_data %>% 
  dplyr::select("HHIDPN", "AYEAR", 
                all_of(vars)) %>% na.omit()

#Variable check-- there's 538 people in the complete data set
colSums(is.na(analytical_sample))
dim(analytical_sample)

#---- **IADLA ----
#Take the IADLA measure closest to ADAMS interview year
analytical_sample %<>% 
  mutate("IADLA" = case_when(AYEAR == 2001 ~ r5iadla_cat, 
                             AYEAR %in% c(2002, 2003) ~ r6iadla_cat, 
                             AYEAR == 2004 ~ r7iadla_cat))
# #Sanity check
# table(analytical_sample$IADLA, useNA = "ifany")

#Get rid of original variables
analytical_sample %<>% dplyr::select(-c(contains("iadla_cat"), "AYEAR")) 

#---- all-way contingency table ----
cross_class <- table(analytical_sample$GENDER, analytical_sample$ETHNIC, 
                     analytical_sample$IADLA) %>% as.data.frame()

#---- plots ----
#Categorical Variables: Sex/Gender, Race/Ethnicity, IADLs
#Continuous Variables: Age, Ed Yrs, MMSE
#Create labeled data
analytical_sample %<>% 
  mutate("GENDER_label" = ifelse(GENDER == 1, "Male", "Female"), 
         "ETHNIC_label" = case_when(ETHNIC == 1 ~ "White", 
                                    ETHNIC == 2 ~ "Black", 
                                    TRUE ~ "Hispanic"), 
         "IADLA_label" = case_when(IADLA == 1 ~ "None", 
                                   IADLA == 2 ~ "One", 
                                   IADLA == 3 ~ "Two", 
                                   IADLA == 4 ~ "Three")) %>% 
  mutate_at(c("GENDER_label", "ETHNIC_label", "IADLA_label"), as.factor) 
analytical_sample$IADLA_label <- fct_relevel(analytical_sample$IADLA_label, 
                                             c("None", "One", "Two", "Three"))
  

#---- *marginal ----
#---- **summary stats ----
sex_gender_plot <- table(analytical_sample$GENDER) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate("labels" = ifelse(Var1 == 1, "Male", "Female"))

race_eth_plot <- table(analytical_sample$ETHNIC) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate("labels" = case_when(Var1 == 1 ~ "White", 
                              Var1 == 2 ~ "Black", 
                              TRUE ~ "Hispanic"))
  
IADL_plot <- table(analytical_sample$IADLA) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate("labels" = case_when(Var1 == 1 ~ "None", 
                              Var1 == 2 ~ "One", 
                              Var1 == 3 ~ "Two", 
                              TRUE ~ "Three"))

Age_plot <- table(exp(analytical_sample$log_AAGE)) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate_at("Var1", as.factor)

Edyrs_plot <- table(analytical_sample$EDYRS) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate_at("Var1", as.factor)

MMSE_plot <- table(analytical_sample$ANMSETOT) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) %>% 
  mutate_at("Var1", as.factor)
  
#---- **categorical plots ----
sex_gender <- 
  ggplot(data = sex_gender_plot) + 
  geom_bar(mapping = aes(x = factor(labels), y = Prop, fill = factor(labels)), 
           stat = "identity") + 
  theme_minimal() + xlab("Sex/Gender") + ylab("Proportion") + 
  theme(legend.position = "none") + ylim(c(0, 1)) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

race_eth <- 
  ggplot(data = race_eth_plot) + 
  geom_bar(mapping = aes(x = factor(labels), y = Prop, fill = factor(labels)), 
           stat = "identity") + ylim(c(0, 1)) +
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

IADLs <- 
  ggplot(data = IADL_plot) + 
  geom_bar(mapping = aes(x = factor(labels), y = Prop, fill = factor(labels)), 
           stat = "identity") + ylim(c(0, 1)) +
  theme_minimal() + xlab("Difficulty with IADLs") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **continuous plots ----
Age <- ggplot(data = Age_plot) + 
  geom_bar(mapping = aes(x = Var1, y = Prop), 
           color = rev(wes_palette("Darjeeling1"))[1],
           fill = rev(wes_palette("Darjeeling1"))[1],
           stat = "identity") + 
  theme_minimal() + xlab("Age") + ylab("Proportion") + 
  theme(legend.position = "none")

Edyrs <- ggplot(data = Edyrs_plot) + 
  geom_bar(mapping = aes(x = Var1, y = Prop), 
           color = rev(wes_palette("Darjeeling1"))[2],
           fill = rev(wes_palette("Darjeeling1"))[2],
           stat = "identity") + 
  theme_minimal() + xlab("Years of Education") + ylab("Proportion") + 
  theme(legend.position = "none")

MMSE <- ggplot(data = MMSE_plot) + 
  geom_bar(mapping = aes(x = Var1, y = Prop), 
           color = rev(wes_palette("Darjeeling1"))[4],
           fill = rev(wes_palette("Darjeeling1"))[4],
           stat = "identity") + 
  theme_minimal() + xlab("MMSE") + ylab("Proportion") + 
  theme(legend.position = "none")

#---- **patchwork plot ----
(((sex_gender + race_eth + IADLs)/Age)/Edyrs)/MMSE

ggsave(filename = "marginal_dists.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
       width = 8, height = 8, units = "in", device = "jpeg")

#---- *2-way categorical ----
#Sex/Gender by Race/Ethnicity
sex_by_race <- 
  ggplot(data = analytical_sample) + 
  geom_bar(mapping = aes(x = factor(GENDER_label), fill = factor(ETHNIC_label)), 
           position = "dodge") + 
  theme_minimal() + xlab("Sex/Gender") + ylab("Count") +
  guides(fill = guide_legend(title = "Race/Ethnicity")) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

sex_by_IADL <- 
  ggplot(data = analytical_sample) + 
  geom_bar(mapping = aes(x = factor(GENDER_label), fill = factor(IADLA_label)), 
           position = "dodge") + 
  theme_minimal() + xlab("Sex/Gender") + ylab("Count") +
  guides(fill = guide_legend(title = "IADL")) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

race_by_IADL <- ggplot(data = analytical_sample) + 
  geom_bar(mapping = aes(x = factor(ETHNIC_label), fill = factor(IADLA_label)), 
           position = "dodge") + 
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Count") +
  guides(fill = guide_legend(title = "IADL")) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#---- **patchwork plot ----
sex_by_race + sex_by_IADL + race_by_IADL

ggsave(filename = "2way_cat_dists.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
       width = 12, height = 3, units = "in", device = "jpeg")

#---- Bayes Stuff ----
#---- **number of runs ----
B = 2

#---- **priors ----
alpha_chain <- matrix(nrow = nrow(cross_class), ncol = B)
alpha_chain[, 1] <- rep(1, nrow(cross_class))

#---- **initiate values ----
pi_chain <- matrix(nrow = nrow(cross_class), ncol = B)
pi_chain[, 1] <- rep(1/nrow(cross_class), nrow(cross_class))

#---- **sampling ----


