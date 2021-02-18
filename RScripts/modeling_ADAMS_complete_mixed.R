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
cross_class_label <- table(analytical_sample$GENDER_label, 
                           analytical_sample$ETHNIC_label, 
                           analytical_sample$IADLA_label) %>% as.data.frame()

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
sex_gender_plot <- table(analytical_sample$GENDER_label) %>% 
  as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) 

race_eth_plot <- table(analytical_sample$ETHNIC_label) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) 

IADL_plot <- table(analytical_sample$IADLA_label) %>% as.data.frame() %>% 
  mutate("Prop" = Freq/nrow(analytical_sample)) 

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
  geom_bar(mapping = aes(x = factor(Var1), y = Prop, fill = factor(Var1)), 
           stat = "identity") + 
  theme_minimal() + xlab("Sex/Gender") + ylab("Proportion") + 
  theme(legend.position = "none") + ylim(c(0, 1)) +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

race_eth <- 
  ggplot(data = race_eth_plot) + 
  geom_bar(mapping = aes(x = factor(Var1), y = Prop, fill = factor(Var1)), 
           stat = "identity") + ylim(c(0, 1)) +
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

IADLs <- 
  ggplot(data = IADL_plot) + 
  geom_bar(mapping = aes(x = factor(Var1), y = Prop, fill = factor(Var1)), 
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

#---- *3-way categorical ----
IADL_by_sex_and_race <- 
  ggplot(data = analytical_sample, aes(x = GENDER_label, y = ETHNIC_label, 
                                       color = IADLA_label)) + 
  geom_point(position = "jitter", alpha = 0.75) + 
  theme_minimal() + ylab("Race/Ethnicity") + xlab("Sex/Gender") + 
  guides(color = guide_legend(title = "IADL")) + 
  scale_color_manual(values = rev(wes_palette("Darjeeling1")))

ggsave(filename = "3way_cat_dists.jpeg", plot = last_plot(), 
       path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
       width = 8, height = 8, units = "in", device = "jpeg")

#---- *cont | categorical ----
for(i in 1:nrow(cross_class_label)){
  if(cross_class_label[i, "Freq"] != 0){
    gender <- cross_class_label[i, "Var1"]
    ethnic <- cross_class_label[i, "Var2"]
    iadl <- cross_class_label[i, "Var3"]
    count <- cross_class_label[i, "Freq"]
    
    data_subset <- analytical_sample %>% 
      filter(GENDER_label == gender & ETHNIC_label == ethnic & 
               IADLA_label == iadl)
    
    #summary stats
    age_plot <- table(exp(data_subset$log_AAGE)) %>% as.data.frame() %>% 
      mutate("Prop" = Freq/nrow(data_subset)) %>% 
      mutate_at("Var1", as.character) %>% 
      mutate_at("Var1", as.numeric)
    missing_ages <- which(!seq(min(exp(analytical_sample$log_AAGE)), 
                               max(exp(analytical_sample$log_AAGE))) %in% 
                            age_plot$Var1) + 69
    age_plot %<>% rbind(as.matrix(cbind(missing_ages, NA, NA)) %>% 
                          set_colnames(c("Var1", "Freq", "Prop"))) %>% 
      mutate_at("Var1", as.factor)
    
    edyrs_plot <- table(data_subset$EDYRS) %>% as.data.frame() %>% 
      mutate("Prop" = Freq/nrow(data_subset)) %>% 
      mutate_at("Var1", as.character) %>% 
      mutate_at("Var1", as.numeric)
    missing_edyrs <- which(!seq(min(analytical_sample$EDYRS), 
                                max(analytical_sample$EDYRS)) %in% 
                             edyrs_plot$Var1) - 1
    edyrs_plot %<>% rbind(as.matrix(cbind(missing_edyrs, NA, NA)) %>% 
                            set_colnames(c("Var1", "Freq", "Prop"))) %>% 
      mutate_at("Var1", as.factor)
    
    mmse_plot <- table(data_subset$ANMSETOT) %>% as.data.frame() %>% 
      mutate("Prop" = Freq/nrow(data_subset)) %>% 
      mutate_at("Var1", as.character) %>% 
      mutate_at("Var1", as.numeric)
    missing_mmse <- which(!seq(min(analytical_sample$ANMSETOT), 
                               max(analytical_sample$ANMSETOT)) %in% 
                            mmse_plot$Var1) - 1
    mmse_plot %<>% rbind(as.matrix(cbind(missing_mmse, NA, NA)) %>% 
                           set_colnames(c("Var1", "Freq", "Prop"))) %>% 
      mutate_at("Var1", as.factor)
    
    #plots
    age <- ggplot(data = age_plot) + 
      geom_bar(mapping = aes(x = Var1, y = Prop), 
               color = rev(wes_palette("Darjeeling1"))[1],
               fill = rev(wes_palette("Darjeeling1"))[1],
               stat = "identity") + 
      theme_minimal() + xlab("Age") + ylab("Proportion") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position = "none") + 
      ggtitle(paste(paste(gender, ethnic, iadl, sep = " | "), "| n =", count)) 
    
    edyrs <- ggplot(data = edyrs_plot) + 
      geom_bar(mapping = aes(x = Var1, y = Prop), 
               color = rev(wes_palette("Darjeeling1"))[2],
               fill = rev(wes_palette("Darjeeling1"))[2],
               stat = "identity") + 
      theme_minimal() + xlab("Years of Education") + ylab("Proportion") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position = "none") + 
      ggtitle(paste(paste(gender, ethnic, iadl, sep = " | "), "| n =", count)) 
    
    mmse <- ggplot(data = mmse_plot) + 
      geom_bar(mapping = aes(x = Var1, y = Prop), 
               color = rev(wes_palette("Darjeeling1"))[4],
               fill = rev(wes_palette("Darjeeling1"))[4],
               stat = "identity") + 
      theme_minimal() + xlab("MMSE") + ylab("Proportion") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      theme(legend.position = "none") + 
      ggtitle(paste(paste(gender, ethnic, iadl, sep = " | "), "| n =", count)) 
    
    #---- **patchwork plot ----
    age + edyrs + mmse
    
    ggsave(filename = paste0("cont_given_cat", i, ".jpeg"), plot = last_plot(), 
           path = "/Users/CrystalShaw/Box/Dissertation/figures/prelim_analyses/", 
           width = 12, height = 3, units = "in", device = "jpeg")
  }
}

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


