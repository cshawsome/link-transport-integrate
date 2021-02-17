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
                                   IADLA == 4 ~ "Three"))

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
  
#---- **categorical plots ----
sex_gender <- 
  ggplot(data = sex_gender_plot) + 
  geom_bar(mapping = aes(x = factor(labels), y = Prop, fill = factor(labels)), 
           stat = "identity") + 
  theme_minimal() + xlab("Sex/Gender") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

race_eth <- 
  ggplot(data = race_eth_plot) + 
  geom_bar(mapping = aes(x = factor(labels), y = Prop, fill = factor(labels)), 
           stat = "identity") + 
  theme_minimal() + xlab("Race/Ethnicity") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

IADLs <- 
  ggplot(data = IADL_plot) + 
  geom_bar(mapping = aes(x = factor(labels), y = Prop, fill = factor(labels)), 
           stat = "identity") + 
  theme_minimal() + xlab("Difficulty with IADLs") + ylab("Proportion") + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rev(wes_palette("Darjeeling1")))

#Make a patchwork plot
sex_gender + race_eth + IADLs


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


