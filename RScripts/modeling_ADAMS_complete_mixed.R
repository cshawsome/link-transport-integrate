#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr")

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


