#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg", "magrittr", "here")

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset_mixed.csv"), 
                         col_types = cols(HHIDPN = col_character()))

#---- **models for priors ----
normal_prior <- readRDS(here::here("priors", "normal_model_25.rds"))
other_prior <- readRDS(here::here("priors", "other_model_25.rds"))
MCI_prior <- readRDS(here::here("priors", "MCI_model_25.rds"))

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- c("AAGE", "ETHNIC_label", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", 
          "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Astroke", "Abmi")

analytical_sample <- ADAMS_subset %>% dplyr::select("HHIDPN", all_of(vars)) %>% 
  mutate("Black" = ifelse(ETHNIC_label == "Black", 1, 0), 
         "Hispanic" = ifelse(ETHNIC_label == "Hispanic", 1, 0))

# #Sanity check
# table(analytical_sample$ETHNIC_label, analytical_sample$Black, useNA = "ifany")
# table(analytical_sample$ETHNIC_label, analytical_sample$Hispanic, 
#       useNA = "ifany")

#---- all-way contingency table ----
cross_class_label <- table(analytical_sample$ETHNIC_label, 
                           analytical_sample$Astroke) %>% as.data.frame()

# #How many are missing from this table?-- only 144! 
# sum(cross_class_label$Freq)

#---- Bayes Stuff ----
#---- **number of runs ----
B = 2

#---- **priors ----


#---- **sampling ----


