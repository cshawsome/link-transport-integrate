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

#---- select variables ----
#based on analysis in priors_latent_classes.R
vars <- c("AAGE", "ETHNIC_label", "ANMSETOT", "ANSER7T", "ANIMMCR", "ANRECYES", 
          "ANWM1TOT", "proxy_cog", "ANDELCOR", "Aiadla", "Astroke", "Abmi")

analytical_sample <- ADAMS_subset %>% 
  dplyr::select("HHIDPN", all_of(vars)) 

#---- all-way contingency table ----
cross_class_label <- table(analytical_sample$ETHNIC_label, 
                           analytical_sample$Astroke) %>% as.data.frame()

# #How many are missing from this table?-- only 144! 
# sum(cross_class_label$Freq)

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


