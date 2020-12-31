#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
#---- **ADAMS ----
ADAMS_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                                "data/cleaned/ADAMS_subset.csv"), 
                         col_types = cols(HHIDPN = col_character()))
  
#---- **RAND ----
RAND_subset <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/RAND_subset.csv"), 
                        col_types = cols(HHIDPN = col_character()))

#---- join data ----
all_data <- left_join(ADAMS_subset, RAND_subset, by = "HHIDPN")

#---- select variables ----
vars <- c("HHIDPN", "GENDER", "ETHNIC", "AAGE_cat", "EDYRS_cat", "ANMSETOT_cat", 
          paste0("r", c(5, 6, 7), "iadla_cat"))

analytical_sample <- all_data %>% dplyr::select(all_of(vars))

#Variable check
colSums(is.na(analytical_sample))

#---- Mixture Model ----
#---- **Step 1: SRS of data ----
samp_size <- floor(0.5*nrow(analytical_sample))
rsamp <- sample_n(analytical_sample, size = samp_size, replace = FALSE)

#---- **Step 2: hyperpriors ----
#---- *** number of latent classes ----
#number of group-level latent classes
group_class_n <- 4
#number of sub-latent classes
sub_class_n <- 5

#---- *** p(group-level latent classes) ----
a_alpha = 1
b_alpha = 0.5

#---- *** p(sub-level latent classes) ----
#From Dunson and Xing (2009)
#These can differ by group-level latent classes if we wish, but we're keeping
# it "simple" for now
a_beta = 0.25
b_beta = 0.25

#Draw priors for one less than the number of sub-level latent classes
v_g <- rgamma(n = (sub_class_n - 1), shape = a_beta, rate = b_beta)
#The last v is fixed at 1
v_g <- c(v_g, 1)

#---- Step 3: sampling ----
#---- **initiate values ----
#Draw priors for one less than the number of group-level latent classes
u_g <- rgamma(n = (group_class_n - 1), shape = a_alpha, rate = b_alpha)
#The last u is fixed at 1
u_g <- c(u_g, 1)

#Draw priors for one less than the number of sub-level latent classes
v_g <- rgamma(n = (sub_class_n - 1), shape = a_beta, rate = b_beta)
#The last v is fixed at 1
v_g <- c(v_g, 1)

#---- **sample beta ----
a_beta_post = a_beta + group_class_n*(sub_class_n - 1)
b_beta_post = 










