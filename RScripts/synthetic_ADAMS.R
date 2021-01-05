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
person_level_vars <- c("GENDER", "ETHNIC", "AAGE_cat", "EDYRS_cat")
measure_level_vars <- c("ANMSETOT_cat", paste0("r", seq(5, 7), "iadla_cat"))

analytical_sample <- all_data %>% 
  dplyr::select("HHIDPN", all_of(person_level_vars), all_of(measure_level_vars))

#Variable check
colSums(is.na(analytical_sample))

#---- Mixture Model ----
#---- **Step 1: SRS of data ----
samp_size <- floor(0.5*nrow(analytical_sample))
rsamp <- sample_n(analytical_sample, size = samp_size, replace = FALSE)

#---- **Step 2: hyperpriors ----
#---- ***number of latent classes ----
#number of group-level latent classes
group_class_n <- 4
#number of sub-latent classes
sub_class_n <- 5

#---- ***assign latent classes ----
#Everyone starts in class 1 at both the person level and measurement level
rsamp[, "dem_group"] <- 1
rsamp[, "measure_group"] <- 1

person_level_vars <- c(person_level_vars, "dem_group")
measure_level_vars <- c(measure_level_vars, "measure_group")

#---- ***person-level latent classes parameter ----
a_alpha = 1
b_alpha = 0.5
alpha <- rgamma(n = 1, shape = a_alpha, rate = b_alpha)

#---- ***measure-level latent classes parameter ----
#From Dunson and Xing (2009)
#These can differ by group-level latent classes if we wish, but we're keeping
# it "simple" for now
a_beta = 0.25
b_beta = 0.25
beta <- rgamma(n = 1, shape = a_beta, rate = b_beta)

#---- Step 3: initiate values ----
#Runs
B = 2

#Draw priors for one less than the number of group-level latent classes
u_g <- rbeta(n = (group_class_n - 1), shape1 = 1, shape2 = alpha)
#The last u is fixed at 1
u_g <- c(u_g, 1)

#Draw priors for one less than the number of sub-level latent classes
v_m <- rbeta(n = (sub_class_n - 1), shape1 = 1, shape2 = beta)
#The last v is fixed at 1
v_m <- c(v_m, 1)

#---- Step 4: initiate chain storage ----
beta_post <- vector(length = B)
alpha_post <- vector(length = B)

#---- Step 5: sampling ----
for(i in 1:B){
  #---- **sample beta ----
  a_beta_post = a_beta + group_class_n*(sub_class_n - 1)
  b_beta_post = b_beta - 4*sum(log(1 - head(v_m, -1)))
  beta_post[i] <- rgamma(n = 1, shape = a_beta_post, rate = b_beta_post)
  
  #---- **sample alpha ----
  a_alpha_post = a_alpha + group_class_n - 1
  b_alpha_post = b_alpha - sum(log(1 - head(u_g, -1)))
  alpha_post[i] = rgamma(n = 1, shape = a_alpha_post, rate = b_alpha_post)
  
  #---- **sample phi ----
  
  
}











