#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "DirichletReg")

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
#store measures and number of levels for the measure variable
measures <- matrix(nrow = 1, ncol = 2) %>% set_colnames(c("MMSE", "IADLA"))
measure_vars <- c("ANMSETOT_cat", paste0("r", seq(5, 7), "iadla_cat"))

analytical_sample <- all_data %>% 
  dplyr::select("HHIDPN", "AYEAR", 
                all_of(person_level_vars), all_of(measure_vars)) %>% na.omit()

#Variable check-- there's 538 people in the complete data set
colSums(is.na(analytical_sample))
dim(analytical_sample)

#---- create measure-level vars ----
timepoints <- 1

#---- **MMSE ----
analytical_sample[, "MMSE_1"] <- analytical_sample[, "ANMSETOT_cat"]
measures[1, "MMSE"] <- max(analytical_sample$MMSE_1, na.rm = TRUE)

#drop original variable
analytical_sample %<>% dplyr::select(-"ANMSETOT_cat")

#---- **IADLA ----
#Take the IADLA measure closest to ADAMS interview year
analytical_sample %<>% 
  mutate("IADLA_1" = case_when(AYEAR == 2001 ~ r5iadla_cat, 
                               AYEAR %in% c(2002, 2003) ~ r6iadla_cat, 
                               AYEAR == 2004 ~ r7iadla_cat))
measures[1, "IADLA"] <- max(analytical_sample$IADLA_1, na.rm = TRUE)

# #Sanity check
# table(analytical_sample$IADLA_1, useNA = "ifany")

#Get rid of original variables
analytical_sample %<>% dplyr::select(-c(contains("iadla_cat"), "AYEAR")) 

#---- Mixture Model ----
#---- **Step 1: SRS of data ----
samp_size <- floor(0.5*nrow(analytical_sample))
rsamp <- sample_n(analytical_sample, size = samp_size, replace = FALSE)

#---- **Step 2: hyperpriors ----
#---- ***number of latent classes ----
#number of person-level latent classes
group_class_n <- 4
#number of assessment-level classes
sub_class_n <- 4

#---- ***assign latent classes ----
#Everyone starts in class 1 at the person level and sub level 
#(normal cog maybe?) 
rsamp[, "group_class"] <- 1
for(j in 1:timepoints){
  rsamp[, paste0("sub_class_", j)] <- 1
}

#Everyone has one set of assessments in this example-- need to generalize later
  #need to sum over the assessments that a person has
rsamp[, "num_assessments"] <- 1
person_level_vars <- c(person_level_vars, "num_assessments")

#---- ***person-level latent classes parameter ----
#From Dunson and Xing (2009)
a_alpha = 0.25
b_alpha = 0.25
alpha <- rgamma(n = 1, shape = a_alpha, rate = b_alpha)

#---- ***sub-level latent classes parameter ----
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
v_gm <- matrix(ncol = sub_class_n, nrow = group_class_n)
for(g in 1:group_class_n){
  v_gm[g, 1:(sub_class_n - 1)] <- 
    rbeta(n = (sub_class_n - 1), shape1 = 1, shape2 = beta)
}

#The last v is fixed at 1 for all groups
v_gm[, ncol(v_gm)] <- 1

#---- Step 4: initiate chain and parameter storage ----
beta_chain <- vector(length = B)
alpha_chain <- vector(length = B)
pi_chain <- matrix(ncol = B, nrow = group_class_n)
omega_gm <- matrix(ncol = sub_class_n, nrow = group_class_n)

# p_M_ij_list <- 
#   lapply(phi_list <- vector(mode = "list", nrow(rsamp)), 
#          function(x) x <- 
#            lapply(x <- vector(mode = "list", 
#                               length(measure_level_vars)), 
#                   function(x) x <- vector(length = sub_class_n)))
# 
phi_list <- 
  lapply(phi_list <- vector(mode = "list", group_class_n),
         function(x) x <- lapply(x <- vector(mode = "list", sub_class_n),
                                 function(x) x <- 
                                   vector(mode = "list", 
                                          length = ncol(measures))))
# 
# lambda_list <- lapply(lambda_list <- vector(mode = "list", group_class_n), 
#                       function(x) x <- vector(mode = "list", 
#                                               length(person_level_vars)))

#---- Step 5: sampling ----
for(b in 1:B){
  #---- **sample beta ----
  a_beta = a_beta + group_class_n*(sub_class_n - 1)
  b_beta = b_beta - sum(log(1 - v_gm[, 1:(sub_class_n - 1)]))
  beta_chain[b] <- rgamma(n = 1, shape = a_beta, rate = b_beta)
  
  #---- **sample alpha ----
  a_alpha = a_alpha + group_class_n - 1
  b_alpha = b_alpha - sum(log(1 - head(u_g, -1)))
  alpha_chain[b] = rgamma(n = 1, shape = a_alpha, rate = b_alpha)
  
  #---- **sample phi ----
  for(k in 1:length(measures)){
    for(g in 1:group_class_n){
      for(m in 1:sub_class_n){
        #these are the a_ks... might want to change this later to be empirical
        # dist
        pars <- rep(1, measures[1, k])
        for(j in 1:timepoints){
          subclass <- paste0("sub_class_", j)
          var <- paste0(colnames(measures)[k], "_", j)
          subset <- rsamp %>% filter(group_class == g & !!sym(subclass) == m)
          if(nrow(subset > 0)){
            counts <- table(subset[, var])
            pars[as.numeric(names(counts))] = pars[as.numeric(names(counts))] + 
              counts
          }
        }
        phi_list[[g]][[m]][[k]] <- rdirichlet(n = 1, alpha = pars)
      }
    }
  }
  
  #---- **sample lambda ----
  for(k in 1:length(person_level_vars)){
    for(g in 1:group_class_n){
      subset <- rsamp %>% filter(dem_group == g)
      for(j in 1:timepoint)
      
      
      cat_count <- max(rsamp[, person_level_vars[k]], na.rm = TRUE)
      if(nrow(subset) > 0){
        #if you use the table function, you might miss some levels
        pars <- vector(length = cat_count)
        for(d in 1:length(pars)){
          pars[d] = sum(subset[, person_level_vars[k]] == d)
        }
      } else{
        pars = rep(0, cat_count)
      }
      lambda_list[[g]][[k]] <- rdirichlet(n = 1, alpha = 1 + pars)
    }
  }
  
  #---- **sample v_gm ----
  for(g in 1:group_class_n){
    people <- which(rsamp[, "dem_group"] == g)
    if(length(people) > 0){
      M_ij_subset <- M_ij[people, ]
      for(m in 1:(sub_class_n - 1)){
        shape1 = 1 + sum(M_ij_subset == m)
        shape2 = beta_chain[b] + sum(M_ij_subset > m)
        v_gm[m, g] <- rbeta(n = 1, shape1 = shape1, shape2 = shape2)
      }
    } else{
      v_gm[1:(sub_class_n - 1), g] <- 
        rbeta(n = (sub_class_n - 1), shape1 = 1, shape2 = beta_chain[b])
    }
  }
  
  #---- ***calculate omega_gm ----
  comp_probs <- 1 - v_gm
  
  omega_gm[1, ] <- v_gm[1, ]
  omega_gm[2, ] <- v_gm[2, ]*comp_probs[1, ]
  
  for(m in 3:sub_class_n){
    colProd <- apply(comp_probs[1:(m - 1), ], 2, prod)
    omega_gm[m, ] <- v_gm[m, ]*colProd
  }
  
  #---- **sample u_g ----
  for(g in 1:(group_class_n - 1)){
    shape1 = as.numeric(sum(1 + table(rsamp$dem_group)[g], na.rm = TRUE))
    shape2 = as.numeric(
      sum(alpha_chain[b] + 
            sum(table(rsamp$dem_group)[(g + 1):group_class_n], na.rm = TRUE)))
    
    u_g[g] <- rbeta(n = 1, shape1 = shape1, shape2 = shape2)
  }
  
  #---- ***calculate pi_g ----
  comp_probs <- 1 - u_g
  
  pi_g[1, b] <- u_g[1]
  pi_g[2, b] <- u_g[2]*comp_probs[1]
  
  for(g in 3:group_class_n){
    pi_g[g, b] <- u_g[g]*prod(comp_probs[1:(g - 1)])
  }
  
  #---- **sample M_ij ----
  for(i in 1:nrow(rsamp)){
    group = as.numeric(rsamp[i, "dem_group"])
    X_ijk = rsamp[i, measure_level_vars]
    for(m in 1:sub_class_n){
      omega_val <- omega_gm[m, group]
      phi_vec <- vector(length = length(measure_level_vars))
      for(k in 1:length(phi_vec)){
        phi_vec[k] <- 
          as.numeric(phi_list[[group]][[m]][[k]])[as.numeric(X_ijk[k])]
      }
      num_prob <- omega_val*prod(phi_vec)
    }
  }
  
}











