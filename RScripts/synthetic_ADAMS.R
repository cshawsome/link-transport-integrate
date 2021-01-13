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
measure_level_vars <- c("ANMSETOT_cat", paste0("r", seq(5, 7), "iadla_cat"))

analytical_sample <- all_data %>% 
  dplyr::select("HHIDPN", 
                all_of(person_level_vars), all_of(measure_level_vars)) %>% 
  na.omit()

#Variable check
colSums(is.na(analytical_sample))
dim(analytical_sample)

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
M_ij <- matrix(1, nrow = nrow(rsamp), ncol = length(measure_level_vars))

#Everyone has the same number of measures/assessments for dementia
rsamp[, "num_assessments"] <- length(measure_level_vars)
person_level_vars <- c(person_level_vars, "num_assessments")

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
v_gm <- matrix(nrow = sub_class_n, ncol = group_class_n)
for(g in 1:group_class_n){
  v_gm[1:(sub_class_n - 1), g] <- 
    rbeta(n = (sub_class_n - 1), shape1 = 1, shape2 = beta)
}

#The last v is fixed at 1 for all groups
v_gm[nrow(v_gm), ] <- 1

#---- Step 4: initiate chain and parameter storage ----
beta_chain <- vector(length = B)
alpha_chain <- vector(length = B)
pi_chain <- matrix(ncol = B, nrow = group_class_n)
omega_gm <- matrix(nrow = sub_class_n, ncol = group_class_n)

p_M_ij_list <- 
  lapply(phi_list <- vector(mode = "list", nrow(rsamp)), 
         function(x) x <- 
           lapply(x <- vector(mode = "list", 
                              length(measure_level_vars)), 
                  function(x) x <- vector(mode = "list", sub_class_n)))

phi_list <- lapply(phi_list <- vector(mode = "list", group_class_n), function(x) 
  x <- lapply(x <- vector(mode = "list", sub_class_n), function(x) 
    x <- vector(mode = "list", length(measure_level_vars))))

lambda_list <- lapply(lambda_list <- vector(mode = "list", group_class_n), 
                      function(x) x <- vector(mode = "list", 
                                              length(person_level_vars)))

#---- Step 5: sampling ----
for(b in 1:B){
  #---- **sample beta ----
  a_beta = a_beta + group_class_n*(sub_class_n - 1)
  b_beta = b_beta - sum(log(1 - v_gm[1:(sub_class_n - 1), ]))
  beta_chain[b] <- rgamma(n = 1, shape = a_beta, rate = b_beta)
  
  #---- **sample alpha ----
  a_alpha = a_alpha + group_class_n - 1
  b_alpha = b_alpha - sum(log(1 - head(u_g, -1)))
  alpha_chain[b] = rgamma(n = 1, shape = a_alpha, rate = b_alpha)
  
  #---- **sample phi ----
  for(k in 1:length(measure_level_vars)){
    for(g in 1:group_class_n){
      for(m in 1:sub_class_n){
        people <- which(M_ij[, k] == m)
        subset <- rsamp[people, ] %>% filter(dem_group == g)
        cat_count <- max(rsamp[, measure_level_vars[k]], na.rm = TRUE)
        if(nrow(subset) > 0){
          #if you use the table function, you might miss some levels
          pars <- vector(length = cat_count)
          for(d in 1:length(pars)){
            pars[d] = sum(subset[, measure_level_vars[k]] == d)
          }
        } else{
          pars = rep(0, cat_count)
        }
        phi_list[[g]][[m]][[k]] <- rdirichlet(n = 1, alpha = 1 + pars)
        }
      }
  }
  
  #---- **sample lambda ----
  for(k in 1:length(person_level_vars)){
    for(g in 1:group_class_n){
      subset <- rsamp %>% filter(dem_group == g)
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
    for(m in 1:(sub_class_n - 1)){
      people <- which(M_ij[, k] == m)
      next_groups <- which(M_ij[, k] > m)
      subset <- rsamp[people, ] %>% filter(dem_group == g)
      shape1 = as.numeric(sum(1 + table(subset$measure_group)[m], na.rm = TRUE))
      shape2 = as.numeric(
        sum(beta_chain[b] + 
              sum(table(subset$measure_group)[(m + 1):sub_class_n], 
                  na.rm = TRUE)))
      v_gm[m, g] <- rbeta(n = 1, shape1 = shape1, shape2 = shape2)
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











