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
#store measures and number of levels for each variable
person_level_vars <- matrix(nrow = 1, ncol = 5) %>% 
  set_colnames(c("GENDER", "ETHNIC", "AAGE_cat", "EDYRS_cat", 
                 "num_assessments"))
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

#---- number of assessments ----
#needs to be generalized later to sum over the number of assessments available
# for a person
analytical_sample %<>% mutate("num_assessments" = 1)

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

#---- ***person-level counts ----
person_level_counts <- analytical_sample %>% 
  dplyr::select(colnames(person_level_vars)) %>% 
  apply(., 2, function(x) max(x, na.rm = TRUE))

person_level_vars[1, names(person_level_counts)] <- person_level_counts

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

phi_list <- 
  lapply(phi_list <- vector(mode = "list", group_class_n),
         function(x) x <- lapply(x <- vector(mode = "list", sub_class_n),
                                 function(x) x <- 
                                   vector(mode = "list", 
                                          length = ncol(measures))))

lambda_list <- lapply(lambda_list <- vector(mode = "list", group_class_n),
                      function(x) x <- vector(mode = "list",
                                              ncol(person_level_vars)))

p_Mij_list <-
  lapply(p_Mij_list <- vector(mode = "list", nrow(rsamp)),
         function(x) x <-
           lapply(x <- vector(mode = "list",
                              length(timepoints)),
                  function(x) x <- vector(length = sub_class_n)))

p_Gi <- matrix(nrow = nrow(rsamp), ncol = group_class_n)

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
  for(k in 1:ncol(measures)){
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
  for(k in 1:ncol(person_level_vars)){
    for(g in 1:group_class_n){
      #these are the a_ks... might want to change this later to be empirical
      # dist
      pars <- rep(1, person_level_vars[1, k])
      var <- paste0(colnames(person_level_vars)[k])
      subset <- rsamp %>% filter(group_class == g)
      if(nrow(subset > 0)){
        counts <- table(subset[, var])
        pars[as.numeric(names(counts))] = pars[as.numeric(names(counts))] + 
          counts
      }
      lambda_list[[g]][[k]] <- rdirichlet(n = 1, alpha = pars)
    }
  }
  
  #---- **sample v_gm ----
  for(g in 1:group_class_n){
    for(m in 1:(sub_class_n - 1)){
      pars <- c(1, beta_chain[b])
      for(j in 1:timepoints){
        subclass <- paste0("sub_class_", j)
        subset <- rsamp %>% filter(group_class == g)
        if(nrow(subset > 0)){
          counts <- table(subset[, subclass])
          pars[1] = pars[1] + sum(counts[which(as.numeric(names(counts)) == m)])
          pars[2] = pars[2] + sum(counts[which(as.numeric(names(counts)) > m)])
        }
      }
      v_gm[g, m] <- rbeta(n = 1, shape1 = pars[1], shape2 = pars[2])
    }
  }
  v_gm[, sub_class_n] <- 1
  
  #---- ***calculate omega_gm ----
  comp_probs <- 1 - v_gm
  
  omega_gm[, 1] <- v_gm[, 1]
  omega_gm[, 2] <- v_gm[, 2]*comp_probs[, 1]
  
  for(m in 3:sub_class_n){
    rowProd <- apply(comp_probs[, 1:(m - 1)], 1, prod)
    omega_gm[, m] <- v_gm[, m]*rowProd
  }
  
  #---- **sample u_g ----
  for(g in 1:(group_class_n - 1)){
    shape1 = as.numeric(sum(1 + table(rsamp$group_class)[g], na.rm = TRUE))
    shape2 = as.numeric(
      sum(alpha_chain[b] + 
            sum(table(rsamp$group_class)[(g + 1):group_class_n], na.rm = TRUE)))
    
    u_g[g] <- rbeta(n = 1, shape1 = shape1, shape2 = shape2)
  }
  
  #---- ***calculate pi_g ----
  comp_probs <- 1 - u_g
  
  pi_chain[1, b] <- u_g[1]
  pi_chain[2, b] <- u_g[2]*comp_probs[1]
  
  for(g in 3:group_class_n){
    pi_chain[g, b] <- u_g[g]*prod(comp_probs[1:(g - 1)])
  }
  
  #---- **sample Mij ----
  for(i in 1:nrow(rsamp)){
    for(j in 1:timepoints){
      for(m in 1:sub_class_n){
        group <- as.numeric(rsamp[i, "group_class"])
        omega_Gim <- omega_gm[group, m]
        num_prod = 1
        
        for(k in 1:ncol(measures)){
          var <- paste0(colnames(measures)[k], "_", j)
          num_prod = 
            prod(num_prod, 
                 phi_list[[group]][[m]][[k]][as.numeric(rsamp[i, var])])
        }
        
        den_sum = 0
        for(s in 1:sub_class_n){
          den_prod = 1
          for(k in 1:ncol(measures)){
            var <- paste0(colnames(measures)[k], "_", j)
            den_prod = 
              prod(den_prod, 
                   phi_list[[group]][[s]][[k]][as.numeric(rsamp[i, var])])
          }
          den_sum = sum(den_sum, omega_gm[group, s]*den_prod)
        }
        
        p_Mij_list[[i]][[j]][[m]] = (omega_Gim*num_prod)/den_sum
      }
      rsamp[i, paste0("sub_class_", j)] <- 
        which(rmultinom(n = 1, size = 1, prob = p_Mij_list[[i]][[j]]) == 1)
    }
  }
  
  #---- **sample Gi ----
  for(i in 1:nrow(rsamp)){
    for(g in 1:group_class_n){
      pi_g <- pi_chain[g, b]
      num_prod_1 = 1
      for(k in 1:ncol(person_level_vars)){
        var <- colnames(person_level_vars)[k]
        lambda_gk <- lambda_list[[g]][[k]][as.numeric(rsamp[i, var])]
        num_prod_2 = 1
        for(j in 1:timepoints){
          num_sum = 0
          for(m in 1:sub_class_n){
            num_prod_3 = 1
            for(k2 in 1:ncol(measures)){
              var2 <- paste0(colnames(measures)[k2], "_", j)
              num_prod_3 = 
                prod(num_prod_3, 
                     phi_list[[g]][[m]][[k2]][as.numeric(rsamp[i, var2])])
            }
            num_sum = sum(num_sum, omega_gm[g, m]*num_prod_3)
          }
          num_prod_2 = prod(num_prod_2, num_sum)
        }
        num_prod_1 = prod(num_prod_1, lambda_gk*num_prod_2)
      }
      
      numerator = pi_g*num_prod_1
      
      dem_sum_1 = 0
      for(f in 1:group_class_n){
        pi_f <- pi_chain[f, b]
        dem_prod_1 = 1
        for(k in 1:ncol(person_level_vars)){
          var <- colnames(person_level_vars)[k]
          lambda_fk <- lambda_list[[f]][[k]][as.numeric(rsamp[i, var])]
          dem_prod_2 = 1
          for(j in 1:timepoints){
            dem_sum_2 = 0
            for(m in 1:sub_class_n){
              dem_prod_3 = 1
              for(k2 in 1:ncol(measures)){
                var2 <- paste0(colnames(measures)[k2], "_", j)
                dem_prod_3 = 
                  prod(dem_prod_3, 
                       phi_list[[f]][[m]][[k2]][as.numeric(rsamp[i, var2])])
              }
              dem_sum_2 = sum(dem_sum_2, omega_gm[f, m]*dem_prod_3)
            }
            dem_prod_2 = prod(dem_prod_2, dem_sum_2)
          }
          dem_prod_1 = prod(dem_prod_1, lambda_fk*dem_prod_2)
        }
        dem_sum_1 = sum(dem_sum_1, pi_f*dem_prod_1)
      }
      p_Gi[i, g] <- numerator/dem_sum_1
    }
  }
  
}












