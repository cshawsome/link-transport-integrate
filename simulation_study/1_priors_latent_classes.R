#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **ADAMS analytic ----
ADAMS_analytic <- 
  read_csv(paste0(path_to_box, "data/ADAMS/cleaned/ADAMS_analytic.csv"))

#---- **variable labels ----
variable_labels <- read_csv(paste0(path_to_box, "data/variable_crosswalk.csv"))

#---- **variable selection results ----
selected_vars <- 
  read_csv(paste0(path_to_box, "analyses/simulation_study/variable_selection/", 
                  "model_coefficients.csv")) 

#---- impute data ----
#---- **pare down list of vars based on ADAMS variable selection ----
remove <- c("ANRECNO", paste0("r", seq(4, 7), "adla"), "avg_proxy_cog_Better", 
            "avg_proxy_cog_Same", "avg_proxy_cog_Worse", "proxy_Spouse", 
            "proxy_Child", "proxy_Other_Relative", "proxy_Other")

ADAMS_subset <- ADAMS_analytic %>% dplyr::select(-all_of(remove))

#---- **check missings ----
needs_imputing <- names(which(colSums(is.na(ADAMS_subset)) != 0))

#---- **define imputation var types ----
hrs_waves <- seq(4, 7)
cog_waves <- seq(5, 7)

sociodem_vars <- c("AAGE", "EDYRS", "Female", "Black", "Hispanic", 
                   "Married/partnered", "Working", "Retired", "Not working")

ADAMS_vars <- c("SELFCOG", "ANMSETOT", "ANSER7T", "ANSCISOR", "ANCACTUS", 
                "ANPRES", "ANAFTOT", "ANBNTTOT", "ANDELCOR", "ANRECYES", 
                "ANRECNO", "ANWM1TOT", "ANWM2TOT", "ANCPTOT", "ANRCPTOT", 
                "ANTMASEC", "ANBWC20", "ANBWC86", "ANIMMCR", "ANSMEM2_Better", 
                "ANSMEM2_Same", "ANSMEM2_Worse", "Dementia", "Other", 
                "MCI")

HRS_vars <- c(paste0("r", cog_waves, "mpart"), paste0("r", hrs_waves, "height"), 
              paste0("r", hrs_waves, "weight"), paste0("r", hrs_waves, "smoken"), 
              paste0("r", hrs_waves, "drinkd"), paste0("r", hrs_waves, "drinkn"), 
              paste0("r", hrs_waves, "hibpe"), paste0("r", hrs_waves, "diabe"), 
              paste0("r", hrs_waves, "hearte"), paste0("r", hrs_waves, "stroke"), 
              paste0("r", hrs_waves, "adla"), paste0("r", hrs_waves, "iadla"), 
              paste0("r", cog_waves, "imrc"), paste0("r", cog_waves, "dlrc"), 
              paste0("r", cog_waves, "ser7"), paste0("r", cog_waves, "bwc20"), 
              paste0("r", seq(5, 6), "bwc86"), paste0("r", cog_waves, "scis"), 
              paste0("r", cog_waves, "cact"), paste0("r", cog_waves, "pres"), 
              paste0("r", cog_waves, "cogtot"), 
              paste0("r", cog_waves, "pstmem_Better"), 
              paste0("r", cog_waves, "pstmem_Same"), 
              paste0("r", cog_waves, "pstmem_Worse"))

not_predictors <- c("HHIDPN", "AYEAR", "White", "ANMSETOT_norm", 
                    "Adem_dx_label_collapsed", "Unimpaired", "Astroke", "Ahibpe", 
                    "Adiabe", "Ahearte", "Abmi", "Aiadla", "Aadla", "Asmoken", 
                    "Adrinkd", "Adrinkn")

# #Sanity check
# union_var_names <- c(sociodem_vars, ADAMS_vars, HRS_vars, not_predictors)
# ncol(ADAMS_analytic) == length(union_var_names)
# 
# colnames(ADAMS_analytic)[which(!colnames(ADAMS_analytic) %in% union_var_names)]
# union_var_names[which(!union_var_names %in% colnames(ADAMS_analytic))]

#---- predictor matrix ----
predict <- 
  matrix(1, nrow = length(needs_imputing), ncol = ncol(ADAMS_analytic)) %>% 
  set_rownames(needs_imputing) %>% set_colnames(colnames(ADAMS_analytic))

#---- **cannot predict themselves ----
predict[needs_imputing, needs_imputing] <- 
  (diag(x = 1, nrow = length(needs_imputing), 
        ncol = length(needs_imputing)) == 0)*1

#---- **non-predictors ----
predict[, not_predictors] <- 0

#---- **ADAMS-ADAMS prediction ----
#use ADAMS variables to predict each other: remove HRS predictors for now
predict[needs_imputing[which(needs_imputing %in% ADAMS_vars)], HRS_vars] <- 0

#---- **HRS-HRS prediction ----
#use HRS variables to predict each other: remove ADAMS predictors 
predict[needs_imputing[which(needs_imputing %in% HRS_vars)], ADAMS_vars] <- 0

#---- **ADAMS extra predictors ----
predict["ANBWC20", paste0("r", cog_waves, "bwc20")] <- 1
predict["ANBWC86", paste0("r", seq(5, 6), "bwc86")] <- 1
predict["ANCACTUS", paste0("r", cog_waves, "cact")] <- 1
predict["SELFCOG", paste0("r", cog_waves, "cogtot")] <- 1
predict["ANDELCOR", paste0("r", cog_waves, "dlrc")] <- 1
predict["ANIMMCR", paste0("r", cog_waves, "imrc")] <- 1
predict["ANPRES", paste0("r", cog_waves, "pres")] <- 1
predict["ANSCISOR", paste0("r", cog_waves, "scis")] <- 1
predict["ANSER7T", paste0("r", cog_waves, "ser7")] <- 1
predict[paste0("ANSMEM2_", c("Better", "Same", "Worse")), 
        paste0("r", cog_waves, "pstmem_Better")] <- 1
predict[paste0("ANSMEM2_", c("Better", "Same", "Worse")), 
        paste0("r", cog_waves, "pstmem_Same")] <- 1
predict[paste0("ANSMEM2_", c("Better", "Same", "Worse")), 
        paste0("r", cog_waves, "pstmem_Worse")] <- 1

#---- **derived vars (do not impute) ----
remove <- c("ANMSETOT_norm", 
            paste0("A", c("stroke", "hibpe", "diabe", "hearte", "bmi", "iadla", 
                          "adla", "smoken", "drinkd", "drinkn"))) 

predict <- predict[!rownames(predict) %in% remove, ]

#---- **sanity check ----
#non-predictors should all be 0
colSums(predict[, not_predictors])

#---- imputation ----
set.seed(20220202)
start <- Sys.time()
fast_impute(predictor_matrix = predict, data = ADAMS_analytic, 
            path_for_output = paste0(path_to_box, "data/ADAMS/cleaned/"),
            method = "PMM", m = 25, maxit = 15)
end <- Sys.time() - start



#---- OLD ----

#---- predictors ----
#unimpaired model predictors
unimpaired_preds <- c("AAGE", "Black", "Hispanic", "ANMSETOT_norm", "ANSER7T", 
                      "ANIMMCR", "ANRECYES", "ANWM1TOT", "proxy_cog")

#other model predictors
other_preds <- c("AAGE", "ANMSETOT_norm", "ANIMMCR", "ANDELCOR")

#mci model predictors
mci_preds <- c("ANMSETOT_norm", "ANIMMCR", "Aiadla", "Astroke", "Abmi")

#---- model ----
bootstrap_models <- function(prop){
  subsample <- sample_frac(ADAMS_subset, size = prop, replace = TRUE)
  
  unimpaired_model <- 
    glm(formula(paste("AUnimpaired ~ ", 
                      paste(unimpaired_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample)
  
  other_model <- 
    glm(formula(paste("AOther ~ ", paste(other_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(AUnimpaired == 0))
  
  mci_model <- 
    glm(formula(paste("AMCI ~ ", paste(mci_preds, collapse = " + "), 
                      collapse = "")), family = "binomial", data = subsample %>% 
          filter(AUnimpaired == 0 & AOther == 0)) 
  
  return(list("unimpaired_betas" = coefficients(unimpaired_model),
              "other_betas" = coefficients(other_model), 
              "mci_betas" = coefficients(mci_model), 
              "unimpaired_cov" = as.vector(vcov(unimpaired_model)), 
              "other_cov" = as.vector(vcov(other_model)), 
              "mci_cov" = as.vector(vcov(mci_model))))
}

bootstrap_runs <- replicate(10000, bootstrap_models(prop = 1), simplify = FALSE)

#---- check distributions ----
for(group in c("unimpaired", "other", "mci")){
  data <- lapply(bootstrap_runs, "[[", paste0(group, "_betas")) %>% 
    do.call(rbind, .) %>% t() %>% as.data.frame() 
  
  for(var in rownames(data)){
    show(hist(as.numeric(data[var, ]), main = paste0(group, " ", var)))
  }
}

#---- format output ----
for(est in c("betas", "cov")){
  for(group in c("unimpaired", "other", "mci")){
    data <- lapply(bootstrap_runs, "[[", paste0(group, "_", est)) %>% 
      do.call(rbind, .) %>% t() %>% as.data.frame() 
    
    if(est == "betas"){
      data %<>% mutate("preds" = c("(Intercept)", get(paste0(group, "_preds"))))
    }
    
    data %>% 
      write_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/priors/", 
                       "ADAMS_test/latent_class_", group, "_", est, ".csv"))
  }
}
