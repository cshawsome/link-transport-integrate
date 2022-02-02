#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
ADAMS_imputed <- 
  readRDS(paste0(path_to_box, "data/ADAMS/cleaned/MI/MI_datasets"))

#---- derive post-imputation variables ----
test <- ADAMS_imputed[[1]]
hrs_waves <- seq(4, 7)

#---- **bmi ----
derive_bmi <- function(data, waves){
  for(wave in waves){
    data %<>% mutate(!!paste0("r", wave, "bmi_derived") := 
                       !!sym(paste0("r", wave, "weight"))/
                       (!!sym(paste0("r", wave, "height"))*
                          !!sym(paste0("r", wave, "height"))))
  }
  return(data)
}

ADAMS_imputed <- lapply(ADAMS_imputed, derive_bmi, hrs_waves)

# #Sanity check
# test <- ADAMS_imputed[[25]]
# View(test[, c("r4weight", "r4height", "r4bmi_derived")])
# head(test[, "r4weight"]/(test[, "r4height"]^2))
# tail(colnames(test))

#---- define variable types ----

cog_waves <- seq(5, 7)

sociodem_vars <- c("AAGE", "EDYRS", "Female", "Black", "Hispanic", 
                   "Married/partnered", "Working", "Retired", "Not working")

ADAMS_vars <- c("SELFCOG", "ANMSETOT", "ANSER7T", "ANSCISOR", "ANCACTUS", 
                "ANPRES", "ANAFTOT", "ANBNTTOT", "ANDELCOR", "ANRECYES", 
                "ANRECNO", "ANWM1TOT", "ANWM2TOT", "ANCPTOT", "ANRCPTOT", 
                "ANTMASEC", "ANBWC20", "ANBWC86", "ANIMMCR", "ANSMEM2_Better", 
                "ANSMEM2_Same", "ANSMEM2_Worse", "avg_proxy_cog_Better", 
                "avg_proxy_cog_Same", "avg_proxy_cog_Worse", "Dementia", "Other", 
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


