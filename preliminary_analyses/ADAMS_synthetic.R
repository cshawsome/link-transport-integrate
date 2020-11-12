#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "stringr", "plyr", "broom", "openxlsx")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data: neuropsych ----
neuropsych_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
neuropsych_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                 "ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych_A <- read_da_dct(neuropsych_data_path_A, 
                                  neuropsych_dict_path_A, HHIDPN = "TRUE")

#---- import data: sociodemographic ----
ADAMS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                  "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                  "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") 

#---- import data: demdx ----
demdx_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1a/adams1ada/ADAMS1AD_R.da")
demdx_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                            "ADAMS/adams1a/adams1asta/ADAMS1AD_R.dct")

ADAMSA_demdx <- read_da_dct(demdx_data_path_A, demdx_dict_path_A, 
                            HHIDPN = "TRUE") 

#---- format data: neuropsych ----
#variables of interest
ADAMS_vars <- c("HHIDPN",
                #item-level
                "SCISOR", "CACTUS", "ANPRES", 
                #total scores
                "ANMSETOT", "SER7T", "ANAFTOT", "CPTOT", "BWC86", "IMMCR", 
                "DELCOR", "RECYES", "RECNO", "WM1TOT", "WM2TOT", "MASEC", 
                "MBSEC", "ANSDMTOT")

not_these <- c("ANRCPTOT")

ADAMSA_assessment <- ADAMS_neuropsych_A %>% 
  dplyr::select(contains(ADAMS_vars)) %>% 
  dplyr::select(-c(contains(not_these)))

#Need to recode 6 = start over to 0
recode_incorrect <- 
  colnames(ADAMSA_assessment)[as.logical(
    str_detect(colnames(ADAMSA_assessment), "BWC"))] 

for(var in recode_incorrect){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] == 6), var] <- 0
}

# #Sanity check
# for(var in colnames(ADAMSA_assessment)){
#   print(var)
#   print(table(ADAMSA_assessment[, var], useNA = "ifany"))
# }

#best trial of repeated trials-- need to account for possible 97 missing value
BWC <- ADAMSA_assessment %>% dplyr::select(contains("BWC")) %>% 
  mutate("ANBWC86" = ANBWC861)
for(i in 1:nrow(BWC)){
  if(BWC[i, "ANBWC86"] == 1 | is.na(BWC[i, "ANBWC86"])) 
    next
  else if(BWC[i, "ANBWC86"] == 0 & BWC[i, "ANBWC862"] == 1)
    BWC[i, "ANBWC86"] = 1
}

ADAMSA_assessment %<>% mutate("BWC" = BWC$ANBWC86)

#Based on data check, I can define the best of the three trials this way
#There are no NAs and there's no 97 after completed prior trials
ADAMSA_assessment %<>% 
  mutate("ANIMMCR" = apply(ADAMSA_assessment %>% 
                             dplyr::select(contains("ANIMMCR")), 1, 
                           max, na.rm = TRUE)) 

# #Sanity check
# View(ADAMSA_assessment %>% dplyr::select(contains("BWC")))
# View(ADAMSA_assessment %>% dplyr::select(contains("ANIMMCR")))

#Drop item-level for these variables
ADAMSA_assessment %<>% dplyr::select(-c("ANBWC861", "ANBWC862", 
                                        paste0("ANIMMCR", seq(1, 3, by = 1))))
#TICS short form total
#Sum Scissor, Cactus, President
TICSTOT <- ADAMSA_assessment %>% 
  dplyr::select("ANSCISOR", "ANCACTUS", "ANPRES")
#Recode 98 = don't know or 99 = refused as 0
TICSTOT[TICSTOT == 98] <- 0
TICSTOT[TICSTOT == 99] <- 0
TICSTOT[, "TOTAL"] = rowSums(TICSTOT)
TICSTOT[TICSTOT$TOTAL > 3, "TOTAL"] <- 97

ADAMSA_assessment %<>% 
  mutate("TICSTOT" = TICSTOT$TOTAL)

# #Sanity Check
# View(ADAMSA_assessment %>%
#        dplyr::select("ANSCISOR", "ANCACTUS", "ANPRES", "TICSTOT"))

#Drop item-level for these variables
ADAMSA_assessment %<>% dplyr::select(-c("ANSCISOR", "ANCACTUS", "ANPRES"))

#---- format data: sociodemographics ----
ADAMS_tracker %<>% 
  dplyr::select("HHIDPN", "AASSESS", "AAGE", "GENDER", "ETHNIC", "EDYRS") %>% 
  mutate("female" = ifelse(GENDER == 1, 0, 1), 
         "ethnic_cat" = case_when(ETHNIC == 1 ~ "Non-hispanic White", 
                                  ETHNIC == 2 ~ "Non-hispanic Black", 
                                  ETHNIC == 3 ~ "Hispanic")) %>% 
  #drop ADAMS variables
  dplyr::select(-c("GENDER", "ETHNIC"))

#---- format data: dem dx ----
ADAMSA_demdx %<>% 
  dplyr::select("HHIDPN", "ADFDX1") %>% 
  mutate("dem_dx" = 
           case_when(ADFDX1 %in% c(1, 2) ~ "Probable/Possible AD", 
                     ADFDX1 %in% c(3, 4) ~ 
                       "Probable/Possible Vascular Dementia", 
                     ADFDX1 %in% 
                       c(5, 8, 14, 23, 24, 25, 26, 27, 21, 28, 29, 30) ~ 
                       "Other",
                     ADFDX1 %in% c(18) ~ "Probable Dementia",
                     ADFDX1 %in% c(10, 13, 15) ~ "Dementia", 
                     ADFDX1 %in% c(20, 22) ~ "MCI", 
                     ADFDX1 == 31 ~ "Normal"), 
         "impaired" = 
           ifelse(dem_dx %in% c("Normal", "Other"), 0, 1)) %>%
  #drop ADAMS variables
  dplyr::select(-c("ADFDX1"))

#---- join all the data ----
ADAMSA <- join_all(list(ADAMS_tracker, ADAMSA_assessment, ADAMSA_demdx), 
                   by = "HHIDPN", type = "left") %>% 
  #filter to those with Wave A assessment
  filter(AASSESS == 1)

#---- look at distributions and models for key covariates ----


#---- generate synthetic data ----

