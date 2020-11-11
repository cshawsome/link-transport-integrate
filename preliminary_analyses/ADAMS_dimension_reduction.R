#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "stringr", "ggcorrplot", "psych", "broom", 
       "openxlsx")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import neuropsych data ----
neuropsych_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                 "ADAMS/adams1a/adams1ada/ADAMS1AN_R.da")
neuropsych_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                 "ADAMS/adams1a/adams1asta/ADAMS1AN_R.dct")

ADAMS_neuropsych_A <- read_da_dct(neuropsych_data_path_A, 
                                  neuropsych_dict_path_A, HHIDPN = "TRUE")

cog_test_labels <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                   "cog_test_meaningful_labels.csv"))
test_labels <- unlist(cog_test_labels$Label)
names(test_labels) <- unlist(cog_test_labels$`Variable Name`)

#---- variables of interest ----
ADAMS_vars <- c("HHIDPN",
                #item-level
                paste0("MSE", seq(1, 22, by = 1)), 
                "SCISOR", "CACTUS", "ANPRES", 
                #total scores
                "ANMSETOT", "SER7T", "ANAFTOT", "CPTOT", "BWC86", "IMMCR", 
                "DELCOR", "RECYES", "RECNO", "WM1TOT", "WM2TOT", "MASEC", 
                "MBSEC", "ANSDMTOT")

not_these <- c("ANMSE11T", "ANRCPTOT")

ADAMSA_assessment <- ADAMS_neuropsych_A %>% 
  dplyr::select(contains(ADAMS_vars)) %>% 
  dplyr::select(-c(contains(not_these))) 

#---- coding correct/incorrect answers ----
#Need to recode 2 = correct (different level) to 1
recode_correct <- c("ANMSE16", "ANMSE17", "ANMSE19", "ANMSE21")

#Need to recode 6 = start over to 0
recode_incorrect <- 
  colnames(ADAMSA_assessment)[as.logical(
    str_detect(colnames(ADAMSA_assessment), "BWC"))] 

for(var in recode_incorrect){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] == 6), var] <- 0
}

for(var in recode_correct){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] == 2), var] <- 1
}

# #Sanity check
# for(var in colnames(ADAMSA_assessment)){
#   print(var)
#   print(table(ADAMSA_assessment[, var], useNA = "ifany"))
# }

#---- code missigness ----
recode_missing_97 <- colnames(ADAMSA_assessment)[as.logical(
  str_detect(colnames(ADAMSA_assessment), "MSE") + 
    str_detect(colnames(ADAMSA_assessment), "BWC") + 
    str_detect(colnames(ADAMSA_assessment), "SER7") + 
    str_detect(colnames(ADAMSA_assessment), "SCISOR") + 
    str_detect(colnames(ADAMSA_assessment), "CACTUS") + 
    str_detect(colnames(ADAMSA_assessment), "PRES") + 
    str_detect(colnames(ADAMSA_assessment), "TOT") + 
    str_detect(colnames(ADAMSA_assessment), "IMMCR") + 
    str_detect(colnames(ADAMSA_assessment), "DELCOR") + 
    str_detect(colnames(ADAMSA_assessment), "REC"))]

recode_missing_995 <- c("ANTMASEC", "ANTMBSEC")

for(var in recode_missing_97){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] >= 97), var] <- NA
}

for(var in recode_missing_995){
  ADAMSA_assessment[which(ADAMSA_assessment[, var] >= 995), var] <- NA
}

# #Sanity check
# for(var in colnames(ADAMSA_assessment)){
#   print(var)
#   print(table(ADAMSA_assessment[, var], useNA = "ifany"))
# }

#---- best trial of repeated trials ----
ADAMSA_assessment %<>% 
  mutate("ANBWC86" = apply(ADAMSA_assessment %>% 
                             dplyr::select(contains("BWC")), 1, 
                           max, na.rm = TRUE), 
         "ANIMMCR" = apply(ADAMSA_assessment %>% 
                             dplyr::select(contains("ANIMMCR")), 1, 
                           max, na.rm = TRUE)) 

ADAMSA_assessment[is.infinite(ADAMSA_assessment[, "ANBWC86"]), 
                  "ANBWC86"] <- NA
ADAMSA_assessment[is.infinite(ADAMSA_assessment[, "ANIMMCR"]), 
                  "ANIMMCR"] <- NA

# #Sanity check
# View(ADAMSA_assessment %>% dplyr::select(contains("BWC")))
# View(ADAMSA_assessment %>% dplyr::select(contains("ANIMMCR")))

#Drop item-level for these variables
ADAMSA_assessment %<>% dplyr::select(-c("ANBWC861", "ANBWC862", 
                                        paste0("ANIMMCR", seq(1, 3, by = 1))))
#---- Create TICS short form total ----
#Sum Scissor, Cactus, President
ADAMSA_assessment %<>% 
  mutate("TICSTOT" = apply(ADAMSA_assessment %>% 
                             dplyr::select("ANSCISOR", "ANCACTUS", "ANPRES"), 1, 
                           function(x) ifelse(sum(is.na(x)) == 3, NA, 
                                              sum(x, na.rm = TRUE))))

# #Sanity Check
# View(ADAMSA_assessment %>% 
#        dplyr::select("ANSCISOR", "ANCACTUS", "ANPRES", "TICSTOT"))

#---- make correlation matrix ----
total_scores <- c("ANMSETOT", "ANSER7T", "TICSTOT", "ANAFTOT", "ANCPTOT", 
                  "ANDCPTOT", "ANDELCOR", "ANRECYES", "ANRECNO", "ANWM1TOT", 
                  "ANWM2TOT", "ANTMASEC", "ANTMBSEC", "ANSDMTOT", "ANBWC86", 
                  "ANIMMCR")

# #Number of people with complete battery (item-level) = 312
# #Number of people with complete battery (total-scores) = 342
# num_complete <- rowSums(is.na(ADAMSA_assessment %>% 
#                                 dplyr::select(all_of(total_scores))))
# table(num_complete, useNA = "ifany")

# #ANMSE 16 and ANMSE 17 have no variance, so need to drop these from the item-
# #level analyses
# apply(na.omit(ADAMSA_assessment[, 2:ncol(ADAMSA_assessment)]), 2, 
#       function(x) sd(x, na.rm = TRUE))

ADAMSA_corr <- 
  cor(ADAMSA_assessment[, which(colnames(ADAMSA_assessment) %in% 
                                  tidyselect::all_of(total_scores))], 
      use = "complete.obs")

colnames(ADAMSA_corr) <- test_labels[colnames(ADAMSA_corr)]
rownames(ADAMSA_corr) <- colnames(ADAMSA_corr)

#Visualize the matrix
ADAMSA_corr_plot <- ggcorrplot(ADAMSA_corr, hc.order = TRUE) +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))

ggsave(paste0("/Users/CrystalShaw/Box/Dissertation/preliminary_analyses/", 
              "ADAMS_PCA/figures/ADAMSA_total_scores_corr.jpeg"),
       plot = ADAMSA_corr_plot, device = "jpeg", width = 5, height = 5,
       units = "in", dpi = 300)

#---- PCA ----
#Checking appropriateness of method
num_complete <- sum(rowSums(is.na(ADAMSA_assessment %>% 
                                dplyr::select(all_of(total_scores)))) == 0)

bartlett_p <- cortest.bartlett(ADAMSA_corr, n = num_complete)

det_corr <- det(ADAMSA_corr)

#First pass PCA-- don't know how many factors we want yet
ADAMSA_PCA_1 <- principal(ADAMSA_corr, nfactors = 10, rotate = "varimax", 
                          n.obs = num_complete)

#Scree plot
jpeg(paste0("/Users/CrystalShaw/Box/Dissertation/preliminary_analyses/", 
            "ADAMS_PCA/figures/PCA_scree.jpeg"), width = 7, height = 5, 
     units = "in", res = 300)
plot(ADAMSA_PCA_1$values, type = "b")
dev.off()

#Second pass PCA-- seems like we want 10 components (88% variance explained)
ADAMSA_PCA <- principal(ADAMSA_corr, nfactors = 10, rotate = "varimax", 
                        n.obs = num_complete)

#---- Reproduced Correlations ----
ADAMSA_factor_residuals <- factor.residuals(ADAMSA_corr, ADAMSA_PCA$loadings)
ADAMSA_factor_residuals <- 
  as.matrix(ADAMSA_factor_residuals[upper.tri(ADAMSA_factor_residuals)])

#Sanity check-- want most residuals to be <0.05 (85% of ours are)
plot(ADAMSA_factor_residuals)
large_resid <- abs(ADAMSA_factor_residuals) > 0.05
sum(large_resid)/nrow(ADAMSA_factor_residuals)

#---- Factor loadings ----
print.psych(ADAMSA_PCA, cut = 0.03, sort = TRUE)

PCA_loadings_10 <- as.data.frame(unclass(ADAMSA_PCA$loadings)) %>% 
  rownames_to_column(var = "test_item") 

PCA_loadings_10 %<>% 
  mutate("Factor_index" = 
           unlist(apply(PCA_loadings_10[, 2:ncol(PCA_loadings_10)], 
                        1, function(x) which(abs(x) == max(abs(x)))))) %>% 
  mutate("Factor" = colnames(PCA_loadings_10)[`Factor_index` + 1])

#View(PCA_loadings_10 %>% dplyr::select("test_item", "Factor"))

#Save matrix of loadings
write_csv(PCA_loadings_10, 
          paste0("/Users/CrystalShaw/Box/Dissertation/preliminary_analyses/",
                 "ADAMS_PCA/tables/PCA_loadings_10.csv"))

#---- create factors ----
#Complete case data
ADAMSA_test_data <- ADAMSA_assessment %>% 
  dplyr::select("HHIDPN", all_of(total_scores)) 

loadings <- as.data.frame(unclass(ADAMSA_PCA$loadings))
factors <- ADAMSA_test_data %>% dplyr::select(all_of(total_scores)) %>% 
  set_colnames(test_labels[total_scores]) %>% 
  #ensure ordering
  dplyr::select(all_of(rownames(loadings)))

factors <- as.matrix(factors) %*% as.matrix(loadings)

ADAMSA_assessment %<>% cbind(factors)

#---- demdx data ----
demdx_data_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "ADAMS/adams1a/adams1ada/ADAMS1AD_R.da")
demdx_dict_path_A <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                            "ADAMS/adams1a/adams1asta/ADAMS1AD_R.dct")

ADAMS_demdx_A <- read_da_dct(demdx_data_path_A, demdx_dict_path_A, 
                             HHIDPN = "TRUE") %>% 
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
           ifelse(dem_dx %in% c("Normal", "Other"), 0, 1))

# #Sanity check
# table(ADAMS_demdx_A$dem_dx, ADAMS_demdx_A$impaired, useNA = "ifany")

#join dementia data
ADAMSA_assessment %<>% left_join(., ADAMS_demdx_A, by = "HHIDPN")

# #data check
# #complete data only
# neuropsych_totals <- ADAMSA_assessment %>%
#   dplyr::select(all_of(total_scores), "impaired") %>% na.omit()
# table(neuropsych_totals$impaired, useNA = "ifany")
# #all data
# neuropsych_totals <- ADAMSA_assessment %>%
#   dplyr::select(all_of(total_scores), "impaired") 
# table(neuropsych_totals$impaired, useNA = "ifany")

#---- sociodemographic data ----
ADAMS_tracker_data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                                  "ADAMS/adams1trk/ADAMS1TRK_R.da")
ADAMS_tracker_dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/",
                                  "ADAMS/adams1trk/adams1trksta/ADAMS1TRK_R.dct")

ADAMS_tracker <- read_da_dct(ADAMS_tracker_data_path, ADAMS_tracker_dict_path, 
                             HHIDPN = "TRUE") %>% 
  dplyr::select("HHIDPN", "AAGE", "GENDER", "ETHNIC", "EDYRS") %>% 
  mutate("female" = ifelse(GENDER == 1, 0, 1), 
         "ethnic_cat" = case_when(ETHNIC == 1 ~ "Non-hispanic White", 
                                  ETHNIC == 2 ~ "Non-hispanic Black", 
                                  ETHNIC == 3 ~ "Hispanic"))

# #Sanity check
# table(ADAMS_tracker$GENDER, ADAMS_tracker$female, useNA = "ifany")
# table(ADAMS_tracker$ETHNIC, ADAMS_tracker$ethnic_cat, useNA = "ifany")

#join sociodemographic data
ADAMSA_assessment %<>% left_join(., ADAMS_tracker, by = "HHIDPN")

#---- univariate models ----
#Running these models in the complete case data only leads to different 
#conclusions
sex <- glm(impaired ~ female, family = binomial(link = "logit"), 
           data = ADAMSA_assessment)
tidy(sex, exponentiate = TRUE, conf.int = TRUE)

age <- glm(impaired ~ AAGE, family = binomial(link = "logit"), 
           data = ADAMSA_assessment)
tidy(age, exponentiate = TRUE, conf.int = TRUE)

ethnicity <- glm(impaired ~ as.factor(ethnic_cat), 
                 family = binomial(link = "logit"), 
                 data = ADAMSA_assessment)
tidy(ethnicity, exponentiate = TRUE, conf.int = TRUE)

education <- glm(impaired ~ EDYRS, 
                 family = binomial(link = "logit"), 
                 data = ADAMSA_assessment)
tidy(education, exponentiate = TRUE, conf.int = TRUE)

#Sanity Check-- missingness in sociodemographic data
ADAMSA_assessment %>% dplyr::select("female", "EDYRS", "ethnic_cat", "AAGE") %>% 
  is.na() %>% colSums()

#---- multivariate model ----
model_complete <- glm(impaired ~ female + AAGE + as.factor(ethnic_cat) + EDYRS + 
                        RC1 + RC2 + RC3 + RC4 + RC5 + RC6 + RC7 + RC8 + RC9 + 
                        RC10, family = binomial(link = "logit"), 
                      data = ADAMSA_assessment)

tidy(model_complete, exponentiate = TRUE, conf.int = TRUE)

model.list <- list("Complete Data" = tidy(model_complete, exponentiate = TRUE, 
                                          conf.int = TRUE))
write.xlsx(model.list, 
           file = paste0("/Users/CrystalShaw/Box/Dissertation/",
                         "preliminary_analyses/ADAMS_PCA/tables/", 
                         "dem_outcome_models.xlsx"))

#---- Covariate Distributions ----
#How different are complete cases from those who have at least one missing 
#assessment
ADAMSA_missing_some <- ADAMSA_assessment %<>%






