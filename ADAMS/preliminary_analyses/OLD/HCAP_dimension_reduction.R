#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "stringr", "ggcorrplot", "psych")

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16da/HC16HP_R.da")
dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16sta/HC16HP_R.dct")

HCAP <- read_da_dct(data_path, dict_path, HHIDPN = "TRUE")
cog_test_labels <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/data/", 
                            "cog_test_meaningful_labels.csv"))

#---- variables of interest ----
HCAP_vars <- c("HHIDPN", "MSE", "TICS", "WLIMM1", "WLIMM2", "WLIMM3", "1066", 
               "SCORE")
not_these <- c("TEST", "INTRO", "H1RMSESCORE", "H1RTICSSCORE", "H1R1066SCORE", 
               "SCORESE", "CESDSCORE", "STSCORE")

HCAP_assessment <- HCAP %>% dplyr::select(contains(HCAP_vars)) %>% 
  dplyr::select(-c(contains(not_these), "H1RMSE11T")) #H1RMSE11T1 needs to stay

#---- coding correct/incorrect answers ----
#Need to recode 2 = correct (different level) to 1
recode_correct <- c("H1RMSE17", "H1RMSE21")

#Need to recode 5 = error to 0
recode_incorrect <- 
  colnames(HCAP_assessment)[as.logical(
    str_detect(colnames(HCAP_assessment), "MSE") + 
      str_detect(colnames(HCAP_assessment), "TICS") + 
      str_detect(colnames(HCAP_assessment), "1066"))] 
recode_incorrect <- recode_incorrect[-which(recode_incorrect == "H1RMSE11T1" | 
                                              recode_incorrect == "H1RMSE13")]
for(var in recode_incorrect){
  HCAP_assessment[which(HCAP_assessment[, var] == 5), var] <- 0
}

for(var in recode_correct){
  HCAP_assessment[which(HCAP_assessment[, var] == 2), var] <- 1
}

# #Sanity check
# for(var in colnames(HCAP_assessment)){
#   print(var)
#   print(table(HCAP_assessment[, var], useNA = "ifany"))
# }

#---- code missigness ----
recode_missing_97 <- colnames(HCAP_assessment)[as.logical(
  str_detect(colnames(HCAP_assessment), "MSE") + 
    str_detect(colnames(HCAP_assessment), "H1RCPIMMSCORE") + 
    str_detect(colnames(HCAP_assessment), "H1RCPDELSCORE"))]

recode_missing_8 <- colnames(HCAP_assessment)[as.logical(
  str_detect(colnames(HCAP_assessment), "TICS") + 
    str_detect(colnames(HCAP_assessment), "1066"))]

recode_missing_996 <- c("H1RNSSCORE", "H1RTMASCORE", "H1RTMBSCORE")
  
for(var in recode_missing_97){
  HCAP_assessment[which(HCAP_assessment[, var] >= 97), var] <- NA
}

for(var in recode_missing_8){
  HCAP_assessment[which(HCAP_assessment[, var] >= 8), var] <- NA
}

for(var in recode_missing_996){
  HCAP_assessment[which(HCAP_assessment[, var] >= 996), var] <- NA
}

# #Sanity check
# for(var in colnames(HCAP_assessment)){
#   print(var)
#   print(table(HCAP_assessment[, var], useNA = "ifany"))
# }

#---- sum scores and averages of repeated trials ----
HCAP_assessment %<>% 
  mutate("H1RMSE12SCORE" = rowSums(HCAP_assessment %>% 
                                     dplyr::select(contains("H1RMSE12")), 
                                   na.rm = TRUE), 
         "H1RWLIMMSCORE" = apply(HCAP_assessment %>% 
                                      dplyr::select(contains("H1RWLIMM")), 1, 
                                 max, na.rm = TRUE)) 
HCAP_assessment[is.infinite(HCAP_assessment[, "H1RWLIMMSCORE"]), 
                "H1RWLIMMSCORE"] <- NA

#Sanity check
View(HCAP_assessment %>% dplyr::select(contains("H1RMSE12")))
View(HCAP_assessment %>% dplyr::select(contains("H1RWLIMM")))
table(HCAP_assessment$H1RWLIMMSCORE, useNA = "ifany")

#Drop item-level for these variables
HCAP_assessment %<>% dplyr::select(-c(paste0("H1RMSE12", LETTERS[seq(1, 5)]), 
                                      paste0("H1RWLIMM", seq(1, 3), "SCORE")))

#---- make correlation matrix ----
# #Number of people with complete battery = 2303
# num_complete <- rowSums(is.na(HCAP_assessment[, 2:ncol(HCAP_assessment)]))
# table(num_complete, useNA = "ifany")

HCAP_corr <- cor(HCAP_assessment[, 2:ncol(HCAP_assessment)], 
                 use = "complete.obs")
colnames(HCAP_corr) <- unlist(cog_test_labels[`Variable Name` == 
                                                colnames(HCAP_corr), "Label"])
rownames(HCAP_corr) <- colnames(HCAP_corr)

#Visualize the matrix
HCAP_corr_plot <- ggcorrplot(HCAP_corr, hc.order = TRUE) + 
  theme(axis.text.x = element_text(size = 5), 
        axis.text.y = element_text(size = 5), 
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6))

ggsave(paste0("/Users/CrystalShaw/Box/Dissertation/preliminary_analyses/",
              "figures/HCAP_corr.jpeg"), 
       plot = HCAP_corr_plot, device = "jpeg", width = 5, height = 5, 
       units = "in", dpi = 300)

#---- PCA ----
#Checking appropriateness of method
num_complete <- sum(rowSums(is.na(
  HCAP_assessment[, 2:ncol(HCAP_assessment)])) == 0)

bartlett_p <- cortest.bartlett(HCAP_corr, n = num_complete)

det_corr <- det(HCAP_corr)

#First pass PCA-- don't know how many factors we want yet
HCAP_PCA_1 <- principal(HCAP_corr, nfactors = 10, rotate = "none", 
                        n.obs = num_complete)

#Scree plot
jpeg(paste0("/Users/CrystalShaw/Box/Dissertation/preliminary_analyses/",
     "figures/PCA_scree.jpeg"), width = 7, height = 5, units = "in", 
     res = 300)
plot(HCAP_PCA_1$values, type = "b")
dev.off()

#Second pass PCA-- seems like we want 4 components
HCAP_PCA <- principal(HCAP_corr, nfactors = 4, rotate = "none", 
                      n.obs = num_complete)

#---- Reproduced Correlations ----
HCAP_factor_residuals <- factor.residuals(HCAP_corr, HCAP_PCA$loadings)
HCAP_factor_residuals <- 
  as.matrix(HCAP_factor_residuals[upper.tri(HCAP_factor_residuals)])

#Sanity check-- want these most residuals to be <0.05 (80% of ours are)
plot(HCAP_factor_residuals)
large_resid <- abs(HCAP_factor_residuals) > 0.05
sum(large_resid)/nrow(HCAP_factor_residuals)

#---- PCA Rotations ----
HCAP_PCA_varimax <- principal(HCAP_corr, nfactors = 32, rotate = "varimax", 
                              n.obs = num_complete)
print.psych(HCAP_PCA_varimax, cut = 0.03, sort = TRUE)

PCA_loadings_32 <- as.data.frame(unclass(HCAP_PCA_varimax$loadings)) %>% 
  rownames_to_column(var = "test_item") 

PCA_loadings_32 %<>% 
  mutate("Factor_index" = 
           unlist(apply(PCA_loadings_32[, 2:ncol(PCA_loadings_32)], 
                              1, function(x) which(abs(x) == max(abs(x)))))) %>% 
  mutate("Factor" = colnames(PCA_loadings_32)[`Factor_index` + 1])

#Save matrix of loadings
write_csv(PCA_loadings_32, 
          paste0("/Users/CrystalShaw/Box/Dissertation/preliminary_analyses/",
                 "tables/PCA_loadings_32.csv"))




