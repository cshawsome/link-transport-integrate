#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "tidyr")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **clean HRS ----
HRS <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv")) %>% 
  dplyr::select(-one_of("GENDER", "RACE", "HISPANIC", "PIWTYPE", "PJ005M1", 
                        "r13drinkd", "r13drinkn", "r13pstmem", "GENDER_label", 
                        "RACE_label", "RACE_White", "RACE_Black", "RACE_Other", 
                        "HISPANIC_indicator", "ETHNIC_label", "PJ005M1_label", 
                        "PJ005M1_collapsed_label", "drinks_per_week", 
                        "drink_cat", "drink_cat_label"))

#---- **analytic HCAP ----
#need HRS bins to match these
HCAP <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv")) 

#---- derive variable bins ----
HRS %<>% 
  mutate(
    #---- **age ----
    "age_cat" = cut(age, breaks = c(70, 85, max(HRS$age)), 
                    include.lowest = TRUE, right = FALSE), 
    #---- **education ----
    "edyrs_cat" = cut(edyrs, 
                      breaks = c(0, 11, 12, max(HRS$edyrs)),
                      include.lowest = TRUE, right = TRUE), 
    #---- **immediate word recall ----
    "immrc_cat" = cut(immrc, breaks = c(0, 6, 8, max(HRS$immrc, na.rm = TRUE)), 
                      include.lowest = TRUE, right = FALSE),
    #---- **delayed word recall ----
    "delrc_cat" = cut(delrc, breaks = c(0, 5, 7, max(HRS$delrc, na.rm = TRUE)), 
                      include.lowest = TRUE, right = FALSE),
    #---- **serial 7s ----
    "ser7_cat" = cut(ser7, breaks = c(0, 4, max(HRS$ser7, na.rm = TRUE)), 
                     include.lowest = TRUE, right = TRUE), 
    #---- **hrs cognition ----
    "hrs_cog_cat" = cut(hrs_cog, 
                        breaks = c(0, 17, 21, 23, 25, 
                                   max(HRS$hrs_cog, na.rm = TRUE)), 
                        include.lowest = TRUE, right = TRUE), 
    #---- **imputed memory ----
    "memimp16_cat" = cut(memimp16, 
                         breaks = c(-2.01, 0.222, 0.697, 1.03, 
                                    max(HRS$memimp16, na.rm = TRUE)), 
                         include.lowest = TRUE, right = TRUE))

#Sanity check bins: these need to match HCAP bins
check_vars <- c("age_cat", "edyrs_cat", "immrc_cat", "delrc_cat", "ser7_cat",
                "hrs_cog_cat", "memimp16_cat")

for(var in check_vars){
  print(table(HCAP[, var]))
  print(table(HRS[, var]))
}

#---- **correct some labels ----
HRS %<>% mutate_at(c("age_cat", "hrs_cog_cat", "memimp16_cat"), as.character)

HRS[which(HRS$age_cat == "[85,107]"), "age_cat"] <- "[85,103]"
HRS[which(HRS$hrs_cog_cat == "(25,35]"), "hrs_cog_cat"] <- "(25,33]"
HRS[which(HRS$memimp16_cat == "(1.03,2.1]"), "memimp16_cat"] <- "(1.03,2.02]"

# #Sanity check bins: do these match now?
# check_vars <- c("age_cat", "hrs_cog_cat", "memimp16_cat")
# 
# for(var in check_vars){
#   print(table(HCAP[, var]))
#   print(table(HRS[, var]))
# }
