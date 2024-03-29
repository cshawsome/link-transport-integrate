#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "here", "magrittr", "tidyr")

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **clean HRS ----
HRS <- 
  read_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_clean.csv"))

#---- **analytic HCAP ----
#need HRS bins to match these
HCAP <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic_for_sim.csv"))

# #Sanity check
# colSums(is.na(HCAP))[colSums(is.na(HCAP)) > 0]

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
    "hrs_cog_superpop_cat" = cut(hrs_cog, 
                                 breaks = c(0, 18, 22, 25, 
                                            max(HRS$hrs_cog, na.rm = TRUE)), 
                                 include.lowest = TRUE, right = TRUE))

#Sanity check bins: these need to match HCAP bins
check_vars <- c("age_cat", "edyrs_cat", "immrc_cat", "delrc_cat", "ser7_cat",
                "hrs_cog_superpop_cat")

for(var in check_vars){
  print(table(HCAP[, var]))
  print(table(HRS[, var]))
}

#---- **correct some labels ----
HRS %<>% 
  mutate_at(c("age_cat", "hrs_cog_superpop_cat"), as.character)

HRS[which(HRS$age_cat == "[85,107]"), "age_cat"] <- "[85,103]"
HRS[which(HRS$hrs_cog_superpop_cat == "(25,35]"), "hrs_cog_superpop_cat"] <- 
  "(25,33]"

#Sanity check bins: do these match now?
check_vars <- c("age_cat", "hrs_cog_superpop_cat")

for(var in check_vars){
  print(table(HCAP[, var]))
  print(table(HRS[, var]))
}

#---- save dataset ----
HRS %>% write_csv(paste0(path_to_box, "data/HRS/cleaned/HRS_analytic.csv"))

