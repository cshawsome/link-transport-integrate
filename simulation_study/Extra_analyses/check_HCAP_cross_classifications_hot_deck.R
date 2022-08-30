#---- Package loading ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse")

options(scipen = 999)

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HCAP data ----
HCAP_analytic <- 
  read_csv(paste0(path_to_box, "data/HCAP/cleaned/HCAP_analytic.csv")) 

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) %>% 
  filter(HCAP %in% colnames(HCAP_analytic)) 

#label data
HCAP_analytic %<>% 
  rename_at(vars(variable_labels$HCAP), ~ variable_labels$data_label) 

#---- create categorical vars ----
HCAP_analytic %<>% 
  mutate(
    #---- **age ----
    "age_cat" = cut(age, breaks = c(min(HCAP_analytic$age), 75, 80, 85, 90, 
                                    max(HCAP_analytic$age)), 
                    include.lowest = TRUE, right = FALSE), 
    #---- **education ----
    "edyrs_cat" = cut(edyrs, breaks = c(0, 0.9, 11, 12, 15, 16, 
                                        max(HCAP_analytic$edyrs)),
                      labels = c("none", "less than HS", "HS", "some college", 
                                 "college", "graduate studies"),
                      include.lowest = TRUE), right = TRUE, 
    #---- **immediate word recall ----
    "immrc_cat" = cut_number(immrc, n = 5), 
    #---- **serial 7s ----
    "ser7_cat" = cut_number(ser7, n = 2), 
    #---- **hrs cognition ----
    "hrs_cog_cat" = cut_number(hrs_cog, n = 5), 
    #---- **adl ----
    "adl_cat" = ifelse(adl == 0, "none", "any"), 
    #---- **iadl ----
    "iadl_cat" = ifelse(iadl == 0, "none", "any"), 
    #---- **bmi ----
    "bmi_cat" = cut_number(bmi, n = 5), 
    #---- **mmse_norm ----
    "mmse_norm_cat" = cut_number(mmse_norm, n = 5),
    #---- **delayed word recall ----
    "delrc_cat" = cut_number(delrc, n = 5),
    #---- **animal naming ----
    "animal_naming_cat" = cut_number(animal_naming, n = 5),
    #---- **word recall yes ----
    "wrc_yes_cat" = case_when(wrc_yes %in% seq(0, 6) ~ "some", 
                              wrc_yes %in% seq(7, 10) ~ "most"), 
    #---- **word recall no ----
    "wrc_no_cat" = case_when(wrc_no %in% seq(0, 6) ~ "some", 
                             wrc_no %in% seq(7, 10) ~ "most"), 
    #---- **immediate story recall ----
    "imm_story_cat" = cut_number(imm_story, n = 5), 
    #---- **delayed story recall ----
    "del_story_cat" = cut_number(del_story, n = 5),
    #---- **immediate constructional praxis ----
    "imm_cp_cat" = cut_number(imm_cp, n = 4), 
    #---- **delayed constructional praxis ----
    "del_cp_cat" = cut_number(del_cp, n = 4),
    #---- **trails A ----
    "trailsA_cat" = cut_number(trailsA, n = 5))

#---- cross tabulate ----
var_select = c("married_partnered", "age_cat", "smoken", "hibpe", "diabe", 
               "hearte", "stroke", "scissor", "cactus", "pres", "female",
               "black", "hispanic", "retired", "not_working", 
               "subj_cog_better", "subj_cog_worse", "moderate_drinking", 
               "heavy_drinking", "edyrs_cat", "immrc_cat", "ser7_cat", "bwc20", 
               "hrs_cog_cat", "adl_cat", "iadl_cat")

var_select_test <- c("married_partnered", "age_cat", "smoken", "hibpe", "diabe", 
                     "hearte", "stroke", "scissor", "cactus", "pres", "female",
                     "black", "hispanic", "retired", "not_working")

cross_tabs <- 
  as.data.frame(table(subset(HCAP_analytic, select = var_select_test)))

