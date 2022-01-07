#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "mice")

options(scipen = 999)

#---- read in data ----
ADAMS_subset <- 
  read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                  "data/ADAMS/cleaned/ADAMS_subset_mixed.csv"), 
           col_types = cols("AYEAR" = col_character(), 
                            "Astroke" = col_character(), 
                            "Ahibpe" = col_character(), 
                            "Adiabe" = col_character(), 
                            "Ahearte" = col_character(), 
                            "Acancre" = col_character(), 
                            "Asmoken" = col_character(), 
                            "White" = col_character(), 
                            "Black" = col_character(), 
                            "Hispanic" = col_character())) %>% 
  mutate_if(is.character, as.factor)

#Set factor levels
sapply(ADAMS_subset, class)
ADAMS_subset %<>% 
  mutate("GENDER_label" = fct_relevel(GENDER_label, "Male"), 
         "ETHNIC_label" = fct_relevel(ETHNIC_label, "White"), 
         "EDYRScat_label" = fct_relevel(EDYRScat_label, "No School Completed", 
                                        "1st-8th Grade", "Some High School", 
                                        "High School Diploma", "Some College"), 
         "AAMARRD_label" = fct_relevel(AAMARRD_label, "Married/partnered"), 
         "AACURRWK_label" = fct_relevel(AACURRWK_label, "Retired"), 
         "DRINKcat_label" = fct_relevel(DRINKcat_label, "No Drinking", 
                                        "Moderate Drinking"))

# #Sanity check
# levels(ADAMS_subset$GENDER_label)
# levels(ADAMS_subset$ETHNIC_label)
# levels(ADAMS_subset$EDYRScat_label)
# levels(ADAMS_subset$AAMARRD_label)
# levels(ADAMS_subset$AACURRWK_label)
# levels(ADAMS_subset$DRINKcat_label)

#---- Z-score ----
ADAMS_subset %<>% mutate_if(is.numeric, scale)