#---- package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "here", "wesanderson", "RColorBrewer")

#---- read in data ----
#---- **ADAMS ----
ADAMS_train <- read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/", 
                               "data/cleaned/ADAMS/ADAMS_train.csv"), 
                        col_types = cols(HHIDPN = col_character())) 

#---- **synthetic ----
for(run in 1:10){
  if(run == 1){
    synthetic_sample <- 
      read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                      "results/ADAMSA/standard_normal/ADAMSA_synthetic_", run, 
                      ".csv")) %>% mutate("sample" = run)
  } else{
    synthetic_sample %<>% 
      rbind(., 
            read_csv(paste0("/Users/CrystalShaw/Box/Dissertation/analyses/", 
                            "results/ADAMSA/standard_normal/ADAMSA_synthetic_", 
                            run, ".csv")) %>% mutate("sample" = run))
  }
}

#---- categorical checks ----
#---- **race/ethnicity x stroke ----
true_counts

#---- continuous checks ----
#---- **density plots ----
#---- **median plots ----
#---- **skew plots ----

#---- impairment classification ----

