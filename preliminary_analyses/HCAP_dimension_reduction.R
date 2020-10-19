#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

#---- source scripts ----
source(paste0("/Users/CrystalShaw/Desktop/Git Repos/useful-scripts/R/", 
              "data_read/read_da_dct.R"))

#---- import data ----
data_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16da/HC16HP_R.da")
dict_path <- paste0("/Users/CrystalShaw/Box/Dissertation/data/HCAP/HC16/", 
                    "HC16sta/HC16HP_R.dct")

HCAP <- read_da_dct(data_path, dict_path, HHIDPN = "TRUE")

#---- variables of interest ----
HCAP_vars <- c("HHIDPN", )

