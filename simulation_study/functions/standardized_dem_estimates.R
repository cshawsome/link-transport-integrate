standardized_dem_estimates <- function(synthetic_data, standard_data){
  
  #---- create age strata ----
  synthetic_data %<>% 
    mutate("age_cat" = 
             cut(age, 
                 breaks = c(70, 75, 80, 85, 90, max(synthetic_data$age)), 
                 include.lowest = TRUE, right = FALSE))
  
  #rename highest strata to match standard data
  synthetic_data %<>% mutate_at("age_cat", as.character)
  synthetic_data[synthetic_data$age_cat == 
                   paste0("[90,", floor(max(synthetic_data$age)),"]"), 
                 "age_cat"] <- "[90,107]"
  
  # #Sanity check
  # table(synthetic_data$age_cat, useNA = "ifany")
  
  #---- standardization tables ----
  agesex_totals <- 
    standard_data %>% group_by(female, age_cat) %>% count() %>% 
    arrange(female) %>% 
    set_colnames(c("female", "age", "standard_count")) %>% 
    dplyr::select("standard_count", everything())
  
  for(race in c("White", "black", "hispanic")){
    assign(paste0(tolower(race), "_dem_risk_table"), 
           synthetic_data %>% filter(!!sym(race) == 1) %>% 
             group_by(female, age_cat) %>% 
             summarise("dem_prop" = mean(Dementia), .groups = "keep") %>% 
             arrange(female) %>% 
             set_colnames(c("female", "age", 
                            paste0(tolower(race), "_dem_risk"))))
  }
  
  #merge tables
  agesex_standardized <- 
    left_join(agesex_totals, white_dem_risk_table, 
              by = c("female", "age")) %>% 
    left_join(., black_dem_risk_table, by = c("female", "age")) %>% 
    left_join(., hispanic_dem_risk_table, by = c("female", "age"))
  
  #expected counts
  for(race in c("white", "black", "hispanic")){
    agesex_standardized %<>% 
      mutate(!!paste0("expected_", race, "_dem_count") := 
               !!sym(paste0(race, "_dem_risk"))*standard_count)
  }
  
  #---- estimates ----
  white_risk <- 
    sum(agesex_standardized$expected_white_dem_count, na.rm = TRUE)/
    nrow(standard_data)
  black_risk <- 
    sum(agesex_standardized$expected_black_dem_count, na.rm = TRUE)/
    nrow(standard_data)
  hispanic_risk <- 
    sum(agesex_standardized$expected_hispanic_dem_count, na.rm = TRUE)/
    nrow(standard_data)
  
  #PR compared to white
  PR_black <- black_risk/white_risk
  PR_hispanic <- hispanic_risk/white_risk
  
  #---- return values ----
  return(c("white_risk" = white_risk, "black_risk" = black_risk, 
           "hispanic_risk" = hispanic_risk, "PR_black" = PR_black, 
           "PR_hispanic" = PR_hispanic))
}
