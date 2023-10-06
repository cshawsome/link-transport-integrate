
#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "here", "ggbreak")

#---- source functions ----
source(here::here("functions", "read_results.R"))

#---- **read in data ----
path_to_box <- "~/Library/CloudStorage/Box-Box/Dissertation/"

#---- ****results ----
results_paths <- 
  list.files(path = paste0(path_to_box, "analyses/simulation_study/results"), 
             full.names = TRUE, pattern = "*ADAMS_prior.csv")

results <- do.call(rbind, lapply(results_paths, read_results)) %>% 
  group_by(dataset_name) %>% slice_head(n = 1000) %>% 
  mutate("HRS_sample_size" = sample_size) %>% 
  mutate("HCAP_sample_size" = HCAP_prop/100*HRS_sample_size) %>% 
  dplyr::select(-one_of("sample_size"))

#---- ******check number of simulation runs ----
#there should be 1000 runs of each scenario
table(results$dataset_name, useNA = "ifany")

#---- ****superpop ----
superpop_impairment_props <- 
  read_csv(paste0(path_to_box, 
                  "data/superpopulations/impairment_class_props.csv"))

superpop_impairment_props$Group <- 
  factor(superpop_impairment_props$Group, 
         levels = c("Unimpaired", "MCI", "Dementia", "Other"))

#---- **proportions ----
results %<>% 
  mutate(true_Unimpaired_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "Unimpaired", "prop"]), 
         true_MCI_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "MCI", "prop"]), 
         true_Dementia_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "Dementia", "prop"]), 
         true_Other_prop = 
           as.numeric(superpop_impairment_props[
             superpop_impairment_props$Group == "Other", "prop"]))


# Figure 5a ---------------------------------------------------------------

truth <- data.frame(
  "class" = factor(c("Unimpaired", "MCI", "Dementia", "Other"), 
                   levels = c("Unimpaired", "MCI", "Dementia", "Other")), 
  prop = c(
    unique(round(na.omit(results$true_Unimpaired_prop), 2)), 
    unique(round(na.omit(results$true_MCI_prop), 2)), 
    unique(round(na.omit(results$true_Dementia_prop), 2)), 
    unique(round(na.omit(results$true_Other_prop), 2))
  )
)

results %>% ungroup() %>% 
  select(HRS_sample_size, HCAP_sample_size, 
         expand.grid(
           c("mean_", "LCI_", "UCI_"), 
           c("Unimpaired", "MCI", "Dementia", "Other")
         ) %>% mutate(var = paste0(Var1, Var2)) %>% 
           pull(var)
  ) %>% 
  mutate(
    across(-c(HRS_sample_size, HCAP_sample_size), 
           function(x) x / HCAP_sample_size, 
           .names = "prop_{.col}")
  ) %>% 
  select(HRS_sample_size, starts_with("prop")) %>% 
  group_by(HRS_sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  ungroup() %>% 
  pivot_longer(
    cols = starts_with("prop"), 
    names_to = c("type", "class"),
    names_pattern = "prop_(.*)_(.*)",
    values_to = "value"
  ) %>% 
  pivot_wider(
    names_from = type,
    values_from = value
  ) %>%
  mutate(
    class = factor(class, levels = c("Unimpaired", "MCI", "Dementia", "Other")), 
    HRS_sample_size = factor(HRS_sample_size)
  ) %>% 
  ggplot(aes(x = HRS_sample_size)) + 
  geom_point(aes(y = mean)) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0.3)) + 
  geom_hline(data = truth, aes(yintercept = prop)) + 
  coord_flip() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0.15, 0.4, by = 0.05), limits = c(0.14, 0.38)) + 
  labs(x = "HRS sample size", y = "Impairment class proportion") + 
  facet_wrap(~class, nrow = 1)


# Figure 5b

results %>% ungroup() %>% 
  select(HRS_sample_size, ends_with("coverage")) %>% 
  group_by(HRS_sample_size) %>% 
  summarise(across(everything(), function(x) sum(x) / 1000 * 100)) %>% 
  ungroup() %>% 
  pivot_longer(
    cols = ends_with("coverage"),
    names_to = "class", 
    values_to = "coverage"
  ) %>% 
  mutate(
    class = str_remove(class, "_coverage") %>% 
      factor(levels = c("Unimpaired", "MCI", "Dementia", "Other")), 
    HRS_sample_size = factor(HRS_sample_size)
  ) %>% 
  ggplot(aes(x = HRS_sample_size, y = coverage, group = class)) + 
  geom_point() + 
  geom_line(aes(linetype = class)) + 
  geom_hline(yintercept = 95, linetype = "dotted") + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(90, 100, by = 2), limits = c(90, 100)) + 
  theme(legend.position = "bottom") + 
  labs(x = "HRS Sample Size", y = "95% interval coverage", 
       linetype = "Impairment Class")


# Figure 6c
truth <- results %>% ungroup() %>% 
  select(HRS_sample_size, starts_with("true_PD")) %>% 
  distinct() %>% 
  pivot_longer(
    cols = c(true_PD_black, true_PD_hispanic), 
    names_to = "race", 
    names_prefix = "true_PD_",
    values_to = "value"
  )

results %>% ungroup() %>% 
  select(HRS_sample_size, contains("BLMM")) %>% 
  select(HRS_sample_size, contains("PD"), -contains("coverage")) %>% 
  group_by(HRS_sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  ungroup() %>% 
  pivot_longer(
    cols = -HRS_sample_size, 
    names_to = c("type", "race"),
    names_pattern = "(.*)_PD_(.*)_BLMM",
    values_to = "value"
  ) %>% 
  pivot_wider(
    names_from = type,
    values_from = value
  ) %>% 
  mutate(HRS_sample_size = factor(HRS_sample_size)) %>% 
  ggplot(aes(x = HRS_sample_size)) + 
  geom_point(aes(y = mean)) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0.3)) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_hline(data = truth, aes(yintercept = value)) +
  coord_flip() + 
  theme_bw() + 
  # scale_y_continuous(breaks = seq(0.15, 0.4, by = 0.05), limits = c(0.14, 0.38)) + 
  labs(x = "HRS sample size", y = "Prevalence Difference (PD)") + 
  facet_wrap(~race, nrow = 1)

# Web figure 3a

results %>% ungroup() %>% 
  select(
    HRS_sample_size, 
    dem_prev_coverage_black_BLMM, 
    dem_prev_coverage_hispanic_BLMM,
    dem_prev_coverage_white_BLMM
  ) %>% 
  group_by(HRS_sample_size) %>% 
  summarise(across(everything(), function(x) sum(x) / 1000 * 100)) %>% 
  ungroup() %>%
  pivot_longer(
    cols = -HRS_sample_size,
    names_to = "race",
    names_pattern = "dem_prev_coverage_(.*)_BLMM",
    values_to = "dem_prev_coverage"
  ) %>% 
  mutate(
    race = str_to_sentence(race) %>% factor(levels = c("White", "Black", "Hispanic")), 
    HRS_sample_size = factor(HRS_sample_size)
  ) %>% 
  ggplot(aes(x = HRS_sample_size, y = dem_prev_coverage, group = race)) + 
  geom_point() + 
  geom_line() + 
  geom_hline(yintercept = 95, linetype = "dotted") + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(80, 100, by = 5), limits = c(80, 100)) + 
  facet_wrap(~race, nrow = 1) + 
  labs(x = "HRS Sample Size", y = "95% interval coverage")


# Web figure 3b
results %>% ungroup() %>% 
  select(HRS_sample_size, 
         starts_with("bias_")) %>% 
  rename(
    bias_Unimpaired_overall = bias_Unimpaired,
    bias_MCI_overall = bias_MCI, 
    bias_Dementia_overall = bias_Dementia, 
    bias_Other_overall = bias_Other
  ) %>% 
  pivot_longer(
    cols = -HRS_sample_size,
    names_pattern = "bias_(.*)_(.*)",
    names_to = c("class", "race"),
    values_to = "percent_bias"
  ) %>%
  mutate(
    race = str_to_sentence(race) %>% factor(levels = c("Overall", "White", "Black", "Hispanic")), 
    HRS_sample_size = factor(HRS_sample_size)
  ) %>% 
  group_by(HRS_sample_size, class, race) %>% 
  summarise(percent_bias = median(percent_bias)) %>% 
  ggplot(aes(x = HRS_sample_size, y = percent_bias, group = class)) + 
  geom_point() + 
  geom_line(aes(linetype = class)) + 
  facet_wrap(~race, nrow = 1)

# Web Figure 3c
# code from section: Appendix Figure XX: RMSE dementia prevalence 
cols_by_race <- expand_grid(c("mean"), 
                            c("dem_prev"), 
                            c("white", "black", "hispanic"), 
                            c("BLMM", "ModHurd", "Hurd", "LKW")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

truth <- 
  data.frame("Race" = c("white", "black", "hispanic"), 
             "true" = c(
               unique(round(na.omit(results$true_dem_prev_white), 2)), 
               unique(round(na.omit(results$true_dem_prev_black), 2)), 
               unique(round(na.omit(results$true_dem_prev_hispanic), 2))))

plot_data <- results %>% ungroup() %>%
  dplyr::select("HRS_sample_size", "HCAP_sample_size",  all_of(cols_by_race)) %>%
  pivot_longer(all_of(cols_by_race),
               names_to = c(".value", "measure1", "measure2", "Race", "Algorithm"), 
               names_sep = "_") %>% left_join(truth) %>%
  mutate("error" = mean - true) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("Race", function(x) str_to_sentence(x)) %>%
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>%
  group_by(HRS_sample_size, Race, Algorithm) %>% 
  summarize_at(c("error", "squared_error"), mean) %>% 
  rename(c("bias" = "error")) %>% 
  mutate("RMSE" = sqrt(squared_error)) %>%
  mutate_at("Algorithm", function(x) 
    factor(x, levels = c("BLMM", "ModHurd", "Hurd", "LKW")))

plot_data %>% 
  filter(Algorithm == "BLMM") %>% 
  ggplot(aes(x = HRS_sample_size, y = RMSE, group = Race)) + 
  geom_line() +
  geom_point(size = 3) + 
  theme_bw() + 
  ylab("RMSE") + xlab("HRS Sample Size") +
  facet_wrap(~ Race, nrow = 1)
  facet_grid(cols = Race) + ylim(0, 0.25) +
  theme(text = element_text(size = 24)) 
