
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
  geom_hline(data = truth, aes(yintercept = prop)) + 
  geom_point(aes(y = mean), size = 1) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0)) + 
  coord_flip() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0.15, 0.4, by = 0.05), limits = c(0.14, 0.38)) + 
  labs(x = "HRS sample size", y = "Impairment class proportion") + 
  facet_wrap(~class, nrow = 1)

ggsave(filename = 
         paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                "figure5a_mean_CI_impairment_class.pdf"), 
       dpi = 300, width = 7, height = 2, units = "in")

# Figure 5b ----

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
  geom_hline(yintercept = 95, linetype = "dotted", color = "grey70") + 
  geom_point() + 
  geom_line(aes(linetype = class)) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(90, 100, by = 2), 
                     limits = c(90, 100), expand = c(0, 0)) + 
  theme(legend.position = "bottom") + 
  labs(x = "HRS Sample Size", y = "95% interval coverage", 
       linetype = "Impairment Class")

ggsave(filename = 
         paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                "figure5b_impairment_class_coverage.pdf"), 
       dpi = 300, width = 6, height = 4, units = "in")

# Figure 6a ----

truth <- results %>% ungroup() %>% 
  select(HRS_sample_size, starts_with("true_dem_prev_")) %>% 
  distinct() %>% 
  pivot_longer(
    cols = -HRS_sample_size, 
    names_to = "race", 
    names_prefix = "true_dem_prev_",
    values_to = "value"
  ) %>% 
  mutate(race = str_to_sentence(race) %>% factor(levels = c("White", "Black", "Hispanic")))
  
results %>% 
  ungroup() %>% 
  select(HRS_sample_size, contains("BLMM")) %>% 
  select(HRS_sample_size, contains("_dem_prev")) %>% 
  group_by(HRS_sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  ungroup() %>% 
  pivot_longer(
    cols = -HRS_sample_size, 
    names_to = c("type", "race"),
    names_pattern = "(.*)_dem_prev_(.*)_BLMM",
    values_to = "value"
  ) %>% 
  pivot_wider(
    names_from = type,
    values_from = value
  ) %>%
  mutate(
    race = str_to_sentence(race) %>% factor(levels = c("White", "Black", "Hispanic")), 
    HRS_sample_size = factor(HRS_sample_size)
  ) %>% 
  ggplot(aes(x = HRS_sample_size)) + 
  geom_hline(data = truth, aes(yintercept = value)) +
  geom_point(aes(y = mean), size = 1) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0)) + 
  coord_flip() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0, 0.45, by = 0.1), limits = c(0, 0.45)) +
  labs(x = "HRS sample size", y = "Prevalence of dementia") + 
  facet_wrap(~race, nrow = 1)

ggsave(filename = paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                         "figure6a_mean_CI_dem_prev_by_race.pdf"), 
       dpi = 300, width = 7, height = 2, units = "in")

# Figure 6b ----
truth <- results %>% ungroup() %>% 
  select(HRS_sample_size, starts_with("true_PR")) %>% 
  distinct() %>% 
  pivot_longer(
    cols = c(true_PR_black, true_PR_hispanic), 
    names_to = "race", 
    names_prefix = "true_PR_",
    values_to = "value"
  )

results %>% ungroup() %>% 
  select(HRS_sample_size, contains("BLMM")) %>% 
  select(HRS_sample_size, contains("_PR_"), -contains("coverage")) %>% 
  group_by(HRS_sample_size) %>% 
  summarise(across(everything(), mean)) %>%
  ungroup() %>% 
  pivot_longer(
    cols = -HRS_sample_size,
    names_to = c("type", "race"),
    names_pattern = "(.*)_PR_(.*)_BLMM",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  ) %>%
  mutate(HRS_sample_size = factor(HRS_sample_size)) %>%
  ggplot(aes(x = HRS_sample_size)) + 
  geom_point(aes(y = mean), size = 1) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey70") + 
  geom_hline(data = truth, aes(yintercept = value)) +
  coord_flip() + 
  theme_bw() + 
  # scale_y_continuous(breaks = seq(0.15, 0.4, by = 0.05), limits = c(0.14, 0.38)) + 
  labs(x = "HRS sample size", y = "Prevalence Ratio (PR)") + 
  facet_wrap(~race, nrow = 1)

ggsave(filename = paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                         "figureXXb_mean_CI_PR.pdf"), 
       dpi = 300, width = 7, height = 2, units = "in")

# Figure 6c ----
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
  geom_point(aes(y = mean), size = 1) + 
  geom_errorbar(aes(ymin = LCI, ymax = UCI, width = 0)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") + 
  geom_hline(data = truth, aes(yintercept = value)) +
  coord_flip() + 
  theme_bw() + 
  # scale_y_continuous(breaks = seq(0.15, 0.4, by = 0.05), limits = c(0.14, 0.38)) + 
  labs(x = "HRS sample size", y = "Prevalence Difference (PD)") + 
  facet_wrap(~race, nrow = 1)

ggsave(filename = paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                         "figureXXc_mean_CI_PD.pdf"), 
       dpi = 300, width = 7, height = 2, units = "in")


# Web figure 3 ----

## 3a percent bias ----
# copied from 7a_make_figures_OLD.R - Figure XX: percent bias overall + race-stratified
by_race_cols <- expand_grid(c("mean", "true"), 
                            c("Unimpaired", "MCI", "Dementia", "Other"), 
                            c("white", "black", "hispanic")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

overall_cols <- paste0(
  expand_grid(c("mean", "true"), 
              c("Unimpaired", "MCI", "Dementia", "Other")) %>% 
    unite("names", everything(), sep = "_") %>% 
    unlist() %>% unname(), 
  "_overall")

cols <- c(overall_cols, by_race_cols)

plot_data <- results %>% ungroup() %>%
  dplyr::select("HRS_sample_size", "HCAP_sample_size",
                #"_overall" suffix is not in the results data
                all_of(str_remove(cols, "_overall"))) %>% 
  #rename the overall cols
  set_colnames(c("HRS_sample_size", "HCAP_sample_size", cols))

for(race in c("white", "black", "hispanic", "overall")){
  plot_data %<>% mutate(!!paste0("mean_total_", race) := 
                          !!sym(paste0("mean_Unimpaired_", race)) + 
                          !!sym(paste0("mean_MCI_", race)) + 
                          !!sym(paste0("mean_Dementia_", race)) + 
                          !!sym(paste0("mean_Other_", race)))
  
  plot_data %<>% mutate(!!paste0("true_total_", race) := 
                          !!sym(paste0("true_Unimpaired_", race)) + 
                          !!sym(paste0("true_MCI_", race)) + 
                          !!sym(paste0("true_Dementia_", race)) + 
                          !!sym(paste0("true_Other_", race)))
}

for(race in c("white", "black", "hispanic", "overall")){
  for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
    for(measure in c("mean", "true")){
      plot_data %<>% 
        mutate(!!paste(measure, class, race, sep = "_") := 
                 !!sym(paste(measure, class, race, sep = "_"))/
                 !!sym(paste(measure, "total", race, sep = "_")))
    }
  }  
}

plot_data %<>%
  pivot_longer(c(all_of(cols)),
               names_to = c(".value", "class", "race"), 
               names_sep = "_") %>%
  mutate("error" = mean - true) %>%
  mutate("percent_error" = error/true*100) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>%
  mutate_at("class", function(x) 
    factor(x, levels = c("Unimpaired", "MCI", "Dementia", "Other"))) %>% 
  mutate_at("race", function(x) 
    factor(x, levels = c("Overall", "White", "Black", "Hispanic"))) %>% 
  group_by(HRS_sample_size, class, race) %>% 
  summarize_at(c("error", "percent_error", "squared_error"), mean) %>% 
  rename(c("bias" = "error", "percent_bias" = "percent_error")) %>% 
  mutate("RMSE" = sqrt(squared_error)) 

ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = percent_bias, group = class)) + 
  geom_line(aes(linetype = class)) + geom_point(size = 1) + 
  geom_hline(yintercept = 0, lty = "dotted", color = "grey70") + theme_bw() + 
  ylab("Percent bias") + 
  facet_grid(cols = vars(race)) + 
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_shape_manual(values = rep(c(19, 15, 17, 18), 3)) +
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  scale_y_continuous(limits = c(-30, 30), breaks = seq(-30, 30, by = 10)) +
  labs(linetype = "Impairment Class", shape = "Impairment Class") +
  theme(legend.position = "bottom")

ggsave(filename = paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                         "webfigure3a_impairment_class_percent_bias.pdf"), 
       dpi = 300, width = 7, height = 2.5, units = "in")

## 3b RMSE ----
ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = RMSE, group = class)) + 
  geom_line(aes(linetype = class)) + geom_point(size = 1) + 
  theme_bw() + ylab("RMSE") + 
  facet_grid(cols = vars(race)) + 
  scale_linetype_manual(values = c(1, 2, 3, 4)) +
  scale_shape_manual(values = rep(c(19, 15, 17, 18), 3)) +
  scale_x_discrete(name = "HRS Sample Size", 
                   breaks = unique(plot_data$HRS_sample_size)) + 
  labs(linetype = "Impairment Class", shape = "Impairment Class") +
  theme(legend.position = "bottom")

ggsave(filename = paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                         "webfigure3b_impairment_class_RMSE.pdf"), 
       dpi = 300, width = 7, height = 2.5, units = "in")

# Web figure 4 ----

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
  geom_hline(yintercept = 95, linetype = "dotted", color = "grey70") + 
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(80, 100, by = 5), 
                     limits = c(80, 100), expand = c(0, 0.06)) + 
  facet_wrap(~race, nrow = 1) + 
  labs(x = "HRS Sample Size", y = "95% interval coverage")

ggsave(filename = 
         paste0(path_to_box, "papers/paper1_bayes_framework/figures/", 
                "webfigure3a_coverage.pdf"), 
       dpi = 300, width = 6, height = 2, units = "in")


