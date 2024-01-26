#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse", "magrittr", "wesanderson", "here", "ggbreak", "ggh4x")

#---- source functions ----
source(here::here("functions", "read_results.R"))

#---- **read in data ----
# path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"
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

# #defining column names
# overall_impairment_cols <-
#   expand_grid(c("mean", "LCI", "UCI"),
#               c("Unimpaired", "MCI", "Dementia", "Other")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# by_race_impairment_cols <-
#   expand_grid(c("mean", "LCI", "UCI"),
#               c("Unimpaired", "MCI", "Dementia", "Other"),
#               c("white", "black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# dem_prev_cols <-
#   expand_grid(c("mean"), "dem_prev", c("white", "black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# PR_cols <-
#   expand_grid(c("mean"), "log_PR", c("black", "hispanic")) %>%
#   unite("names", everything(), sep = "_") %>% unlist() %>% unname()
# 
# pt_est_cols_cores <- c(overall_impairment_cols, by_race_impairment_cols, 
#                        dem_prev_cols, PR_cols)
# 
# #selecting correct columns
# calc_props <- 
#   pt_est_cols_cores[-c(which(str_detect(pt_est_cols_cores, "dem_prev")), 
#                        which(str_detect(pt_est_cols_cores, "log_PR")))]
# 
# results[, paste0(calc_props, "_prop")] <- 
#   results[, calc_props]/results$HCAP_sample_size

#---- Figure XX: dementia prevalence ----
#---- **plot data ----
truth <- 
  data.frame("Race" = c("White", "Black", "Hispanic"), 
             "prev" = c(
               unique(round(na.omit(results$true_dem_prev_white), 2)), 
               unique(round(na.omit(results$true_dem_prev_black), 2)), 
               unique(round(na.omit(results$true_dem_prev_hispanic), 2)))) %>% 
  mutate_at("Race", 
            function(x) factor(x, levels = c("White", "Black", "Hispanic")))

#columns to select
colnames <- expand_grid(c("mean", "LCI", "UCI"), 
                        c("dem_prev"), 
                        c("white", "black", "hispanic"), 
                        c("BLMM", "LKW", "Hurd", "ModHurd")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% 
  group_by(HRS_sample_size) %>% 
  summarise_at(colnames, mean) %>% 
  pivot_longer(all_of(colnames),
               names_to = c(".value", "Race"), 
               names_pattern = "(.*?)_(.*)") %>% 
  mutate_at("Race", function(x) str_remove(x, "dem_prev_")) %>% 
  separate(Race, into = c("Race", "Algorithm"), sep = "_") %>%
  mutate_at("Race", str_to_sentence) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>% 
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  )

#---- **plot ----
ptshp <- c("LTI" = 21, "ModHurd" = 24, "Hurd" = 23, "LKW" = 22)

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = Algorithm, color = Algorithm)) + 
  geom_vline(aes(xintercept = prev), data = truth) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0, linewidth = 0.8, 
                position = position_dodge(-0.8)) +
  geom_point(size = 3, position = position_dodge(-0.8), fill = "white") + 
  facet_grid(cols = vars(Race)) + theme_bw() + 
  xlab("Prevalence of dementia") + ylab('"HRS" sample size') + 
  scale_x_continuous(limits = c(0, 0.45), breaks = seq(0.0, 0.45, by = 0.1)) +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24), 
        panel.grid = element_blank(),
        legend.key.width = unit(2, "line"))  

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXX_mean_CI_dem_prev.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = Algorithm)) + 
  geom_vline(aes(xintercept = prev), data = truth) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI, linetype = Algorithm), width = 0, linewidth = 0.8,
                position = position_dodge(-0.8)) +
  geom_point(size = 3, position = position_dodge(-0.8), fill = "white") + 
  facet_grid(cols = vars(Race)) + theme_bw() + 
  xlab("Prevalence of dementia") + ylab('"HRS" sample size') + 
  scale_x_continuous(limits = c(0, 0.45), breaks = seq(0.0, 0.45, by = 0.1)) +
  scale_shape_manual(values = ptshp) + 
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
  theme(text = element_text(size = 24), 
        panel.grid = element_blank(),
        legend.key.width = unit(2, "line")) 

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXX_mean_CI_dem_prev_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

#---- Figure XX: 95% CI coverage dementia prevalence ----
#---- **plot data ----
#columns to select
colnames <- expand_grid(c("dem_prev_coverage"), 
                        c("white", "black", "hispanic"), 
                        c("BLMM", "LKW", "Hurd", "ModHurd")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% 
  group_by(HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(colnames, mean) %>% 
  pivot_longer(all_of(colnames), names_to = c(".value", "Race"), 
               names_pattern = "(.*?)_(.*)") %>%
  mutate_at("Race", function(x) str_remove(x, "prev_coverage_")) %>% 
  separate(Race, into = c("Race", "Algorithm"), sep = "_") %>%
  mutate_at("Race", str_to_sentence) %>%
  mutate_at("Race", function(x) 
    factor(x, levels = c("White", "Black", "Hispanic"))) %>%
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  ) %>% 
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("coverage_percent" = dem*100)

#---- **plot ----

ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = coverage_percent, group = Algorithm, 
           color = Algorithm, shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  facet_grid(cols = vars(Race)) +
  geom_hline(yintercept = 95, lty = "dotted", linewidth = 0.8) +
  theme_bw() + ylab("95% interval coverage") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24), 
        legend.key.width = unit(2, "line"))

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXX_dem_prev_coverage.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in") 


ggplot(data = plot_data, 
       aes(x = HRS_sample_size, y = coverage_percent, group = Algorithm, 
           # color = Algorithm, 
           shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  facet_grid(cols = vars(Race)) +
  geom_hline(yintercept = 95, lty = "dotted", linewidth = 0.8) +
  theme_bw() + ylab("95% interval coverage") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))     

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXX_dem_prev_coverage_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in") 

#---- Figure XXa-b: PR + PD dementia ----
#---- **plot data ----
truth <- 
  data.frame("Comparison" = c("Black vs. White", "Hispanic vs. White"), 
             "PR" = c(unique(na.omit(round(results$true_dem_prev_black/
                                             results$true_dem_prev_white, 2))), 
                      unique(na.omit(round(results$true_dem_prev_hispanic/
                                             results$true_dem_prev_white, 2)))), 
             "PD" = c(unique(na.omit(round(results$true_dem_prev_black -
                                             results$true_dem_prev_white, 2))), 
                      unique(na.omit(round(results$true_dem_prev_hispanic -
                                             results$true_dem_prev_white, 2))))) %>% 
  mutate_at("Comparison", 
            function(x) 
              factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  pivot_longer(c("PR", "PD"), names_to = "measure") %>% 
  mutate("null" = ifelse(measure == "PR", 1, 0))

#columns to select
PR_colnames <- expand_grid(c("mean", "LCI", "UCI"), "PR",
                           c("black", "hispanic"), 
                           c("BLMM", "LKW", "Hurd", "ModHurd")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

PD_colnames <- expand_grid(c("mean", "LCI", "UCI"), "PD",
                           c("black", "hispanic"), 
                           c("BLMM", "LKW", "Hurd", "ModHurd")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

#---- **PR plot ----
plot_data <- results %>% 
  group_by(HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(PR_colnames, function(x) exp(mean(log(x)))) %>% 
  pivot_longer(all_of(PR_colnames),
               names_to = c(".value", "measure", "Race", "Algorithm"), 
               names_sep = "_") %>%
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  ) 

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = Algorithm, color = Algorithm)) + 
  geom_vline(aes(xintercept = value), 
             data = truth %>% filter(measure == "PR")) +
  geom_vline(aes(xintercept = null), 
             data = truth %>% filter(measure == "PR"), lty = "dashed") +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0, linewidth = 0.8, 
                position = position_dodge(-0.8)) +
  geom_point(size = 3, position = position_dodge(-0.8), fill = "white") + 
  scale_shape_manual(values = ptshp) + 
  facet_grid(cols = vars(Comparison)) + theme_bw() +
  xlab("Prevalence Ratio (PR)") + ylab('"HRS" sample size') + 
  theme(text = element_text(size = 24), 
        panel.grid = element_blank(), 
        legend.key.width = unit(2, "line"))  

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXa_PR.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = Algorithm)) + 
  geom_vline(aes(xintercept = value), 
             data = truth %>% filter(measure == "PR")) +
  geom_vline(aes(xintercept = null), 
             data = truth %>% filter(measure == "PR"), lty = "dashed") +
  geom_errorbar(aes(xmin = LCI, xmax = UCI, linetype = Algorithm), 
                width = 0, linewidth = 0.8, 
                position = position_dodge(-0.8)) +
  geom_point(size = 3, position = position_dodge(-0.8), fill = "white") + 
  scale_shape_manual(values = ptshp) + 
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
  facet_grid(cols = vars(Comparison)) + theme_bw() +
  xlab("Prevalence Ratio (PR)") + ylab('"HRS" sample size') + 
  theme(text = element_text(size = 24), 
        panel.grid = element_blank(), 
        legend.key.width = unit(2, "line"))  

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXa_PR_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")


#---- **PD plot ----
plot_data <- results %>% 
  group_by(HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(PD_colnames, mean) %>% 
  pivot_longer(all_of(PD_colnames),
               names_to = c(".value", "measure", "Race", "Algorithm"), 
               names_sep = "_") %>%
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  ) 

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = Algorithm, color = Algorithm)) + 
  geom_vline(aes(xintercept = value), 
             data = truth %>% filter(measure == "PD")) +
  geom_vline(aes(xintercept = null), 
             data = truth %>% filter(measure == "PD"), lty = "dashed") +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0, linewidth = 0.8, 
                position = position_dodge(-0.8)) +
  geom_point(size = 3, position = position_dodge(-0.8), fill = "white") + 
  scale_shape_manual(values = ptshp) + 
  facet_grid(cols = vars(Comparison)) + theme_bw() +
  xlab("Prevalence Difference (PD)") + ylab('"HRS" sample size') + 
  theme(text = element_text(size = 24),
        panel.grid = element_blank(), 
        legend.key.width = unit(2, "line")) 

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXb_PD.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

ggplot(data = plot_data, 
       aes(x = mean, y = HRS_sample_size, shape = Algorithm)) + 
  geom_vline(aes(xintercept = value), 
             data = truth %>% filter(measure == "PD")) +
  geom_vline(aes(xintercept = null), 
             data = truth %>% filter(measure == "PD"), lty = "dashed") +
  geom_errorbar(aes(xmin = LCI, xmax = UCI, linetype = Algorithm), 
                width = 0, linewidth = 0.8, 
                position = position_dodge(-0.8)) +
  geom_point(size = 3, position = position_dodge(-0.8), fill = "white") + 
  scale_shape_manual(values = ptshp) + 
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
  facet_grid(cols = vars(Comparison)) + theme_bw() +
  xlab("Prevalence Difference (PD)") + ylab('"HRS" sample size') + 
  theme(text = element_text(size = 24),
        panel.grid = element_blank(), 
        legend.key.width = unit(2, "line")) 

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXb_PD_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

#---- Figure XXa-b: 95% CI coverage PR + PD ----
#---- **plot data ----
#columns to select
colnames <- expand_grid(c("PR", "PD"), "coverage",
                        c("black", "hispanic"), 
                        c("BLMM", "LKW", "Hurd", "ModHurd")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

plot_data <- results %>% 
  group_by(HRS_sample_size, HCAP_sample_size) %>% 
  summarise_at(colnames, mean) %>% 
  pivot_longer(all_of(colnames), 
               names_to = c("measure", ".value", "Race", "Algorithm"), 
               names_sep = "_") %>%
  mutate("Comparison" = case_when(Race == "black" ~ "Black vs. White", 
                                  Race == "hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate_at("HRS_sample_size", as.factor) %>% 
  mutate("coverage_percent" = coverage*100) %>%
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  ) 

#---- **PR plot ----

ggplot(data = plot_data %>% filter(measure == "PR"), 
       aes(x = HRS_sample_size, y = coverage_percent, group = Algorithm, 
           shape = Algorithm, color = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  facet_grid(cols = vars(Comparison)) + 
  geom_hline(yintercept = 95, lty = "dotted", linewidth = 0.8) + 
  theme_bw() + ylab("95% interval coverage\n(PR)") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))  

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXa_PR_coverage.pdf"), 
       dpi = 300, width = 13.25, height = 4, units = "in")

ggplot(data = plot_data %>% filter(measure == "PR"), 
       aes(x = HRS_sample_size, y = coverage_percent, group = Algorithm, 
           shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  facet_grid(cols = vars(Comparison)) + 
  geom_hline(yintercept = 95, lty = "dotted", linewidth = 0.8) + 
  theme_bw() + ylab("95% interval coverage\n(PR)") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))  

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXa_PR_coverage_bw.pdf"), 
       dpi = 300, width = 13.25, height = 4, units = "in")

#---- **PD plot ----
ggplot(data = plot_data %>% filter(measure == "PD"), 
       aes(x = HRS_sample_size, y = coverage_percent, group = Algorithm, 
           shape = Algorithm, color = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  facet_grid(cols = vars(Comparison)) +
  geom_hline(yintercept = 95, lty = "dotted", linewidth = 0.8) + 
  #scale_y_continuous(limits = c(0, 102), breaks = c(0, seq(80, 100, by = 5))) + 
  #scale_y_cut(breaks = c(79), space = 0.2, which = c(1, 2), scales = c(5, 0.5)) +
  theme_bw() + ylab("95% interval coverage\n(PD)") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))   

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXb_PD_coverage.pdf"), 
       dpi = 300, width = 13.25, height = 4, units = "in")

ggplot(data = plot_data %>% filter(measure == "PD"), 
       aes(x = HRS_sample_size, y = coverage_percent, group = Algorithm, 
           shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  facet_grid(cols = vars(Comparison)) +
  geom_hline(yintercept = 95, lty = "dotted", linewidth = 0.8) + 
  #scale_y_continuous(limits = c(0, 102), breaks = c(0, seq(80, 100, by = 5))) + 
  #scale_y_cut(breaks = c(79), space = 0.2, which = c(1, 2), scales = c(5, 0.5)) +
  theme_bw() + ylab("95% interval coverage\n(PD)") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))   

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "figureXXb_PD_coverage_bw.pdf"), 
       dpi = 300, width = 13.25, height = 4, units = "in")

#---- Appendix Figure XX: bias dementia prevalence ----
#---- **plot data ----
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
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  )

#---- **plot ----

ggplot(data = plot_data, aes(x = HRS_sample_size, y = bias, group = Algorithm, 
                             shape = Algorithm, color = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  geom_hline(yintercept = 0, lty = "dotted", linewidth = 0.8) +
  theme_bw() + ylab("Bias") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  facet_grid(cols = vars(Race)) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line")) 

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_dem_prev_bias.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

ggplot(data = plot_data, aes(x = HRS_sample_size, y = bias, group = Algorithm, 
                             shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  geom_hline(yintercept = 0, lty = "dotted", linewidth = 0.8) +
  theme_bw() + ylab("Bias") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  facet_grid(cols = vars(Race)) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line")) 

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_dem_prev_bias_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4, units = "in")

#---- Appendix Figure XX: RMSE dementia prevalence ----
ggplot(data = plot_data, aes(x = HRS_sample_size, y = RMSE, group = Algorithm, 
                             shape = Algorithm, color = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  theme_bw() + ylab("RMSE") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  facet_grid(cols = vars(Race)) + ylim(0, 0.25) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))      

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_dem_prev_rmse.jpeg"), 
       dpi = 300, width = 13.5, height = 4, units = "in") 

ggplot(data = plot_data, aes(x = HRS_sample_size, y = RMSE, group = Algorithm, 
                             shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  theme_bw() + ylab("RMSE") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  facet_grid(cols = vars(Race)) + ylim(0, 0.25) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"))      

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_dem_prev_rmse_bw.jpeg"), 
       dpi = 300, width = 13.5, height = 4, units = "in") 

#---- Appendix Figure XX: PR + PD bias ----
#---- **plot data ----
cols_by_race <- expand_grid(c("mean"), 
                            c("PR", "PD"), 
                            c("black", "hispanic"), 
                            c("BLMM", "ModHurd", "Hurd", "LKW")) %>% 
  unite("names", everything(), sep = "_") %>% unlist() %>% unname()

truth <- 
  data.frame("Comparison" = c("Black vs. White", "Hispanic vs. White"),
             "race" = c("black", "hispanic"),
             "PR" = c(unique(na.omit(round(results$true_dem_prev_black/
                                             results$true_dem_prev_white, 2))), 
                      unique(na.omit(round(results$true_dem_prev_hispanic/
                                             results$true_dem_prev_white, 2)))), 
             "PD" = c(unique(na.omit(round(results$true_dem_prev_black -
                                             results$true_dem_prev_white, 2))), 
                      unique(na.omit(round(results$true_dem_prev_hispanic -
                                             results$true_dem_prev_white, 2))))) %>% 
  pivot_longer(c("PR", "PD")) %>% rename("true" = "value")

# JZ: the left_join(truth) line gives the following warning: 
# Joining with `by = join_by(race)`
# Warning in View :
#   Detected an unexpected many-to-many relationship between `x` and `y`.
# truth should be joined to the results by both race and PR/PD, not race alone
truth <- truth %>% rename(measure = name)
# After updating the code, the figure is different compared to before. 
# I saved my updated figures named "appendix_figureXX_PR_PD_bias_JZ.jpeg" and
# "appendix_figureXX_PR_PD_rmse_JZ.jpeg" in the same figures folder. 

plot_data <- results %>% ungroup() %>%
  dplyr::select("HRS_sample_size", "HCAP_sample_size",  all_of(cols_by_race)) %>%
  pivot_longer(all_of(cols_by_race),
               names_to = c(".value", "measure", "race", "Algorithm"), 
               names_sep = "_") %>% 
  # left_join(truth) %>% # this gives the warning above
  left_join(truth, by = c("race", "measure")) %>% # my update here
  mutate("error" = mean - true) %>%
  mutate("squared_error" = error^2) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  mutate_at("race", function(x) str_to_sentence(x)) %>% 
  mutate("Comparison" = case_when(race == "Black" ~ "Black vs. White", 
                                  race == "Hispanic" ~ "Hispanic vs. White")) %>%
  mutate_at("Comparison", function(x) 
    factor(x, levels = c("Black vs. White", "Hispanic vs. White"))) %>%
  mutate_at("HRS_sample_size", as.factor) %>%
  group_by(HRS_sample_size, Comparison, measure, Algorithm) %>% 
  summarize_at(c("error", "squared_error"), mean) %>% 
  rename(c("bias" = "error")) %>% 
  mutate("RMSE" = sqrt(squared_error)) %>% 
  mutate_at("measure", function(x) 
    factor(x, levels = c("PR", "PD"))) %>%
  mutate(
    Algorithm = ifelse(Algorithm == "BLMM", "LTI", Algorithm) %>% 
      factor(levels = c("LTI", "ModHurd", "Hurd", "LKW"))
  ) 

#---- **plot ----

ggplot(data = plot_data, aes(x = HRS_sample_size, y = bias, group = Algorithm, 
                             shape = Algorithm, color = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  geom_hline(yintercept = 0, lty = "dotted", linewidth = 0.8) + 
  theme_bw() + ylab("Bias") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  ylim(c(-2.5, 3)) +
  # facet_grid(rows = vars(measure), cols = vars(Comparison)) +
  facet_grid(measure ~ Comparison, scales = "free_y") +
  facetted_pos_scales(
    y = list(
      measure == "PR" ~ scale_y_continuous(limits = c(-2, 3)),
      measure == "PD" ~ scale_y_continuous(limits = c(-0.08, 0.08),
                                           breaks = c(-0.08, -0.04, 0, 0.04, 0.08))
    )
  ) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"),
        panel.spacing = unit(1, "lines"))      

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_PR_PD_bias.pdf"), 
       dpi = 300, width = 13.5, height = 4.5, units = "in") 

ggplot(data = plot_data, aes(x = HRS_sample_size, y = bias, group = Algorithm, 
                             shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  geom_hline(yintercept = 0, lty = "dotted", linewidth = 0.8) + 
  theme_bw() + ylab("Bias") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  ylim(c(-2.5, 3)) +
  # facet_grid(rows = vars(measure), cols = vars(Comparison)) +
  facet_grid(measure ~ Comparison, scales = "free_y") +
  facetted_pos_scales(
    y = list(
      measure == "PR" ~ scale_y_continuous(limits = c(-2, 3)),
      measure == "PD" ~ scale_y_continuous(limits = c(-0.08, 0.08),
                                           breaks = c(-0.08, -0.04, 0, 0.04, 0.08))
    )
  ) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"),
        panel.spacing = unit(1, "lines"))      

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_PR_PD_bias_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4.5, units = "in") 

#---- Appendix Figure XX: PR + PD RMSE ----
ggplot(data = plot_data, aes(x = HRS_sample_size, y = RMSE, group = Algorithm, 
                             color = Algorithm, shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  theme_bw() + ylab("RMSE") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  ylim(c(0, 3)) +
  facet_grid(measure ~ Comparison, scales = "free_y") +
  facetted_pos_scales(
    y = list(
      measure == "PR" ~ scale_y_continuous(limits = c(0, 3)),
      measure == "PD" ~ scale_y_continuous(limits = c(0, 0.08), breaks = c(0, 0.04, 0.08))
    )
  ) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"),
        panel.spacing = unit(1, "lines"))      

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_PR_PD_rmse.pdf"), 
       dpi = 300, width = 13.5, height = 4.5, units = "in") 

ggplot(data = plot_data, aes(x = HRS_sample_size, y = RMSE, group = Algorithm, 
                             shape = Algorithm)) + 
  geom_line(aes(linetype = Algorithm), linewidth = 0.8) + 
  geom_point(size = 3, fill = "white") + 
  theme_bw() + ylab("RMSE") + xlab('"HRS" Sample Size') +
  scale_shape_manual(values = ptshp) + 
  ylim(c(0, 3)) +
  facet_grid(measure ~ Comparison, scales = "free_y") +
  facetted_pos_scales(
    y = list(
      measure == "PR" ~ scale_y_continuous(limits = c(0, 3)),
      measure == "PD" ~ scale_y_continuous(limits = c(0, 0.08), breaks = c(0, 0.04, 0.08))
    )
  ) +
  theme(text = element_text(size = 24),
        legend.key.width = unit(2, "line"),
        panel.spacing = unit(1, "lines"))      

ggsave(filename = paste0(path_to_box, 
                         "papers/paper2_model_comparison_ADAMS_prior/figures/", 
                         "appendix_figureXX_PR_PD_rmse_bw.pdf"), 
       dpi = 300, width = 13.5, height = 4.5, units = "in") 
