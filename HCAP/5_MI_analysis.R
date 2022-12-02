#---- Package loading + options ----
if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("tidyverse")

#---- source functions ----
source(here::here("HCAP", "functions", "analysis_function.R"))
source(here::here("HCAP", "functions", "generate_synthetic_function.R"))

#---- read in data ----
path_to_box <- "/Users/crystalshaw/Library/CloudStorage/Box-Box/Dissertation/"

#---- **HCAP data ----
HCAP_imputed <- 
  readRDS(paste0(path_to_box, "analyses/HCAP/HCAP_MI_hotdeck")) %>% 
  lapply(., function(x) x %<>% mutate("(Intercept)" = 1))

#---- **variable labels ----
variable_labels <- 
  read_csv(paste0(path_to_box, "data/variable_crosswalk.csv")) 

#---- **cell ID key ----
cell_ID_key <- read_csv(paste0(path_to_box, "data/cell_ID_key.csv")) %>% 
  mutate_all(as.character)

#---- **impairment class color palette ----
color_palette <- read_csv(paste0(path_to_box, "data/color_palette.csv")) 

#---- **contrasts matrix ----
A <- read_csv(paste0(path_to_box, "data/contrasts_matrix.csv")) %>% as.matrix()

#---- define vars ----
#---- **selected variables ----
selected_vars <- 
  read_csv(paste0(path_to_box, "data/variable_selection/", 
                  "model_coefficients.csv"))

#categorical vars (notation from Schafer 1997)
W <- c("black", "hispanic", "stroke")

#continuous vars (notation from Schafer 1997)
Z <- selected_vars[str_detect(selected_vars$data_label, "_Z"), 
                   "data_label"] %>% unlist() %>% as.vector()

#---- run analysis ----
set.seed(20221125)

MI_results <- lapply(HCAP_imputed, function(x) 
  analysis_function(warm_up = 100, starting_props = c(0.25, 0.25, 0.25, 0.25),
                    dataset_to_copy = x,
                    orig_means = 
                      x %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
                      colMeans() %>% t() %>% as.data.frame(),
                    orig_sds = 
                      x %>% dplyr::select(all_of(str_remove(Z, "_Z"))) %>%
                      apply(., 2, sd) %>% t() %>% as.data.frame(),
                    calibration_sample = FALSE, calibration_prop = NA,
                    calibration_sample_name = NA,
                    path_to_data = paste0(path_to_box,"data/"),
                    categorical_vars = W, continuous_vars = Z, id_var = "HHIDPN",
                    variable_labels = variable_labels, cell_ID_key = cell_ID_key,
                    color_palette = color_palette, contrasts_matrix = A,
                    kappa_0_mat = 
                      read_csv(paste0(path_to_box, 
                                      "data/tuning/kappa_0_matrix_HCAP.csv")),
                    nu_0_mat = 
                      read_csv(paste0(path_to_box, 
                                      "data/tuning/nu_0_matrix_HCAP.csv")),
                    num_synthetic = 1000))

results <- do.call(rbind, MI_results) %>% as.data.frame()

#save results
write_csv(results, 
          paste0(path_to_box, "analyses/HCAP/MI_prediction_results.csv"))

#---- compare results with HCAP ----
#---- **proportions of cognitive impairment classes ----
#---- ****read in results ----
results <- 
  read_csv(paste0(path_to_box, "analyses/HCAP/MI_prediction_results.csv"))

#---- ****Rubin's Rules ----
#---- ******calculate within-SE ----
for(class in c("Unimpaired", "MCI", "Dementia", "Other")){
  results %<>%
    mutate(!!sym(paste0("SE_", class)) :=
             rowMeans(cbind(
               abs(!!sym(paste0("mean_", class)) -
                     !!sym(paste0("LCI_", class))),
               abs(!!sym(paste0("mean_", class)) -
                     !!sym(paste0("UCI_", class)))))/1.96)
}

within_vars <- 
  colMeans(
    (results[, paste0("SE_", c("Unimpaired", "MCI", "Dementia", "Other"))])^2)

between_vars <- 
  apply(results[, paste0("mean_", c("Unimpaired", "MCI", "Dementia", "Other"))], 
        2, var)

total_vars <- within_vars + between_vars + between_vars/nrow(results)

#---- ****plot data ----
plot_data <- 
  data.frame("class" = c("Unimpaired", "MCI", "Dementia", "Other"), 
             "mean" = NA, "LCI" = NA, "UCI" = NA, "Color" = NA)

plot_data[, "mean"] <- colMeans(results[, paste0("mean_", plot_data$class)])
plot_data[, "LCI"] <- plot_data$mean - 1.96*sqrt(total_vars)
plot_data[, "UCI"] <- plot_data$mean + 1.96*sqrt(total_vars)
plot_data[, "Color"] <- color_palette$Color

#scale to proportions
plot_data[, c("mean", "LCI", "UCI")] <- 
  plot_data[, c("mean", "LCI", "UCI")]/nrow(HCAP_imputed[[1]])

#study
plot_data %<>% mutate("Algorithm" = "Bayesian Latent Class\nMixture Model")

#HCAP data 
HCAP_results <- 
  data.frame("class" = c("Unimpaired", "MCI", "Dementia", "Other"), 
             "mean" = c(NA, 0.231028037, .1375700935, NA), 
             "LCI" = c(NA, .1846242991, .1034616822, NA), 
             "UCI" = c(NA, .2638429907, .1733345794, NA), 
             "Color" = plot_data$Color, 
             "Algorithm" = "HCAP")

plot_data %<>% rbind(HCAP_results)

#---- ****plot ----
plot_data$class <- 
  factor(plot_data$class, 
         levels = rev(c("Unimpaired", "MCI", "Dementia", "Other")))

ggplot(data = plot_data, 
       aes(y = class, x = mean, color = class, shape = Algorithm)) + theme_bw() + 
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.10) + 
  geom_point(size = 3) +
  xlab("") + ylab("Prevalence") + theme(legend.position = "none") + 
  scale_color_manual(values = unique(plot_data$Color[order(plot_data$class)])) + 
  theme(text = element_text(size = 12)) + 
  guides(shape = guide_legend(title = "Algorithm"))
  

ggsave(filename = 
         paste0(path_to_box, "figures/chapter_6/", 
                "figure6.4_algorithmic_dementia_classification.jpeg"), 
       dpi = 300, width = 5, height = 3, units = "in")

#---- Rubin's rules function ----
rubin_rules <- function(data, var, impairment_class){
  pt_est <- mean(unlist(results[, paste0(impairment_class, "_logOR_", var)]))
  within_var <- mean(unlist(results[, paste0(impairment_class, "_SE_", var)]^2))
  between_var <- var(results[, paste0(impairment_class, "_logOR_", var)])
  total_var <- within_var + between_var + between_var/nrow(data)
  
  return(paste0(var, " ", impairment_class, ": ", 
                round(exp(pt_est), 2), 
                " (", round(exp(pt_est - 1.96*sqrt(total_var)), 2), ", " 
                , round(exp(pt_est + 1.96*sqrt(total_var)), 2), ")"))
}

#---- **increased dementia risk with 5-year age increase ----
rubin_rules(results, "age5", "dementia")

#---- **increased mci risk with 5-year age increase ----
rubin_rules(results, "age5", "mci")

