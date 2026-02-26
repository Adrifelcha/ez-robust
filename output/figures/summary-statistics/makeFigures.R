######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()

# Establish output directory
output_dir <- here("output", "figures", "summary-statistics")


######################################################
# Plot summary statistics for full simulation study
######################################################
label <- "full"
main_dir <- here("output", "RData", "full-simulation")

plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "bound_mean", point_cex = 0.5, custom_title_label = label)
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "drift_mean", point_cex = 0.5, custom_title_label = label)
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "bound_mean", point_cex = 0.5, custom_title_label = label,
                     third_param = "drift_mean", third_param_low = -1, third_param_high = 1)
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "drift_mean", point_cex = 0.5, custom_title_label = label,
                     third_param = "bound_mean", third_param_low = 3.5, third_param_high = 4)

plot_RTdiff_by_param_ratio(main_dir, output_dir, point_alpha = 1, point_cex = 0.5, custom_title_label = label)


######################################################
# Plot summary statistics for full simulation study
######################################################
label <- "restricted"
main_dir <- here("output", "RData", "restricted-simulation")

plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "bound_mean", point_cex = 0.5, custom_title_label = label)
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "drift_mean", point_cex = 0.5, custom_title_label = label)
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "bound_mean", point_cex = 0.5, custom_title_label = label,
                     third_param = "drift_mean", third_param_low = -1, third_param_high = 1)
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.7, x_param = "drift_mean", point_cex = 0.5, custom_title_label = label,
                     third_param = "bound_mean", third_param_low = 3.5, third_param_high = 4)

plot_RTdiff_by_param_ratio(main_dir, output_dir, point_alpha = 1, point_cex = 0.5, custom_title_label = label)