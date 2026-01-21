library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


main_dir <- here("output", "RData_simStudy_results")
output_dir <- here("output", "figures_summStats")

plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.5, true_param = "bound_mean")
plot_RTdiff_by_param(main_dir, output_dir, point_alpha = 0.5, true_param = "drift_mean")

