library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


main_dir <- here("output", "simStudy_results")
output_dir <- here("output", "figures_AUC")

# Plot with conditions on the x-axis
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "condition")

# Plot with beta levels on the x-axis
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "beta") 
