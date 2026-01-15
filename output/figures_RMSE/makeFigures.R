library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


main_dir <- here("output", "RData_simStudy_results")
output_dir <- here("output", "figures_RMSE")

# Parameter to plot (can be: "bound_mean", "drift_mean", "nondt_mean", "betaweight")
parameter <- "betaweight"

# This function generates three PDFs: RMSE, Bias, and Variance
# All three plots use the same data computation (called once)
plot_RMSEgrid(main_dir = main_dir, output_dir = output_dir, parameter = parameter, plot_by = "beta") 

# Optional: Plot with conditions on the x-axis instead
#plot_RMSEgrid(main_dir = main_dir, output_dir = output_dir, parameter = parameter, plot_by = "condition")