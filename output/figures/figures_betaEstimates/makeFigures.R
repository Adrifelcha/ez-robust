library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()

main_dir <- here("output", "RData_simStudy_results")
output_dir <- here("output", "figures_betaEstimates")

# Generate grid plot showing mean posterior beta estimates across conditions
# with error bars showing +/- 1 SD
plot_betaEstimateGrid(main_dir = main_dir, output_dir = output_dir, y_axis_ticks = 5)
