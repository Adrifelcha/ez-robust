library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

# Define directories
seed_dir <- here("demos", "simulation-study", "samples")
output_dir <- here("output", "RData_simStudy_results")
#
# Run the memory-safe processing function
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)

source(here("output", "figures_AUC", "makeFigures.R"))
source(here("output", "figures_BayesFactors", "makeFigures.R"))
source(here("output", "figures_ROC", "makeFigures.R"))