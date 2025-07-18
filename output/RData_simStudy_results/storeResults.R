library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

# Define directories
seed_dir <- here("demos", "simulation-study", "samples-full_Jul8")
output_dir <- here("output", "RData_simStudy_results")
#
# Run the memory-safe processing function
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)
