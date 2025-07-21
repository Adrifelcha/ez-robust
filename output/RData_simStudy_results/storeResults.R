rm(list = ls())
library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

# Define directories
seed_dir <- here("demos", "simulation-study", "samples")
output_dir <- here("output", "RData_simStudy_results")
#
# Run the memory-safe processing function
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)

load(here("demos", "simulation-study", "samples", "seed-433.RData"))
# load(here("demos", "simulation-study", "samples", "seed-433.RData"))