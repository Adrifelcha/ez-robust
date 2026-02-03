######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

######################################
# Full simulation study
######################################
# Define directories
seed_dir <- here("demos", "simStudy_full", "samples")
output_dir <- here("output", "RData", "full-simulation")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)

######################################
# Restricted simulation study
######################################
# Define directories
seed_dir <- here("demos", "simStudy_restricted", "samples")
output_dir <- here("output", "RData", "restricted-simulation")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)

