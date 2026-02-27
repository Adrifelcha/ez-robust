######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

###################################################################################
###  F U L L     S I M U L A T I O N     S T U D Y  ###############################
###################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define directories
seed_dir <- here("demos", "simStudy_full", "01-wide-parameters", "samples")
output_dir <- here("output", "RData", "full-simulation", "wide-parameters")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Restricted parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define directories
seed_dir <- here("demos", "simStudy_full", "02-restricted-parameters", "samples")
output_dir <- here("output", "RData", "full-simulation", "restricted-parameters")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)

###################################################################################
###  C E L L     S I M U L A T I O N     S T U D Y  ###############################
###################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# High drift and high bound
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define directories
seed_dir <- here("demos", "simStudy_cells", "highDrift-highBound", "samples")
output_dir <- here("output", "RData", "cell-simulation", "highDrift-highBound")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Low drift and low bound
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define directories
seed_dir <- here("demos", "simStudy_cells", "lowDrift-lowBound", "samples")
output_dir <- here("output", "RData", "cell-simulation", "lowDrift-lowBound")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# High drift and low bound
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define directories
seed_dir <- here("demos", "simStudy_cells", "highDrift-lowBound", "samples")
output_dir <- here("output", "RData", "cell-simulation", "highDrift-lowBound")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Low drift and high bound
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seed_dir <- here("demos", "simStudy_cells", "lowDrift-highBound", "samples")
output_dir <- here("output", "RData", "cell-simulation", "lowDrift-highBound")
# Process and compress RData files for this simulation study
process_sim_data_by_cell(seed_dir = seed_dir, output_dir = output_dir)