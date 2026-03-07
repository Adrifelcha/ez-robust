######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()

########################################################################################
# F U L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
# Establish output directory
output_dir <- here("output", "figures", "AUC-grid", "full-simulation")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "wideRange"
main_dir <- here("output", "RData", "full-simulation", "wide-parameters")

plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label) 
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Restricted range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "restrictedRange"
main_dir <- here("output", "RData", "full-simulation", "restricted-parameters")

plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label) 
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)

########################################################################################
# C E L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "AUC-grid", "cell-simulation")
main_dir <- here("output", "RData", "cell-simulation")
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "beta")


# Low drift and low bound
cell_dir <- here("output", "RData", "cell-simulation", "lowDrift-lowBound")
plot_AUCgrid(main_dir = cell_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label)
plot_AUCgrid(main_dir = cell_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)

# High drift and high bound
label <- "highDrift-highBound"
cell_dir <- here("output", "RData", "cell-simulation", "highDrift-highBound")
plot_AUCgrid(main_dir = cell_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label)
plot_AUCgrid(main_dir = cell_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)

# High drift and low bound
