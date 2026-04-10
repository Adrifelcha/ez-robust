######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


########################################################################################
# C E L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "AUC-grid", "cell-simulation")
main_dir <- here("output", "RData", "cell-simulation")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AUC nested by parameter (new function)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Default settings: x_param = "drift", T = 40, beta = 0.2
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean")
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", t_level_select = 160)

# Alternate test: x-axis is boundary parameter
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean")
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", t_level_select = 160)

# Plot nested simulation study results with T = 40 and T = 160
plot_AUC_nested(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", run_diagnose = FALSE)

# Simplified 2 x 2 nested layout at one trial level T
plot_AUC_nested_2x2(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", run_diagnose = FALSE)


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

plot_AUCgrid_full(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label) 
plot_AUCgrid_full(main_dir = main_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Restricted range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "restrictedRange"
main_dir <- here("output", "RData", "full-simulation", "restricted-parameters")

plot_AUCgrid_full(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label) 
plot_AUCgrid_full(main_dir = main_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)


