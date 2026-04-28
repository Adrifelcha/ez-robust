######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


########################################################################################
# FOLLOW-UP   S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "AUC-grid", "cell-simulation")
main_dir <- here("output", "RData", "cell-simulation")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Side-by-side: AUC as a function of true parameter values (panels distinguish between low and high)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mu_drift on the x-axis, beta = 0.2
# T = 40 (default)
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean")
# T = 160
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", t_level_select = 160)

# Mu_bound on the x-axis, beta = 0.2
# T = 40 (default)
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean")
# T = 160
plot_AUC_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", t_level_select = 160)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 x 2 AUC as a function of betaweight parameter for one T level, across high/low bound and low/high drift.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_AUC_nested_2x2(main_dir = main_dir, output_dir = output_dir, run_diagnose = FALSE)


########################################################################################
# F U L L     S I M U L A T I O N     S T U D Y  #######################################
#########################################################################################
# These simulations include 5 levels of trial size (T) and sample size (P)
# Establish output directory
output_dir <- here("output", "figures", "AUC-grid", "full-simulation")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "wideRange"
main_dir <- here("output", "RData", "full-simulation", "wide-parameters")

plot_AUCgrid_full(main_dir = main_dir, output_dir = output_dir, custom_title_label = label) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Restricted range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "restrictedRange"
main_dir <- here("output", "RData", "full-simulation", "restricted-parameters")

plot_AUCgrid_full(main_dir = main_dir, output_dir = output_dir, custom_title_label = label) 


