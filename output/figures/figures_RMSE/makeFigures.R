######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


########################################################################################
# C E L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "figures_RMSE", "cell-simulation")
main_dir <- here("output", "RData", "cell-simulation")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RMSE nested by parameter (RMSE, bias, variance vs binned population x-axis)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parameter: which estimand (columns in simStudy_Beta): "betaweight", "drift_mean", "bound_mean", "nondt_mean".
# Rows still slice by simulation betaweight (beta_levels_select). Y-axis label uses beta, mu[nu], mu[alpha], or mu[tau].
# Default: x-axis drift (mu_nu), two columns = low vs high population bound; beta_levels_select default c(0, 0.2, 0.4).
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "betaweight")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "drift_mean")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "bound_mean")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "nondt_mean")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "betaweight", t_level_select = 160)
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "drift_mean", t_level_select = 160)
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "bound_mean", t_level_select = 160)
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", parameter = "nondt_mean", t_level_select = 160)

# Alternate: x-axis boundary (mu_alpha), columns = low vs high drift
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "betaweight")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "drift_mean")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "bound_mean")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "nondt_mean")
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "betaweight", t_level_select = 160)
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "drift_mean", t_level_select = 160)
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "bound_mean", t_level_select = 160)
plot_RMSE_nested_by_param(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", parameter = "nondt_mean", t_level_select = 160)


########################################################################################
# F U L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
# Parameter to plot (can be: "bound_mean", "drift_mean", "nondt_mean", "betaweight")
parameter <- "betaweight"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "wideRange"
main_dir <- here("output", "RData", "full-simulation", "wide-parameters")
output_dir <- here("output", "figures", "figures_RMSE", "full-simulation", label)

plot_RMSEgrid(main_dir = main_dir, output_dir = output_dir, parameter = parameter, plot_by = "beta")
# plot_RMSEgrid(main_dir = main_dir, output_dir = output_dir, parameter = parameter, plot_by = "condition")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Restricted range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "restrictedRange"
main_dir <- here("output", "RData", "full-simulation", "restricted-parameters")
output_dir <- here("output", "figures", "figures_RMSE", "full-simulation", label)

plot_RMSEgrid(main_dir = main_dir, output_dir = output_dir, parameter = parameter, plot_by = "beta")
# plot_RMSEgrid(main_dir = main_dir, output_dir = output_dir, parameter = parameter, plot_by = "condition")

