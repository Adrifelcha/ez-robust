######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


########################################################################################
# C E L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "RMSE-grid", "cell-simulation")
main_dir <- here("output", "RData", "cell-simulation")

plot_RMSE_meanDrift_beta_sideBySide(main_dir = main_dir, output_dir = output_dir, annotate_bias_highbound_segment = FALSE)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combined 3x4 grid: drift (cols 1-2) and betaweight (cols 3-4), low vs high bound per pair.
# Writes three PDFs (RMSE, Bias, Variance). Empty panels were fixed by assigning rbind to
# panel_raw_drift[[...]] / panel_raw_beta[[...]] directly (not via a temporary list handle).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_bias_nested_drift_betaweight_by_bound(main_dir = main_dir, output_dir = output_dir)
plot_bias_nested_drift_betaweight_by_bound(main_dir = main_dir, output_dir = output_dir, t_level_select = 160)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# One row, two columns: split_by = low vs high drift OR bound; x-axis = binned drift or bound;
# all beta_levels_select overlaid per panel. Three PDFs (RMSE, Bias, Variance).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_RMSE_nested_onerow_by_split(
  main_dir = main_dir, output_dir = output_dir,
  split_by = "bound_mean", x_axis_param = "drift_mean", parameter = "betaweight"
)

plot_RMSE_nested_onerow_by_split(
  main_dir = main_dir, output_dir = output_dir,
  split_by = "bound_mean", x_axis_param = "drift_mean", parameter = "drift_mean"
)


plot_RMSE_nested_onerow_by_split(
  main_dir = main_dir, output_dir = output_dir, t_level_select = 160,
  split_by = "bound_mean", x_axis_param = "drift_mean", parameter = "betaweight"
)
# plot_RMSE_nested_onerow_by_split(main_dir, output_dir, split_by = "drift_mean", x_axis_param = "bound_mean", parameter = "drift_mean")


########################################################################################
# F U L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "wideRange"
main_dir <- here("output", "RData", "full-simulation", "wide-parameters")
output_dir <- here("output", "figures", "RMSE-grid", "full-simulation", label)

plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "betaweight")
plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "drift_mean")
plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "bound_mean")
plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "nondt_mean")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Wide range of parameter values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label <- "restrictedRange"
main_dir <- here("output", "RData", "full-simulation", "restricted-parameters")
output_dir <- here("output", "figures", "RMSE-grid", "full-simulation", label)

plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "betaweight")
plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "drift_mean")
plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "bound_mean")
plot_RMSE_fullGrid(main_dir = main_dir, output_dir = output_dir, parameter = "nondt_mean")