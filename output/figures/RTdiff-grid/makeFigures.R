######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()


########################################################################################
# C E L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "RTdiff-grid", "cell-simulation")
main_dir <- here("output", "RData", "cell-simulation")

######################################
# Contamination effect nested plot
######################################


# Beta-difference of RT differences (diffdiff)
plot_RTdiffdiff_nested(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", beta_alt_select = 0.2)
# Alternate param (optional)
plot_RTdiffdiff_nested(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", beta_alt_select = 0.2)



# Default: x_param = "drift_mean", T = 40, beta = 0
plot_RTcontamination_nested(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean")
# Optional IQR version
plot_RTcontamination_nested(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", show_IQrange = TRUE)


# Default test: x_param = "drift", T = 40, beta = 0
plot_RTdiff_nested(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", show_IQrange = FALSE)
plot_RTdiff_nested(main_dir = main_dir, output_dir = output_dir, x_param = "drift_mean", show_IQrange = TRUE)
# Alternate x-axis parameter test
#plot_RTdiff_nested(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", show_IQrange = FALSE)
#plot_RTdiff_nested(main_dir = main_dir, output_dir = output_dir, x_param = "bound_mean", show_IQrange = TRUE)
# Alternate beta level test (closest available beta is used if exact match is absent)
#plot_RTdiff_nested(main_dir = main_dir, output_dir = output_dir, beta_level_select = 0.2)

######################################
# Predicted-mean RT differences (EZ)
######################################
# Default: x_param = "drift_mean", T = 40, beta = 0
plot_RTdiff_predNested(main_dir = main_dir,output_dir = output_dir,x_param = "drift_mean")
plot_RTdiff_predNested(main_dir = main_dir,output_dir = output_dir,x_param = "drift_mean", show_IQrange = TRUE)


# Alternate x-axis parameter
plot_RTdiff_predNested(
  main_dir = main_dir,
  output_dir = output_dir,
  x_param = "bound_mean"
)

# Alternate beta level
plot_RTdiff_predNested(
  main_dir = main_dir,
  output_dir = output_dir,
  x_param = "drift_mean",
  beta_level_select = 0.2
)

########################################################################################
# F U L L     S I M U L A T I O N     S T U D Y  #######################################
########################################################################################
output_dir <- here("output", "figures", "RTdiff-grid", "full-simulation")

label <- "full"
main_dir <- here("output", "RData", "full-simulation", "wide-parameters")
plot_RTdiff_full(main_dir, output_dir, point_alpha = 0.7, x_param = "bound_mean", point_cex = 0.5, custom_title_label = label)
plot_RTdiff_full(main_dir, output_dir, point_alpha = 0.7, x_param = "drift_mean", point_cex = 0.5, custom_title_label = label)

label <- "restricted"
main_dir <- here("output", "RData", "full-simulation", "restricted-parameters")
plot_RTdiff_full(main_dir, output_dir, point_alpha = 0.7, x_param = "bound_mean", point_cex = 0.5, custom_title_label = label)
plot_RTdiff_full(main_dir, output_dir, point_alpha = 0.7, x_param = "drift_mean", point_cex = 0.5, custom_title_label = label)