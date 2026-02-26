######################################
# Load libraries and functions
######################################
library(here)
source(here("src", "load_allFunctions.R"))

load_allCustomFunctions()

# Establish output directory
output_dir <- here("output", "figures", "AUC-grid")

###########################################
# Plot AUC grid for full simulation study
###########################################
label <- "full"
main_dir <- here("output", "RData", "full-simulation")

plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label) 
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)

##################################################
# Plot AUC grid for restricted simulation study
##################################################
label <- "restricted"
main_dir <- here("output", "RData", "restricted-simulation")

plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "beta", custom_title_label = label) 
plot_AUCgrid(main_dir = main_dir, output_dir = output_dir, plot_by = "condition", custom_title_label = label)