library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

# Define the four subfolders
subfolders <- c("EZ_contaminated", "EZ_clean", "EZRobust_clean", "EZRobust_contaminated")

# Loop through each subfolder
for (subfolder in subfolders) {
    # Get the full path to the subfolder
    folder_path <- here("output", "RData_simStudy_results", subfolder)
    # Create output directory based on subfolder
    output_dir <- here("output", "figures_BayesFactors", subfolder)
    # Create the grid of plots
    y_range <- switch(subfolder, "EZ_clean" = c(-2.5,4),
                                 "EZ_contaminated" = c(-3,4),
                                 "EZRobust_clean" = c(-2.5,4),
                                 "EZRobust_contaminated" = c(-2.5,4))
    cat("Creating grid of plots for", subfolder, "\n")
    plot_grid(subfolder_path = folder_path, output_dir = output_dir,
              plot_type = "BF")
    
    # Get all RData files in the folder
    rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)
    
    cat("Creating individual plots for", subfolder, "\n")
    # Loop through each RData file
    for (resultsFile in rdata_files) {
          # Extract filename for output naming
          filename <- basename(resultsFile)
          # Print progress
          cat("Processing:", filename, "\n")
          # Create the plot
          plot_BF_jitter(resultsFile, show_x_axis = TRUE,
                        show_y_axis = TRUE, output_dir = output_dir)
    }
}