library(here)

# Define the four subfolders
subfolders <- c("EZ_clean", "EZ_contaminated", "EZRobust_clean", "EZRobust_contaminated")

# Loop through each subfolder
for (subfolder in subfolders) {
  # Get the full path to the subfolder
  folder_path <- here("output", "RData_simStudy_results", subfolder)
  # Create output directory based on subfolder
  output_dir <- here("output", "figures_BetaRecovery", subfolder)
  # Create the grid of plots
  y_range <- c(-1.2,1.2)
  
  cat("Creating grid of plots for", subfolder, "\n")
  plot_grid(subfolder_path = folder_path, output_dir = output_dir,
          plot_type = "Betas")

  # Get all RData files in the folder
  rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)
  
  # Loop through each RData file
  for (resultsFile in rdata_files) {
    # Extract filename for output naming
    filename <- basename(resultsFile)
    
    # Print progress
    cat("Processing:", filename, "\n")
    
    # Create the plot
    plot_cellBetaEstimates(resultsFile, show_frequency_bars = TRUE,
                          show_x_axis = TRUE, output_dir = output_dir)
    
  }
}
