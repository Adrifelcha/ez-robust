library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

# Define the four subfolders
subfolders <- c("EZ_clean", "EZ_contaminated", "EZRobust_clean", "EZRobust_contaminated")

# Loop through each subfolder
for(subfolder in subfolders){
    # Get the full path to the subfolder
    folder_path <- here("output", "simStudy_results", subfolder)   
    # Create output directory based on subfolder
    output_dir <- here("output", "figures_ROC", subfolder)
    # Create the grid of plots
    cat("Creating grid of plots for", subfolder, "\n")
    plot_grid(subfolder_path = folder_path, output_dir = output_dir,
              plot_type = "ROC")

    # Get all RData files in the folder
    rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)
    cat("Creating individual plots for", subfolder, "\n")
    # Loop through each RData file
      for(resultsFile in rdata_files){
          # Extract filename for output naming
          filename <- basename(resultsFile)          
          # Create output directory based on subfolder
          output_dir <- here("output", "figures_ROC", subfolder)
          # Print progress
          cat("Processing:", filename, "\n")
          # Create the ROC plot
          plot_cellROC(resultsFile, epsilon = 0.05, levelsM = 1000,
                      show_legend = TRUE, output_dir = output_dir)
    }
}

cat("All ROC plots have been generated!\n")