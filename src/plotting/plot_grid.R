##########################################################################
# This custom function creates a grid of plots showing the results
# across all simulation conditions (trial and participant size.)
##########################################################################
# Input:
# subfolder_path: The path to the folder containing the simulation results
#                 .RData files.
# output_dir: The path to the folder where the grid of plots will be saved.
# plot_type: The type of plot to generate (ROC, BF, or Betas).
##########################################################################

plot_grid <- function(subfolder_path, output_dir, plot_type = NULL, y_range = NULL) {
    # List all .RData files in the folder
      all_files <- list.files(subfolder_path, pattern = "\\.RData$", full.names = TRUE)
    # Filter and keep only files that match the expected naming scheme
      rdata_files <- all_files[grepl("_P\\d+T\\d+_", basename(all_files))]
    
    if(length(rdata_files) == 0){
      warning(paste("No .RData files with pattern _PxxTxx_ found in subfolder_path provided."))
      return(invisible(NULL))
    }

    # Get the P and T values from the filenames
    filenames <- basename(rdata_files)
    p_values <- as.numeric(sub(".*_P(\\d+)T.*", "\\1", filenames))
    t_values <- as.numeric(sub(".*T(\\d+)_.*", "\\1", filenames))
    # Get the unique P and T values, sorted ascending
    p_levels <- sort(unique(p_values))
    t_levels <- sort(unique(t_values))

    # Setup the plot layout with minimal space between plots
    par(mfrow = c(length(t_levels), length(p_levels)),
        oma = c(2, 1, 3, 3), # outer margins c(bottom, left, top, right)
        mai = c(0.15, 0.15, 0, 0)) # inner margins

    # Loop through each T level (rows)
    for(t_level in t_levels){
        # Loop through each P level (columns)
        for(p_level in p_levels){
            # Find the file for the current P and T level
            file_index <- which(p_values == p_level & t_values == t_level)
            
            # If the file exists, create the plot
            if(length(file_index) > 0){
                resultsFile <- rdata_files[file_index]

                cat("Plotting cell for", resultsFile, "\n")
                # Show X axis only for the bottom row
                if(plot_type == "ROC") {
                  plot_cellROC(resultsFile, bty = "o",
                               show_x_axis = (t_level == max(t_levels)),
                               show_y_axis = (p_level == min(p_levels)),
                               show_legend = (t_level == min(t_levels)) && (p_level == max(p_levels)))
                }else{if(plot_type == "BF") {
                        plot_BF_jitter(resultsFile,
                                        show_x_axis = (t_level == max(t_levels)),
                                        show_y_axis = (p_level == min(p_levels)),
                                        bty = "o", line = -0.75, y_range = y_range)
                      }else{if(plot_type == "Betas") {
                              plot_cellBetaEstimates(resultsFile,
                                                      show_x_axis = (t_level == max(t_levels)),
                                                      show_y_axis = (p_level == min(p_levels)),
                                                      bty = "o", y_range = y_range)
                            }else{
                              stop("Invalid plot_type. Please choose from 'ROC', 'BF', or 'Betas'.")
                            }
                      }
              }
            }else{
                    # If no file, plot an empty space
                    plot.new()
            }

            # Add Participant (column) labels on the top of the first row
            if (t_level == min(t_levels)) {
              mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 1.75, font = 2, las = 1)
            }

            # Add Trial (row) labels on the right of the last column, rotated
            if (p_level == max(p_levels)) {
              mtext(paste("T =", t_level), side = 4, line = 1.2, cex = 1.75, font = 2, las = 0)
            }
        }
    }

    # Save the plot
    output_filename <- paste0(basename(subfolder_path), "_", plot_type, "_grid.pdf")
    dev.copy(pdf, file.path(output_dir, output_filename), width = 12, height = 12)
    dev.off()
}


plot_grid(subfolder_path = folder_path, output_dir = output_dir,
          plot_type = "Betas", y_range = y_range)