# Main function to create a grid plot showing meanRT - medianRT difference vs a chosen true parameter
# Creates a 5 (trial levels) x 4 (beta groups and data types) grid
plot_RTdiff_by_drift <- function(main_dir, output_dir, y_range = NULL, x_range = NULL, 
                                  point_alpha = 1, true_param = "bound_mean") {
  
  # Create output directory if it doesn't exist  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define conditions (clean vs contaminated for each beta group)
  clean_conditions <- c("EZ_clean", "EZRobust_clean")
  contaminated_conditions <- c("EZ_contaminated", "EZRobust_contaminated")
  
  # Infer T levels from the first subfolder
  first_condition_path <- file.path(main_dir, clean_conditions[1])
  all_files <- list.files(first_condition_path, pattern = "\\.RData$", full.names = TRUE)
  rdata_files <- all_files[grepl("_P\\d+T\\d+_", basename(all_files))]
  
  # Get T values (trial levels) - will be the rows
  filenames <- basename(rdata_files)
  t_values <- as.numeric(sub(".*T(\\d+)_.*", "\\1", filenames))
  t_levels <- sort(unique(t_values))
  
  # Get P values (participant levels) - will be merged across all P
  p_values <- as.numeric(sub(".*_P(\\d+)T.*", "\\1", filenames))
  p_levels <- sort(unique(p_values))
  
  cat("\n============================================================\n")
  cat("Creating RT difference vs", true_param, "grid plot\n")
  cat("Grid dimensions:", length(t_levels), "rows x 4 columns\n")
  cat("Trial levels:", paste(t_levels, collapse = ", "), "\n")  
  cat("============================================================\n\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Collect all data across P levels for each T level
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Collecting data for all cells...\n")
  total_t_levels <- length(t_levels)
  total_p_levels <- length(p_levels)
  total_iterations <- total_t_levels * total_p_levels
  
  # Initialize overall progress
  overall_progress_bar <- txtProgressBar(min = 0, max = total_iterations, 
                                         style = 3, width = 50, char = "=")
  iteration_count <- 0
  
  all_plot_data <- list()
  
  for (t_idx in seq_along(t_levels)) {
    t_level <- t_levels[t_idx]
    cell_key <- paste("T", t_level, sep = "_")
    all_plot_data[[cell_key]] <- list()
    
    cat("\n  Processing trial level", t_idx, "of", total_t_levels, ": T =", t_level, "\n")
    
    # Columns 1-2: Beta = 0 (Clean, Contaminated)
    # Columns 3-4: Beta = 0.4 (Clean, Contaminated)
    
    # Initialize data collectors for this T level
    beta0_clean_data <- list()
    beta0_contaminated_data <- list()
    beta04_clean_data <- list()
    beta04_contaminated_data <- list()
    
    # Loop through all P levels to merge data
    for (p_idx in seq_along(p_levels)) {
      p_level <- p_levels[p_idx]
      iteration_count <- iteration_count + 1
      
      # Update progress bar
      setTxtProgressBar(overall_progress_bar, iteration_count)
      
      # Load Beta = 0 data (Clean conditions)
      for (condition in clean_conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path) && file.exists(file_path)) {
          load(file_path)
          
          # Check if summStats exists and has the right structure
          if (is.null(simStudy_Beta$summStats) || !is.matrix(simStudy_Beta$summStats)) {
            next  # Skip this file if summStats is not available
          }
          
          # Filter for beta = 0
          beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0)
          if (length(beta_indices) > 0) {
            # Extract meanRT, medianRT, and chosen true parameter for beta = 0
            # Ensure we extract numeric values - handle list matrices
            meanRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "meanRT"])
            medianRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "medianRT"])
            meanRT_vals <- as.numeric(meanRT_vals)
            medianRT_vals <- as.numeric(medianRT_vals)
            # 'true' is stored as a list-matrix; unlist before coercion
            param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
            
            # Calculate RT difference
            rt_diff <- meanRT_vals - medianRT_vals
            
            # Store (only if we have valid data)
            valid <- !is.na(rt_diff) & !is.na(param_vals)
            if (any(valid)) {
              beta0_clean_data[[length(beta0_clean_data) + 1]] <- data.frame(
                rt_diff = rt_diff[valid],
                param_value = param_vals[valid],
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
      
      # Load Beta = 0 data (Contaminated conditions)
      for (condition in contaminated_conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path) && file.exists(file_path)) {
          load(file_path)
          
          # Check if summStats exists and has the right structure
          if (is.null(simStudy_Beta$summStats) || !is.matrix(simStudy_Beta$summStats)) {
            next  # Skip this file if summStats is not available
          }
          
          # Filter for beta = 0
          beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0)
          if (length(beta_indices) > 0) {
            # Ensure we extract numeric values - handle list matrices
            meanRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "meanRT"])
            medianRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "medianRT"])
            meanRT_vals <- as.numeric(meanRT_vals)
            medianRT_vals <- as.numeric(medianRT_vals)
            # 'true' is stored as a list-matrix; unlist before coercion
            param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
            
            rt_diff <- meanRT_vals - medianRT_vals
            valid <- !is.na(rt_diff) & !is.na(param_vals)
            if (any(valid)) {
              beta0_contaminated_data[[length(beta0_contaminated_data) + 1]] <- data.frame(
                rt_diff = rt_diff[valid],
                param_value = param_vals[valid],
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
      
      # Load Beta = 0.4 data (Clean conditions)
      # NOTE: Only use condition 1 data (where beta applies, X == 1)
      for (condition in clean_conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path) && file.exists(file_path)) {
          load(file_path)
          
          # Check if summStats exists and has the right structure
          if (is.null(simStudy_Beta$summStats) || !is.matrix(simStudy_Beta$summStats)) {
            next  # Skip this file if summStats is not available
          }
          
          # Filter for beta = 0.4 AND X = 1 (condition 1) in a single step
          beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0.4)
          if (length(beta_indices) > 0) {
            # Further filter for X = 1 if X column exists in summStats
            if ("X" %in% colnames(simStudy_Beta$summStats)) {
              x_vals <- unlist(simStudy_Beta$summStats[beta_indices, "X"])
              x_vals <- as.numeric(x_vals)
              # Keep only indices where X == 1
              beta_indices <- beta_indices[which(x_vals == 1)]
              if (length(beta_indices) == 0) {
                next  # Skip if no condition 1 data
              }
            }
            
            # Extract values only for the filtered indices
            meanRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "meanRT"])
            medianRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "medianRT"])
            meanRT_vals <- as.numeric(meanRT_vals)
            medianRT_vals <- as.numeric(medianRT_vals)
            # 'true' is stored as a list-matrix; unlist before coercion
            param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
            
            rt_diff <- meanRT_vals - medianRT_vals
            valid <- !is.na(rt_diff) & !is.na(param_vals)
            if (any(valid)) {
              beta04_clean_data[[length(beta04_clean_data) + 1]] <- data.frame(
                rt_diff = rt_diff[valid],
                param_value = param_vals[valid],
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
      
      # Load Beta = 0.4 data (Contaminated conditions)
      # NOTE: Only use condition 1 data (where beta applies, X == 1)
      for (condition in contaminated_conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path) && file.exists(file_path)) {
          load(file_path)
          
          # Filter for beta = 0.4 AND X = 1 (condition 1) in a single step
          beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0.4)
          if (length(beta_indices) > 0) {
            # Further filter for X = 1 if X column exists in summStats
            if ("X" %in% colnames(simStudy_Beta$summStats)) {
              x_vals <- unlist(simStudy_Beta$summStats[beta_indices, "X"])
              x_vals <- as.numeric(x_vals)
              # Keep only indices where X == 1
              beta_indices <- beta_indices[which(x_vals == 1)]
              if (length(beta_indices) == 0) {
                next  # Skip if no condition 1 data
              }
            }
            
            # Extract values only for the filtered indices
            meanRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "meanRT"])
            medianRT_vals <- unlist(simStudy_Beta$summStats[beta_indices, "medianRT"])
            meanRT_vals <- as.numeric(meanRT_vals)
            medianRT_vals <- as.numeric(medianRT_vals)
            # 'true' is stored as a list-matrix; unlist before coercion
            param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
            
            rt_diff <- meanRT_vals - medianRT_vals
            valid <- !is.na(rt_diff) & !is.na(param_vals)
            if (any(valid)) {
              beta04_contaminated_data[[length(beta04_contaminated_data) + 1]] <- data.frame(
                rt_diff = rt_diff[valid],
                param_value = param_vals[valid],
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    } # end P loop
    
    # Combine all data for this T level
    cat("    Combining data for T =", t_level, "...\n")
    if (length(beta0_clean_data) > 0) {
      all_plot_data[[cell_key]]$beta0_clean <- do.call(rbind, beta0_clean_data)
    } else {
      all_plot_data[[cell_key]]$beta0_clean <- data.frame(rt_diff = numeric(0), param_value = numeric(0))
    }
    
    if (length(beta0_contaminated_data) > 0) {
      all_plot_data[[cell_key]]$beta0_contaminated <- do.call(rbind, beta0_contaminated_data)
    } else {
      all_plot_data[[cell_key]]$beta0_contaminated <- data.frame(rt_diff = numeric(0), param_value = numeric(0))
    }
    
    if (length(beta04_clean_data) > 0) {
      all_plot_data[[cell_key]]$beta04_clean <- do.call(rbind, beta04_clean_data)
    } else {
      all_plot_data[[cell_key]]$beta04_clean <- data.frame(rt_diff = numeric(0), param_value = numeric(0))
    }
    
    if (length(beta04_contaminated_data) > 0) {
      all_plot_data[[cell_key]]$beta04_contaminated <- do.call(rbind, beta04_contaminated_data)
    } else {
      all_plot_data[[cell_key]]$beta04_contaminated <- data.frame(rt_diff = numeric(0), param_value = numeric(0))
    }
    
    cat("    ✓ Completed T =", t_level, "\n")
  }
  
  # Close the overall progress bar
  close(overall_progress_bar)
  cat("\nData collection complete!\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Determine axis ranges if not provided
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (is.null(y_range)) {
    cat("\nDetermining y-axis range...\n")
    all_rt_diff <- numeric(0)
    for (cell_key in names(all_plot_data)) {
      for (data_name in names(all_plot_data[[cell_key]])) {
        if (nrow(all_plot_data[[cell_key]][[data_name]]) > 0) {
          all_rt_diff <- c(all_rt_diff, all_plot_data[[cell_key]][[data_name]]$rt_diff)
        }
      }
    }
    if (length(all_rt_diff) > 0) {
      y_range <- range(all_rt_diff, na.rm = TRUE)
      y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
    } else {
      y_range <- c(-0.1, 0.1)
    }
    cat("Y-axis range:", y_range[1], "to", y_range[2], "\n")
  }
  
  if (is.null(x_range)) {
    cat("Determining x-axis range...\n")
    all_drift <- numeric(0)
    for (cell_key in names(all_plot_data)) {
      for (data_name in names(all_plot_data[[cell_key]])) {
        if (nrow(all_plot_data[[cell_key]][[data_name]]) > 0) {
          all_drift <- c(all_drift, all_plot_data[[cell_key]][[data_name]]$param_value)
        }
      }
    }
    if (length(all_drift) > 0) {
      x_range <- range(all_drift, na.rm = TRUE)
      x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)
    } else {
      x_range <- c(-3, 3)
    }
    cat("X-axis range:", x_range[1], "to", x_range[2], "\n\n")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create the plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("\n")
  cat("============================================================\n")
  cat("Creating plot...\n")
  cat("============================================================\n")
  
  output_filename <- "RTdiff_by_drift_grid.pdf"
  pdf(file.path(output_dir, output_filename), width = 16, height = 14)
  
  # Setup plot layout: 5 rows x 5 columns (4 plots + 1 gap column between groups)
  # Use layout() to create custom spacing with gap between columns 2 and 3
  n_rows <- length(t_levels)
  
  # Create layout matrix with 5 columns: plot, plot, gap, plot, plot
  layout_matrix <- matrix(c(1, 2, 0, 3, 4,
                            5, 6, 0, 7, 8,
                            9, 10, 0, 11, 12,
                            13, 14, 0, 15, 16,
                            17, 18, 0, 19, 20), 
                          nrow = n_rows, ncol = 5, byrow = TRUE)
  
  # Use layout with custom widths: equal for plots, 0.3 for gap
  layout(layout_matrix, widths = c(1, 1, 0.5, 1, 1))
  
  # Set margins
  par(oma = c(6, 7, 3, 3), mar = c(1, 1.5, 0.5, 0.5))
  
  # Plot each cell (rows: T levels, columns: 4 data groups)
  for (t_level in rev(t_levels)) {  # Reverse to plot high to low
    cell_key <- paste("T", t_level, sep = "_")
    
    # Column 1: Beta = 0, Clean
    plot_data <- all_plot_data[[cell_key]]$beta0_clean
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha, 
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = TRUE)
    
    # Column 2: Beta = 0, Contaminated  
    plot_data <- all_plot_data[[cell_key]]$beta0_contaminated
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha,
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = FALSE)
    
    # Column 3: Gap (handled by layout with width 0.5)
    
    # Column 4: Beta = 0.4, Clean
    plot_data <- all_plot_data[[cell_key]]$beta04_clean
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha,
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = FALSE)
    
    # Column 5: Beta = 0.4, Contaminated
    plot_data <- all_plot_data[[cell_key]]$beta04_contaminated
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha,
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = FALSE)
  }
  
  # Add column labels
  mtext(expression(paste("MeanRT - MedianRT")), side = 2, line = 3, cex = 2.5, outer = TRUE)
  mtext(paste("True parameter:", true_param), side = 1, line = 4.5, cex = 2.3, outer = TRUE)
  
  # Add group labels (Beta = 0 and Beta = 0.4)
  mtext(expression(paste(beta, " = 0")), side = 3, line = 0.5, at = 0.25, cex = 2, outer = TRUE)
  mtext(expression(paste(beta, " = 0.4")), side = 3, line = 0.5, at = 0.75, cex = 2, outer = TRUE)
  
  # Add data type labels
  line_topMargin_2 <- 0
  mtext("Clean", side = 3, line = line_topMargin_2, at = 0.125, cex = 1.5, outer = TRUE)
  mtext("Contaminated", side = 3, line = line_topMargin_2, at = 0.375, cex = 1.5, outer = TRUE)
  mtext("Clean", side = 3, line = line_topMargin_2, at = 0.625, cex = 1.5, outer = TRUE)
  mtext("Contaminated", side = 3, line = line_topMargin_2, at = 0.875, cex = 1.5, outer = TRUE)
  
  dev.off()
  
  cat("\n")
  cat("============================================================\n")
  cat("Plot saved to:", file.path(output_dir, output_filename), "\n")
  cat("Plotting complete!\n")
  cat("============================================================\n\n")
}

###################################################################
# Helper function to plot a single cell with scatter points
###################################################################
plot_cell_scatter <- function(plot_data, x_range, y_range, point_alpha, 
                              show_x_axis = FALSE, show_y_axis = FALSE) {
  
  # Create empty plot
  plot(NA, NA, xlim = x_range, ylim = y_range, 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  
  # Add points if data exists
  if (nrow(plot_data) > 0) {
    # Use adjustcolor or rgb for transparency
    point_color <- rgb(0, 0, 0, alpha = point_alpha)  # Black with transparency
    
    points(plot_data$param_value, plot_data$rt_diff, 
           col = point_color, pch = 16, cex = 0.5)
  }
  
  # Add axes if needed
  if (show_x_axis) {
    x_at <- pretty(x_range, n = 5)
    axis(1, at = x_at, labels = x_at, cex.axis = 2)
  }
  
  if (show_y_axis) {
    y_at <- pretty(y_range, n = 5)
    axis(2, at = y_at, labels = y_at, cex.axis = 2, las = 1)
  }
  
  # Add horizontal line at y = 0
  abline(h = 0, lty = 2, col = "gray50", lwd = 1)
}
