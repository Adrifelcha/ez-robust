##########################################################################
# This custom function creates a grid of plots showing the 
# meanRT - medianRT difference across values of a chosen true parameter
##########################################################################
# Input:
# main_dir: The path to the folder containing the simulation results
# output_dir: The path to the folder where the grid of plots will be saved.
# y_range: The range of the y-axis (meanRT - medianRT difference).
# x_range: The range of the x-axis (true parameter value).
# point_alpha: The transparency of the points.
# true_param: The true parameter to plot against (bound_mean, drift_mean, nondt_mean, betaweight).
##########################################################################

plot_RTdiff_by_param <- function(main_dir, output_dir, y_range = NULL, x_range = NULL,
                                 point_alpha = 1, true_param = "bound_mean", point_cex = 0.5) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create output directory, if it doesn't exist
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Infer conditions from folder names in main_dir
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  all_folders <- list.dirs(main_dir, full.names = FALSE, recursive = FALSE)
  # Separate into clean and contaminated based on folder name suffix
  clean_conditions <- sort(all_folders[grepl("_clean$", all_folders)])
  contaminated_conditions <- sort(all_folders[grepl("_contaminated$", all_folders)])
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Load first RData file to get simulation settings
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  first_condition_path <- file.path(main_dir, clean_conditions[1])
  all_files <- list.files(first_condition_path, pattern = "\\.RData$", full.names = TRUE)
  first_file <- all_files[1]
  load(first_file)
  p_levels <- sort(simStudy_Beta$settings_summary$participant_levels)
  t_levels <- sort(simStudy_Beta$settings_summary$trial_levels)
  total_t_levels <- length(t_levels)
  total_p_levels <- length(p_levels)
  total_iterations <- total_t_levels * total_p_levels  

  cat("\n============================================================\n")
  cat("Creating RT difference vs", true_param, "grid plot\n")
  cat("Grid dimensions:", length(t_levels), "rows x 4 columns\n")
  cat("Trial levels:", paste(t_levels, collapse = ", "), "\n")  
  cat("============================================================\n\n")
  cat("Collecting data for all cells...\n")
    
  # Initialize overall progress
  overall_progress_bar <- txtProgressBar(min = 0, max = total_iterations, 
                                         style = 3, width = 50, char = "=")
  iteration_count <- 0  
  all_plot_data <- list()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Collect data for all cells
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # We loop over levels of trials-per-condition (T)
  for(t_idx in seq_along(t_levels)) {
      # Identify this trial level
      t_level <- t_levels[t_idx]
      cell_key <- paste("T", t_level, sep = "_")

      # Initialize storage for this T level
      all_plot_data[[cell_key]] <- list()      
      beta0_clean_data <- list()
      beta0_contaminated_data <- list()
      beta04_clean_data <- list()
      beta04_contaminated_data <- list()

      cat("\n  Processing trial level", t_idx, "of", total_t_levels, ": T =", t_level, "\n")
      
      # Loop through all P levels to merge data
      for(p_idx in seq_along(p_levels)) {
          # Identify this participant level
          p_level <- p_levels[p_idx]          
          
          # Update progress bar
          iteration_count <- iteration_count + 1
          setTxtProgressBar(overall_progress_bar, iteration_count)
          
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Clean conditions
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
          # Load data from the first "Clean" condition
          # (The only difference between EZ and Robust lies on the estimates)
          # (both conditions compute standard and robust statistics from the same data)
          condition = clean_conditions[1]
          # Identify the RData file for this condition and participant level
          pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
          file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]

          # If the file exists...
          if(!is.na(file_path) && file.exists(file_path)) {
              # Load the data
              load(file_path)              
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              # Beta = 0 (Column 1)
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
              beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0)
              # Check that this beta value was used in the simulation
              if (length(beta_indices) > 0) {
                # Extract meanRT, medianRT, and compute the difference
                meanRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "meanRT"]))
                medianRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "medianRT"]))            
                rt_diff <- meanRT_vals - medianRT_vals
                # Extract the true parameter value for this beta value
                param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
                  # Each population-level parameter repeats across participants per condition
                  param_vals <- rep(param_vals, each = p_level*2)
                            
                # Store the data in the beta0_clean_data list                                          
                aqui0_clean <- length(beta0_clean_data) + 1
                beta0_clean_data[[aqui0_clean]] <- data.frame(rt_diff = rt_diff,
                                                      param_value = param_vals,
                                                      stringsAsFactors = FALSE)
              }
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              # Beta = 0.4 (Column 3)
              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
              # Look for simulations with beta = 0.4              
              beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0.4)
              # Check 
              if (length(beta_indices) > 0) {
                # Extract the X values for this beta value
                x_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "X"]))                  
                # Extract the meanRT and medianRT values for this beta value
                meanRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "meanRT"]))
                medianRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "medianRT"]))
                rt_diff <- meanRT_vals - medianRT_vals
                # Extract the true parameter value for this beta value
                param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
                  # Each population-level parameter repeats across participants per condition
                  param_vals <- rep(param_vals, each = p_level*2)
                
                # Keep only the data from condition 1 (X == 1)
                keep <- which(x_vals == 1)
                rt_diff <- rt_diff[keep]
                param_vals <- param_vals[keep]
                               
                # Store the data in the beta04_clean_data list
                aqui04_clean <- length(beta04_clean_data) + 1
                beta04_clean_data[[aqui04_clean]] <- data.frame(rt_diff = rt_diff,
                                                     param_value = param_vals,
                                                     stringsAsFactors = FALSE)                                  
              }
          } # Close clean condition
          
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          # Contaminated conditions
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
          # Load data from the first "Contaminated" condition
          # (The only difference between EZ and Robust lies on the estimates)
          # (both conditions compute standard and robust statistics from the same data)
          condition = contaminated_conditions[1]
          # Identify the RData file for this condition and participant level
          pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
          file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
          
          # If the file exists...
          if(!is.na(file_path) && file.exists(file_path)) {
            # Load the data
            load(file_path)          

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Beta = 0 (Column 2)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
            # Look for simulations with beta = 0
            beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0)
            # Check that this beta value was used in the simulation
            if(length(beta_indices) > 0) {
              # Extract meanRT, medianRT, and compute the difference              
              meanRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "meanRT"]))
              medianRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "medianRT"]))
              rt_diff <- meanRT_vals - medianRT_vals
              # Extract the true parameter value for this beta value
              param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
                # Each population-level parameter repeats across participants per condition
                param_vals <- rep(param_vals, each = p_level*2)
              
              # Store the data in the beta0_contaminated_data list
              aqui0_cont <- length(beta0_contaminated_data) + 1
              beta0_contaminated_data[[aqui0_cont]] <- data.frame(rt_diff = rt_diff,
                                                           param_value = param_vals,
                                                           stringsAsFactors = FALSE)
            }
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Beta = 0.4 (Column 4)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
            # Look for simulations with beta = 0.4              
            beta_indices <- which(simStudy_Beta$true[, "betaweight"] == 0.4)
            # Check 
            if (length(beta_indices) > 0) {
              # Extract the X values for this beta value
              x_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "X"]))                  
              # Extract the meanRT and medianRT values for this beta value
              meanRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "meanRT"]))
              medianRT_vals <- as.numeric(unlist(simStudy_Beta$summStats[beta_indices, "medianRT"]))
              rt_diff <- meanRT_vals - medianRT_vals
              # Extract the true parameter value for this beta value
              param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, true_param]))
                # Each population-level parameter repeats across participants per condition
                param_vals <- rep(param_vals, each = p_level*2)
              
              # Keep only the data from condition 1 (X == 1)
              keep <- which(x_vals == 1)
              rt_diff <- rt_diff[keep]
              param_vals <- param_vals[keep]
                              
              # Store the data in the beta04_clean_data list
              aqui04_cont <- length(beta04_contaminated_data) + 1
              beta04_contaminated_data[[aqui04_cont]] <- data.frame(rt_diff = rt_diff,
                                                              param_value = param_vals,
                                                              stringsAsFactors = FALSE)                                  
            }
          } # Close contaminated condition
      } # end P loop
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Combine all participant size (P) levels for this T level
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      cat("    Combining data for T =", t_level, "...\n")      
      # Initialize empty data frame to use if either beta = 0 or beta = 0.4 was not used
      empty_data <- data.frame(rt_diff = numeric(0), param_value = numeric(0))
      # Check if beta = 0 was used     
      if (length(beta0_clean_data) > 0) {
        all_plot_data[[cell_key]]$beta0_clean <- do.call(rbind, beta0_clean_data)
        all_plot_data[[cell_key]]$beta0_contaminated <- do.call(rbind, beta0_contaminated_data)
      } else {       
        all_plot_data[[cell_key]]$beta0_clean <- empty_data
        all_plot_data[[cell_key]]$beta0_contaminated <- empty_data
      }      
      # Check if beta = 0.4 was used     
      if (length(beta04_clean_data) > 0) {
        all_plot_data[[cell_key]]$beta04_clean <- do.call(rbind, beta04_clean_data)
        all_plot_data[[cell_key]]$beta04_contaminated <- do.call(rbind, beta04_contaminated_data)
      } else {
        all_plot_data[[cell_key]]$beta04_clean <- empty_data
        all_plot_data[[cell_key]]$beta04_contaminated <- empty_data
      }
      cat("    ✓ Completed T =", t_level, "\n")
  }
  
  # Close the overall progress bar
  close(overall_progress_bar)
  cat("\nData collection complete!\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define plotting space
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Identify the range of the y-axis (RT difference)
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
    y_range <- range(all_rt_diff, na.rm = TRUE)
    y_range <- c(y_range[1] - diff(y_range) * 0.05, y_range[2] + diff(y_range) * 0.05)
    cat("Y-axis range:", y_range[1], "to", y_range[2], "\n")
  }
  
  # Identify the range of the x-axis (true parameter value)
  if (is.null(x_range)) {
    cat("Determining x-axis range...\n")
    all_param_values <- numeric(0)
    for (cell_key in names(all_plot_data)) {
      for (data_name in names(all_plot_data[[cell_key]])) {
        if (nrow(all_plot_data[[cell_key]][[data_name]]) > 0) {
          all_param_values <- c(all_param_values, all_plot_data[[cell_key]][[data_name]]$param_value)
        }
      }
    }    
    x_range <- range(all_param_values, na.rm = TRUE)
    x_range <- c(x_range[1] - diff(x_range) * 0.05, x_range[2] + diff(x_range) * 0.05)    
    cat("X-axis range:", x_range[1], "to", x_range[2], "\n\n")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create the plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("\n")
  cat("============================================================\n")
  cat("Creating plot...\n")
  cat("============================================================\n")
  
  output_filename <- paste0("RTdiff_by_", true_param, "_grid.pdf")
  pdf(file.path(output_dir, output_filename), width = 16, height = 14)
  
  # Setup plot layout: 5 rows x 4 columns (4 plots + 1 gap column between groups)  
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
  for(t_level in rev(t_levels)) {  # Reverse to plot high to low
    cell_key <- paste("T", t_level, sep = "_")
    
    # Column 1: Beta = 0, Clean
    plot_data <- all_plot_data[[cell_key]]$beta0_clean
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha, 
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = TRUE, point_cex = point_cex)
    
    # Column 2: Beta = 0, Contaminated  
    plot_data <- all_plot_data[[cell_key]]$beta0_contaminated
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha,
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = FALSE, point_cex = point_cex)
            
    # Column 3: Beta = 0.4, Clean
    plot_data <- all_plot_data[[cell_key]]$beta04_clean
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha,
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = FALSE, point_cex = point_cex)
    
    # Column 4: Beta = 0.4, Contaminated
    plot_data <- all_plot_data[[cell_key]]$beta04_contaminated
    plot_cell_scatter(plot_data, x_range, y_range, point_alpha,
                     show_x_axis = (t_level == min(t_levels)), 
                     show_y_axis = FALSE, point_cex = point_cex)
  }
  
  # Add column labels
  mtext(expression(paste("MeanRT - MedianRT")), side = 2, line = 3, cex = 2.7, outer = TRUE)
  mtext(x_axis_label(true_param), side = 1, line = 4.5, cex = x_axis_label_cex(true_param), outer = TRUE)
  
  # Add group labels (Beta = 0 and Beta = 0.4)
  mtext(expression(paste(beta, " = 0.0")), side = 3, line = 0.5, at = 0.25, cex = 2, outer = TRUE)
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
                              show_x_axis = FALSE, show_y_axis = FALSE, point_cex = 0.5) {
  
  # Create empty plot
  plot(NA, NA, xlim = x_range, ylim = y_range, 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  
  # Add points if data exists
  if (nrow(plot_data) > 0) {
    # Use adjustcolor or rgb for transparency
    point_color <- rgb(0, 0, 0, alpha = point_alpha)  # Black with transparency
    
    points(plot_data$param_value, plot_data$rt_diff, 
           col = point_color, pch = 16, cex = point_cex)
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


# Helper function to define the x-axis label to be printed on the margin
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_axis_label <- function(true_param){
  if(true_param == "bound_mean"){
    return(expression(paste(mu[alpha])))
  } else if(true_param == "nondt_mean"){
    return(expression(paste(mu[tau[0]])))
  } else if(true_param == "drift_mean"){
    return(expression(paste("Population-level intercept (", mu[nu], ")")))
  } else {
    return(expression(paste("True parameter:", true_param)))
  }
}

# Helper function to define the size of the margin text for the x-axis label
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_axis_label_cex <- function(true_param){
  if(true_param == "drift_mean"){  return(2.7)
  } else {  return(3.5) }
}