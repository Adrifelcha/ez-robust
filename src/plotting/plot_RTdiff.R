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
# x_param: The true parameter to plot against (bound_mean, drift_mean, nondt_mean, betaweight).
##########################################################################
plot_RTdiff_full <- function(main_dir, output_dir, x_param = "bound_mean", custom_title_label = NULL,
                                 third_param = NULL, third_param_low = NULL, third_param_high = NULL,
                                 y_range = NULL, x_range = NULL, point_alpha = 1, point_cex = 0.5, colored = FALSE) {
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create output directory, if it doesn't exist
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
           
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Infer conditions from folder names in main_dir
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    condition_info <- read_conditions(main_dir = main_dir)
    all_folders <- condition_info$conditions
    # Separate into clean and contaminated based on folder name suffix
    clean_condition <- sort(all_folders[grepl("_clean$", all_folders)])[1]
    contaminated_condition <- sort(all_folders[grepl("_contaminated$", all_folders)])[1]
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load first RData file to get simulation settings
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    first_condition_path <- file.path(main_dir, clean_condition)
    all_files <- list.files(first_condition_path, pattern = "\\.RData$", full.names = TRUE)
    first_file <- all_files[1]
    load(first_file)
    # Identify trial levels
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    t_levels <- sort(simStudy_Beta$settings_summary$trial_levels)
    total_t_levels <- length(t_levels)
    # Identify participant levels
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    p_levels <- sort(simStudy_Beta$settings_summary$participant_levels)
    # The number of participants don't change the summary statistic difference
    if(is.null(third_param)){
      # By default, we will only use the first participant level
      p_levels <- p_levels[1]
    } else {
      # But we'll add the second participant level to the list if
      # we are conditioning on a third parameter (which reduces observations)
      p_levels <- p_levels[1:2]
    }  
    total_p_levels <- length(p_levels)
    total_iterations <- total_t_levels * total_p_levels  

    cat("\n============================================================\n")
    cat("Creating RT difference vs", x_param, "grid plot\n")
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

            # Identify the RData file for this condition and participant level
            pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
            file_path <- list.files(file.path(main_dir, clean_condition), pattern = pattern, full.names = TRUE)[1]

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
                  param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, x_param]))
                    # Each population-level parameter repeats across participants per condition
                    param_vals <- rep(param_vals, each = p_level*2)
                  
                  # Filter by third_param if specified
                  if (!is.null(third_param)) {
                    # Extract the values for the third parameter for this beta value
                    third_param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, third_param]))
                      # Each population-level parameter repeats across participants per condition
                      third_param_vals <- rep(third_param_vals, each = p_level*2)
                    keep <- third_param_vals >= third_param_low & third_param_vals <= third_param_high
                    rt_diff <- rt_diff[keep]
                    param_vals <- param_vals[keep]
                  }
                              
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
                  param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, x_param]))
                    # Each population-level parameter repeats across participants per condition
                    param_vals <- rep(param_vals, each = p_level*2)
                  
                  # Keep only the data from condition 1 (X == 1)
                  keep <- which(x_vals == 1)
                  rt_diff <- rt_diff[keep]
                  param_vals <- param_vals[keep]
                  
                  # Filter by third_param if specified
                  if (!is.null(third_param)) {
                    third_param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, third_param]))
                    third_param_vals <- rep(third_param_vals, each = p_level*2)
                    third_param_vals <- third_param_vals[keep]  # Apply same X==1 filter
                    keep_third <- third_param_vals >= third_param_low & third_param_vals <= third_param_high
                    rt_diff <- rt_diff[keep_third]
                    param_vals <- param_vals[keep_third]
                  }
                                
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
            
            # Identify the RData file for this condition and participant level
            pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
            file_path <- list.files(file.path(main_dir, contaminated_condition), pattern = pattern, full.names = TRUE)[1]
            
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
                param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, x_param]))
                  # Each population-level parameter repeats across participants per condition
                  param_vals <- rep(param_vals, each = p_level*2)
                
                # Filter by third_param if specified
                if (!is.null(third_param)) {
                  third_param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, third_param]))
                  third_param_vals <- rep(third_param_vals, each = p_level*2)
                  keep <- third_param_vals >= third_param_low & third_param_vals <= third_param_high
                  rt_diff <- rt_diff[keep]
                  param_vals <- param_vals[keep]
                }
                
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
                param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, x_param]))
                  # Each population-level parameter repeats across participants per condition
                  param_vals <- rep(param_vals, each = p_level*2)
                
                # Keep only the data from condition 1 (X == 1)
                keep <- which(x_vals == 1)
                rt_diff <- rt_diff[keep]
                param_vals <- param_vals[keep]
                
                # Filter by third_param if specified
                if (!is.null(third_param)) {
                  third_param_vals <- as.numeric(unlist(simStudy_Beta$true[beta_indices, third_param]))
                  third_param_vals <- rep(third_param_vals, each = p_level*2)
                  third_param_vals <- third_param_vals[keep]  # Apply same X==1 filter
                  keep_third <- third_param_vals >= third_param_low & third_param_vals <= third_param_high
                  rt_diff <- rt_diff[keep_third]
                  param_vals <- param_vals[keep_third]
                }
                                
                # Store the data in the beta04_contaminated_data list
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
    cat("Rendering panels...\n")
    
    output_filename <- file_name(x_param, third_param, custom_title_label)
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
    # Use layout with custom widths: equal for plots, narrow gap in the middle
    layout_widths <- c(1, 1, 0.3, 1, 1)
    layout(layout_matrix, widths = layout_widths)  
    
    # Positions (in outer coordinates) for row labels on the far right
    row_label_pos <- rev(seq(from = 1 / (2 * n_rows),
                            to   = 1 - 1 / (2 * n_rows),
                            length.out = n_rows))
    # Set margins (increase top outer margin to give space for beta/group labels)
    par(oma = c(6, 7, 5, 5), mar = c(1, 1.5, 0.5, 0.5))

    # Compute normalized x-positions for the centers of the four plot columns
    total_width <- sum(layout_widths)
    col1_center <- 0.5 / total_width          # center of first column
    col2_center <- 1.5 / total_width          # center of second column
    col3_center <- 2.8 / total_width          # center of third column (after gap)
    col4_center <- 3.8 / total_width          # center of fourth column
    
    # Plot each cell (rows: T levels, columns: 4 data groups)
    render_total <- length(t_levels) * 4
    render_progress_bar <- txtProgressBar(min = 0, max = render_total,
                                          style = 3, width = 50, char = "=")
    render_step <- 0
    row_index <- 0
    for(t_level in rev(t_levels)) {  # Reverse to plot high to low
      row_index <- row_index + 1
      cell_key <- paste("T", t_level, sep = "_")
      
      # Column 1: Beta = 0, Clean
      plot_data <- all_plot_data[[cell_key]]$beta0_clean
      make_scatterplot_full(plot_data, x_range, y_range, point_alpha, show_x_axis = (t_level == min(t_levels)), 
                      show_y_axis = TRUE, point_cex = point_cex)
      render_step <- render_step + 1
      setTxtProgressBar(render_progress_bar, render_step)
      
      # Column 2: Beta = 0, Contaminated  
      plot_data <- all_plot_data[[cell_key]]$beta0_contaminated
      make_scatterplot_full(plot_data, x_range, y_range, point_alpha, show_x_axis = (t_level == min(t_levels)), 
                        show_y_axis = FALSE, point_cex = point_cex)
      render_step <- render_step + 1
      setTxtProgressBar(render_progress_bar, render_step)
              
      # Column 3: Beta = 0.4, Clean
      plot_data <- all_plot_data[[cell_key]]$beta04_clean
      make_scatterplot_full(plot_data, x_range, y_range, point_alpha, show_x_axis = (t_level == min(t_levels)), 
                      show_y_axis = FALSE, point_cex = point_cex)
      render_step <- render_step + 1
      setTxtProgressBar(render_progress_bar, render_step)
      
      # Column 4: Beta = 0.4, Contaminated
      plot_data <- all_plot_data[[cell_key]]$beta04_contaminated
      make_scatterplot_full(plot_data, x_range, y_range, point_alpha, show_x_axis = (t_level == min(t_levels)), 
                      show_y_axis = FALSE, point_cex = point_cex)
      render_step <- render_step + 1
      setTxtProgressBar(render_progress_bar, render_step)
      
      # Add T-level label on the far right for this row
      mtext(text = paste("T =", t_level),
            side = 4, line = 2, at = row_label_pos[row_index],
            cex = 2.5, outer = TRUE)
    }
    close(render_progress_bar)
    cat("\nPanel rendering complete.\n")
    
    # Add column labels
    mtext(expression(paste("MeanRT - MedianRT")), side = 2, line = 3, cex = 2.7, outer = TRUE)
    mtext(x_axis_label(x_param), side = 1, line = 4.5, cex = x_axis_label_cex(x_param), outer = TRUE)
    
    # Add group labels (Beta = 0 and Beta = 0.4), centered over the two-column groups
    mtext(expression(paste(beta, " = 0.0")), side = 3, line = 1.5,
          at = (col1_center + col2_center) / 2, cex = 2, outer = TRUE)
    mtext(expression(paste(beta, " = 0.4")), side = 3, line = 1.5,
          at = (col3_center + col4_center) / 2, cex = 2, outer = TRUE)
    
    # Add data type labels, centered over each individual column
    line_topMargin_2 <- 0.3
    mtext("Clean", side = 3, line = line_topMargin_2, at = col1_center, cex = 1.5, outer = TRUE)
    mtext("Contaminated", side = 3, line = line_topMargin_2, at = col2_center, cex = 1.5, outer = TRUE)
    mtext("Clean", side = 3, line = line_topMargin_2, at = col3_center, cex = 1.5, outer = TRUE)
    mtext("Contaminated", side = 3, line = line_topMargin_2, at = col4_center, cex = 1.5, outer = TRUE)
    
    dev.off()
    
    cat("\n")    
    cat("Plot saved to:", file.path(output_dir, output_filename), "\n")        
}

########################################################################################################
# Helper function to plot single-cell scatterplots for the full simulation study
########################################################################################################
make_scatterplot_full <- function(plot_data, x_range, y_range, point_alpha, 
                              show_x_axis = FALSE, show_y_axis = FALSE, point_cex = 0.5) {
    # Create empty plot
    plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
    
    # Add points to the scatter plot
    # ------------------------------------------------------------------
    x_vals <- plot_data$param_value
    y_vals <- plot_data$rt_diff
    
    # Use adjustcolor or rgb for transparency for points
    point_color <- rgb(0, 0, 0, alpha = point_alpha)  # Black with opacity point_alpha  
    points(x_vals, y_vals, col = point_color, pch = 16, cex = point_cex)
    
    # Add a thick colored line showing a moving average over x_vals
    # ------------------------------------------------------------------
    # Sort by x to get a sensible curve
    ord <- order(x_vals)
    x_sorted <- x_vals[ord]
    y_sorted <- y_vals[ord]
    n <- length(x_sorted)
      
    # Define window width as a proportion of the x-range
    x_span <- max(x_sorted) - min(x_sorted)  
    window_width <- x_span * 0.10  # 10% of the x-range
      
    smoothed_y <- rep(NA, n)
    for (i in 1:n){
      in_window <- (x_sorted >= x_sorted[i] - window_width/2) &
                    (x_sorted <= x_sorted[i] + window_width/2)
      smoothed_y[i] <- mean(y_sorted[in_window])    
    }
    
    lines(x_sorted, smoothed_y, col = "darkred", lwd = 3)
        
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
x_axis_label <- function(x_param){
  if(x_param == "bound_mean"){
    return(expression(paste(mu[alpha])))
  } else if(x_param == "nondt_mean"){
    return(expression(paste(mu[tau[0]])))
  } else if(x_param == "drift_mean"){
    return(expression(paste("Population-level intercept (", mu[nu], ")")))
  } else {
    return(expression(paste("True parameter:", x_param)))
  }
}

# Helper function to define the size of the margin text for the x-axis label
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_axis_label_cex <- function(x_param){
  if(x_param == "drift_mean"){  return(2.7)
  } else {  return(3.5) }
}

# Helper function to lower and upper bounds for the background box per x parameter
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
box_bounds <- function(x_param){
  if(x_param == "drift_mean"){  return(c(-1,1))
  }else{if(x_param == "bound_mean"){  return(c(3.5,4.5))
  }else{if(x_param == "nondt_mean"){  return(c(0.1,0.4))
  }else{return(NULL)}}}
}


file_name <- function(x_param, third_param = NULL, custom_title_label = NULL){
  if(is.null(third_param)){
    if(is.null(custom_title_label)){
      return(paste0("RTdiff_by_", x_param, ".pdf"))
    }else{
      return(paste0("RTdiff_by_", x_param, "_", custom_title_label, ".pdf"))
    }
  }else{
    if(is.null(custom_title_label)){
      return(paste0("RTdiff_by_", x_param, "_fixed_", third_param, ".pdf"))
    }else{
      return(paste0("RTdiff_by_", x_param, "_fixed_", third_param, "_", custom_title_label, ".pdf"))
    }
  }
}

##########################################################################
# This custom function creates a nested-by-parameter plot showing the 
# absolute meanRT - medianRT difference across values of a chosen true parameter
##########################################################################
# Input:
# main_dir: The path to the folder containing the simulation results
# output_dir: The path to the folder where the grid of plots will be saved.
# y_range: The range of the y-axis (absolute meanRT - medianRT difference).
# x_range: The range of the x-axis (true parameter value).
# x_param: The true parameter to plot against (drift_mean, bound_mean).
# t_level_select: The T level to plot (e.g., 40).
# beta_level_select: The beta level to plot (e.g., 0.0).
# custom_title_label: A custom title label for the plot.
# show_IQrange: Whether to show the IQ range (e.g., TRUE).
##########################################################################
plot_RTdiff_nested <- function(main_dir, output_dir, t_level_select = 40,
                               beta_level_select = 0, x_param = "drift_mean",
                               custom_title_label = NULL, show_IQrange = FALSE) {
  if (!x_param %in% c("drift_mean", "bound_mean")) {
    stop("Invalid x_param. Use 'drift_mean' or 'bound_mean'.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  condition_info <- read_conditions(main_dir = main_dir)
  if (!identical(condition_info$layout, "nested")) {
    stop("plot_RTdiff_nested_by_param() requires nested simulation folders.")
  }
  
  conditions <- condition_info$conditions
  clean_condition <- sort(conditions[grepl("_clean$", conditions)])[1]
  contaminated_condition <- sort(conditions[grepl("_contaminated$", conditions)])[1]
  if (is.na(clean_condition) || is.na(contaminated_condition)) {
    stop("Could not infer clean/contaminated condition folders.")
  }
  
  if (identical(x_param, "drift_mean")) {
    x_col <- "drift_mean"
    fixed_col <- "bound_mean"
    fixed_levels <- c("lowBound", "highBound")
    panel_titles <- c(expression(paste("Low ", mu[alpha])),
                      expression(paste("High ", mu[alpha])))
    x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
  } else {
    x_col <- "bound_mean"
    fixed_col <- "drift_mean"
    fixed_levels <- c("lowDrift", "highDrift")
    panel_titles <- c(expression(paste("Low ", mu[nu])),
                      expression(paste("High ", mu[nu])))
    x_axis_label <- expression(paste("True population mean boundary separation (", mu[alpha], ")"))
  }
  
  bin_spec <- define_param_bins(x_param)
  x_ticks <- bin_spec$ticks
  
  files_clean <- c()
  files_cont <- c()
  for (combo in condition_info$parameter_cells) {
    combo_dir <- file.path(main_dir, combo)
    pattern <- paste0("sim_P[0-9]+T", t_level_select, "_.*\\.RData$")
    files_clean <- c(files_clean, list.files(file.path(combo_dir, clean_condition), pattern = pattern, full.names = TRUE))
    files_cont <- c(files_cont, list.files(file.path(combo_dir, contaminated_condition), pattern = pattern, full.names = TRUE))
  }
  all_files <- c(files_clean, files_cont)
  progress_bar <- txtProgressBar(min = 0, max = max(1, length(all_files)), style = 3, width = 50, char = "=")
  
  panel_df <- data.frame(
    panel = character(0),
    condition = character(0),
    x_level = character(0),
    x_value = numeric(0),
    median = numeric(0),
    q1 = numeric(0),
    q3 = numeric(0),
    stringsAsFactors = FALSE
  )
  
  step <- 0
  if (length(files_clean) > 0) {
    for (file_path in files_clean) {
      out <- extract_rt_param_by_bins(
        file_path = file_path,
        beta_level_select = beta_level_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec,
        show_IQrange = show_IQrange
      )
      if (!is.null(out) && nrow(out) > 0) {
        out$condition <- "clean"
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  if (length(files_cont) > 0) {
    for (file_path in files_cont) {
      out <- extract_rt_param_by_bins(
        file_path = file_path,
        beta_level_select = beta_level_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec,
        show_IQrange = show_IQrange
      )
      if (!is.null(out) && nrow(out) > 0) {
        out$condition <- "contaminated"
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  close(progress_bar)
  
  if (nrow(panel_df) == 0) {
    stop("No nested RT data found for selected T/beta.")
  }
  
  if (isTRUE(show_IQrange)) {
    agg <- aggregate(cbind(median, q1, q3) ~ panel + condition + x_level + x_value,
                     data = panel_df, FUN = median, na.rm = TRUE)
    y_min <- min(agg$q1, agg$median, agg$q3, na.rm = TRUE)
    y_max <- max(agg$q1, agg$median, agg$q3, na.rm = TRUE)
  } else {
    agg <- aggregate(median ~ panel + condition + x_level + x_value,
                     data = panel_df, FUN = median, na.rm = TRUE)
    agg$q1 <- NA_real_
    agg$q3 <- NA_real_
    y_min <- min(agg$median, na.rm = TRUE)
    y_max <- max(agg$median, na.rm = TRUE)
  }
  y_range <- c(y_min, y_max)
  y_pad <- if (diff(y_range) > 0) diff(y_range) * 0.08 else 0.05
  y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  y_ticks <- pretty(y_range, n = 5)
  
  metric_tag <- if (isTRUE(show_IQrange)) "IQrange" else "median"
  beta_tag <- gsub("\\.", "", as.character(beta_level_select))
  beta_tag <- sub("^0+", "", beta_tag)
  if (nchar(beta_tag) == 0) beta_tag <- "0"
  x_param_tag <- sub("_mean$", "", x_param)
  output_filename <- paste0(
    metric_tag, "_by_", x_param_tag, "_T", t_level_select,
    "_beta", beta_tag,
    ifelse(is.null(custom_title_label), "", paste0("_", custom_title_label)),
    ".pdf"
  )
  output_path <- file.path(output_dir, output_filename)
  
  pdf(output_path, width = 10, height = 5.8)
  par(mfrow = c(1, 2), oma = c(1.8, 2.6, 0.9, 0.1), mar = c(2.4, 2.4, 1.4, 0.15), cex = 1.0)
  
  panel_left <- agg[agg$panel == fixed_levels[1], ]
  panel_right <- agg[agg$panel == fixed_levels[2], ]
  # Reduce only inner margins between panels, preserving plot region size.
  par(mar = c(2.6, 2.5, 1.0, 0.2))
  make_RTplot_nestedcell(panel_left, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                     show_y_axis = TRUE, show_legend = FALSE, show_IQrange = show_IQrange)
  par(mar = c(2.6, 0.6, 1.0, 0.6))
  make_RTplot_nestedcell(panel_right, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                     show_y_axis = FALSE, show_legend = TRUE, show_IQrange = show_IQrange)
  
  if (isTRUE(show_IQrange)) {
    mtext(expression(paste("|MeanRT - MedianRT|")), side = 2, line = 0.4, outer = TRUE, cex = 2.2, font = 1)
  } else {
    mtext(expression(paste("Median of |MeanRT - MedianRT|")), side = 2, line = 0.4, outer = TRUE, cex = 2.2, font = 1)
  }
  mtext(x_axis_label, side = 1, line = 0.8, outer = TRUE, cex = 2.2, font = 1)
  mtext(panel_titles[1], side = 3, line = -1.2, at = 0.25, outer = TRUE, cex = 2.2, font = 2)
  mtext(panel_titles[2], side = 3, line = -1.2, at = 0.75, outer = TRUE, cex = 2.2, font = 2)
  dev.off()
  
  cat("RT nested-by-parameter plot saved to:", output_path, "\n")
}

##########################################################################
# Alternative: beta-difference of RT differences (meanRT - medianRT)
##########################################################################
# For each participant, compute raw (meanRT - medianRT) at beta = 0 and at a
# non-zero beta (default 2). Plot, per bin/panel, the median of:
#    (meanRT - medianRT)|beta_alt  -  (meanRT - medianRT)|beta0
# Two lines per panel: clean vs contaminated datasets.
plot_RTdiffdiff_nested <- function(main_dir, output_dir, t_level_select = 40,
                                   beta_alt_select = 0.2, x_param = "drift_mean",
                                   custom_title_label = NULL) {
  if (!x_param %in% c("drift_mean", "bound_mean")) {
    stop("Invalid x_param. Use 'drift_mean' or 'bound_mean'.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  condition_info <- read_conditions(main_dir = main_dir)
  if (!identical(condition_info$layout, "nested")) {
    stop("plot_RTdiffdiff_nested() requires nested simulation folders.")
  }
  conditions <- condition_info$conditions
  clean_condition <- sort(conditions[grepl("_clean$", conditions)])[1]
  contaminated_condition <- sort(conditions[grepl("_contaminated$", conditions)])[1]
  if (is.na(clean_condition) || is.na(contaminated_condition)) {
    stop("Could not infer clean/contaminated condition folders.")
  }
  if (identical(x_param, "drift_mean")) {
    x_col <- "drift_mean"
    fixed_col <- "bound_mean"
    fixed_levels <- c("lowBound", "highBound")
    panel_titles <- c(expression(paste("Low ", mu[alpha])),
                      expression(paste("High ", mu[alpha])))
    x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
  } else {
    x_col <- "bound_mean"
    fixed_col <- "drift_mean"
    fixed_levels <- c("lowDrift", "highDrift")
    panel_titles <- c(expression(paste("Low ", mu[nu])),
                      expression(paste("High ", mu[nu])))
    x_axis_label <- expression(paste("True population mean boundary separation (", mu[alpha], ")"))
  }
  bin_spec <- define_param_bins(x_param)
  x_ticks <- bin_spec$ticks
  files_clean <- c()
  files_cont <- c()
  for (combo in condition_info$parameter_cells) {
    combo_dir <- file.path(main_dir, combo)
    pattern <- paste0("sim_P[0-9]+T", t_level_select, "_.*\\.RData$")
    files_clean <- c(files_clean, list.files(file.path(combo_dir, clean_condition), pattern = pattern, full.names = TRUE))
    files_cont <- c(files_cont, list.files(file.path(combo_dir, contaminated_condition), pattern = pattern, full.names = TRUE))
  }
  all_files <- c(files_clean, files_cont)
  progress_bar <- txtProgressBar(min = 0, max = max(1, length(all_files)), style = 3, width = 50, char = "=")
  panel_df <- data.frame(
    panel = character(0),
    condition = character(0),
    x_level = character(0),
    x_value = numeric(0),
    median = numeric(0),
    stringsAsFactors = FALSE
  )
  step <- 0
  if (length(files_clean) > 0) {
    for (file_path in files_clean) {
      out <- extract_rt_diffdiff_by_bins(
        file_path = file_path,
        beta0 = 0,
        beta_alt = beta_alt_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec
      )
      if (!is.null(out) && nrow(out) > 0) {
        out$condition <- "clean"
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  if (length(files_cont) > 0) {
    for (file_path in files_cont) {
      out <- extract_rt_diffdiff_by_bins(
        file_path = file_path,
        beta0 = 0,
        beta_alt = beta_alt_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec
      )
      if (!is.null(out) && nrow(out) > 0) {
        out$condition <- "contaminated"
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  close(progress_bar)
  if (nrow(panel_df) == 0) {
    stop("No nested RT data found for selected T / beta levels.")
  }
  agg <- aggregate(median ~ panel + condition + x_level + x_value,
                   data = panel_df, FUN = median, na.rm = TRUE)
  y_min <- min(agg$median, na.rm = TRUE)
  y_max <- max(agg$median, na.rm = TRUE)
  y_range <- c(y_min, y_max)
  y_pad <- if (diff(y_range) > 0) diff(y_range) * 0.08 else 0.05
  y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  y_ticks <- pretty(y_range, n = 5)
  beta_tag <- paste0("0_vs_", gsub("\\.", "", as.character(beta_alt_select)))
  x_param_tag <- sub("_mean$", "", x_param)
  output_filename <- paste0(
    "diffdiff_by_", x_param_tag, "_T", t_level_select,
    "_beta", beta_tag,
    ifelse(is.null(custom_title_label), "", paste0("_", custom_title_label)),
    ".pdf"
  )
  output_path <- file.path(output_dir, output_filename)
  pdf(output_path, width = 10, height = 5.8)
  par(mfrow = c(1, 2), oma = c(1.8, 2.6, 0.9, 0.1), mar = c(2.4, 2.4, 1.4, 0.15), cex = 1.0)
  panel_left <- agg[agg$panel == fixed_levels[1], ]
  panel_right <- agg[agg$panel == fixed_levels[2], ]
  par(mar = c(2.6, 2.5, 1.0, 0.2))
  make_RTdiffdiff_nestedcell(panel_left, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                             show_y_axis = TRUE, show_legend = FALSE)
  par(mar = c(2.6, 0.6, 1.0, 0.6))
  make_RTdiffdiff_nestedcell(panel_right, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                             show_y_axis = FALSE, show_legend = TRUE)
  mtext(expression(paste("Median of [(MeanRT - MedianRT)|", beta[alt], " - (MeanRT - MedianRT)|", beta[0], "]")),
        side = 2, line = 0.4, outer = TRUE, cex = 1.8, font = 1)
  mtext(x_axis_label, side = 1, line = 0.8, outer = TRUE, cex = 2.0, font = 1)
  mtext(panel_titles[1], side = 3, line = -1.2, at = 0.25, outer = TRUE, cex = 2.0, font = 2)
  mtext(panel_titles[2], side = 3, line = -1.2, at = 0.75, outer = TRUE, cex = 2.0, font = 2)
  dev.off()
  cat("RT diffdiff nested-by-parameter plot saved to:", output_path, "\n")
}

# Helper: compute per-file bin medians of (diff_betaAlt - diff_beta0), raw diffs
extract_rt_diffdiff_by_bins <- function(file_path, beta0, beta_alt, x_col, fixed_col, fixed_levels, bin_spec) {
  e <- new.env(parent = emptyenv())
  load(file_path, envir = e)
  if (!exists("simStudy_Beta", envir = e, inherits = FALSE)) {
    return(NULL)
  }
  sim_beta <- get("simStudy_Beta", envir = e)
  true_df <- sim_beta$true
  summ_df <- sim_beta$summStats
  betas <- as.numeric(true_df[, "betaweight"])
  idx0 <- which(abs(betas - beta0) < 1e-8)
  idxA <- which(abs(betas - beta_alt) < 1e-8)
  # Strict: require exact beta levels; if absent in this file, skip it.
  if (length(idx0) == 0 || length(idxA) == 0) return(NULL)
  # Unpaired medians per bin at beta0 and beta_alt, then difference of medians.
  mean0 <- as.numeric(unlist(summ_df[idx0, "meanRT"]))
  med0 <- as.numeric(unlist(summ_df[idx0, "medianRT"]))
  diff0 <- mean0 - med0
  x0 <- as.numeric(true_df[idx0, x_col])
  fixed0 <- as.numeric(true_df[idx0, fixed_col])
  meanA <- as.numeric(unlist(summ_df[idxA, "meanRT"]))
  medA <- as.numeric(unlist(summ_df[idxA, "medianRT"]))
  diffA <- meanA - medA
  xA <- as.numeric(true_df[idxA, x_col])
  fixedA <- as.numeric(true_df[idxA, fixed_col])
  if (length(diff0) == 0 || length(diffA) == 0) return(NULL)
  fixed_ranges <- rt_low_high_ranges(fixed_col)
  # Helper to compute median per bin for a given beta set
  medians_by_bin <- function(x_vals, fixed_vals, diffs) {
    panels <- list()
    for (panel_name in fixed_levels) {
      in_panel <- if (panel_name == fixed_levels[1]) {
        fixed_vals >= fixed_ranges$low[1] & fixed_vals <= fixed_ranges$low[2]
      } else {
        fixed_vals >= fixed_ranges$high[1] & fixed_vals <= fixed_ranges$high[2]
      }
      if (!any(in_panel, na.rm = TRUE)) next
      px <- x_vals[in_panel]
      pd <- diffs[in_panel]
      x_ranges <- rt_low_high_ranges(x_col)
      x_subgroup <- ifelse(
        px >= x_ranges$low[1] & px <= x_ranges$low[2],
        "low",
        ifelse(px >= x_ranges$high[1] & px <= x_ranges$high[2], "high", NA_character_)
      )
      panel_out <- data.frame(panel = character(0), x_level = character(0), bin = integer(0),
                              x_value = numeric(0), median = numeric(0), stringsAsFactors = FALSE)
      for (x_level in c("low", "high")) {
        idx_level <- which(x_subgroup == x_level)
        if (length(idx_level) == 0) next
        breaks <- bin_spec$breaks[[x_level]]
        vals_x <- px[idx_level]
        vals_d <- pd[idx_level]
        bin_id <- cut(vals_x, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
        for (b in seq_len(length(breaks) - 1)) {
          in_bin <- which(bin_id == b)
          if (length(in_bin) == 0) next
          med_stat <- median(vals_d[in_bin], na.rm = TRUE)
          x_center <- mean(c(breaks[b], breaks[b + 1]))
          panel_out <- rbind(panel_out, data.frame(
            panel = panel_name, x_level = x_level, bin = b,
            x_value = x_center, median = med_stat, stringsAsFactors = FALSE
          ))
        }
      }
      panels[[panel_name]] <- panel_out
    }
    do.call(rbind, panels)
  }
  by0 <- medians_by_bin(x0, fixed0, diff0)
  byA <- medians_by_bin(xA, fixedA, diffA)
  if (is.null(by0) || is.null(byA) || nrow(by0) == 0 || nrow(byA) == 0) return(NULL)
  merged <- merge(byA, by0, by = c("panel", "x_level", "bin", "x_value"), suffixes = c("_alt", "_0"))
  if (nrow(merged) == 0) return(NULL)
  out <- data.frame(
    panel = merged$panel,
    x_level = merged$x_level,
    x_value = merged$x_value,
    median = merged$median_alt - merged$median_0,
    stringsAsFactors = FALSE
  )
  out
}

make_RTdiffdiff_nestedcell <- function(panel_df, x_ticks, y_ticks, y_range,
                                      show_y_axis = TRUE, show_legend = FALSE) {
  x_range <- range(x_ticks, na.rm = TRUE)
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE), cex.axis = 1.35)
  if (show_y_axis) {
    axis(2, at = y_ticks, las = 1, cex.axis = 1.7)
  }
  line_cols <- c(clean = "#1b9e77", contaminated = "#d95f02")
  point_pch <- c(clean = 16, contaminated = 15)
  for (cond in c("clean", "contaminated")) {
    dd <- panel_df[panel_df$condition == cond, ]
    if (nrow(dd) == 0) next
    for (lev in unique(dd$x_level)) {
      dd_lev <- dd[dd$x_level == lev, ]
      ord <- order(dd_lev$x_value)
      x_ord <- dd_lev$x_value[ord]
      med_ord <- dd_lev$median[ord]
      lines(x_ord, med_ord, col = line_cols[[cond]], lwd = 2.6, lty = 1)
      points(x_ord, med_ord, col = line_cols[[cond]], pch = point_pch[[cond]], cex = 1.0)
    }
  }
  abline(h = 0, lty = 3, col = "gray60", lwd = 1)
  if (show_legend) {
    legend("topright",
           legend = c("Clean dataset", "Contaminated dataset"),
           col = c(line_cols[["clean"]], line_cols[["contaminated"]]),
           lty = c(1, 1),
           pch = c(point_pch[["clean"]], point_pch[["contaminated"]]),
           lwd = c(2.6, 2.6),
           bty = "n", cex = 1.1)
  }
}

##########################################################################
# Alternative nested plot: contamination effect on RT summaries
##########################################################################
# This function mirrors plot_RTdiff_nested() layout, but compares clean vs
# contaminated datasets directly within matched simulation files.
# It plots:
#   - |meanRT_clean - meanRT_contaminated|
#   - |medianRT_clean - medianRT_contaminated|
# with lines distinguished by summary statistic.
plot_RTcontamination_nested <- function(main_dir, output_dir, t_level_select = 40,
                                        beta_level_select = 0, x_param = "drift_mean",
                                        custom_title_label = NULL, show_IQrange = FALSE) {
  if (!x_param %in% c("drift_mean", "bound_mean")) {
    stop("Invalid x_param. Use 'drift_mean' or 'bound_mean'.")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  condition_info <- read_conditions(main_dir = main_dir)
  if (!identical(condition_info$layout, "nested")) {
    stop("plot_RTcontamination_nested() requires nested simulation folders.")
  }
  conditions <- condition_info$conditions
  clean_condition <- sort(conditions[grepl("_clean$", conditions)])[1]
  contaminated_condition <- sort(conditions[grepl("_contaminated$", conditions)])[1]
  if (is.na(clean_condition) || is.na(contaminated_condition)) {
    stop("Could not infer clean/contaminated condition folders.")
  }
  
  if (identical(x_param, "drift_mean")) {
    x_col <- "drift_mean"
    fixed_col <- "bound_mean"
    fixed_levels <- c("lowBound", "highBound")
    panel_titles <- c(expression(paste("Low ", mu[alpha])),
                      expression(paste("High ", mu[alpha])))
    x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
  } else {
    x_col <- "bound_mean"
    fixed_col <- "drift_mean"
    fixed_levels <- c("lowDrift", "highDrift")
    panel_titles <- c(expression(paste("Low ", mu[nu])),
                      expression(paste("High ", mu[nu])))
    x_axis_label <- expression(paste("True population mean boundary separation (", mu[alpha], ")"))
  }
  
  bin_spec <- define_param_bins(x_param)
  x_ticks <- bin_spec$ticks
  
  paired_files <- list()
  normalize_pair_key <- function(path_vec) {
    nm <- basename(path_vec)
    nm <- tolower(nm)
    nm <- gsub("contaminated", "", nm, fixed = TRUE)
    nm <- gsub("clean", "", nm, fixed = TRUE)
    nm <- gsub("__+", "_", nm)
    nm <- gsub("_\\.", ".", nm)
    nm
  }
  for (combo in condition_info$parameter_cells) {
    combo_dir <- file.path(main_dir, combo)
    pattern <- paste0("sim_P[0-9]+T", t_level_select, "_.*\\.RData$")
    clean_files <- list.files(file.path(combo_dir, clean_condition), pattern = pattern, full.names = TRUE)
    cont_files <- list.files(file.path(combo_dir, contaminated_condition), pattern = pattern, full.names = TRUE)
    if (length(clean_files) == 0 || length(cont_files) == 0) next
    clean_keys <- normalize_pair_key(clean_files)
    cont_keys <- normalize_pair_key(cont_files)
    clean_map <- split(clean_files, clean_keys)
    cont_map <- split(cont_files, cont_keys)
    shared_keys <- intersect(names(clean_map), names(cont_map))
    if (length(shared_keys) > 0) {
      for (k in shared_keys) {
        n_pair <- min(length(clean_map[[k]]), length(cont_map[[k]]))
        if (n_pair == 0) next
        for (j in seq_len(n_pair)) {
          paired_files[[length(paired_files) + 1]] <- list(clean = clean_map[[k]][j], contaminated = cont_map[[k]][j])
        }
      }
    } else {
      # Fallback: pair in sorted order when filenames differ systematically by condition tags.
      clean_sorted <- sort(clean_files)
      cont_sorted <- sort(cont_files)
      n_pair <- min(length(clean_sorted), length(cont_sorted))
      if (n_pair > 0) {
        for (j in seq_len(n_pair)) {
          paired_files[[length(paired_files) + 1]] <- list(clean = clean_sorted[j], contaminated = cont_sorted[j])
        }
      }
    }
  }
  
  progress_bar <- txtProgressBar(min = 0, max = max(1, length(paired_files)), style = 3, width = 50, char = "=")
  panel_df <- data.frame(
    panel = character(0),
    stat = character(0),
    x_level = character(0),
    x_value = numeric(0),
    median = numeric(0),
    q1 = numeric(0),
    q3 = numeric(0),
    stringsAsFactors = FALSE
  )
  
  step <- 0
  if (length(paired_files) > 0) {
    for (pair in paired_files) {
      out <- extract_rt_contam_by_bins(
        file_clean = pair$clean,
        file_contaminated = pair$contaminated,
        beta_level_select = beta_level_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec,
        show_IQrange = show_IQrange
      )
      if (!is.null(out) && nrow(out) > 0) {
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  close(progress_bar)
  
  if (nrow(panel_df) == 0) stop("No matched clean/contaminated nested RT data found for selected T/beta.")
  
  if (isTRUE(show_IQrange)) {
    agg <- aggregate(cbind(median, q1, q3) ~ panel + stat + x_level + x_value,
                     data = panel_df, FUN = median, na.rm = TRUE)
    y_min <- min(agg$q1, agg$median, agg$q3, na.rm = TRUE)
    y_max <- max(agg$q1, agg$median, agg$q3, na.rm = TRUE)
  } else {
    agg <- aggregate(median ~ panel + stat + x_level + x_value,
                     data = panel_df, FUN = median, na.rm = TRUE)
    agg$q1 <- NA_real_
    agg$q3 <- NA_real_
    y_min <- min(agg$median, na.rm = TRUE)
    y_max <- max(agg$median, na.rm = TRUE)
  }
  
  y_range <- c(y_min, y_max)
  y_pad <- if (diff(y_range) > 0) diff(y_range) * 0.08 else 0.05
  y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  y_ticks <- pretty(y_range, n = 5)
  
  metric_tag <- if (isTRUE(show_IQrange)) "IQrange_contamination" else "median_contamination"
  beta_tag <- gsub("\\.", "", as.character(beta_level_select))
  beta_tag <- sub("^0+", "", beta_tag); if (nchar(beta_tag) == 0) beta_tag <- "0"
  x_param_tag <- sub("_mean$", "", x_param)
  output_filename <- paste0(
    metric_tag, "_by_", x_param_tag, "_T", t_level_select,
    "_beta", beta_tag,
    ifelse(is.null(custom_title_label), "", paste0("_", custom_title_label)),
    ".pdf"
  )
  output_path <- file.path(output_dir, output_filename)
  
  pdf(output_path, width = 10, height = 5.8)
  par(mfrow = c(1, 2), oma = c(1.8, 2.6, 0.9, 0.1), mar = c(2.4, 2.4, 1.4, 0.15), cex = 1.0)
  panel_left <- agg[agg$panel == fixed_levels[1], ]
  panel_right <- agg[agg$panel == fixed_levels[2], ]
  par(mar = c(2.6, 2.5, 1.0, 0.2))
  make_RTcontam_nestedcell(panel_left, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                           show_y_axis = TRUE, show_legend = FALSE, show_IQrange = show_IQrange)
  par(mar = c(2.6, 0.6, 1.0, 0.6))
  make_RTcontam_nestedcell(panel_right, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                           show_y_axis = FALSE, show_legend = TRUE, show_IQrange = show_IQrange)
  
  if (isTRUE(show_IQrange)) {
    mtext(expression(paste("|Clean - Contaminated|")), side = 2, line = 0.4, outer = TRUE, cex = 2.2, font = 1)
  } else {
    mtext(expression(paste("Median of |Clean - Contaminated|")), side = 2, line = 0.4, outer = TRUE, cex = 2.2, font = 1)
  }
  mtext(x_axis_label, side = 1, line = 0.8, outer = TRUE, cex = 2.2, font = 1)
  mtext(panel_titles[1], side = 3, line = -1.2, at = 0.25, outer = TRUE, cex = 2.2, font = 2)
  mtext(panel_titles[2], side = 3, line = -1.2, at = 0.75, outer = TRUE, cex = 2.2, font = 2)
  dev.off()
  cat("RT contamination nested-by-parameter plot saved to:", output_path, "\n")
}

# Helper: compute contamination differences by shared bins in a clean/contaminated file pair.
extract_rt_contam_by_bins <- function(file_clean, file_contaminated, beta_level_select,
                                      x_col, fixed_col, fixed_levels, bin_spec, show_IQrange = FALSE) {
  align_to_target_length <- function(x, target_n) {
    if (target_n <= 0 || length(x) == 0) return(numeric(0))
    if (length(x) == target_n) return(x)
    if (length(x) == 1) return(rep(x, target_n))
    rep_factor <- target_n / length(x)
    if (is.finite(rep_factor) && rep_factor >= 1 && abs(rep_factor - round(rep_factor)) < 1e-8) {
      return(rep(x, each = as.integer(round(rep_factor))))
    }
    x[seq_len(min(length(x), target_n))]
  }
  
  e_clean <- new.env(parent = emptyenv())
  e_cont <- new.env(parent = emptyenv())
  load(file_clean, envir = e_clean)
  load(file_contaminated, envir = e_cont)
  if (!exists("simStudy_Beta", envir = e_clean, inherits = FALSE)) return(NULL)
  if (!exists("simStudy_Beta", envir = e_cont, inherits = FALSE)) return(NULL)
  sim_clean <- get("simStudy_Beta", envir = e_clean)
  sim_cont <- get("simStudy_Beta", envir = e_cont)
  
  true_clean <- sim_clean$true
  true_cont <- sim_cont$true
  summ_clean <- sim_clean$summStats
  summ_cont <- sim_cont$summStats
  if (is.null(true_clean) || is.null(true_cont) || is.null(summ_clean) || is.null(summ_cont)) return(NULL)
  
  beta_clean <- as.numeric(true_clean[, "betaweight"])
  beta_cont <- as.numeric(true_cont[, "betaweight"])
  idx_clean <- which(abs(beta_clean - beta_level_select) < 1e-8)
  idx_cont <- which(abs(beta_cont - beta_level_select) < 1e-8)
  if (length(idx_clean) == 0 || length(idx_cont) == 0) {
    ub_clean <- sort(unique(beta_clean)); ub_clean <- ub_clean[!is.na(ub_clean)]
    ub_cont <- sort(unique(beta_cont)); ub_cont <- ub_cont[!is.na(ub_cont)]
    shared_beta <- intersect(ub_clean, ub_cont)
    if (length(shared_beta) == 0) return(NULL)
    beta_use <- shared_beta[which.min(abs(shared_beta - beta_level_select))]
    idx_clean <- which(abs(beta_clean - beta_use) < 1e-8)
    idx_cont <- which(abs(beta_cont - beta_use) < 1e-8)
  }
  if (length(idx_clean) == 0 || length(idx_cont) == 0) return(NULL)
  
  mean_clean <- as.numeric(unlist(summ_clean[idx_clean, "meanRT"]))
  mean_cont <- as.numeric(unlist(summ_cont[idx_cont, "meanRT"]))
  median_clean <- as.numeric(unlist(summ_clean[idx_clean, "medianRT"]))
  median_cont <- as.numeric(unlist(summ_cont[idx_cont, "medianRT"]))
  
  n <- min(length(mean_clean), length(mean_cont), length(median_clean), length(median_cont))
  if (n == 0) return(NULL)
  mean_diff <- abs(mean_clean[seq_len(n)] - mean_cont[seq_len(n)])
  median_diff <- abs(median_clean[seq_len(n)] - median_cont[seq_len(n)])
  
  x_vals_seed <- as.numeric(true_clean[idx_clean, x_col])
  fixed_vals_seed <- as.numeric(true_clean[idx_clean, fixed_col])
  x_vals <- align_to_target_length(x_vals_seed, n)
  fixed_vals <- align_to_target_length(fixed_vals_seed, n)
  n2 <- min(n, length(x_vals), length(fixed_vals))
  if (n2 == 0) return(NULL)
  x_vals <- x_vals[seq_len(n2)]
  fixed_vals <- fixed_vals[seq_len(n2)]
  mean_diff <- mean_diff[seq_len(n2)]
  median_diff <- median_diff[seq_len(n2)]
  
  fixed_ranges <- rt_low_high_ranges(fixed_col)
  fixed_groups <- ifelse(
    fixed_vals >= fixed_ranges$low[1] & fixed_vals <= fixed_ranges$low[2],
    fixed_levels[1],
    ifelse(
      fixed_vals >= fixed_ranges$high[1] & fixed_vals <= fixed_ranges$high[2],
      fixed_levels[2],
      NA_character_
    )
  )
  
  out <- data.frame(
    panel = character(0),
    stat = character(0),
    x_level = character(0),
    x_value = numeric(0),
    median = numeric(0),
    q1 = numeric(0),
    q3 = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (panel_name in fixed_levels) {
    in_panel <- fixed_groups == panel_name
    if (!any(in_panel, na.rm = TRUE)) next
    panel_x <- x_vals[in_panel]
    panel_mean <- mean_diff[in_panel]
    panel_median <- median_diff[in_panel]
    
    x_ranges <- rt_low_high_ranges(x_col)
    x_subgroup <- ifelse(
      panel_x >= x_ranges$low[1] & panel_x <= x_ranges$low[2],
      "low",
      ifelse(panel_x >= x_ranges$high[1] & panel_x <= x_ranges$high[2], "high", NA_character_)
    )
    
    for (x_level in c("low", "high")) {
      idx_level <- which(x_subgroup == x_level)
      if (length(idx_level) == 0) next
      breaks <- bin_spec$breaks[[x_level]]
      vals_x <- panel_x[idx_level]
      vals_mean <- panel_mean[idx_level]
      vals_median <- panel_median[idx_level]
      bin_id <- cut(vals_x, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
      for (b in seq_len(length(breaks) - 1)) {
        in_bin <- which(bin_id == b)
        if (length(in_bin) == 0) next
        med_mean <- median(vals_mean[in_bin], na.rm = TRUE)
        med_median <- median(vals_median[in_bin], na.rm = TRUE)
        q1_mean <- if (isTRUE(show_IQrange)) as.numeric(quantile(vals_mean[in_bin], probs = 0.25, na.rm = TRUE, type = 7)) else NA_real_
        q3_mean <- if (isTRUE(show_IQrange)) as.numeric(quantile(vals_mean[in_bin], probs = 0.75, na.rm = TRUE, type = 7)) else NA_real_
        q1_median <- if (isTRUE(show_IQrange)) as.numeric(quantile(vals_median[in_bin], probs = 0.25, na.rm = TRUE, type = 7)) else NA_real_
        q3_median <- if (isTRUE(show_IQrange)) as.numeric(quantile(vals_median[in_bin], probs = 0.75, na.rm = TRUE, type = 7)) else NA_real_
        
        out <- rbind(out, data.frame(
          panel = panel_name, stat = "meanRT", x_level = x_level,
          x_value = mean(c(breaks[b], breaks[b + 1])),
          median = med_mean, q1 = q1_mean, q3 = q3_mean,
          stringsAsFactors = FALSE
        ))
        out <- rbind(out, data.frame(
          panel = panel_name, stat = "medianRT", x_level = x_level,
          x_value = mean(c(breaks[b], breaks[b + 1])),
          median = med_median, q1 = q1_median, q3 = q3_median,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  out
}

make_RTcontam_nestedcell <- function(panel_df, x_ticks, y_ticks, y_range,
                                     show_y_axis = TRUE, show_legend = FALSE, show_IQrange = FALSE) {
  x_range <- range(x_ticks, na.rm = TRUE)
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE), cex.axis = 1.35)
  if (show_y_axis) axis(2, at = y_ticks, las = 1, cex.axis = 1.7)
  
  line_cols <- c(meanRT = "#1f78b4", medianRT = "#984ea3")
  fill_cols <- c(
    meanRT = rgb(31/255, 120/255, 180/255, alpha = 0.18),
    medianRT = rgb(152/255, 78/255, 163/255, alpha = 0.18)
  )
  point_pch <- c(meanRT = 16, medianRT = 16)
  
  for (st in c("meanRT", "medianRT")) {
    dd <- panel_df[panel_df$stat == st, ]
    if (nrow(dd) == 0) next
    for (lev in unique(dd$x_level)) {
      dd_lev <- dd[dd$x_level == lev, ]
      ord <- order(dd_lev$x_value)
      x_ord <- dd_lev$x_value[ord]
      med_ord <- dd_lev$median[ord]
      if (isTRUE(show_IQrange)) {
        q1_ord <- dd_lev$q1[ord]
        q3_ord <- dd_lev$q3[ord]
        valid_iqr <- !(is.na(q1_ord) | is.na(q3_ord))
        if (sum(valid_iqr) >= 2) {
          x_poly <- c(x_ord[valid_iqr], rev(x_ord[valid_iqr]))
          y_poly <- c(q1_ord[valid_iqr], rev(q3_ord[valid_iqr]))
          polygon(x_poly, y_poly, border = NA, col = fill_cols[[st]])
          lines(x_ord[valid_iqr], q1_ord[valid_iqr], col = line_cols[[st]], lwd = 1.8, lty = 1)
          lines(x_ord[valid_iqr], q3_ord[valid_iqr], col = line_cols[[st]], lwd = 1.8, lty = 1)
        }
      }
      lines(x_ord, med_ord, col = line_cols[[st]], lwd = 2.6, lty = 1)
      points(x_ord, med_ord, col = line_cols[[st]], pch = point_pch[[st]], cex = 1.0)
    }
  }
  abline(h = 0, lty = 3, col = "gray60", lwd = 1)
  
  if (show_legend) {
    if (isTRUE(show_IQrange)) {
      legend("topright",
             legend = c("MeanRT clean-contaminated difference", "MedianRT clean-contaminated difference", "Q1 and Q3 bounds"),
             col = c(line_cols[["meanRT"]], line_cols[["medianRT"]], "gray30"),
             lty = c(1, 1, 1),
             pch = c(point_pch[["meanRT"]], point_pch[["medianRT"]], NA),
             lwd = c(2.6, 2.6, 1.8),
             bty = "n", cex = 1.0)
    } else {
      legend("topright",
             legend = c("MeanRT clean-contaminated difference", "MedianRT clean-contaminated difference"),
             col = c(line_cols[["meanRT"]], line_cols[["medianRT"]]),
             lty = c(1, 1),
             pch = c(point_pch[["meanRT"]], point_pch[["medianRT"]]),
             lwd = c(2.6, 2.6),
             bty = "n", cex = 1.0)
    }
  }
}

########################################################################################################
# Helper function to plot the absolute meanRT - medianRT difference across values of a chosen true parameter
########################################################################################################
# Low/high true-parameter ranges from nested simulation design.
rt_low_high_ranges <- function(param_col) {
  if (identical(param_col, "drift_mean")) {
    return(list(low = c(0, 1), high = c(2, 3)))
  }
  if (identical(param_col, "bound_mean")) {
    return(list(low = c(2, 2.5), high = c(3.5, 4)))
  }
  stop("Unsupported parameter for low/high range mapping.")
}

# Define the bins for the x-axis based on the true parameter
define_param_bins <- function(x_param) {
  if(x_param %in% c("drift", "drift_mean")){
      low_breaks <- seq(0, 1, length.out = 5)
      high_breaks <- seq(2, 3, length.out = 5)
      return(list(breaks = list(low = low_breaks, high = high_breaks),
                  ticks = sort(unique(c(low_breaks, high_breaks)))))
  }
  if (!x_param %in% c("bound", "bound_mean")) {
    stop("Invalid x_param in define_param_bins(). Use 'drift_mean' or 'bound_mean'.")
  }
  low_breaks <- seq(2, 2.5, length.out = 3)
  high_breaks <- seq(3.5, 4, length.out = 3)
  list(breaks = list(low = low_breaks, high = high_breaks),
       ticks = sort(unique(c(low_breaks, high_breaks))))
}

extract_rt_param_by_bins <- function(file_path, beta_level_select, x_col, fixed_col, fixed_levels, bin_spec, show_IQrange = FALSE) {
  e <- new.env(parent = emptyenv())
  load(file_path, envir = e)
  if (!exists("simStudy_Beta", envir = e, inherits = FALSE)) {
    return(NULL)
  }
  sim_beta <- get("simStudy_Beta", envir = e)
  true_df <- sim_beta$true
  summ_df <- sim_beta$summStats
  
  beta_vals <- as.numeric(true_df[, "betaweight"])
  idx <- which(abs(beta_vals - beta_level_select) < 1e-8)
  if (length(idx) == 0) {
    unique_betas <- sort(unique(beta_vals))
    unique_betas <- unique_betas[!is.na(unique_betas)]
    if (length(unique_betas) == 0) {
      return(NULL)
    }
    beta_use <- unique_betas[which.min(abs(unique_betas - beta_level_select))]
    idx <- which(abs(beta_vals - beta_use) < 1e-8)
  }
  if (length(idx) == 0) {
    return(NULL)
  }
  
  mean_vals <- as.numeric(unlist(summ_df[idx, "meanRT"]))
  median_vals <- as.numeric(unlist(summ_df[idx, "medianRT"]))
  rt_diff_raw <- mean_vals - median_vals
  
  x_vals_seed <- as.numeric(true_df[idx, x_col])
  fixed_vals_seed <- as.numeric(true_df[idx, fixed_col])
  if (length(x_vals_seed) == 0 || length(rt_diff_raw) == 0) {
    return(NULL)
  }
  
  rep_factor <- length(rt_diff_raw) / length(x_vals_seed)
  if (is.finite(rep_factor) && rep_factor >= 1 && abs(rep_factor - round(rep_factor)) < 1e-8) {
    rep_factor <- as.integer(round(rep_factor))
    x_vals <- rep(x_vals_seed, each = rep_factor)
    fixed_vals <- rep(fixed_vals_seed, each = rep_factor)
  } else {
    n <- min(length(rt_diff_raw), length(x_vals_seed), length(fixed_vals_seed))
    rt_diff_raw <- rt_diff_raw[seq_len(n)]
    x_vals <- x_vals_seed[seq_len(n)]
    fixed_vals <- fixed_vals_seed[seq_len(n)]
  }
  
  fixed_ranges <- rt_low_high_ranges(fixed_col)
  fixed_groups <- ifelse(
    fixed_vals >= fixed_ranges$low[1] & fixed_vals <= fixed_ranges$low[2],
    fixed_levels[1],
    ifelse(
      fixed_vals >= fixed_ranges$high[1] & fixed_vals <= fixed_ranges$high[2],
      fixed_levels[2],
      NA_character_
    )
  )
  
  out <- data.frame(
    panel = character(0),
    x_level = character(0),
    x_value = numeric(0),
    median = numeric(0),
    q1 = numeric(0),
    q3 = numeric(0),
    stringsAsFactors = FALSE
  )
  for (panel_name in fixed_levels) {
    in_panel <- fixed_groups == panel_name
    if (!any(in_panel, na.rm = TRUE)) next
    panel_x <- x_vals[in_panel]
    panel_rt <- rt_diff_raw[in_panel]
    
    x_ranges <- rt_low_high_ranges(x_col)
    x_subgroup <- ifelse(
      panel_x >= x_ranges$low[1] & panel_x <= x_ranges$low[2],
      "low",
      ifelse(
        panel_x >= x_ranges$high[1] & panel_x <= x_ranges$high[2],
        "high",
        NA_character_
      )
    )
    for (x_level in c("low", "high")) {
      idx_level <- x_subgroup == x_level
      if (!any(idx_level, na.rm = TRUE)) next
      breaks <- bin_spec$breaks[[x_level]]
      vals_x <- panel_x[idx_level]
      vals_rt <- panel_rt[idx_level]
      bin_id <- cut(vals_x, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
      for (b in seq_len(length(breaks) - 1)) {
        in_bin <- which(bin_id == b)
        if (length(in_bin) == 0) next
        abs_vals <- abs(vals_rt[in_bin])
        med_stat <- median(abs_vals, na.rm = TRUE)
        q1_stat <- if (isTRUE(show_IQrange)) as.numeric(quantile(abs_vals, probs = 0.25, na.rm = TRUE, type = 7)) else NA_real_
        q3_stat <- if (isTRUE(show_IQrange)) as.numeric(quantile(abs_vals, probs = 0.75, na.rm = TRUE, type = 7)) else NA_real_
        out <- rbind(out, data.frame(
          panel = panel_name,
          x_level = x_level,
          x_value = mean(c(breaks[b], breaks[b + 1])),
          median = med_stat,
          q1 = q1_stat,
          q3 = q3_stat,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  out
}

make_RTplot_nestedcell <- function(panel_df, x_ticks, y_ticks, y_range,
                               show_y_axis = TRUE, show_legend = FALSE, show_IQrange = FALSE) {
  x_range <- range(x_ticks, na.rm = TRUE)
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE), cex.axis = 1.35)
  if (show_y_axis) {
    axis(2, at = y_ticks, las = 1, cex.axis = 1.7)
  }
  
  line_cols <- c(clean = "#1b9e77", contaminated = "#d95f02")
  fill_cols <- c(
    clean = rgb(27/255, 158/255, 119/255, alpha = 0.22),
    contaminated = rgb(217/255, 95/255, 2/255, alpha = 0.22)
  )
  point_pch <- c(clean = 16, contaminated = 15)
  for (cond in c("clean", "contaminated")) {
    dd <- panel_df[panel_df$condition == cond, ]
    if (nrow(dd) == 0) next
    for (lev in unique(dd$x_level)) {
      dd_lev <- dd[dd$x_level == lev, ]
      ord <- order(dd_lev$x_value)
      x_ord <- dd_lev$x_value[ord]
      med_ord <- dd_lev$median[ord]
      # Median absolute difference: thick solid
      # Optional IQR ribbon and quartile lines
      if (isTRUE(show_IQrange)) {
        q1_ord <- dd_lev$q1[ord]
        q3_ord <- dd_lev$q3[ord]
        valid_iqr <- !(is.na(q1_ord) | is.na(q3_ord))
        if (sum(valid_iqr) >= 2) {
          x_poly <- c(x_ord[valid_iqr], rev(x_ord[valid_iqr]))
          y_poly <- c(q1_ord[valid_iqr], rev(q3_ord[valid_iqr]))
          polygon(x_poly, y_poly, border = NA, col = fill_cols[[cond]])
          lines(x_ord[valid_iqr], q1_ord[valid_iqr], col = line_cols[[cond]], lwd = 1.8, lty = 1)
          lines(x_ord[valid_iqr], q3_ord[valid_iqr], col = line_cols[[cond]], lwd = 1.8, lty = 1)
        }
      }
      lines(x_ord, med_ord, col = line_cols[[cond]], lwd = 2.6, lty = 1)
      points(x_ord, med_ord, col = line_cols[[cond]], pch = point_pch[[cond]], cex = 1.0)
    }
  }
  abline(h = 0, lty = 3, col = "gray60", lwd = 1)
  
  if (show_legend) {
    if (isTRUE(show_IQrange)) {
      legend_entries <- c("Clean median difference", "Contaminated median difference", "Q1 and Q3 bounds")
      lg <- legend("topright",
                   legend = legend_entries,
                   col = c(line_cols[["clean"]], line_cols[["contaminated"]], rgb(1, 1, 1, 0)),
                   lty = c(1, 1, 1),
                   pch = c(point_pch[["clean"]], point_pch[["contaminated"]], NA),
                   lwd = c(2.6, 2.6, 1.8),
                   bty = "n", cex = 1.1, plot = FALSE)
      legend("topright",
             legend = legend_entries,
             col = c(line_cols[["clean"]], line_cols[["contaminated"]], rgb(1, 1, 1, 0)),
             lty = c(1, 1, 1),
             pch = c(point_pch[["clean"]], point_pch[["contaminated"]], NA),
             lwd = c(2.6, 2.6, 1.8),
             bty = "n", cex = 1.1)
      y_mid <- lg$text$y[3]
      dy <- diff(par("usr")[3:4]) * 0.006
      x_start <- lg$rect$left + lg$rect$w * 0.02
      x_end <- lg$rect$left + lg$rect$w * 0.16
      # Draw a light underlay then two parallel colored quartile cue lines.
      segments(x_start, y_mid, x_end, y_mid, col = "white", lwd = 2.6, xpd = NA)
      segments(x_start, y_mid + dy, x_end, y_mid + dy, col = line_cols[["clean"]], lwd = 1.8, lty = 1, xpd = NA)
      segments(x_start, y_mid - dy, x_end, y_mid - dy, col = line_cols[["contaminated"]], lwd = 1.8, lty = 1, xpd = NA)
    } else {
      legend("topright",
             legend = c("Clean median", "Contaminated median"),
             col = c(line_cols[["clean"]], line_cols[["contaminated"]]),
             lty = c(1, 1),
             pch = c(point_pch[["clean"]], point_pch[["contaminated"]]),
             lwd = c(2.6, 2.6),
             bty = "n", cex = 1.1)
    }
  }
}




# #######################################################################################################
# Third custom function: RT differences vs parameter using EZ theoretical mean (M^{pred})
# #######################################################################################################
# This plot mirrors plot_RTdiff_nested layout (1x2 panels for low/high of the non-selected parameter).
# Y-axis shows the median absolute difference between theoretical mean RT (EZ forward equations) and
# empirical meanRT or medianRT, producing four lines per panel:
#   - Clean meanRT, Contaminated meanRT, Clean medianRT, Contaminated medianRT
# Inputs:
#   - main_dir: nested simulation directory (low/high combinations)
#   - output_dir: where to save the PDF
#   - t_level_select: T level to use (default 40)
#   - beta_level_select: beta level to use (default 0)
#   - x_param: "drift_mean" or "bound_mean" (default "drift_mean")
#   - custom_title_label: optional suffix for filename
plot_RTdiff_predNested <- function(main_dir, output_dir, t_level_select = 40,
                                   beta_level_select = 0, x_param = "drift_mean",
                                   custom_title_label = NULL, show_IQrange = FALSE) {
  if (!x_param %in% c("drift_mean", "bound_mean")) {
    stop("Invalid x_param. Use 'drift_mean' or 'bound_mean'.")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  condition_info <- read_conditions(main_dir = main_dir)
  if (!identical(condition_info$layout, "nested")) {
    stop("plot_RTdiff_predNested requires nested simulation folders.")
  }
  conditions <- condition_info$conditions
  clean_condition <- sort(conditions[grepl("_clean$", conditions)])[1]
  contaminated_condition <- sort(conditions[grepl("_contaminated$", conditions)])[1]
  if (is.na(clean_condition) || is.na(contaminated_condition)) {
    stop("Could not infer clean/contaminated condition folders.")
  }
  if (identical(x_param, "drift_mean")) {
    x_col <- "drift_mean"
    fixed_col <- "bound_mean"
    fixed_levels <- c("lowBound", "highBound")
    panel_titles <- c(expression(paste("Low ", mu[alpha])),
                      expression(paste("High ", mu[alpha])))
    x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
  } else {
    x_col <- "bound_mean"
    fixed_col <- "drift_mean"
    fixed_levels <- c("lowDrift", "highDrift")
    panel_titles <- c(expression(paste("Low ", mu[nu])),
                      expression(paste("High ", mu[nu])))
    x_axis_label <- expression(paste("True population mean boundary separation (", mu[alpha], ")"))
  }
  bin_spec <- define_param_bins(x_param)
  x_ticks <- bin_spec$ticks
  files_clean <- c()
  files_cont <- c()
  for (combo in condition_info$parameter_cells) {
    combo_dir <- file.path(main_dir, combo)
    pattern <- paste0("sim_P[0-9]+T", t_level_select, "_.*\\.RData$")
    files_clean <- c(files_clean, list.files(file.path(combo_dir, clean_condition), pattern = pattern, full.names = TRUE))
    files_cont <- c(files_cont, list.files(file.path(combo_dir, contaminated_condition), pattern = pattern, full.names = TRUE))
  }
  all_files <- c(files_clean, files_cont)
  progress_bar <- txtProgressBar(min = 0, max = max(1, length(all_files)), style = 3, width = 50, char = "=")
  panel_df <- data.frame(
    panel = character(0),
    condition = character(0),
    stat = character(0),              # "meanRT" or "medianRT"
    x_level = character(0),
    x_value = numeric(0),
    diff_val = numeric(0),
    stringsAsFactors = FALSE
  )
  step <- 0
  if (length(files_clean) > 0) {
    for (file_path in files_clean) {
      out <- extract_pred_diffs_by_bins(
        file_path = file_path,
        beta_level_select = beta_level_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec,
        show_IQrange = show_IQrange
      )
      if (!is.null(out) && nrow(out) > 0) {
        out$condition <- "clean"
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  if (length(files_cont) > 0) {
    for (file_path in files_cont) {
      out <- extract_pred_diffs_by_bins(
        file_path = file_path,
        beta_level_select = beta_level_select,
        x_col = x_col,
        fixed_col = fixed_col,
        fixed_levels = fixed_levels,
        bin_spec = bin_spec,
        show_IQrange = show_IQrange
      )
      if (!is.null(out) && nrow(out) > 0) {
        out$condition <- "contaminated"
        panel_df <- rbind(panel_df, out)
      }
      step <- step + 1
      setTxtProgressBar(progress_bar, step)
    }
  }
  close(progress_bar)
  if (nrow(panel_df) == 0) stop("No nested RT data found for selected T/beta.")
  # Aggregate across files: median per panel/condition/stat/bin
  if (isTRUE(show_IQrange)) {
    agg <- aggregate(cbind(diff_val, q1, q3) ~ panel + condition + stat + x_level + x_value,
                     data = panel_df, FUN = median, na.rm = TRUE)
    y_min <- min(c(agg$q1, agg$diff_val), na.rm = TRUE)
    y_max <- max(c(agg$q3, agg$diff_val), na.rm = TRUE)
  } else {
    agg <- aggregate(diff_val ~ panel + condition + stat + x_level + x_value,
                     data = panel_df, FUN = median, na.rm = TRUE)
    y_min <- min(agg$diff_val, na.rm = TRUE)
    y_max <- max(agg$diff_val, na.rm = TRUE)
  }
  y_range <- c(y_min, y_max)
  y_pad <- if (diff(y_range) > 0) diff(y_range) * 0.08 else 0.05
  y_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  y_ticks <- pretty(y_range, n = 5)
  beta_tag <- gsub("\\.", "", as.character(beta_level_select))
  beta_tag <- sub("^0+", "", beta_tag); if (nchar(beta_tag) == 0) beta_tag <- "0"
  x_param_tag <- sub("_mean$", "", x_param)
  metric_tag <- ifelse(isTRUE(show_IQrange), "predDiffIQrange", "predDiffMedian")
  output_filename <- paste0(
    metric_tag, "_by_", x_param_tag, "_T", t_level_select,
    "_beta", beta_tag,
    ifelse(is.null(custom_title_label), "", paste0("_", custom_title_label)),
    ".pdf"
  )
  output_path <- file.path(output_dir, output_filename)
  pdf(output_path, width = 10, height = 5.8)
  par(mfrow = c(1, 2), oma = c(1.8, 2.6, 0.9, 0.1), mar = c(2.4, 2.4, 1.4, 0.15), cex = 1.0)
  for (i in seq_along(fixed_levels)) {
    panel_key <- fixed_levels[i]
    show_y_axis <- (i == 1)
    show_legend <- (i == 1)
    dd <- agg[agg$panel == panel_key, ]
    make_RTpred_nestedcell(dd, x_ticks = x_ticks, y_ticks = y_ticks, y_range = y_range,
                           show_y_axis = show_y_axis, show_legend = show_legend,
                           show_IQrange = show_IQrange)
  }
  y_label <- if (isTRUE(show_IQrange)) {
    expression(paste("|", M^{pred}, " - empirical|"))
  } else {
    expression(paste("Median of |", M^{pred}, " - empirical|"))
  }
  mtext(y_label, side = 2, line = 0.4, outer = TRUE, cex = 2.1, font = 1)
  mtext(x_axis_label, side = 1, line = 0.8, outer = TRUE, cex = 2.1, font = 1)
  mtext(panel_titles[1], side = 3, line = -1.2, at = 0.25, outer = TRUE, cex = 2.1, font = 2)
  mtext(panel_titles[2], side = 3, line = -1.2, at = 0.75, outer = TRUE, cex = 2.1, font = 2)
  dev.off()
  cat("RT predicted-mean difference plot saved to:", output_path, "\n")
}

# Compute per-file bin medians of |M_pred - meanRT| and |M_pred - medianRT|
extract_pred_diffs_by_bins <- function(file_path, beta_level_select, x_col, fixed_col, fixed_levels, bin_spec,
                                       show_IQrange = FALSE) {
  extract_true_param_vector <- function(df, row_idx, param_base) {
    if (length(row_idx) == 0) return(numeric(0))
    # 1) Direct column (scalar per row or list-column)
    if (param_base %in% colnames(df)) {
      x <- df[row_idx, param_base, drop = TRUE]
      vals <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
      if (length(vals) > 0 && any(!is.na(vals))) return(vals)
    }
    # 2) Expanded participant columns, e.g., drift1, drift_1, drift[1]
    patt <- paste0("^", param_base, "(\\[?[0-9]+\\]?|_[0-9]+)$")
    cols <- grep(patt, colnames(df), value = TRUE)
    if (length(cols) > 0) {
      vals <- suppressWarnings(as.numeric(unlist(df[row_idx, cols, drop = FALSE], use.names = FALSE)))
      if (length(vals) > 0) return(vals)
    }
    numeric(0)
  }
  
  align_to_target_length <- function(x, target_n) {
    if (target_n <= 0 || length(x) == 0) return(numeric(0))
    if (length(x) == target_n) return(x)
    if (length(x) == 1) return(rep(x, target_n))
    rep_factor <- target_n / length(x)
    if (is.finite(rep_factor) && rep_factor >= 1 && abs(rep_factor - round(rep_factor)) < 1e-8) {
      return(rep(x, each = as.integer(round(rep_factor))))
    }
    x[seq_len(min(length(x), target_n))]
  }
  
  e <- new.env(parent = emptyenv())
  load(file_path, envir = e)
  if (!exists("simStudy_Beta", envir = e, inherits = FALSE)) return(NULL)
  sim_beta <- get("simStudy_Beta", envir = e)
  true_df <- sim_beta$true
  summ_df <- sim_beta$summStats
  beta_vals <- as.numeric(true_df[, "betaweight"])
  idx <- which(abs(beta_vals - beta_level_select) < 1e-8)
  if (length(idx) == 0) {
    unique_betas <- sort(unique(beta_vals)); unique_betas <- unique_betas[!is.na(unique_betas)]
    if (length(unique_betas) == 0) return(NULL)
    beta_use <- unique_betas[which.min(abs(unique_betas - beta_level_select))]
    idx <- which(abs(beta_vals - beta_use) < 1e-8)
  }
  if (length(idx) == 0) return(NULL)
  mean_vals <- as.numeric(unlist(summ_df[idx, "meanRT"]))
  median_vals <- as.numeric(unlist(summ_df[idx, "medianRT"]))
  # Use participant-level true parameters when available.
  # Fallback to population means only if individual columns are absent.
  v_all <- extract_true_param_vector(true_df, idx, "drift")
  a_all <- extract_true_param_vector(true_df, idx, "bound")
  t0_all <- extract_true_param_vector(true_df, idx, "nondt")
  if (length(v_all) == 0) v_all <- suppressWarnings(as.numeric(unlist(true_df[idx, "drift_mean"], use.names = FALSE)))
  if (length(a_all) == 0) a_all <- suppressWarnings(as.numeric(unlist(true_df[idx, "bound_mean"], use.names = FALSE)))
  if (length(t0_all) == 0) t0_all <- suppressWarnings(as.numeric(unlist(true_df[idx, "nondt_mean"], use.names = FALSE)))
  target_n <- min(length(mean_vals), length(median_vals))
  v_all <- align_to_target_length(v_all, target_n)
  a_all <- align_to_target_length(a_all, target_n)
  t0_all <- align_to_target_length(t0_all, target_n)
  target_n <- min(target_n, length(v_all), length(a_all), length(t0_all))
  if (target_n == 0) return(NULL)
  mean_vals <- mean_vals[seq_len(target_n)]
  median_vals <- median_vals[seq_len(target_n)]
  v_all <- v_all[seq_len(target_n)]
  a_all <- a_all[seq_len(target_n)]
  t0_all <- t0_all[seq_len(target_n)]
  # Theoretical predicted mean via EZ forward equations
  s_bv <- exp(-a_all * v_all)
  m_pred <- t0_all + (a_all / (2 * v_all)) * ((1 - s_bv) / (1 + s_bv))
  # X and fixed values for binning/panels
  x_vals_seed <- as.numeric(true_df[idx, x_col])
  fixed_vals_seed <- as.numeric(true_df[idx, fixed_col])
  x_vals <- align_to_target_length(x_vals_seed, length(m_pred))
  fixed_vals <- align_to_target_length(fixed_vals_seed, length(m_pred))
  n2 <- min(length(x_vals), length(fixed_vals), length(m_pred))
  if (n2 == 0) return(NULL)
  x_vals <- x_vals[seq_len(n2)]
  fixed_vals <- fixed_vals[seq_len(n2)]
  m_pred <- m_pred[seq_len(n2)]
  mean_vals <- mean_vals[seq_len(n2)]
  median_vals <- median_vals[seq_len(n2)]
  fixed_ranges <- rt_low_high_ranges(fixed_col)
  fixed_groups <- ifelse(
    fixed_vals >= fixed_ranges$low[1] & fixed_vals <= fixed_ranges$low[2],
    fixed_levels[1],
    ifelse(
      fixed_vals >= fixed_ranges$high[1] & fixed_vals <= fixed_ranges$high[2],
      fixed_levels[2],
      NA_character_
    )
  )
  out <- data.frame(
    panel = character(0),
    stat = character(0),
    x_level = character(0),
    x_value = numeric(0),
    diff_val = numeric(0),
    q1 = numeric(0),
    q3 = numeric(0),
    stringsAsFactors = FALSE
  )
  for (panel_name in fixed_levels) {
    in_panel <- fixed_groups == panel_name
    if (!any(in_panel, na.rm = TRUE)) next
    panel_x <- x_vals[in_panel]
    pred_panel <- m_pred[in_panel]
    mean_panel <- mean_vals[in_panel]
    median_panel <- median_vals[in_panel]
    x_ranges <- rt_low_high_ranges(x_col)
    x_subgroup <- ifelse(
      panel_x >= x_ranges$low[1] & panel_x <= x_ranges$low[2],
      "low",
      ifelse(
        panel_x >= x_ranges$high[1] & panel_x <= x_ranges$high[2],
        "high",
        NA_character_
      )
    )
    for (xl in c("low", "high")) {
      idx_level <- which(x_subgroup == xl)
      if (length(idx_level) == 0) next
      vals_x <- panel_x[idx_level]
      d_mean <- abs(pred_panel[idx_level] - mean_panel[idx_level])
      d_median <- abs(pred_panel[idx_level] - median_panel[idx_level])
      breaks <- bin_spec$breaks[[xl]]
      bin_id <- cut(vals_x, breaks = breaks, include.lowest = TRUE, right = TRUE, labels = FALSE)
      for (b in seq_len(length(breaks) - 1)) {
        in_bin <- which(bin_id == b)
        if (length(in_bin) == 0) next
        q1_mean <- if (isTRUE(show_IQrange)) as.numeric(quantile(d_mean[in_bin], probs = 0.25, na.rm = TRUE)) else NA_real_
        q3_mean <- if (isTRUE(show_IQrange)) as.numeric(quantile(d_mean[in_bin], probs = 0.75, na.rm = TRUE)) else NA_real_
        q1_median <- if (isTRUE(show_IQrange)) as.numeric(quantile(d_median[in_bin], probs = 0.25, na.rm = TRUE)) else NA_real_
        q3_median <- if (isTRUE(show_IQrange)) as.numeric(quantile(d_median[in_bin], probs = 0.75, na.rm = TRUE)) else NA_real_
        out <- rbind(out, data.frame(
          panel = panel_name,
          stat = "meanRT",
          x_level = xl,
          x_value = mean(c(breaks[b], breaks[b + 1])),
          diff_val = median(d_mean[in_bin], na.rm = TRUE),
          q1 = q1_mean,
          q3 = q3_mean,
          stringsAsFactors = FALSE
        ))
        out <- rbind(out, data.frame(
          panel = panel_name,
          stat = "medianRT",
          x_level = xl,
          x_value = mean(c(breaks[b], breaks[b + 1])),
          diff_val = median(d_median[in_bin], na.rm = TRUE),
          q1 = q1_median,
          q3 = q3_median,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  out
}

# Plot a nested cell with four lines: clean/contaminated x meanRT/medianRT diffs
make_RTpred_nestedcell <- function(panel_df, x_ticks, y_ticks, y_range,
                                  show_y_axis = TRUE, show_legend = FALSE, show_IQrange = FALSE) {
  x_range <- range(x_ticks, na.rm = TRUE)
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE), cex.axis = 1.35)
  if (show_y_axis) axis(2, at = y_ticks, las = 1, cex.axis = 1.7)
  # Colors: keep dataset-coded color, stat-coded lty
  line_cols <- c(clean = "#1b9e77", contaminated = "#d95f02")
  fill_cols <- c(
    clean = rgb(27/255, 158/255, 119/255, alpha = 0.15),
    contaminated = rgb(217/255, 95/255, 2/255, alpha = 0.15)
  )
  lty_by_stat <- c(meanRT = 1, medianRT = 2)  # solid for meanRT, dashed for medianRT
  pch_by_combo <- list(
    clean_meanRT = 16, contaminated_meanRT = 15,
    clean_medianRT = 1, contaminated_medianRT = 0
  )
  combos <- list(
    list(condition = "clean", stat = "meanRT"),
    list(condition = "contaminated", stat = "meanRT"),
    list(condition = "clean", stat = "medianRT"),
    list(condition = "contaminated", stat = "medianRT")
  )
  for (cmb in combos) {
    cond <- cmb$condition; stat <- cmb$stat
    dd <- panel_df[panel_df$condition == cond & panel_df$stat == stat, ]
    if (nrow(dd) == 0) next
    for (lev in unique(dd$x_level)) {
      dd_lev <- dd[dd$x_level == lev, ]
      ord <- order(dd_lev$x_value)
      x_ord <- dd_lev$x_value[ord]
      y_ord <- dd_lev$diff_val[ord]
      if (isTRUE(show_IQrange) && all(c("q1", "q3") %in% names(dd_lev))) {
        q1_ord <- dd_lev$q1[ord]
        q3_ord <- dd_lev$q3[ord]
        valid_iqr <- !(is.na(q1_ord) | is.na(q3_ord))
        if (sum(valid_iqr) >= 2) {
          x_poly <- c(x_ord[valid_iqr], rev(x_ord[valid_iqr]))
          y_poly <- c(q1_ord[valid_iqr], rev(q3_ord[valid_iqr]))
          polygon(x_poly, y_poly, border = NA, col = fill_cols[[cond]])
          lines(x_ord[valid_iqr], q1_ord[valid_iqr], col = line_cols[[cond]], lwd = 1.5, lty = lty_by_stat[[stat]])
          lines(x_ord[valid_iqr], q3_ord[valid_iqr], col = line_cols[[cond]], lwd = 1.5, lty = lty_by_stat[[stat]])
        }
      }
      lines(x_ord, y_ord,
            col = line_cols[[cond]], lwd = 2.2, lty = lty_by_stat[[stat]])
      points(x_ord, y_ord,
             col = line_cols[[cond]],
             pch = pch_by_combo[[paste(cond, stat, sep = "_")]], cex = 1.0)
    }
  }
  if (show_legend) {
    if (isTRUE(show_IQrange)) {
      legend("topright",
             legend = c("Clean meanRT (median)", "Contaminated meanRT (median)",
                        "Clean medianRT (median)", "Contaminated medianRT (median)",
                        "Q1/Q3 bounds and IQR ribbon"),
             col = c(line_cols[["clean"]], line_cols[["contaminated"]],
                     line_cols[["clean"]], line_cols[["contaminated"]], "gray30"),
             lty = c(lty_by_stat[["meanRT"]], lty_by_stat[["meanRT"]],
                     lty_by_stat[["medianRT"]], lty_by_stat[["medianRT"]], 1),
             pch = c(pch_by_combo$clean_meanRT, pch_by_combo$contaminated_meanRT,
                     pch_by_combo$clean_medianRT, pch_by_combo$contaminated_medianRT, NA),
             lwd = c(2.2, 2.2, 2.2, 2.2, 1.5), bty = "n", cex = 1.0)
    } else {
      legend("topright",
             legend = c("Clean meanRT", "Contaminated meanRT", "Clean medianRT", "Contaminated medianRT"),
             col = c(line_cols[["clean"]], line_cols[["contaminated"]],
                     line_cols[["clean"]], line_cols[["contaminated"]]),
             lty = c(lty_by_stat[["meanRT"]], lty_by_stat[["meanRT"]],
                     lty_by_stat[["medianRT"]], lty_by_stat[["medianRT"]]),
             pch = c(pch_by_combo$clean_meanRT, pch_by_combo$contaminated_meanRT,
                     pch_by_combo$clean_medianRT, pch_by_combo$contaminated_medianRT),
             lwd = 2, bty = "n", cex = 1.1)
    }
  }
}