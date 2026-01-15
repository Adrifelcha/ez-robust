# Main function to create the grid of RMSE, MSE, Bias, and Variance plots
# Generates four PDF files: one for RMSE, one for MSE, one for Bias, and one for Variance
plot_RMSEgrid <- function(main_dir, output_dir, parameter = "bound_mean", plot_by = "condition", 
                         y_range_rmse = NULL, y_range_mse = NULL, y_range_bias = NULL, y_range_variance = NULL) {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Validate inputs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Validate plot_by argument
  if (!plot_by %in% c("condition", "beta")) {
    stop("Invalid 'plot_by' argument. Choose 'condition' or 'beta'.")
  }  
  # Validate parameter argument
  valid_params <- c("bound_mean", "drift_mean", "nondt_mean", "betaweight")
  if (!parameter %in% valid_params) {
    stop(paste("Invalid 'parameter' argument. Choose one of:", paste(valid_params, collapse = ", ")))
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create output directory, if it doesn't exist
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define simulation conditions and labels
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  conditions <- c("EZ_contaminated", "EZRobust_contaminated", "EZ_clean", "EZRobust_clean")
  condition_labels <- c("EZ x Contaminated", "Robust x Contaminated", "EZ x Clean","Robust x Clean")
  
  # Infer P and T levels from the first subfolder
  first_condition_path <- file.path(main_dir, conditions[1])
  all_files <- list.files(first_condition_path, pattern = "\\.RData$", full.names = TRUE)
  rdata_files <- all_files[grepl("_P\\d+T\\d+_", basename(all_files))]
  
  # Get P and T values
  filenames <- basename(rdata_files)
  p_values <- as.numeric(sub(".*_P(\\d+)T.*", "\\1", filenames))
  t_values <- as.numeric(sub(".*T(\\d+)_.*", "\\1", filenames))
  p_levels <- sort(unique(p_values))
  t_levels <- sort(unique(t_values))
  
  # Print start message
  cat("\n")
  cat("============================================================\n")
  cat("Creating RMSE, MSE, Bias, and Variance grid plots for parameter:", parameter, "\n")
  cat("Plot by:", plot_by, "\n")
  cat("Grid dimensions:", length(p_levels), "x", length(t_levels), "(P x T)\n")
  cat("============================================================\n\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute RMSE, MSE, Bias, and Variance for all cells
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Computing metrics for all cells...\n")
  all_rmse_data <- list()
  all_mse_data <- list()
  all_bias_data <- list()
  all_variance_data <- list()
  total_cells <- length(t_levels) * length(p_levels)
  current_cell <- 0
  
  for (t_level in t_levels) {
    for (p_level in p_levels) {
      current_cell <- current_cell + 1
      cell_key <- paste(p_level, t_level, sep = "_")
      cat("Computing cell", current_cell, "of", total_cells, ": P =", p_level, ", T =", t_level, "\n")
      
      
      rmse_data_list <- list()
      mse_data_list <- list()
      bias_data_list <- list()
      variance_data_list <- list()
      
      for (condition in conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path) && file.exists(file_path)) {
          rmse_data <- get_cellRMSE(resultsFile = file_path, parameter = parameter)
          for (beta_level_char in names(rmse_data$rmse_by_beta)) {
            rmse_val <- rmse_data$rmse_by_beta[[beta_level_char]]
            mse_val <- rmse_data$mse_by_beta[[beta_level_char]]
            bias_val <- rmse_data$bias_by_beta[[beta_level_char]]
            variance_val <- rmse_data$variance_by_beta[[beta_level_char]]
            
            if (!is.na(rmse_val)) {
              rmse_data_list[[length(rmse_data_list) + 1]] <- data.frame(
                condition = condition, beta = as.numeric(beta_level_char), rmse = rmse_val, stringsAsFactors = FALSE)
            }
            if (!is.na(mse_val)) {
              mse_data_list[[length(mse_data_list) + 1]] <- data.frame(
                condition = condition, beta = as.numeric(beta_level_char), mse = mse_val, stringsAsFactors = FALSE)
            }
            if (!is.na(bias_val)) {
              bias_data_list[[length(bias_data_list) + 1]] <- data.frame(
                condition = condition, beta = as.numeric(beta_level_char), bias = bias_val, stringsAsFactors = FALSE)
            }
            if (!is.na(variance_val)) {
              variance_data_list[[length(variance_data_list) + 1]] <- data.frame(
                condition = condition, beta = as.numeric(beta_level_char), variance = variance_val, stringsAsFactors = FALSE)
            }
          }
        }
      }
      
      # Store the data for this cell
      if (length(rmse_data_list) > 0) {
        all_rmse_data[[cell_key]] <- do.call(rbind, rmse_data_list)
      } else {
        all_rmse_data[[cell_key]] <- NULL
      }
      if (length(mse_data_list) > 0) {
        all_mse_data[[cell_key]] <- do.call(rbind, mse_data_list)
      } else {
        all_mse_data[[cell_key]] <- NULL
      }
      if (length(bias_data_list) > 0) {
        all_bias_data[[cell_key]] <- do.call(rbind, bias_data_list)
      } else {
        all_bias_data[[cell_key]] <- NULL
      }
      if (length(variance_data_list) > 0) {
        all_variance_data[[cell_key]] <- do.call(rbind, variance_data_list)
      } else {
        all_variance_data[[cell_key]] <- NULL
      }
    }
  }
  
  # Determine y-axis ranges if not provided
  if (is.null(y_range_rmse)) {
    cat("\nDetermining RMSE y-axis range...\n")
    max_rmse <- 0
    for (cell_key in names(all_rmse_data)) {
      if (!is.null(all_rmse_data[[cell_key]])) {
        max_rmse <- max(max_rmse, all_rmse_data[[cell_key]]$rmse, na.rm = TRUE)
      }
    }
    y_range_rmse <- c(0, max_rmse * 1.1)
    cat("RMSE y-axis range:", y_range_rmse[1], "to", y_range_rmse[2], "\n")
  }
  
  if (is.null(y_range_mse)) {
    cat("Determining MSE y-axis range...\n")
    max_mse <- 0
    for (cell_key in names(all_mse_data)) {
      if (!is.null(all_mse_data[[cell_key]])) {
        max_mse <- max(max_mse, all_mse_data[[cell_key]]$mse, na.rm = TRUE)
      }
    }
    y_range_mse <- c(0, max_mse * 1.1)
    cat("MSE y-axis range:", y_range_mse[1], "to", y_range_mse[2], "\n")
  }
  
  if (is.null(y_range_bias)) {
    cat("Determining Bias y-axis range...\n")
    min_bias <- Inf
    max_bias <- -Inf
    for (cell_key in names(all_bias_data)) {
      if (!is.null(all_bias_data[[cell_key]])) {
        min_bias <- min(min_bias, all_bias_data[[cell_key]]$bias, na.rm = TRUE)
        max_bias <- max(max_bias, all_bias_data[[cell_key]]$bias, na.rm = TRUE)
      }
    }
    range_bias <- max_bias - min_bias
    y_range_bias <- c(min_bias - range_bias * 0.1, max_bias + range_bias * 0.1)
    cat("Bias y-axis range:", y_range_bias[1], "to", y_range_bias[2], "\n")
  }
  
  if (is.null(y_range_variance)) {
    cat("Determining Variance y-axis range...\n")
    max_variance <- 0
    for (cell_key in names(all_variance_data)) {
      if (!is.null(all_variance_data[[cell_key]])) {
        max_variance <- max(max_variance, all_variance_data[[cell_key]]$variance, na.rm = TRUE)
      }
    }
    y_range_variance <- c(0, max_variance * 1.1)
    cat("Variance y-axis range:", y_range_variance[1], "to", y_range_variance[2], "\n\n")
  }
  
  # Generate RMSE plot
  cat("Creating RMSE plot...\n")
  create_metric_plot(all_rmse_data, "rmse", "RMSE", parameter, plot_by, output_dir, 
                     conditions, condition_labels, p_levels, t_levels, y_range_rmse)
  
  # Generate MSE plot
  cat("\nCreating MSE plot...\n")
  create_metric_plot(all_mse_data, "mse", "MSE", parameter, plot_by, output_dir, 
                     conditions, condition_labels, p_levels, t_levels, y_range_mse)
  
  # Generate Bias plot
  cat("\nCreating Bias plot...\n")
  create_metric_plot(all_bias_data, "bias", "Bias", parameter, plot_by, output_dir, 
                     conditions, condition_labels, p_levels, t_levels, y_range_bias, 
                     add_zero_line = TRUE)
  
  # Generate Variance plot
  cat("\nCreating Variance plot...\n")
  create_metric_plot(all_variance_data, "variance", "Variance", parameter, plot_by, output_dir, 
                     conditions, condition_labels, p_levels, t_levels, y_range_variance)
  
  cat("\nAll plots completed!\n")
}

###################################################################
################# H E L P E R   F U N C T I O N S ####################
###################################################################

# Helper function to create a plot for any metric (RMSE, Bias, or Variance)
create_metric_plot <- function(all_data, metric_name, metric_label, parameter, plot_by, output_dir,
                               conditions, condition_labels, p_levels, t_levels, y_range, 
                               add_zero_line = FALSE) {
  
  # Define PDF output file
  output_filename <- paste0(metric_label, "_grid_", parameter, "_by_", plot_by, ".pdf")
  pdf(file.path(output_dir, output_filename), width = 12, height = 12)
  
  # Setup plot layout
  par(mfrow = c(length(t_levels), length(p_levels)),
      oma = c(7, 7, 3, 3), # bottom, left, top, right
      mar = c(1, 1.5, 0, 0))
  
  # Loop through each T level (rows) from high to low for plotting
  for (t_level in rev(t_levels)) {
    for (p_level in p_levels) {
      cell_key <- paste(p_level, t_level, sep = "_")
      metric_df <- all_data[[cell_key]]
      
      # Plotting logic for one cell
      if (!is.null(metric_df) && nrow(metric_df) > 0) {
        show_x_axis <- (t_level == min(t_levels))
        show_y_axis <- (p_level == min(p_levels))
        
        if (plot_by == "condition") {
          show_legend <- (t_level == max(t_levels)) && (p_level == max(p_levels))
          plot_cell_by_condition_metric(metric_df, metric_name, conditions, condition_labels, 
                                       y_range, show_x_axis, show_y_axis, show_legend, add_zero_line)
        } else {
          show_legend <- (t_level == max(t_levels)) && (p_level == min(p_levels))
          highlight_cell <- ifelse(as.numeric(p_level) * as.numeric(t_level) == 6400, 
                                   TRUE, FALSE)
          plot_cell_by_beta_metric(metric_df, metric_name, conditions, condition_labels, 
                                   y_range, show_x_axis, show_y_axis, show_legend,
                                   highlight_cell, add_zero_line)
        }
        
      } else {
        plot.new() # Draw an empty plot if no data
      }
      
      # Add cell labels
      if (t_level == max(t_levels)) {
        mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 2.5, font = 2)
      }
      if (p_level == max(p_levels)) {
        mtext(paste("T =", t_level), side = 4, line = 1.85, cex = 2.5, font = 2, las = 0)
      }
    }
  }
  
  # Add common outer labels
  mtext(paste(metric_label, "(", parameter, ")"), side = 2, line = 3.8, cex = 2.5, outer = TRUE)
  if (plot_by == "beta") {
    mtext(expression(paste("Effect size (", beta, ")")),
          side = 1, line = 5.7, cex = 2.5, outer = TRUE)
  } else {
    mtext("Condition", side = 1, line = 5.7, cex = 2.5, outer = TRUE)
  }

  dev.off()
  cat(metric_label, "grid plot saved to:", file.path(output_dir, output_filename), "\n")
}

###################################################################
################# P L O T     T Y P E  ################################
################## F U N C T I O N S ###################################
###################################################################

# BY CONDITION: The factor levels on the X axis are the conditions
plot_cell_by_condition_metric <- function(metric_df, metric_name, conditions, condition_labels, 
                                          y_range, show_x_axis, show_y_axis, show_legend, add_zero_line = FALSE) {
  # Create an empty plot
  plot(NA, NA, xlim = c(1, 4), ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty="o")
  
  # Add horizontal line at y=0 for bias plots
  if (add_zero_line) {
    abline(h = 0, lty = 2, col = "gray50", lwd = 1.5)
  }
  
  if(show_y_axis){ 
    y_at <- seq(y_range[1], y_range[2], length.out = 6)
    axis(2, at = y_at, round(y_at, 2), las = 1) 
  }
  
  x_at <- seq(1.2,3.7,length.out = length(condition_labels))
  
  if (show_x_axis) {
    axis(1, at = x_at, labels = FALSE) # Remove axis labels to add them manually
    split_labels <- strsplit(condition_labels, " x ")
    # Manually add labels on two lines
    for (i in seq_along(split_labels)) {
      mtext(split_labels[[i]][1], side = 1, line = 1.3, at = x_at[i], cex = 1.05, font = 2)
      mtext(split_labels[[i]][2], side = 1, line = 2.7, at = x_at[i], cex = 1.05, font = 2)
    }
  }
  
  beta_levels <- sort(unique(metric_df$beta))
  beta_colors <- c("#2ea02d", "#de8520", "#d62728", "#9467bd")
  
  # Plot lines and points for each beta level
  for (i in seq_along(beta_levels)) {
    beta_val <- beta_levels[i]
    subset_df <- metric_df[metric_df$beta == beta_val, ]
    ordered_subset <- subset_df[match(conditions, subset_df$condition), ]
    
    lines(x_at, ordered_subset[[metric_name]], col = beta_colors[i], lwd = 2)
    points(x_at, ordered_subset[[metric_name]], col = beta_colors[i], pch = 19, cex = 1.2)
  }
  
  # Add legend to the top-right plot
  if (show_legend) {
    legend("topright",
           legend = sapply(beta_levels, function(x) as.expression(bquote(beta == .(x)))),
           col = beta_colors,
           lwd = 2, pch = 19, bty = "n", cex = 1.2)
  }
}

# Helper function to plot a single grid cell with Beta levels on the X-axis
plot_cell_by_beta_metric <- function(metric_df, metric_name, conditions, condition_labels, y_range, 
                                     show_x_axis, show_y_axis, show_legend, highlight_cell = FALSE, 
                                     highlight_color = "white", add_zero_line = FALSE) {
  
  beta_levels <- sort(unique(metric_df$beta))
  offset <- 0.01
  xlim <- c(min(beta_levels) - offset, max(beta_levels) + offset)
  # Create an empty plot
  plot(NA, NA, xlim = xlim, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty="o")

  if(highlight_cell){
    # Add colored background using polygon
    polygon(x = c(xlim[1], xlim[2], xlim[2], xlim[1]), 
            y = c(y_range[1], y_range[1], y_range[2], y_range[2]), 
            col = highlight_color, border = NA)
  }

  # Add horizontal line at y=0 for bias plots
  if (add_zero_line) {
    abline(h = 0, lty = 2, col = "gray50", lwd = 1.5)
  }

  if (show_y_axis) {
    y_at <- seq(y_range[1], y_range[2], length.out = 6)
    axis(2, at = y_at, round(y_at, 2), las = 1, cex.axis = 2.5) 
  }

  if (show_x_axis) {
    labels <- beta_levels
    axis(1, at = beta_levels, labels = labels, line = 1, cex.axis = 2.5)
  }

  
  condition_colors <- c("#d3540b", "#160f0fea", "#47D647", "#E982FF")
  widths <- c(5,5,4,3)
  styles <- c(2,1,4,3)
  points <- c(19,17,15,18)

  # Plot lines and points for each condition
  for (i in seq_along(conditions)) {
    condition_val <- conditions[i]
    subset_df <- metric_df[metric_df$condition == condition_val, ]
    ordered_subset <- subset_df[order(subset_df$beta), ]
    
    lines(ordered_subset$beta, ordered_subset[[metric_name]], col = condition_colors[i], 
          lwd = widths[i], lty = styles[i])
    points(ordered_subset$beta, ordered_subset[[metric_name]], col = condition_colors[i], 
           pch = points[i], cex = 2)
  }
  
  # Add legend to the top-right plot
  if (show_legend) {
    legend("topleft",
           legend = condition_labels,
           col = condition_colors,
           lwd = 2, pch = 19, bty = "n", cex = 1.7)
  }
}
