# Main function to create a grid of plots showing mean posterior beta estimates
# across simulation conditions with error bars showing +/- 1 SD
plot_betaEstimateGrid <- function(main_dir, output_dir, y_range = NULL, y_axis_ticks = 6) {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Validate inputs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (!dir.exists(main_dir)) {
    stop(paste("Main directory does not exist:", main_dir))
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create output directory, if it doesn't exist
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define simulation conditions and labels
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  conditions <- c("EZ_contaminated", "EZRobust_contaminated", "EZ_clean", "EZRobust_clean")
  condition_labels <- c("EZ x Contaminated", "Robust x Contaminated", "EZ x Clean", "Robust x Clean")
  
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
  cat("Creating beta estimate grid plots\n")
  cat("Grid dimensions:", length(p_levels), "x", length(t_levels), "(P x T)\n")
  cat("============================================================\n\n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract beta estimates for all cells
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Extracting beta estimates for all cells...\n")
  all_beta_data <- list()
  total_cells <- length(t_levels) * length(p_levels)
  current_cell <- 0
  
  for (t_level in t_levels) {
    for (p_level in p_levels) {
      current_cell <- current_cell + 1
      cell_key <- paste(p_level, t_level, sep = "_")
      cat("Processing cell", current_cell, "of", total_cells, ": P =", p_level, ", T =", t_level, "\n")
      
      # Storage for this cell: list indexed by condition, each containing estimates by beta level
      cell_data <- list()
      
      # Load data from all conditions
      for (condition in conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        
        if (!is.na(file_path) && file.exists(file_path)) {
          # Load the data (simStudy_Beta is loaded from the RData file)
          load(file_path)
          
          # Extract true beta levels and estimates
          true_betas <- as.vector(unlist(simStudy_Beta$true[, "betaweight"]))
          beta_estimates <- as.vector(unlist(simStudy_Beta$estimates[, "betaweight"]))
          
          # Get unique beta levels
          beta_levels <- sort(unique(true_betas))
          
          # Store estimates by beta level for this condition
          condition_estimates <- list()
          for (beta_val in beta_levels) {
            beta_indices <- which(true_betas == beta_val)
            valid_indices <- !is.na(beta_estimates[beta_indices])
            if (sum(valid_indices) > 0) {
              condition_estimates[[as.character(beta_val)]] <- beta_estimates[beta_indices][valid_indices]
            } else {
              condition_estimates[[as.character(beta_val)]] <- numeric(0)
            }
          }
          cell_data[[condition]] <- condition_estimates
        }
      }
      
      # Compute summary statistics across conditions for each beta level
      beta_summary_list <- list()
      if (length(cell_data) > 0) {
        # Get all beta levels from first condition (should be same across conditions)
        first_condition <- names(cell_data)[1]
        beta_levels <- sort(as.numeric(names(cell_data[[first_condition]])))
        
        for (beta_val in beta_levels) {
          beta_char <- as.character(beta_val)
          
          # Compute mean and SD for each condition at this beta level
          condition_means <- numeric(length(conditions))
          condition_sds <- numeric(length(conditions))
          condition_counts <- numeric(length(conditions))
          
          for (i in seq_along(conditions)) {
            condition <- conditions[i]
            if (condition %in% names(cell_data) && beta_char %in% names(cell_data[[condition]])) {
              estimates <- cell_data[[condition]][[beta_char]]
              if (length(estimates) > 0) {
                condition_means[i] <- mean(estimates)
                condition_sds[i] <- sd(estimates)
                condition_counts[i] <- length(estimates)
              } else {
                condition_means[i] <- NA
                condition_sds[i] <- NA
                condition_counts[i] <- 0
              }
            } else {
              condition_means[i] <- NA
              condition_sds[i] <- NA
              condition_counts[i] <- 0
            }
          }
          
          # Store summary
          beta_summary_list[[beta_char]] <- data.frame(
            beta = beta_val,
            condition = conditions,
            condition_mean = condition_means,
            condition_sd = condition_sds,
            condition_count = condition_counts,
            stringsAsFactors = FALSE
          )
        }
        
        # Combine all beta summaries for this cell
        if (length(beta_summary_list) > 0) {
          all_beta_data[[cell_key]] <- do.call(rbind, beta_summary_list)
        } else {
          all_beta_data[[cell_key]] <- NULL
        }
      } else {
        all_beta_data[[cell_key]] <- NULL
      }
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Determine y-axis range if not provided
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (is.null(y_range)) {
    cat("\nDetermining y-axis range...\n")
    min_val <- Inf
    max_val <- -Inf
    for (cell_key in names(all_beta_data)) {
      if (!is.null(all_beta_data[[cell_key]])) {
        cell_df <- all_beta_data[[cell_key]]
        # Consider both means and means +/- SD for range
        means_plus_sd <- cell_df$condition_mean + cell_df$condition_sd
        means_minus_sd <- cell_df$condition_mean - cell_df$condition_sd
        min_val <- min(min_val, means_minus_sd, na.rm = TRUE)
        max_val <- max(max_val, means_plus_sd, na.rm = TRUE)
      }
    }
    range_val <- max_val - min_val
    y_range <- c(min_val - range_val * 0.05, max_val + range_val * 0.2)
    cat("Y-axis range:", y_range[1], "to", y_range[2], "\n\n")
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Generate plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Creating beta estimate grid plot...\n")
  output_filename <- "betaEstimate_grid.pdf"
  pdf(file.path(output_dir, output_filename), width = 12, height = 12)
  
  # Setup plot layout
  par(mfrow = c(length(t_levels), length(p_levels)),
      oma = c(7, 8, 3, 3), # bottom, left, top, right
      mar = c(1, 1.5, 0, 0))
  
  # Loop through each T level (rows) from high to low for plotting
  for (t_level in rev(t_levels)) {
    for (p_level in p_levels) {
      cell_key <- paste(p_level, t_level, sep = "_")
      cell_df <- all_beta_data[[cell_key]]
      
      # Plotting logic for one cell
      if (!is.null(cell_df) && nrow(cell_df) > 0) {
        show_x_axis <- (t_level == min(t_levels))
        show_y_axis <- (p_level == min(p_levels))
        show_legend <- (t_level == max(t_levels)) && (p_level == max(p_levels))
        
        plot_cell_beta_estimates(cell_df, conditions, condition_labels, y_range, 
                                 show_x_axis, show_y_axis, show_legend, y_axis_ticks)
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
  mtext(expression(paste("Mean posterior estimate of ", beta)), side = 2, line = 4, cex = 3, outer = TRUE)
  mtext(expression(paste("True effect size (", beta, ")")),
        side = 1, line = 5.7, cex = 2.5, outer = TRUE)
  
  dev.off()
  cat("Beta estimate grid plot saved to:", file.path(output_dir, output_filename), "\n")
  cat("\nAll plots completed!\n")
}

###################################################################
################# H E L P E R   F U N C T I O N ####################
###################################################################

# Helper function to plot beta estimates for a single cell
plot_cell_beta_estimates <- function(cell_df, conditions, condition_labels, y_range,
                                     show_x_axis, show_y_axis, show_legend, y_axis_ticks) {
  
  # Get unique beta levels
  beta_levels <- sort(unique(cell_df$beta))
  offset <- 0.01
  xlim <- c(min(beta_levels) - offset, max(beta_levels) + offset)
  
  # Create an empty plot
  plot(NA, NA, xlim = xlim, ylim = y_range, xlab = "", ylab = "", 
       xaxt = "n", yaxt = "n", bty = "o")
  
  # Add diagonal reference line (y = x)
  abline(a = 0, b = 1, lty = 2, col = "gray70", lwd = 1.5)
  
  # Setup axes
  if (show_y_axis) {
    # Use pretty() to generate nicely rounded, evenly spaced tick marks
    y_at <- pretty(y_range, n = y_axis_ticks)
    # Ensure we don't go outside the range
    y_at <- y_at[y_at >= y_range[1] & y_at <= y_range[2]]
    axis(2, at = y_at, labels = y_at, las = 1, cex.axis = 2.5)
  }
  
  if (show_x_axis) {
    # Draw ticks at the bottom margin of the plot
    axis(1, at = beta_levels, labels = FALSE, line = 0, cex.axis = 2.5)
    
    # Calculate a small offset for close labels (in data coordinates)
    x_range <- xlim[2] - xlim[1]
    label_offset <- x_range * 0.015
    
    # Manually place labels with custom positioning to avoid overlap
    for (i in seq_along(beta_levels)) {
      beta_val <- beta_levels[i]
      
      # Determine label position: 0 and 0.4 centered, 0.1 left, 0.2 right
      if (i == 1 || beta_val == max(beta_levels)) {
        # First or last: use exact position
        x_pos <- beta_val
      } else if (i == 2) {
        # Second value (0.1): offset to the left
        x_pos <- beta_val - label_offset
      } else {
        # Third value (0.2): offset to the right
        x_pos <- beta_val + label_offset
      }
      
      # Draw label using mtext (in data coordinates when outer = FALSE)
      # Use cex = 2.0 to match cex.axis = 2.5 (axis labels use different scaling)
      mtext(text = beta_val, side = 1, line = 1.5, at = x_pos, cex = 1.7, outer = FALSE)
    }
  }
  
  # Define colors, line widths, line types, and point styles for each condition
  condition_colors <- c("#d3540b", "#160f0fea", "#47D647", "#E982FF")
  widths <- c(5, 5, 4, 3)
  styles <- c(2, 1, 4, 3)
  points <- c(19, 17, 15, 18)
  
  # Calculate offset for error bars (so they don't overlap)
  n_conditions <- length(conditions)
  bar_width_offset <- (max(beta_levels) - min(beta_levels)) * 0.015  # Small offset
  offsets <- seq(-(n_conditions - 1) / 2, (n_conditions - 1) / 2, length.out = n_conditions) * bar_width_offset
  
  # Plot lines, points, and error bars for each condition
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    condition_df <- cell_df[cell_df$condition == condition & !is.na(cell_df$condition_mean), ]
    
    if (nrow(condition_df) > 0) {
      # Sort by beta level
      condition_df <- condition_df[order(condition_df$beta), ]
      
      # X positions with offset for error bars
      x_positions <- condition_df$beta + offsets[i]
      
      # Plot line connecting means
      lines(x_positions, condition_df$condition_mean, 
            col = condition_colors[i], lwd = widths[i], lty = styles[i])
      
      # Plot points at means
      points(x_positions, condition_df$condition_mean, 
             col = condition_colors[i], pch = points[i], cex = 2)
      
      # Add error bars (+/- 1 SD)
      for (j in seq_len(nrow(condition_df))) {
        mean_val <- condition_df$condition_mean[j]
        sd_val <- condition_df$condition_sd[j]
        
        if (!is.na(sd_val) && sd_val > 0) {
          # Upper whisker
          segments(x_positions[j], mean_val, 
                   x_positions[j], mean_val + sd_val,
                   col = condition_colors[i], lwd = 2)
          # Lower whisker
          segments(x_positions[j], mean_val, 
                   x_positions[j], mean_val - sd_val,
                   col = condition_colors[i], lwd = 2)
          # Horizontal caps
          cap_width <- bar_width_offset * 0.3
          segments(x_positions[j] - cap_width, mean_val + sd_val,
                   x_positions[j] + cap_width, mean_val + sd_val,
                   col = condition_colors[i], lwd = 2)
          segments(x_positions[j] - cap_width, mean_val - sd_val,
                   x_positions[j] + cap_width, mean_val - sd_val,
                   col = condition_colors[i], lwd = 2)
        }
      }
    }
  }
  
  # Add custom legend (manually positioned to reduce empty space at top)
  if (show_legend) {
    # Calculate positions within plot area (closer to top)
    x_left <- xlim[1] + (xlim[2] - xlim[1]) * 0.02  # Left margin
    y_top <- y_range[2] - (y_range[2] - y_range[1]) * 0.02  # Close to top
    y_spacing <- (y_range[2] - y_range[1]) * 0.12  # Vertical spacing between items (increased)
    
    # Draw legend items for each condition
    for (i in seq_along(conditions)) {
      y_pos <- y_top - (i - 1) * y_spacing
      
      # Draw line segment
      line_length <- (xlim[2] - xlim[1]) * 0.08
      segments(x_left*0.8, y_pos, x_left + line_length, y_pos,
               col = condition_colors[i], lwd = 2, lty = styles[i])
      
      # Draw point symbol
      points(x_left + line_length * 0.5, y_pos,
             col = condition_colors[i], pch = points[i], cex = 2)
      
      # Add text label (in black for better readability)
      text(x_left + line_length * 1.2, y_pos, condition_labels[i],
           col = "black", cex = 1.7, adj = 0)
    }
  }
}