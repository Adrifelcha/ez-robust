# Main function to create the grid of AUC plots
plot_AUCgrid <- function(main_dir, output_dir, plot_by = "condition", y_range = NULL) {
  
  # Validate plot_by argument
  if (!plot_by %in% c("condition", "beta")) {
    stop("Invalid 'plot_by' argument. Choose 'condition' or 'beta'.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define conditions and labels
  conditions <- c("EZ_clean", "EZ_contaminated", "EZRobust_clean", "EZRobust_contaminated")
  condition_labels <- c("EZ x Clean", "EZ x Outliers", "Robust x Clean", "Robust x Outliers")
  
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
  
  # Determine y-axis range if not provided
  if (is.null(y_range)) {
    # (calculation logic remains the same)
    min_auc <- 1.0
    for (t_level in t_levels) {
      for (p_level in p_levels) {
        for (condition in conditions) {
          pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
          file_path_list <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)
          if (length(file_path_list) > 0) {
            roc_data <- get_cellROCs(resultsFile = file_path_list[1])
            for (beta_level_char in names(roc_data$tpr_list)) {
              auc <- compute_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
              if (auc < min_auc) min_auc <- auc
            }
          }
        }
      }
    }
    y_range <- c(0.49, 1.0)
  }
  
  # Define PDF output file
  output_filename <- paste0("AUC_grid_by_", plot_by, ".pdf")
  pdf(file.path(output_dir, output_filename), width = 12, height = 12)
  
  # Setup plot layout
  par(mfrow = c(length(t_levels), length(p_levels)),
      oma = c(3.5, 5.5, 3, 3), # bottom, left, top, right
      mar = c(1, 1.5, 0, 0))
  
  # Loop through each T level (rows) from high to low
  for (t_level in rev(t_levels)) {
    for (p_level in p_levels) {
      
      # Collect AUC data for the current cell
      auc_data_list <- list()
      for (condition in conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path)) {
          roc_data <- get_cellROCs(resultsFile = file_path)
          for (beta_level_char in names(roc_data$tpr_list)) {
            auc <- compute_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
            auc_data_list[[length(auc_data_list) + 1]] <- data.frame(
              condition = condition, beta = as.numeric(beta_level_char), auc = auc, stringsAsFactors = FALSE)
          }
        }
      }
      
      # Plotting logic for one cell
      if (length(auc_data_list) > 0) {
        auc_df <- do.call(rbind, auc_data_list)
        show_x_axis <- (t_level == min(t_levels))
        show_y_axis <- (p_level == min(p_levels))
        show_legend <- (t_level == max(t_levels)) && (p_level == max(p_levels))
        
        if (plot_by == "condition") {
          plot_cell_by_condition(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend)
        } else {
          highlight_cell <- ifelse(as.numeric(p_level) * as.numeric(t_level) == 6400, 
                                   TRUE, FALSE)
          plot_cell_by_beta(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend,
                            highlight_cell)
        }
        
      } else {
        plot.new() # Draw an empty plot if no data
      }
      
      # Add cell labels
      if (t_level == max(t_levels)) {
        mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 2, font = 2)
      }
      if (p_level == max(p_levels)) {
        mtext(paste("T =", t_level), side = 4, line = 1.5, cex = 2, font = 2, las = 0)
      }
    }
  }
  
  # Add common outer labels
  mtext("Area Under Curve (AUC)", side = 2, line = 2.8, cex = 2, outer = TRUE)

  dev.off()
  cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
}

#######################################################################
################# P L O T     T Y P E  ################################
################## F U N C T I O N S ###################################
#######################################################################

# BY CONDITION: The factor levels on the X axis are the conditions
plot_cell_by_condition <- function(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend) {
  # Create an empty plot
  plot(NA, NA, xlim = c(1, 4), ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty="o")
  
  if(show_y_axis){ 
    y_at <- seq(y_range[1], y_range[2], length.out = 6)
    axis(2, at = y_at, round(y_at,1), las = 1) 
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
  
  beta_levels <- sort(unique(auc_df$beta))
  beta_colors <- c("#2ea02d", "#de8520", "#d62728", "#9467bd")
  
  # Plot lines and points for each beta level
  for (i in seq_along(beta_levels)) {
    beta_val <- beta_levels[i]
    subset_df <- auc_df[auc_df$beta == beta_val, ]
    ordered_subset <- subset_df[match(conditions, subset_df$condition), ]
    
    lines(x_at, ordered_subset$auc, col = beta_colors[i], lwd = 2)
    points(x_at, ordered_subset$auc, col = beta_colors[i], pch = 19, cex = 1.2)
  }
  
  # Add legend to the top-right plot
  if (show_legend) {
    legend("bottomright",
           legend = sapply(beta_levels, function(x) as.expression(bquote(beta == .(x)))),
           col = beta_colors,
           lwd = 2, pch = 19, bty = "n", cex = 1.2)
  }
}

# Helper function to plot a single grid cell with Beta levels on the X-axis
plot_cell_by_beta <- function(auc_df, conditions, condition_labels, y_range, show_x_axis, 
                              show_y_axis, show_legend, highlight_cell = FALSE, highlight_color = "#dedb9c") {
  
  beta_levels <- sort(unique(auc_df$beta))
  offset <- 0.05
  xlim <- c(min(beta_levels) - offset, max(beta_levels) + offset)
  # Create an empty plot
  plot(NA, NA, xlim = xlim, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty="o")

  if(highlight_cell){
    # Add colored background using polygon
    polygon(x = c(xlim[1], xlim[2], xlim[2], xlim[1]), 
            y = c(y_range[1], y_range[1], y_range[2], y_range[2]), 
            col = highlight_color, border = NA)
  }

  if (show_y_axis) {
    y_at <- seq(y_range[1], y_range[2], length.out = 6)
    axis(2, at = y_at, round(y_at,1), las = 1, cex.axis = 1.5) 
  }

  if (show_x_axis) {
    labels <- sapply(beta_levels, function(x) as.expression(bquote(beta == .(x))))
    axis(1, at = beta_levels, labels = labels, line = 1, cex.axis = 1.5)
  }

  
  condition_colors <- c("#129412", "#d3540b", "#de77f3", "#160f0fea")
  widths <- c(5,5,4,3)
  styles <- c(1,3,2,4)
  points <- c(19,17,15,18)

  # Plot lines and points for each condition
  for (i in seq_along(conditions)) {
    condition_val <- conditions[i]
    subset_df <- auc_df[auc_df$condition == condition_val, ]
    ordered_subset <- subset_df[order(subset_df$beta), ]
    
    lines(ordered_subset$beta, ordered_subset$auc, col = condition_colors[i], lwd = widths[i], lty = styles[i])
    points(ordered_subset$beta, ordered_subset$auc, col = condition_colors[i], pch = points[i], cex = 2)
  }
  
  # Add legend to the top-right plot
  if (show_legend) {
    legend("bottomright",
           legend = condition_labels,
           col = condition_colors,
           lwd = 2, pch = 19, bty = "n", cex = 1.5)
  }
}