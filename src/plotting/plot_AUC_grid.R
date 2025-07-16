plot_AUCgrid_byCondition <- function(main_dir, output_dir, y_range = NULL) {
  
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
  
  # Filter and keep only files that match the expected naming scheme
  rdata_files <- all_files[grepl("_P\\d+T\\d+_", basename(all_files))]
    
  # Get the P and T values from the filenames
  filenames <- basename(rdata_files)
  p_values <- as.numeric(sub(".*_P(\\d+)T.*", "\\1", filenames))
  t_values <- as.numeric(sub(".*T(\\d+)_.*", "\\1", filenames))
  # Get the unique P and T values, sorted ascending
  p_levels <- sort(unique(p_values))
  t_levels <- sort(unique(t_values))
  
  # Determine overall y-axis range if not provided
  if (is.null(y_range)) {
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
    y_range <- c(min_auc * 0.98, 1.0)
  }
  
  # Define PDF output file
  output_filename <- paste0("AUC_grid.pdf")
  pdf(file.path(output_dir, output_filename), width = 12, height = 12)
  
  # Setup plot layout
  par(mfrow = c(length(t_levels), length(p_levels)),
      oma = c(4, 4, 3, 3), # bottom, left, top, right
      mar = c(1, 1, 0, 0))

  # Loop through each T level (rows) from high to low
  for (t_level in rev(t_levels)) {
    # Loop through each P level (columns)
    for (p_level in p_levels) {
      
      auc_data_list <- list()
      
      # Collect AUC data for the current cell from all conditions
      for (condition in conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path_list <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)
        
        if (length(file_path_list) > 0) {
          file_path <- file_path_list[1]
          roc_data <- get_cellROCs(resultsFile = file_path)
          
          for (beta_level_char in names(roc_data$tpr_list)) {
            tpr <- roc_data$tpr_list[[beta_level_char]]
            auc <- compute_AUC(roc_data$fpr, tpr)
            
            auc_data_list[[length(auc_data_list) + 1]] <- data.frame(
              condition = condition,
              beta = as.numeric(beta_level_char),
              auc = auc,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      
      x_at <- seq(1.1,3.7,length.out = length(condition_labels))
      # Plotting logic for one cell
      if (length(auc_data_list) > 0) {
        auc_df <- do.call(rbind, auc_data_list)
        
        show_x_axis <- (t_level == min(t_levels))
        show_y_axis <- (p_level == min(p_levels))
        
        # Create an empty plot
        plot(NA, NA, xlim = c(1, 4), ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty="o")
        
        if (show_y_axis) {
          axis(2, las = 1)
        }
        
        if (show_x_axis) {
          # Remove axis labels to add them manually
          axis(1, at = x_at, labels = FALSE)
          
          # Split labels for multi-line display
          split_labels <- strsplit(condition_labels, " x ")
          
          # Manually add labels on two lines using mtext
          # Adjust 'line' to position text, 'cex' for font size.
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
        if ((t_level == max(t_levels)) && (p_level == max(p_levels))) {
          legend("bottomright",
                 legend = sapply(beta_levels, function(x) as.expression(bquote(beta == .(x)))),
                 col = beta_colors,
                 lwd = 2,
                 pch = 19,
                 bty = "n",
                 cex = 1.2)
        }
        
      } else {
        plot.new() # Draw an empty plot if no data
      }
      
      # Add Participant (column) labels on top of the first row
      if (t_level == max(t_levels)) {
        mtext(paste("P =", p_level), side = 3, line = 1, cex = 1.5, font = 2)
      }
      
      # Add Trial (row) labels on the right of the last column
      if (p_level == max(p_levels)) {
        mtext(paste("T =", t_level), side = 4, line = 1, cex = 1.5, font = 2, las = 0)
      }
    }
  }
  
  # Add common labels
  mtext("Area Under Curve (AUC)", side = 2, line = 2, cex = 1.75, outer = TRUE)
  
  dev.off()
  cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
}

#Example usage (commented out)# 
 main_dir <- here("output", "simStudy_results")
 output_dir <- here("output", "figures_AUC_grid")
# 
#plot_AUCgrid_byCondition(main_dir = main_dir, output_dir = output_dir) 