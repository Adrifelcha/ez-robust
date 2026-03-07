# Main function to create the grid of AUC plots
plot_AUCgrid <- function(main_dir, output_dir, plot_by = "condition", y_range = NULL,
                         custom_title_label = NULL) {
  
  # Validate plot_by argument
  if (!plot_by %in% c("condition", "beta")) {
    stop("Invalid 'plot_by' argument. Choose 'condition' or 'beta'.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define conditions and labels
  conditions <- c("EZ_contaminated", "EZRobust_contaminated", "EZ_clean", "EZRobust_clean")
  condition_labels <- c("EZ x Contaminated", "Robust x Contaminated", "EZ x Clean","Robust x Clean")
  
  # Check if condition folders exist directly under main_dir
  condition_paths <- file.path(main_dir, conditions)
  all_conditions_exist <- all(dir.exists(condition_paths))
  
  if (all_conditions_exist) {
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
                    auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                    if (auc < min_auc) min_auc <- auc
                  }
                }
              }
            }
          }
          y_range <- c(0.49, 1.0)
        }
        
        # Define PDF output file name
        if(is.null(custom_title_label)){
          output_filename <- paste0("AUC_by_", plot_by, ".pdf")
        } else {
          output_filename <- paste0("AUC_by_", plot_by, "_", custom_title_label, ".pdf")
        }
        pdf(file.path(output_dir, output_filename), width = 12, height = 12)
        
        # Setup plot layout
        par(mfrow = c(length(t_levels), length(p_levels)),
            oma = c(7, 7, 3, 3), # bottom, left, top, right
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
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
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
              
              if (plot_by == "condition") {
                show_legend <- (t_level == max(t_levels)) && (p_level == max(p_levels))
                plot_cell_by_condition(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend)
              } else {
                show_legend <- (t_level == max(t_levels)) && (p_level == min(p_levels))
                #show_legend <- (as.numeric(t_level) == 80) && (p_level == min(p_levels))
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
              mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 2.5, font = 2)
            }
            if (p_level == max(p_levels)) {
              mtext(paste("T =", t_level), side = 4, line = 1.85, cex = 2.5, font = 2, las = 0)
            }
          }
        }
        
        # Add common outer labels
        mtext("Area Under Curve (AUC)", side = 2, line = 3.8, cex = 2.5, outer = TRUE)
        mtext(expression(paste("Effect size (", beta, ")")),
              side = 1, line = 5.7, cex = 2.5, outer = TRUE)

          dev.off()
          cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
  } else {
    # Handle case when condition folders are nested under combination folders
    # Define combination folders: drift-boundary combinations
    combinations <- c("lowDrift-lowBound", "lowDrift-highBound", 
                      "highDrift-lowBound", "highDrift-highBound")
    
    # Fixed parameters for this structure
    p_level <- 160
    t_levels <- c(40, 160)
    
    # Determine y-axis range if not provided
    if (is.null(y_range)) {
      min_auc <- 1.0
      for (t_level in t_levels) {
        for (combination in combinations) {
          combination_path <- file.path(main_dir, combination)
          for (condition in conditions) {
            condition_path <- file.path(combination_path, condition)
            if (dir.exists(condition_path)) {
              pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
              file_path_list <- list.files(condition_path, pattern = pattern, full.names = TRUE)
              if (length(file_path_list) > 0) {
                roc_data <- get_cellROCs(resultsFile = file_path_list[1])
                for (beta_level_char in names(roc_data$tpr_list)) {
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                  if (auc < min_auc) min_auc <- auc
                }
              }
            }
          }
        }
      }
      y_range <- c(0.49, 1.0)
    }
    
    # Define PDF output file name
    if(is.null(custom_title_label)){
      output_filename <- paste0("AUC_by_", plot_by, ".pdf")
    } else {
      output_filename <- paste0("AUC_by_", plot_by, "_", custom_title_label, ".pdf")
    }
    pdf(file.path(output_dir, output_filename), width = 12, height = 12)
    
    # Setup layout: 4 rows (2 boundary levels x 2 T levels) x 2 columns (drift levels)
    # Row 1: lowBound, lowDrift (T=40, T=160)
    # Row 2: lowBound, highDrift (T=40, T=160)
    # Row 3: highBound, lowDrift (T=40, T=160)
    # Row 4: highBound, highDrift (T=40, T=160)
    # Layout matrix: each main cell has 2 subplots (T=40 top, T=160 bottom)
    layout_matrix <- matrix(c(1, 3,    # lowBound-lowDrift T=40, lowBound-highDrift T=40
                               2, 4,    # lowBound-lowDrift T=160, lowBound-highDrift T=160
                               5, 7,    # highBound-lowDrift T=40, highBound-highDrift T=40
                               6, 8),   # highBound-lowDrift T=160, highBound-highDrift T=160
                             nrow = 4, ncol = 2, byrow = TRUE)
    layout(layout_matrix, widths = c(1, 1), heights = c(1, 1, 1, 1))
    
    # Set margins with space for labels
    par(oma = c(7, 7, 5, 3), mar = c(2, 2, 2, 1))
    
    # Define boundary and drift levels for ordering
    boundary_levels <- c("lowBound", "highBound")
    drift_levels <- c("lowDrift", "highDrift")
    
    # Loop order must match layout: boundary -> drift -> T level
    plot_idx <- 0
    for (boundary_idx in seq_along(boundary_levels)) {
      boundary_level <- boundary_levels[boundary_idx]
      for (drift_idx in seq_along(drift_levels)) {
        drift_level <- drift_levels[drift_idx]
        combination <- paste0(drift_level, "-", boundary_level)
        combination_path <- file.path(main_dir, combination)
        
        for (t_idx in seq_along(t_levels)) {
          t_level <- t_levels[t_idx]
          plot_idx <- plot_idx + 1
          
          # Collect AUC data for the current cell
          auc_data_list <- list()
          for (condition in conditions) {
            condition_path <- file.path(combination_path, condition)
            if (dir.exists(condition_path)) {
              pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
              file_path <- list.files(condition_path, pattern = pattern, full.names = TRUE)[1]
              if (!is.na(file_path) && file.exists(file_path)) {
                roc_data <- get_cellROCs(resultsFile = file_path)
                for (beta_level_char in names(roc_data$tpr_list)) {
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                  auc_data_list[[length(auc_data_list) + 1]] <- data.frame(
                    condition = condition, beta = as.numeric(beta_level_char), auc = auc, stringsAsFactors = FALSE)
                }
              }
            }
          }
          
          # Determine which axes to show
          show_x_axis <- TRUE  # Always show X axis
          show_y_axis <- (t_level == max(t_levels)) && (drift_idx == 1)  # Only T=160 plots in left column
          
          # Plotting logic for one cell
          if (length(auc_data_list) > 0) {
            auc_df <- do.call(rbind, auc_data_list)
            
            # Determine legend position
            if (plot_by == "condition") {
              show_legend <- (boundary_idx == 1) && (drift_idx == 2) && (t_level == max(t_levels))
              plot_cell_by_condition(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend)
            } else {
              show_legend <- (boundary_idx == 1) && (drift_idx == 1) && (t_level == max(t_levels))
              highlight_cell <- FALSE
              plot_cell_by_beta(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend,
                                highlight_cell)
            }
          } else {
            plot.new() # Draw an empty plot if no data
          }
          
          # Add T-level label on each subplot (top-left corner)
          mtext(paste("T =", t_level), side = 3, line = 0.5, adj = 0, cex = 1.8, col = "black", font = 2)
        }
      }
    }
    
    # Add labels for the four main cells (boundary and drift combinations)
    # Boundary labels span horizontally across both columns
    # Low boundary (top half, rows 1-2)
    mtext("Low boundary", side = 3, line = 2.5, at = 0.5, cex = 2.5, outer = TRUE, font = 2)
    # High boundary (bottom half, rows 3-4)
    mtext("High boundary", side = 1, line = 1, at = 0.5, cex = 2.5, outer = TRUE, font = 2)
    
    # Drift labels span vertically across all rows
    # Low drift (left column)
    mtext("Low drift", side = 2, line = 1, at = 0.5, cex = 2.5, outer = TRUE, font = 2)
    # High drift (right column)
    mtext("High drift", side = 4, line = 1, at = 0.5, cex = 2.5, outer = TRUE, font = 2)
    
    # Add common outer labels
    mtext("Area Under Curve (AUC)", side = 2, line = 3.8, cex = 2.5, outer = TRUE)
    mtext(expression(paste("Effect size (", beta, ")")),
          side = 1, line = 5.7, cex = 2.5, outer = TRUE)
    
    dev.off()
    cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
  }
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
                              show_y_axis, show_legend, highlight_cell = FALSE, highlight_color = "white") {
  
  beta_levels <- sort(unique(auc_df$beta))
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

  if (show_y_axis) {
    y_at <- seq(y_range[1], y_range[2], length.out = 6)
    axis(2, at = y_at, round(y_at,1), las = 1, cex.axis = 2.5) 
  }

  if (show_x_axis) {
    #labels <- sapply(beta_levels, function(x) as.expression(bquote(beta == .(x))))
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
    subset_df <- auc_df[auc_df$condition == condition_val, ]
    ordered_subset <- subset_df[order(subset_df$beta), ]
    
    lines(ordered_subset$beta, ordered_subset$auc, col = condition_colors[i], lwd = widths[i], lty = styles[i])
    points(ordered_subset$beta, ordered_subset$auc, col = condition_colors[i], pch = points[i], cex = 2)
  }
  
  # Add legend to the top-right plot
  if (show_legend) {
    legend("topleft",
           legend = condition_labels,
           col = condition_colors,
           lwd = 2, pch = 19, bty = "n", cex = 1.7)
  }
}