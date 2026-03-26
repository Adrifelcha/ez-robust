###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
########################### A U C   `F U L L   G R I D   ######################
###############################################################################
# This function plots the AUC grid for a full simulation study.
# Rows are the participant levels, columns are the trial levels.
# Within each cell, we plot the AUC (y-axis) for each beta level (x-axis).
# Lines differentiate between comparison conditions (EZ/Robust x Clean/Contaminated)
plot_AUCgrid_full <- function(main_dir, output_dir, plot_by = "condition", y_range = NULL,
                              custom_title_label = NULL, run_diagnose = FALSE) {  
      # Read in the comparison conditions labels
      condition_info <- read_conditions(main_dir = main_dir)  
      conditions <- condition_info$conditions
      condition_labels <- condition_info$condition_labels

      # Read in the participant and trial levels
      size_info <- read_size(parent_dirs = file.path(main_dir, conditions))
      p_levels <- size_info$p_levels
      t_levels <- size_info$t_levels

      # Initiate prograss bar if diagnostics are not run
      progress_bar <- NULL
      progress_step <- 0
      if(!run_diagnose){
          n_cells <- length(t_levels) * length(p_levels)
          n_checks <- n_cells * length(conditions)
          progress_max <- n_checks + if (is.null(y_range)) n_checks else 0
          progress_bar <- txtProgressBar(min = 0, max = progress_max, style = 3)
      }
      if (is.null(y_range)) {
        auc_values <- numeric(0)
        for (t_level in t_levels) {
          for (p_level in p_levels) {
            for (condition in conditions) {
              pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
              file_path_list <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)
              if (length(file_path_list) > 0) {
                roc_data <- get_cellROCs(resultsFile = file_path_list[1], run_diagnose = run_diagnose)
                for (beta_level_char in names(roc_data$tpr_list)) {
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                  auc_values <- c(auc_values, auc)
                }
              }
              if (!is.null(progress_bar)) {
                progress_step <- progress_step + 1
                setTxtProgressBar(progress_bar, progress_step)
              }
            }
          }
        }
        if (length(auc_values) > 0) {
          auc_min <- min(auc_values, na.rm = TRUE)
          auc_max <- max(auc_values, na.rm = TRUE)
          auc_pad <- max((auc_max - auc_min) * 0.05, 0.01)
          y_range <- c(max(0, auc_min - auc_pad), min(1, auc_max + auc_pad))
        } else {
          y_range <- c(0.49, 1.0)
        }
      }
      output_filename <- if (is.null(custom_title_label)) {
        paste0("AUC_by_", plot_by, ".pdf")
      } else {
        paste0("AUC_by_", plot_by, "_", custom_title_label, ".pdf")
      }
      pdf(file.path(output_dir, output_filename), width = 12, height = 12)
      par(mfrow = c(length(t_levels), length(p_levels)), oma = c(7, 7, 3, 3), mar = c(1, 1.5, 0, 0))
      for(t_level in rev(t_levels)) {
          for(p_level in p_levels) {
              auc_data_list <- list()
              for(condition in conditions){
                  pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
                  file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
                  if(!is.na(file_path)) {
                      roc_data <- get_cellROCs(resultsFile = file_path, run_diagnose = run_diagnose)
                      for(beta_level_char in names(roc_data$tpr_list)){
                          auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                          auc_data_list[[length(auc_data_list) + 1]] <- data.frame(condition = condition, 
                                                                                  beta = as.numeric(beta_level_char), 
                                                                                  auc = auc, stringsAsFactors = FALSE)
                      }
                  }
                  if(!is.null(progress_bar)){
                    progress_step <- progress_step + 1
                    setTxtProgressBar(progress_bar, progress_step)
                  }
              }
              if(length(auc_data_list) > 0){
                auc_df <- do.call(rbind, auc_data_list)
                show_x_axis <- (t_level == min(t_levels))
                show_y_axis <- (p_level == min(p_levels))
                if(plot_by == "condition"){
                  show_legend <- (t_level == max(t_levels)) && (p_level == max(p_levels))
                  plot_cell_by_condition(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend)
                }else{
                  show_legend <- (t_level == max(t_levels)) && (p_level == min(p_levels))
                  highlight_cell <- ifelse(as.numeric(p_level) * as.numeric(t_level) == 6400, TRUE, FALSE)
                  plot_cell_by_beta(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend,
                                    highlight_cell)
                }
              }else{      plot.new()          }

              if(t_level == max(t_levels)){
                mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 2.5, font = 2)
              }
              if(p_level == max(p_levels)){
                mtext(paste("T =", t_level), side = 4, line = 1.85, cex = 2.5, font = 2, las = 0)
              }
          }
      }
      mtext("Area Under Curve (AUC)", side = 2, line = 3.8, cex = 2.5, outer = TRUE)
      mtext(expression(paste("Effect size (", beta, ")")), side = 1, line = 5.7, cex = 2.5, outer = TRUE)
      
      if(!is.null(progress_bar)){    close(progress_bar)   }
      
      dev.off()
      cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
}


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
########################### A U C   `N E S T E D   G R I D   ##################
###############################################################################
# This function plots the AUC grid for a nested simulation study.
# Top two rows are low drift, bottom two rows are high drift.
# Left two columns are low bound, right two columns are high bound.
# For each combination of Drift x Boundary, we have two T levels (40 and 160).
# Within each cell, we plot the AUC (y-axis) for each beta level (x-axis).
# Lines differentiate between comparison conditions (EZ/Robust x Clean/Contaminated)
plot_AUC_nested <- function(main_dir, output_dir, plot_by = "condition", y_range = NULL,
                            custom_title_label = NULL, run_diagnose = FALSE) {
  condition_info <- read_conditions(main_dir = main_dir)
  combinations <- condition_info$parameter_cells  
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
  size_info <- read_size(parent_dirs = nested_paths)
  
  p_level <- size_info$p_levels[1]
  t_levels <- size_info$t_levels
  progress_bar <- NULL
  progress_step <- 0
  if (!run_diagnose) {
    n_cells <- length(t_levels) * length(combinations)
    n_checks <- n_cells * length(conditions)
    progress_max <- n_checks + if (is.null(y_range)) n_checks else 0
    progress_bar <- txtProgressBar(min = 0, max = progress_max, style = 3)
  }
  if (is.null(y_range)) {
    auc_values <- numeric(0)
    for (t_level in t_levels) {
      for (combination in combinations) {
        combination_path <- file.path(main_dir, combination)
        for (condition in conditions) {
          condition_path <- file.path(combination_path, condition)
          if (dir.exists(condition_path)) {
            pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
            file_path_list <- list.files(condition_path, pattern = pattern, full.names = TRUE)
            if (length(file_path_list) > 0) {
              roc_data <- get_cellROCs(resultsFile = file_path_list[1], run_diagnose = run_diagnose)
              for (beta_level_char in names(roc_data$tpr_list)) {
                auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                auc_values <- c(auc_values, auc)
              }
            }
          }
          if (!is.null(progress_bar)) {
            progress_step <- progress_step + 1
            setTxtProgressBar(progress_bar, progress_step)
          }
        }
      }
    }
    if (length(auc_values) > 0) {
      auc_min <- min(auc_values, na.rm = TRUE)
      auc_max <- max(auc_values, na.rm = TRUE)
      auc_pad <- max((auc_max - auc_min) * 0.05, 0.01)
      y_range <- c(max(0, auc_min - auc_pad), min(1, auc_max + auc_pad))
    } else {
      y_range <- c(0.49, 1.0)
    }
  }
  output_filename <- if (is.null(custom_title_label)) {
    paste0("AUC_by_", plot_by, ".pdf")
  } else {
    paste0("AUC_by_", plot_by, "_", custom_title_label, ".pdf")
  }
  pdf(file.path(output_dir, output_filename), width = 11, height = 11)
  layout_matrix <- matrix(c(1, 0, 3, 2, 0, 4, 0, 0, 0, 5, 0, 7, 6, 0, 8), nrow = 5, ncol = 3, byrow = TRUE)
  layout(layout_matrix, widths = c(0.68, 0.01, 0.68), heights = c(0.85, 0.85, 0.05, 0.85, 0.85))
  par(cex = 0.85, oma = c(2.0, 3.4, 1.8, 0.8), mar = c(0.3, 2.0, 0.1, 0.3))
  boundary_levels <- c("lowBound", "highBound")
  drift_levels <- c("lowDrift", "highDrift")
  for (drift_idx in seq_along(drift_levels)) {
    drift_level <- drift_levels[drift_idx]
    for (boundary_idx in seq_along(boundary_levels)) {
      boundary_level <- boundary_levels[boundary_idx]
      combination <- paste0(drift_level, "-", boundary_level)
      combination_path <- file.path(main_dir, combination)
      for (t_idx in seq_along(t_levels)) {
        t_level <- t_levels[t_idx]
        auc_data_list <- list()
        for (condition in conditions) {
          condition_path <- file.path(combination_path, condition)
          if (dir.exists(condition_path)) {
            pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
            file_path <- list.files(condition_path, pattern = pattern, full.names = TRUE)[1]
            if (!is.na(file_path) && file.exists(file_path)) {
              roc_data <- get_cellROCs(resultsFile = file_path, run_diagnose = run_diagnose)
              for (beta_level_char in names(roc_data$tpr_list)) {
                auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                auc_data_list[[length(auc_data_list) + 1]] <- data.frame(
                  condition = condition, beta = as.numeric(beta_level_char), auc = auc, stringsAsFactors = FALSE
                )
              }
            }
          }
          if (!is.null(progress_bar)) {
            progress_step <- progress_step + 1
            setTxtProgressBar(progress_bar, progress_step)
          }
        }
        show_x_axis <- (t_level == max(t_levels)) && (drift_idx == 2)
        show_y_axis <- (boundary_idx == 1)
        if (length(auc_data_list) > 0) {
          auc_df <- do.call(rbind, auc_data_list)
          if (plot_by == "condition") {
            show_legend <- (boundary_idx == 1) && (drift_idx == 2) && (t_level == max(t_levels))
            plot_cell_by_condition(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend)
          } else {
            show_legend <- (boundary_idx == 2) && (drift_idx == 1) && (t_level == min(t_levels))
            plot_cell_by_beta(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend,
                              highlight_cell = FALSE, legend_position = "bottomright")
          }
        } else {
          plot.new()
        }
        usr <- par("usr")
        text(x = usr[1] + 0.03 * diff(usr[1:2]), y = usr[4] - 0.04 * diff(usr[3:4]),
             labels = paste("T =", t_level), adj = c(0, 1), cex = 2.2, font = 2, col = "black")
      }
    }
  }
  tryCatch({
    mtext("Area Under Curve (AUC)", side = 2, line = 3.8, cex = 2.5, outer = TRUE)
    mtext(expression(paste("Effect size (", beta, ")")), side = 1, line = 5.7, cex = 2.5, outer = TRUE)
    mtext("low bound", side = 3, line = 0.2, at = 0.25, cex = 1.8, outer = TRUE, font = 2)
    mtext("high bound", side = 3, line = 0.2, at = 0.75, cex = 1.8, outer = TRUE, font = 2)
    mtext("low drift", side = 2, line = 1.3, at = 0.75, cex = 1.8, outer = TRUE, font = 2)
    mtext("high drift", side = 2, line = 1.3, at = 0.25, cex = 1.8, outer = TRUE, font = 2)
  }, error = function(e) {
    message("Skipping outer labels due to graphics device constraints: ", e$message)
  })
  if (!is.null(progress_bar)) {
    close(progress_bar)
  }
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
    axis(2, at = y_at, round(y_at,1), las = 1, cex.axis = 1.0) 
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
                              show_y_axis, show_legend, highlight_cell = FALSE, highlight_color = "white",
                              legend_position = "topleft") {
  
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
    axis(2, at = y_at, round(y_at,1), las = 1, cex.axis = 1.0) 
  }

  if (show_x_axis) {
    #labels <- sapply(beta_levels, function(x) as.expression(bquote(beta == .(x))))
    labels <- beta_levels
    axis(1, at = beta_levels, labels = labels, line = 0.3, cex.axis = 1.3)
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
    legend(legend_position,
           legend = condition_labels,
           col = condition_colors,
           lwd = 2, pch = 19, bty = "n", cex = 1.7)
  }
}