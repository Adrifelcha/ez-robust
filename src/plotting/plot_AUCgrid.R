# Main function to create the grid of AUC plots
plot_AUCgrid <- function(main_dir, output_dir, plot_by = "condition", y_range = NULL,
                         custom_title_label = NULL, run_diagnose = TRUE) {
  
  # Validate plot_by argument
  if (!plot_by %in% c("condition", "beta")) {
    stop("Invalid 'plot_by' argument. Choose 'condition' or 'beta'.")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Infer simulation layout, condition labels, and (if needed) parameter-cell folders
  condition_info <- read_conditions(main_dir = main_dir)
  if (!is.null(condition_info$message)) {
    message(condition_info$message)
  }
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  
  all_conditions_exist <- identical(condition_info$layout, "direct")
  
  if (all_conditions_exist) {
        # Infer P and T levels from all condition subfolders
        size_info <- read_size(parent_dirs = file.path(main_dir, conditions))
        p_levels <- size_info$p_levels
        t_levels <- size_info$t_levels
        progress_bar <- NULL
        progress_step <- 0
        if (!run_diagnose) {
          n_cells <- length(t_levels) * length(p_levels)
          n_checks <- n_cells * length(conditions)
          progress_max <- n_checks + if (is.null(y_range)) n_checks else 0
          progress_bar <- txtProgressBar(min = 0, max = progress_max, style = 3)
        }
        
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
                  roc_data <- get_cellROCs(resultsFile = file_path_list[1], run_diagnose = run_diagnose)
                  for (beta_level_char in names(roc_data$tpr_list)) {
                    auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                    if (auc < min_auc) min_auc <- auc
                  }
                }
                if (!is.null(progress_bar)) {
                  progress_step <- progress_step + 1
                  setTxtProgressBar(progress_bar, progress_step)
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
                roc_data <- get_cellROCs(resultsFile = file_path, run_diagnose = run_diagnose)
                for (beta_level_char in names(roc_data$tpr_list)) {
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                  auc_data_list[[length(auc_data_list) + 1]] <- data.frame(
                    condition = condition, beta = as.numeric(beta_level_char), auc = auc, stringsAsFactors = FALSE)
                }
              }
              if (!is.null(progress_bar)) {
                progress_step <- progress_step + 1
                setTxtProgressBar(progress_bar, progress_step)
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

          if (!is.null(progress_bar)) {
            close(progress_bar)
          }
          dev.off()
          cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
  } else {
    # Handle case when condition folders are nested under combination folders
    # Use inferred parameter-cell folders for nested structure
    combinations <- condition_info$parameter_cells
    if (length(combinations) == 0) {
      stop("No nested combination folders with valid condition subfolders were found.")
    }    
    
    # Reuse inferred conditions from nested combination folders
    conditions <- condition_info$conditions
    condition_labels <- condition_info$condition_labels
    
    # Infer size levels from all combination-condition folders
    nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
    size_info <- read_size(parent_dirs = nested_paths)
    if (length(size_info$p_levels) != 1) {
      stop("Nested AUC grid layout expects exactly one participant level.")
    }
    if (length(size_info$t_levels) != 2) {
      stop("Nested AUC grid layout expects exactly two trial levels.")
    }
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
                roc_data <- get_cellROCs(resultsFile = file_path_list[1], run_diagnose = run_diagnose)
                for (beta_level_char in names(roc_data$tpr_list)) {
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                  if (auc < min_auc) min_auc <- auc
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
      y_range <- c(0.49, 1.0)
    }
    
    # Define PDF output file name
    if(is.null(custom_title_label)){
      output_filename <- paste0("AUC_by_", plot_by, ".pdf")
    } else {
      output_filename <- paste0("AUC_by_", plot_by, "_", custom_title_label, ".pdf")
    }
    # Use a standard device size and scale plot elements within it
    pdf(file.path(output_dir, output_filename), width = 11, height = 11)
    
    # Setup layout: 4 rows (2 boundary levels x 2 T levels) x 2 columns (drift levels)
    # Loop order: drift -> boundary -> T level
    # Actual layout:
    #   Row 1: lowDrift-lowBound T=40 (col 1), lowDrift-highBound T=40 (col 2)
    #   Row 2: lowDrift-lowBound T=160 (col 1), lowDrift-highBound T=160 (col 2)
    #   Row 3: highDrift-lowBound T=40 (col 1), highDrift-highBound T=40 (col 2)
    #   Row 4: highDrift-lowBound T=160 (col 1), highDrift-highBound T=160 (col 2)
    # Supercells (2x2 grid):
    #   Top-left: lowDrift-lowBound (rows 1-2, col 1)
    #   Top-right: lowDrift-highBound (rows 1-2, col 2)
    #   Bottom-left: highDrift-lowBound (rows 3-4, col 1)
    #   Bottom-right: highDrift-highBound (rows 3-4, col 2)
    # Each supercell contains 2 panels: T=40 (top) and T=160 (bottom)
    # Layout with explicit spacer row/column to separate drift blocks and bound columns
    # 5x3 grid where row 3 and column 2 are blank spacers (0)
    layout_matrix <- matrix(c(
      1, 0, 3,   # Row 1: lowDrift-lowBound T=40, lowDrift-highBound T=40
      2, 0, 4,   # Row 2: lowDrift-lowBound T=160, lowDrift-highBound T=160
      0, 0, 0,   # Spacer between lowDrift and highDrift
      5, 0, 7,   # Row 4: highDrift-lowBound T=40, highDrift-highBound T=40
      6, 0, 8    # Row 5: highDrift-lowBound T=160, highDrift-highBound T=160
    ), nrow = 5, ncol = 3, byrow = TRUE)
    layout(layout_matrix, widths = c(0.8, 0.1, 0.8), heights = c(0.85, 0.85, 0.2, 0.85, 0.85))
    
    # Keep margins compact, but reserve more room for y-axis and outer labels
    par(cex = 0.85, oma = c(2.0, 2.8, 1.8, 0.8), mar = c(1.0, 2.2, 0.5, 0.3))
    
    # Define boundary and drift levels for ordering
    boundary_levels <- c("lowBound", "highBound")
    drift_levels <- c("lowDrift", "highDrift")
    
    # Loop order must match layout: drift -> boundary -> T level
    plot_idx <- 0
    for (drift_idx in seq_along(drift_levels)) {
      drift_level <- drift_levels[drift_idx]
      for (boundary_idx in seq_along(boundary_levels)) {
        boundary_level <- boundary_levels[boundary_idx]
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
                roc_data <- get_cellROCs(resultsFile = file_path, run_diagnose = run_diagnose)
                for (beta_level_char in names(roc_data$tpr_list)) {
                  auc <- get_AUC(roc_data$fpr, roc_data$tpr_list[[beta_level_char]])
                  auc_data_list[[length(auc_data_list) + 1]] <- data.frame(
                    condition = condition, beta = as.numeric(beta_level_char), auc = auc, stringsAsFactors = FALSE)
                }
              }
            }
            if (!is.null(progress_bar)) {
              progress_step <- progress_step + 1
              setTxtProgressBar(progress_bar, progress_step)
            }
          }
          
          # Determine which axes to show
          show_x_axis <- (t_level == max(t_levels))  # Show X axis labels only on bottom row
          show_y_axis <- (boundary_idx == 1)  # Label y-axis only on the first column
          
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
          
          # Add T-level label inside each panel (top-left corner)
          usr <- par("usr")
          text(x = usr[1] + 0.03 * diff(usr[1:2]),
               y = usr[4] - 0.04 * diff(usr[3:4]),
               labels = paste("T =", t_level),
               adj = c(0, 1), cex = 1.3, font = 2, col = "black")
        }
      }
    }
    
    # Draw borders around the four supercells (2x2 grid of drift x boundary)
    # Each supercell spans 2 rows x 1 column
    # Get outer boundaries of each supercell region
    
    # Save current par settings
    old_par <- par(no.readonly = TRUE)
    # Border/overlay geometry can fail on tight devices; don't abort full plot.
    tryCatch({
    
    # Get coordinates from plots at the edges of each supercell
    # Layout is 5 rows × 3 columns (with spacer row/column)
    # Supercell 1 (lowDrift-lowBound): rows 1-2, col 1 - Top-left
    # Top edge: plot 1 (row 1, col 1)
    par(mfg = c(1, 1, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc1_left <- grconvertX(0, from = "npc", to = "ndc")
    sc1_top <- grconvertY(1, from = "npc", to = "ndc")
    sc1_right <- grconvertX(1, from = "npc", to = "ndc")
    
    # Bottom edge: plot 2 (row 2, col 1)
    par(mfg = c(2, 1, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc1_bottom <- grconvertY(0, from = "npc", to = "ndc")
    
    # Supercell 2 (lowDrift-highBound): rows 1-2, col 2 - Top-right
    # Top edge: plot 3 (row 1, col 2)
    par(mfg = c(1, 3, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc2_left <- grconvertX(0, from = "npc", to = "ndc")
    sc2_top <- grconvertY(1, from = "npc", to = "ndc")
    sc2_right <- grconvertX(1, from = "npc", to = "ndc")
    
    # Bottom edge: plot 4 (row 2, col 2)
    par(mfg = c(2, 3, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc2_bottom <- grconvertY(0, from = "npc", to = "ndc")
    
    # Supercell 3 (highDrift-lowBound): rows 3-4, col 1 - Bottom-left
    # Top edge: plot 5 (row 3, col 1)
    par(mfg = c(4, 1, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc3_left <- grconvertX(0, from = "npc", to = "ndc")
    sc3_top <- grconvertY(1, from = "npc", to = "ndc")
    sc3_right <- grconvertX(1, from = "npc", to = "ndc")
    
    # Bottom edge: plot 6 (row 4, col 1)
    par(mfg = c(5, 1, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc3_bottom <- grconvertY(0, from = "npc", to = "ndc")
    
    # Supercell 4 (highDrift-highBound): rows 3-4, col 2 - Bottom-right
    # Top edge: plot 7 (row 3, col 2)
    par(mfg = c(4, 3, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc4_left <- grconvertX(0, from = "npc", to = "ndc")
    sc4_top <- grconvertY(1, from = "npc", to = "ndc")
    sc4_right <- grconvertX(1, from = "npc", to = "ndc")
    
    # Bottom edge: plot 8 (row 4, col 2)
    par(mfg = c(5, 3, 5, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sc4_bottom <- grconvertY(0, from = "npc", to = "ndc")
    
    # Calculate padding to place borders outside panels
    padding <- 0.015  # Padding in normalized device coordinates
    
    # Draw borders in device coordinates
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xpd = NA)
    plot.new()
    
    # Draw borders around each supercell with padding (as lines, not filled rectangles)
    # Top-left: lowDrift-lowBound
    # Top border
    segments(x0 = sc1_left - padding, y0 = sc1_top + padding, 
             x1 = sc1_right + padding, y1 = sc1_top + padding,
             col = "black", lwd = 3)
    # Bottom border
    segments(x0 = sc1_left - padding, y0 = sc1_bottom - padding, 
             x1 = sc1_right + padding, y1 = sc1_bottom - padding,
             col = "black", lwd = 3)
    # Left border
    segments(x0 = sc1_left - padding, y0 = sc1_bottom - padding, 
             x1 = sc1_left - padding, y1 = sc1_top + padding,
             col = "black", lwd = 3)
    # Right border
    segments(x0 = sc1_right + padding, y0 = sc1_bottom - padding, 
             x1 = sc1_right + padding, y1 = sc1_top + padding,
             col = "black", lwd = 3)
    
    # Top-right: highDrift-lowBound
    # Top border
    segments(x0 = sc2_left - padding, y0 = sc2_top + padding, 
             x1 = sc2_right + padding, y1 = sc2_top + padding,
             col = "black", lwd = 3)
    # Bottom border
    segments(x0 = sc2_left - padding, y0 = sc2_bottom - padding, 
             x1 = sc2_right + padding, y1 = sc2_bottom - padding,
             col = "black", lwd = 3)
    # Left border
    segments(x0 = sc2_left - padding, y0 = sc2_bottom - padding, 
             x1 = sc2_left - padding, y1 = sc2_top + padding,
             col = "black", lwd = 3)
    # Right border
    segments(x0 = sc2_right + padding, y0 = sc2_bottom - padding, 
             x1 = sc2_right + padding, y1 = sc2_top + padding,
             col = "black", lwd = 3)
    
    # Bottom-left: lowDrift-highBound
    # Top border
    segments(x0 = sc3_left - padding, y0 = sc3_top + padding, 
             x1 = sc3_right + padding, y1 = sc3_top + padding,
             col = "black", lwd = 3)
    # Bottom border
    segments(x0 = sc3_left - padding, y0 = sc3_bottom - padding, 
             x1 = sc3_right + padding, y1 = sc3_bottom - padding,
             col = "black", lwd = 3)
    # Left border
    segments(x0 = sc3_left - padding, y0 = sc3_bottom - padding, 
             x1 = sc3_left - padding, y1 = sc3_top + padding,
             col = "black", lwd = 3)
    # Right border
    segments(x0 = sc3_right + padding, y0 = sc3_bottom - padding, 
             x1 = sc3_right + padding, y1 = sc3_top + padding,
             col = "black", lwd = 3)
    
    # Bottom-right: highDrift-highBound
    # Top border
    segments(x0 = sc4_left - padding, y0 = sc4_top + padding, 
             x1 = sc4_right + padding, y1 = sc4_top + padding,
             col = "black", lwd = 3)
    # Bottom border
    segments(x0 = sc4_left - padding, y0 = sc4_bottom - padding, 
             x1 = sc4_right + padding, y1 = sc4_bottom - padding,
             col = "black", lwd = 3)
    # Left border
    segments(x0 = sc4_left - padding, y0 = sc4_bottom - padding, 
             x1 = sc4_left - padding, y1 = sc4_top + padding,
             col = "black", lwd = 3)
    # Right border
    segments(x0 = sc4_right + padding, y0 = sc4_bottom - padding, 
             x1 = sc4_right + padding, y1 = sc4_top + padding,
             col = "black", lwd = 3)
    
    # Restore par for labels
    par(old_par)
    
    # Add labels for supercells - bound labels on top margin, drift labels on left margin
    # Top margin: Bound labels (one per column)
    # Convert supercell center positions to normalized figure coordinates
    # Low drift: left supercell column (col 1 in layout) - position at horizontal center of left column
    low_drift_center_x <- (sc1_left + sc1_right) / 2
    low_drift_x_nfc <- grconvertX(low_drift_center_x, from = "ndc", to = "nfc")
    mtext("low bound", side = 3, line = 1.1, at = low_drift_x_nfc, 
          cex = 1.8, outer = TRUE, font = 2)
    
    # High drift: right supercell column (col 2 in layout) - position at horizontal center of right column
    high_drift_center_x <- (sc2_left + sc2_right) / 2
    high_drift_x_nfc <- grconvertX(high_drift_center_x, from = "ndc", to = "nfc")
    mtext("high bound", side = 3, line = 1.1, at = high_drift_x_nfc, 
          cex = 1.8, outer = TRUE, font = 2)
    
    # Left margin: Drift labels (one per supercell row)
    # Low drift: top supercell row (rows 1-2 in layout) - position at vertical center of top row
    low_bound_center_y <- (sc1_top + sc1_bottom) / 2
    low_bound_y_nfc <- grconvertY(low_bound_center_y, from = "ndc", to = "nfc")
    mtext("low drift", side = 2, line = 0.6, at = low_bound_y_nfc, 
          cex = 1.8, outer = TRUE, font = 2)
    
    # High drift: bottom supercell row (rows 3-4 in layout) - position at vertical center of bottom row
    high_bound_center_y <- (sc3_top + sc3_bottom) / 2
    high_bound_y_nfc <- grconvertY(high_bound_center_y, from = "ndc", to = "nfc")
    mtext("high drift", side = 2, line = 0.6, at = high_bound_y_nfc, 
          cex = 1.8, outer = TRUE, font = 2)
    }, error = function(e) {
      par(old_par)
      message("Skipping supercell border overlay due to device size constraints: ", e$message)
    })
    
    # Add common outer labels (guard against invalid graphics state on tight devices)
    tryCatch({
      mtext("Area Under Curve (AUC)", side = 2, line = 3.8, cex = 2.5, outer = TRUE)
      mtext(expression(paste("Effect size (", beta, ")")),
            side = 1, line = 5.7, cex = 2.5, outer = TRUE)
    }, error = function(e) {
      message("Skipping outer labels due to graphics device constraints: ", e$message)
    })
    
    if (!is.null(progress_bar)) {
      close(progress_bar)
    }
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