###############################################################################
#################  M A I N   F U N C T I O N S   ##############################
###############################################################################
# These functions generate pdf figures showing a grid of AUC plots
###############################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A U C   F U L L   G R I D   #################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function plots the AUC grid for a full simulation study.
# Rows are the participant levels, columns are the trial levels.
# Within each cell, we plot the AUC (y-axis) for each beta level (x-axis).
# Lines differentiate between comparison conditions (EZ/Robust x Clean/Contaminated)
plot_AUC_fullGrid <- function(main_dir, output_dir, y_range = NULL,
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
        paste0("AUC_by_beta.pdf")
      } else {
        paste0("AUC_by_beta_", custom_title_label, ".pdf")
      }
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      pdf(file.path(output_dir, output_filename), width = 12, height = 12)
      par(mfrow = c(length(t_levels), length(p_levels)), oma = c(6, 6, 3, 3), mar = c(1, 1.5, 0, 0))
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
                show_legend <- (t_level == max(t_levels)) && (p_level == min(p_levels))
                highlight_cell <- ifelse(as.numeric(p_level) * as.numeric(t_level) == 6400, TRUE, FALSE)
                plotCell_AUC_by_beta(auc_df, conditions, condition_labels, y_range, show_x_axis, show_y_axis, show_legend,
                                  highlight_cell, label_lines_length = 3.4, legend_cex = 1.5,
                                  label_point_cex = c(1.6, 2.1, 2.2, 2.8), y_cex = 2, x_cex = 2.1)
                
              }else{      plot.new()          }

              if(t_level == max(t_levels)){
                mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 2.5, font = 2)
              }
              if(p_level == max(p_levels)){
                mtext(paste("T =", t_level), side = 4, line = 1.85, cex = 2.5, font = 2, las = 0)
              }
          }
      }
      mtext("Area Under Curve (AUC)", side = 2, line = 3, cex = 2.5, outer = TRUE)
      mtext(expression(paste("True effect size (", beta, ")")), side = 1, line = 4.5, cex = 2.5, outer = TRUE)
      
      if(!is.null(progress_bar)){    close(progress_bar)   }
      
      dev.off()
      cat("AUC grid plot saved to:", file.path(output_dir, output_filename), "\n")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A U C   N E S T E D   2 x 2  ###############################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function plots a simplified 2x2 AUC grid for a nested simulation study
# at a single selected trial level T (default 40).
# Layout:
# - Columns: lowBound (left), highBound (right)
# - Rows:    lowDrift (top), highDrift (bottom)
# Within each cell, we plot AUC (y) against beta (x). Lines differentiate
# comparison conditions (EZ/Robust x Clean/Contaminated).
plot_AUC_2x2Grid <- function(main_dir, output_dir, y_range = NULL,
                             custom_title_label = NULL, run_diagnose = FALSE, t_level_select = 40) {
  condition_info <- read_conditions(main_dir = main_dir)
  combinations <- condition_info$parameter_cells
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
  size_info <- read_size(parent_dirs = nested_paths)
  
  # Infer P and T; select one T
  p_level <- size_info$p_levels[1]
  t_levels <- size_info$t_levels
  if (!(t_level_select %in% t_levels)) {
    warning(paste0("Requested T = ", t_level_select, " not found; using T = ", t_levels[1], " instead."))
    t_level <- t_levels[1]
  } else {
    t_level <- t_level_select
  }
  
  # Optional progress bar
  progress_bar <- NULL
  progress_step <- 0
  if (!run_diagnose) {
    # 4 cells x number of conditions (for data scan) + same again if y_range needs computing
    n_cells <- 4
    n_checks <- n_cells * length(conditions)
    progress_max <- n_checks + if (is.null(y_range)) n_checks else 0
    progress_bar <- txtProgressBar(min = 0, max = progress_max, style = 3)
  }
  
  # Compute y_range across the four panels at the selected T
  if(is.null(y_range)){
        auc_values <- numeric(0)
        for(combination in combinations){
            combination_path <- file.path(main_dir, combination)
            for(condition in conditions){
                condition_path <- file.path(combination_path, condition)

                if(dir.exists(condition_path)){
                    pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
                    file_path_list <- list.files(condition_path, pattern = pattern, full.names = TRUE)
                    if(length(file_path_list) > 0){
                        roc_data <- get_cellROCs(resultsFile = file_path_list[1], run_diagnose = run_diagnose)
                        for(beta_level_char in names(roc_data$tpr_list)){
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

        if(length(auc_values) > 0){
            auc_min <- min(auc_values, na.rm = TRUE)
            auc_max <- max(auc_values, na.rm = TRUE)
            auc_pad <- max((auc_max - auc_min) * 0.05, 0.01)
            y_range <- c(max(0, auc_min - auc_pad), min(1, auc_max + auc_pad))
        } else {
            y_range <- c(0.49, 1.0)
        }
  }
  
  # Filenames
  suffix <- if (is.null(custom_title_label)) "" else paste0("_", custom_title_label)
  base_name <- paste0("AUC_2x2_by_beta_T", t_level, suffix)
  pdf_filename <- file.path(output_dir, paste0(base_name, ".pdf"))
  
  # Make sure output dir exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  pdf(pdf_filename, width = 10, height = 8)
  # 2x2 without spacer columns/rows
  par(cex = 1.0, oma = c(4.8, 4.2, 2, 4.4), mar = c(0.6, 1.0, 0.6, 0.2))
  layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE),
         widths = c(1, 1), heights = c(1, 1))

  boundary_levels <- c("lowBound", "highBound")
  drift_levels <- c("lowDrift", "highDrift")

  for (drift_idx in seq_along(drift_levels)) {
    drift_level <- drift_levels[drift_idx]
    for (boundary_idx in seq_along(boundary_levels)) {
      boundary_level <- boundary_levels[boundary_idx]
      combination <- paste0(drift_level, "-", boundary_level)
      combination_path <- file.path(main_dir, combination)

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

      show_x_axis <- (drift_idx == 2) # bottom row only
      show_y_axis <- (boundary_idx == 1) # first column only
      if (length(auc_data_list) > 0) {
        auc_df <- do.call(rbind, auc_data_list)
        # Legend on top-right panel (bottom-right may be busy with x-axis)
        show_legend <- (drift_idx == 1) && (boundary_idx == 2)
        plotCell_AUC_by_beta(auc_df, conditions, condition_labels, y_range,
                          show_x_axis, show_y_axis, show_legend,
                          highlight_cell = FALSE, legend_position = "bottomright",
                          legend_cex = 1.4, y_cex = 1.5, x_cex = 1.7)

      } else {  plot.new()   }
    }
  }

  # Outer labels
  tryCatch({
    mtext("Area Under Curve (AUC)", side = 2, line = 2.3, cex = 2.2, outer = TRUE, font = 1)
    mtext(expression(paste("True effect size (", beta, ")")), side = 1, line = 3.8, cex = 2.2, outer = TRUE, font = 1)
    mtext(expression(paste("Low ", mu[alpha])), side = 3, line = -0.4, at = 0.25, cex = 2, outer = TRUE, font = 2)
    mtext(expression(paste("High ", mu[alpha])), side = 3, line = -0.4, at = 0.75, cex = 2, outer = TRUE, font = 2)
    mtext(expression(atop("Low", mu[nu])), side = 4, line = 0.2, at = 0.75, cex = 2, outer = TRUE, font = 2, las = 1)
    mtext(expression(atop("High", mu[nu])), side = 4, line = 0.2, at = 0.25, cex = 2, outer = TRUE, font = 2, las = 1)
  }, error = function(e) {
    message("Skipping outer labels due to graphics device constraints: ", e$message)
  })

  if (!is.null(progress_bar)) {
    close(progress_bar)
  }

  dev.off()

  cat("AUC 2x2 grid saved to:", pdf_filename, "\n")
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A U C   B Y   T R U E   M E A N   P A R A M E T E R    ##################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function creates a 1x2 nested AUC plot at one selected T and beta level.
# X-axis is one selected parameter (drift or bound), binned into:
# - 5 bins for low values
# - 5 bins for high values
# Y-axis is AUC. Curves are comparison conditions.
plot_AUC_byTrueMean <- function(main_dir, output_dir, t_level_select = 40, beta_level_select = 0.2,
                                     x_param = "drift_mean", custom_title_label = NULL, run_diagnose = FALSE) {
  if (!x_param %in% c("drift_mean", "bound_mean")) {
    stop("Invalid x_param. Use 'drift_mean' or 'bound_mean'.")
  }
  
  condition_info <- read_conditions(main_dir = main_dir)
  combinations <- condition_info$parameter_cells
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
  size_info <- read_size(parent_dirs = nested_paths)
  p_level <- size_info$p_levels[1]
  t_levels <- size_info$t_levels
  t_level <- if (t_level_select %in% t_levels) t_level_select else t_levels[1]
  
  if (x_param == "drift_mean") {
    x_levels <- c("lowDrift", "highDrift")
    fixed_levels <- c("lowBound", "highBound")
    x_col <- "drift_mean"
    panel_titles <- c(expression(paste("Low ", mu[alpha])), expression(paste("High ", mu[alpha])))
    x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
  } else {
    x_levels <- c("lowBound", "highBound")
    fixed_levels <- c("lowDrift", "highDrift")
    x_col <- "bound_mean"
    panel_titles <- c(expression(paste("Low ", mu[nu])), expression(paste("High ", mu[nu])))
    x_axis_label <- expression(paste("True population mean boundary separation (", mu[alpha], ")"))
  }
  x_param_short <- if (identical(x_param, "drift_mean")) "drift" else "bound"
  bin_spec <- get_param_bin_spec(x_param_short)
  
  progress_bar <- NULL
  progress_step <- 0
  if (!run_diagnose) {
    n_files <- 0
    for (fixed_level in fixed_levels) {
      for (x_level in x_levels) {
        combo <- combinations[grepl(fixed_level, combinations) & grepl(x_level, combinations)][1]
        for (condition in conditions) {
          condition_path <- file.path(main_dir, combo, condition)
          pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
          n_files <- n_files + length(list.files(condition_path, pattern = pattern, full.names = TRUE))
        }
      }
    }
    progress_bar <- txtProgressBar(min = 0, max = max(1, n_files), style = 3)
  }
  
  # Gather panel data (raw points first)
  panel_data_raw <- list()
  all_auc <- numeric(0)
  all_x <- numeric(0)
  for (fixed_level in fixed_levels) {
    panel_key <- fixed_level
    panel_data_raw[[panel_key]] <- data.frame(
      condition = character(0), x_value = numeric(0), auc = numeric(0), x_level = character(0),
      stringsAsFactors = FALSE
    )
    
    for (x_level in x_levels) {
      combo <- combinations[grepl(fixed_level, combinations) & grepl(x_level, combinations)][1]
      if (is.na(combo)) next
      
      for (condition in conditions) {
        condition_path <- file.path(main_dir, combo, condition)
        if (!dir.exists(condition_path)) next
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_list <- list.files(condition_path, pattern = pattern, full.names = TRUE)
        if (length(file_list) == 0) next
        
        for (file_path in file_list) {
          out <- extract_auc_param_by_bins(
            file_path = file_path,
            beta_level_select = beta_level_select,
            x_col = x_col,
            x_level = x_level,
            bin_spec = bin_spec,
            run_diagnose = run_diagnose
          )
          if (!is.null(out) && nrow(out) > 0) {
            out$condition <- condition
            panel_data_raw[[panel_key]] <- rbind(panel_data_raw[[panel_key]], out[, c("condition", "x_value", "auc", "x_level", "bin_id")])
            all_auc <- c(all_auc, out$auc)
            all_x <- c(all_x, out$x_value)
          }
          if (!is.null(progress_bar)) {
            progress_step <- progress_step + 1
            setTxtProgressBar(progress_bar, progress_step)
          }
        }
      }
    }
  }
  
  if (!is.null(progress_bar)) close(progress_bar)
  
  if (length(all_auc) == 0 || length(all_x) == 0) {
    stop("No valid data points were found for the selected T and beta level.")
  }
  
  y_pad <- max((max(all_auc) - min(all_auc)) * 0.08, 0.01)
  y_range <- c(max(0, min(all_auc) - y_pad), 1)
  x_pad <- (max(all_x) - min(all_x)) * 0.05
  x_range <- c(min(all_x) - x_pad, max(all_x) + x_pad)
  x_ticks <- sort(unique(c(bin_spec$breaks[[x_levels[1]]], bin_spec$breaks[[x_levels[2]]])))
  y_ticks <- pretty(y_range, n = 6)
  
  # Aggregate repeated files within each fixed panel: one point per condition x bin.
  panel_data <- list()
  for (fixed_level in fixed_levels) {
    dd <- panel_data_raw[[fixed_level]]
    if (nrow(dd) == 0) {
      panel_data[[fixed_level]] <- dd
    } else {
      agg <- aggregate(auc ~ condition + x_level + bin_id + x_value, data = dd, FUN = mean)
      panel_data[[fixed_level]] <- agg[order(agg$x_value), ]
    }
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  suffix <- if (is.null(custom_title_label)) "" else paste0("_", custom_title_label)
  output_filename <- paste0("AUC_by_", x_param, "_T", t_level, suffix, ".pdf")
  output_path <- file.path(output_dir, output_filename)
  
  pdf(output_path, width = 10.6, height = 6.4)
  par(mfrow = c(1, 2), oma = c(2.1, 3.2, 0.9, 0.2), mar = c(2.4, 2.4, 1.4, 0.15), cex = 1.0)
  
  for (i in seq_along(fixed_levels)) {
    panel_key <- fixed_levels[i]
    show_y_axis <- (i == 1)
    show_legend <- (i == 1)
    plotCell_AUC_by_TrueMean(panel_df = panel_data[[panel_key]], conditions = conditions,
                        condition_labels = condition_labels, x_range = x_range, y_range = y_range,
                        x_ticks = x_ticks, y_ticks = y_ticks,
                        show_y_axis = show_y_axis, show_legend = show_legend, panel_title = NULL,
                        legend_position = "bottomleft")
  }
  
  # Top panel descriptors as true outer labels (prevents clipping inside panels).
  mtext(panel_titles[1], side = 3, line = -1.2, at = 0.25, outer = TRUE, cex = 2.2, font = 2)
  mtext(panel_titles[2], side = 3, line = -1.2, at = 0.75, outer = TRUE, cex = 2.2, font = 2)
  mtext("Area Under Curve (AUC)", side = 2, line = 1.5, outer = TRUE, cex = 2.2, font = 1)
  mtext(x_axis_label, side = 1, line = 1.3, outer = TRUE, cex = 2.2, font = 1)
  dev.off()
  
  cat("AUC nested-by-parameter plot saved to:", output_path, "\n")
}


###############################################################################
#################  N E S T E D   F U N C T I O N S   ##########################
###############################################################################
# These functions generate each cell on any of the AUC figures
###############################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting AUC as a function of the true beta level across comparison conditions 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotCell_AUC_by_beta <- function(auc_df, conditions, condition_labels, y_range, show_x_axis, 
                              show_y_axis, show_legend, highlight_cell = FALSE, highlight_color = "white",
                              legend_position = "topleft", legend_cex = 1.7, y_cex = 1.0, x_cex = 1.3,
                              label_lines_length = 4, label_point_cex = c(1.6, 2, 2, 2.3)) {  
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
      axis(2, at = y_at, round(y_at,1), las = 1, cex.axis = y_cex) 
    }

    if (show_x_axis) {
      #labels <- sapply(beta_levels, function(x) as.expression(bquote(beta == .(x))))
      labels <- beta_levels
      axis(1, at = beta_levels, labels = labels, line = 0.3, cex.axis = x_cex)
    }

    
    condition_colors <- c("#d3540b", "#160f0fea", "#47D647", "#E982FF")
    widths <- c(4,4,4,5)
    styles <- c(2,1,4,3)
    points <- c(19,17,15,18)
    point_cex <- c(1.9, 2.2, 2.2, 2.5)

    # Plot lines and points for each condition
    for (i in seq_along(conditions)) {
      condition_val <- conditions[i]
      subset_df <- auc_df[auc_df$condition == condition_val, ]
      ordered_subset <- subset_df[order(subset_df$beta), ]
      
      lines(ordered_subset$beta, ordered_subset$auc, col = condition_colors[i], lwd = widths[i], lty = styles[i])
      points(ordered_subset$beta, ordered_subset$auc, col = condition_colors[i], pch = points[i], cex = point_cex[i])
    }
    
    # Add legend to the top-right plot
    if(show_legend){
        n_leg <- length(condition_labels)
        legend(legend_position, legend = condition_labels,
               col = condition_colors[seq_len(n_leg)],
               pch = points[seq_len(n_leg)],
               pt.cex = label_point_cex[seq_len(n_leg)],
               lwd = widths[seq_len(n_leg)], lty = styles[seq_len(n_leg)], seg.len = label_lines_length,
               bty = "n", cex = legend_cex)
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting AUC as a function of the true mean parameter level across comparison conditions 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotCell_AUC_by_TrueMean <- function(panel_df, conditions, condition_labels, x_range, y_range, x_ticks, y_ticks,
                                show_y_axis = TRUE, show_legend = FALSE, panel_title = NULL,
                                legend_position = "bottomright", legend_cex = 1.4) {
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  axis(1, at = x_ticks, cex.axis = 1.3)
  if (show_y_axis) axis(2, at = y_ticks, las = 1, cex.axis = 1.3)
  if (!is.null(panel_title)) title(main = panel_title, cex.main = 1.8, font.main = 2, line = 1.2)
  
  condition_colors <- c("#d3540b", "#160f0f", "#47D647", "#E982FF")
  widths <- c(4, 4, 5, 5)
  styles <- c(2, 1, 4, 3)
  points <- c(19, 17, 15, 18)
  # One cex per condition (index i); a vector here is recycled across bins within points()
  point_cex <- c(1.9, 2.2, 2.2, 2.5)
  
  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    dd <- panel_df[panel_df$condition == cond, ]
    if (nrow(dd) == 0) next
    
    # Connect only within the same low/high selected-parameter level.
    for (lev in unique(dd$x_level)) {
      dd_lev <- dd[dd$x_level == lev, ]
      ord <- order(dd_lev$x_value)
      lines(dd_lev$x_value[ord], dd_lev$auc[ord], col = condition_colors[i], lwd = widths[i], lty = styles[i])
      points(dd_lev$x_value[ord], dd_lev$auc[ord], col = condition_colors[i], pch = points[i],
             cex = point_cex[i])
    }
  }
  
  if (show_legend) {
    n_leg <- length(condition_labels)
    legend(legend_position,
           legend = condition_labels,
           col = condition_colors[seq_len(n_leg)],
           pch = points[seq_len(n_leg)],
           pt.cex = c(1.7, 2, 2, 2.2)[seq_len(n_leg)],
           lwd = widths[seq_len(n_leg)], lty = styles[seq_len(n_leg)], seg.len = 4,
           bty = "n", cex = legend_cex)
  }
}



###############################################################################
#################  A U X I L I A R Y   F U N C T I O N S   ####################
###############################################################################
# These auxiliary functions are used to process the data for the AUC figures
###############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Defining fixed parameter bins based on simulation design.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get_param_bin_spec <- function(x_param) {
  if (x_param %in% c("drift", "drift_mean")) {
    low_breaks <- seq(0, 1, length.out = 5)    # 4 bins
    high_breaks <- seq(2, 3, length.out = 5)   # 4 bins
    return(list(
      breaks = list(lowDrift = low_breaks, highDrift = high_breaks),
      centers = list(
        lowDrift = (low_breaks[-1] + low_breaks[-length(low_breaks)]) / 2,
        highDrift = (high_breaks[-1] + high_breaks[-length(high_breaks)]) / 2
      )
    ))
  }
  if (!x_param %in% c("bound", "bound_mean")) {
    stop("Invalid x_param in get_param_bin_spec(). Use 'drift_mean' or 'bound_mean'.")
  }
  low_breaks <- seq(2, 2.5, length.out = 3)    # 2 bins
  high_breaks <- seq(3.5, 4, length.out = 3)   # 2 bins
  list(
    breaks = list(lowBound = low_breaks, highBound = high_breaks),
    centers = list(
      lowBound = (low_breaks[-1] + low_breaks[-length(low_breaks)]) / 2,
      highBound = (high_breaks[-1] + high_breaks[-length(high_breaks)]) / 2
    )
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extracting multiple (x_bin_center, auc) points from one file at selected beta level.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_auc_param_by_bins <- function(file_path, beta_level_select, x_col, x_level, bin_spec, run_diagnose = FALSE) {
  e <- new.env(parent = emptyenv())
  load(file_path, envir = e)
  if (!exists("simStudy_Beta", envir = e, inherits = FALSE)) return(NULL)
  
  sim_beta <- get("simStudy_Beta", envir = e)
  true_betas <- as.numeric(sim_beta$true[, "betaweight"])
  x_vals <- as.numeric(sim_beta$true[, x_col])
  bf_results <- compute_BF_ROPE(true_betas = true_betas, chains = sim_beta$beta_chains, run_diagnose = run_diagnose)
  bf_df <- bf_results$bayes_factors
  
  # Keep aligned lengths defensively.
  n <- min(nrow(bf_df), length(true_betas), length(x_vals))
  if (n == 0) return(NULL)
  bf_df <- bf_df[seq_len(n), , drop = FALSE]
  true_betas <- true_betas[seq_len(n)]
  x_vals <- x_vals[seq_len(n)]
  
  breaks <- bin_spec$breaks[[x_level]]
  centers <- bin_spec$centers[[x_level]]
  if (is.null(breaks) || is.null(centers)) return(NULL)
  
  out <- data.frame(x_value = numeric(0), auc = numeric(0), x_level = character(0), bin_id = integer(0),
                    stringsAsFactors = FALSE)
  thresholds <- seq(-10, 10, length.out = 500)
  thresholds[1] <- -Inf
  thresholds[length(thresholds)] <- Inf
  
  for (b in seq_len(length(breaks) - 1)) {
    in_bin <- if (b < (length(breaks) - 1)) {
      x_vals >= breaks[b] & x_vals < breaks[b + 1]
    } else {
      x_vals >= breaks[b] & x_vals <= breaks[b + 1]
    }
    
    idx_null <- in_bin & (abs(true_betas - 0) < 1e-8)
    idx_alt <- in_bin & (abs(true_betas - beta_level_select) < 1e-8)
    # Enforce strict beta levels: if alt not present in this bin, skip.
    if (!any(idx_null) || !any(idx_alt)) next
    
    null_log_bf <- bf_df$log_bayes_factor[idx_null]
    alt_log_bf <- bf_df$log_bayes_factor[idx_alt]
    fpr <- vapply(thresholds, function(m) mean(null_log_bf > m, na.rm = TRUE), numeric(1))
    tpr <- vapply(thresholds, function(m) mean(alt_log_bf > m, na.rm = TRUE), numeric(1))
    auc <- get_AUC(fpr, tpr)
    
    out <- rbind(out, data.frame(
      x_value = centers[b],
      auc = auc,
      x_level = x_level,
      bin_id = b,
      stringsAsFactors = FALSE
    ))
  }
  
  out
}
