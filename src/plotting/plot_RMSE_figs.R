###############################################################################
#################  M A I N   F U N C T I O N S   ##############################
###############################################################################
# These functions generate pdf figures showing a grid of RMSE diagnostics
# including Bias and Variance
###############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R M S E   F U L L   G R I D   ##############################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function plots the RMSE grid for a full simulation study.
# It generates three versions of the grid: RMSE, Bias, and Variance.
# Rows are the participant levels, columns are the trial levels.
# Within each cell, we plot the metric (y-axis) for each beta level (x-axis).
# Lines differentiate between comparison conditions (EZ/Robust x Clean/Contaminated)
plot_RMSE_fullGrid <- function(main_dir, output_dir, parameter = "betaweight", highlight_cell = FALSE,
                         y_range_rmse = NULL, y_range_mse = NULL, y_range_bias = NULL, y_range_variance = NULL) {
  
  # Validate inputs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  valid_params <- c("bound_mean", "drift_mean", "nondt_mean", "betaweight")
  if (!parameter %in% valid_params) {
    stop(paste("Invalid 'parameter' argument. Choose one of:", paste(valid_params, collapse = ", ")))
  }

  # Create output directory, if it doesn't exist
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
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
  
  all_rmse_data <- list()
  all_bias_data <- list()
  all_variance_data <- list()
  total_cells <- length(t_levels) * length(p_levels)

  cat("Computing metrics for all simulation cells...\n")
  pb_metrics <- txtProgressBar(min = 0, max = total_cells, style = 3, char = "=")
  current_cell <- 0

  for (t_level in t_levels) {
    for (p_level in p_levels) {
      current_cell <- current_cell + 1
      cell_key <- paste(p_level, t_level, sep = "_")
      setTxtProgressBar(pb_metrics, current_cell)

      rmse_data_list <- list()
      bias_data_list <- list()
      variance_data_list <- list()
      
      for (condition in conditions) {
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_path <- list.files(file.path(main_dir, condition), pattern = pattern, full.names = TRUE)[1]
        if (!is.na(file_path) && file.exists(file_path)) {
          rmse_data <- get_cellRMSE(resultsFile = file_path, parameter = parameter)
          for (beta_level_char in names(rmse_data$rmse_by_beta)) {
            rmse_val <- rmse_data$rmse_by_beta[[beta_level_char]]            
            bias_val <- rmse_data$bias_by_beta[[beta_level_char]]
            variance_val <- rmse_data$variance_by_beta[[beta_level_char]]
            
            if (!is.na(rmse_val)) {
              rmse_data_list[[length(rmse_data_list) + 1]] <- data.frame(
                condition = condition, beta = as.numeric(beta_level_char), rmse = rmse_val, stringsAsFactors = FALSE)
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
  close(pb_metrics)
  cat("\n")

  all_by_metric <- list(rmse = all_rmse_data, bias = all_bias_data, variance = all_variance_data)
  y_ranges <- list(rmse = y_range_rmse, bias = y_range_bias, variance = y_range_variance)
  metric_keys <- c("rmse", "bias", "variance")
  metric_labels <- c(rmse = "RMSE", bias = "Bias", variance = "Variance")

  # Y-axis limits when not supplied (RMSE/Variance: [0, max*1.1]; Bias: padded min-max)
  n_y_auto <- sum(vapply(metric_keys, function(m) is.null(y_ranges[[m]]), logical(1)))
  if (n_y_auto > 0) {
    cat("Processing plotting space\n")
    pb_y <- txtProgressBar(min = 0, max = n_y_auto, style = 3, char = "=")
    y_auto_step <- 0
  }
  for (m in metric_keys) {
    if (!is.null(y_ranges[[m]])) next
    vals <- unlist(
      lapply(all_by_metric[[m]], function(d) if (!is.null(d) && m %in% names(d)) d[[m]] else NULL),
      use.names = FALSE
    )
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) {
      y_ranges[[m]] <- if (m == "bias") c(-0.1, 0.1) else c(0, 0.01)
    } else if (m == "bias") {
      min_b <- min(vals, na.rm = TRUE)
      max_b <- max(vals, na.rm = TRUE)
      span <- max_b - min_b
      y_ranges[[m]] <- c(min_b - span * 0.1, max_b + span * 0.1)
    } else {
      y_ranges[[m]] <- c(0, max(vals, na.rm = TRUE) * 1.1)
    }
    y_auto_step <- y_auto_step + 1
    setTxtProgressBar(pb_y, y_auto_step)
  }
  if (n_y_auto > 0) {
    close(pb_y)
    cat("\n")
  }

  for (m in metric_keys) {
    lab <- metric_labels[[m]]

    add_zero_line <- identical(m, "bias")
    output_filename <- paste0(lab, "_", parameter, ".pdf")
    pdf(file.path(output_dir, output_filename), width = 12, height = 12)

    par(mfrow = c(length(t_levels), length(p_levels)),
        oma = c(7, 9, 3.5, 3.5),
        mar = c(1, 1.5, 0, 0))

    for (t_level in rev(t_levels)) {
      for (p_level in p_levels) {
        cell_key <- paste(p_level, t_level, sep = "_")
        metric_df <- all_by_metric[[m]][[cell_key]]

        if (!is.null(metric_df) && nrow(metric_df) > 0) {
          show_x_axis <- (t_level == min(t_levels))
          show_y_axis <- (p_level == min(p_levels))          
          t_index <- which(t_levels == t_level)
          p_index <- which(p_levels == p_level)

          if(highlight_cell){
            highlight_this_cell <- (t_index == p_index)
          }else{
            highlight_this_cell <- FALSE
          }
          
          legend_at = "bottomleft"
          if(m == "rmse"){
            show_legend <- (t_level == min(t_levels)) && (p_level == min(p_levels))
          }else{if(m == "bias"){
            show_legend <- (t_level == min(t_levels)) && (p_level == max(p_levels))
          }else{
            show_legend <- (t_level == max(t_levels)) && (p_level == max(p_levels))
            legend_at = "topleft"
          }}
          
          plotCell_RMSEmetric_by_beta(
            metric_df, m, conditions, condition_labels,
            y_ranges[[m]], show_x_axis, show_y_axis, show_legend,
            highlight_this_cell, add_zero_line, highlight_color = "#d3b672",
            legend_at = legend_at
          )
        } else {
          plot.new()
        }

        if (t_level == max(t_levels)) {
          mtext(paste("P =", p_level), side = 3, line = 0.5, cex = 2.5, font = 2)
        }
        if (p_level == max(p_levels)) {
          mtext(paste("T =", t_level), side = 4, line = 1.85, cex = 2.5, font = 2, las = 0)
        }
      }
    }

    mtext(plot_RMSEgrid_outer_ylab(m, parameter), side = 2, line = 4, cex = 3, outer = TRUE)
    mtext(expression(paste("True effect size (", beta, ")")), side = 1, line = 5.7, cex = 2.5, outer = TRUE)

    dev.off()
    cat(lab, "grid plot saved to:", file.path(output_dir, output_filename), "\n")
  }
}

###############################################################################
#################  N E S T E D   F U N C T I O N S   ##########################
###############################################################################
# These functions generate each cell on any of the RMSE figures
###############################################################################


# One P x T cell of the full RMSE grid: beta on x-axis, metric on y-axis, one line per condition.
plotCell_RMSEmetric_by_beta <- function(metric_df, metric_name, conditions, condition_labels, y_range,
                                        show_x_axis, show_y_axis, show_legend, highlight_cell = FALSE,
                                        highlight_color = "white", add_zero_line = FALSE, legend_at = "bottomleft") {
  
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
    axis(2, at = y_at, round(y_at, 2), las = 1, cex.axis = 2) 
  }

  if (show_x_axis) {
    labels <- beta_levels
    axis(1, at = beta_levels, labels = labels, line = 1, cex.axis = 2.1)
    axis(1, at = beta_levels[3], labels = labels[3], line = 1, cex.axis = 2.1)
  }

  
  condition_colors <- c("#d3540b", "#160f0fea", "#47D647", "#E982FF")
  widths <- c(4,4,4,5)
  styles <- c(2,1,4,3)
  point_pch <- c(19, 17, 15, 18)
  point_cex <- c(2, 2.2, 2.2, 2.75)

  # Plot lines and points for each condition
  for (i in seq_along(conditions)) {
    condition_val <- conditions[i]
    subset_df <- metric_df[metric_df$condition == condition_val, ]
    ordered_subset <- subset_df[order(subset_df$beta), ]
    
    lines(ordered_subset$beta, ordered_subset[[metric_name]], col = condition_colors[i], 
          lwd = widths[i], lty = styles[i])
    points(ordered_subset$beta, ordered_subset[[metric_name]], col = condition_colors[i], 
           pch = point_pch[i], cex = point_cex[i])
  }
  
  # Add legend to the top-right plot
  if (show_legend) {
    legend_point_cex <- c(1.8, 2.1, 2.2, 2.85)
    lines_order <- c(3, 4, 1, 2)
    legend(legend_at, legend = condition_labels[lines_order], col = condition_colors[lines_order], 
           lty = styles[lines_order], lwd = widths[lines_order]-2, pch = point_pch[lines_order], bty = "n", 
           cex = 1.4, xpd = TRUE, pt.cex = legend_point_cex[lines_order], seg.len = 4)
  }
}

###############################################################################
#################  A U X I L I A R Y   F U N C T I O N S   ####################
###############################################################################
# These auxiliary functions are used to process the data for the RMSE figures
###############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate a suitable outer y-axis label for the RMSE figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_RMSEgrid_outer_ylab <- function(metric_name, parameter) {
  key <- paste(metric_name, parameter, sep = "@")
  lab <- switch(key,
    "rmse@betaweight" = bquote("RMSE for" ~ beta),
    "rmse@drift_mean" = bquote("RMSE for" ~ mu[nu]),
    "rmse@bound_mean" = bquote("RMSE for" ~ mu[alpha]),
    "rmse@nondt_mean" = bquote("RMSE for" ~ mu[tau]),
    "bias@betaweight" = bquote("Estimation bias for" ~ beta),
    "bias@drift_mean" = bquote("Estimation bias for" ~ mu[nu]),
    "bias@bound_mean" = bquote("Estimation bias for" ~ mu[alpha]),
    "bias@nondt_mean" = bquote("Estimation bias for" ~ mu[tau]),
    "variance@betaweight" = bquote("Estimation variance for" ~ beta),
    "variance@drift_mean" = bquote("Estimation variance for" ~ mu[nu]),
    "variance@bound_mean" = bquote("Estimation variance for" ~ mu[alpha]),
    "variance@nondt_mean" = bquote("Estimation variance for" ~ mu[tau]),
    bquote("Metric for" ~ beta)
  )
  lab
}


normalize_nested_panel_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  df$beta_level <- as.numeric(as.character(df$beta_level))
  df$condition <- as.character(df$condition)
  df
}






plot_RMSE_nested_by_param <- function(main_dir, output_dir, t_level_select = 40,
                                      beta_levels_select = c(0, 0.2, 0.4),
                                      x_param = "drift_mean",
                                      parameter = "betaweight",
                                      custom_title_label = NULL,
                                      y_range_rmse = NULL, y_range_bias = NULL, y_range_variance = NULL) {
  valid_params <- c("betaweight", "drift_mean", "bound_mean", "nondt_mean")
  if (!parameter %in% valid_params) {
    stop(paste("Invalid parameter. Choose one of:", paste(valid_params, collapse = ", ")))
  }
  if (!x_param %in% c("drift_mean", "bound_mean")) {
    stop("Invalid x_param. Use \"drift_mean\" or \"bound_mean\".")
  }

  condition_info <- read_conditions(main_dir = main_dir)
  if (condition_info$layout != "nested") {
    stop("plot_RMSE_nested_by_param() requires a nested simulation layout (parameter-cell folders).")
  }
  combinations <- condition_info$parameter_cells
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
  size_info <- read_size(parent_dirs = nested_paths)
  p_level <- size_info$p_levels[1]
  t_levels <- size_info$t_levels
  t_level <- if (t_level_select %in% t_levels) t_level_select else t_levels[1]

  if (identical(x_param, "drift_mean")) {
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
  bin_spec <- get_param_bin_spec_nested(x_param_short)

  progress_bar <- NULL
  progress_step <- 0
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
  n_files <- n_files * length(beta_levels_select)
  progress_bar <- txtProgressBar(min = 0, max = max(1, n_files), style = 3)

  panel_data_raw <- list()
  all_rmse <- numeric(0)
  all_bias <- numeric(0)
  all_var <- numeric(0)
  all_x <- numeric(0)

  for (fixed_level in fixed_levels) {
    panel_key <- fixed_level
    panel_data_raw[[panel_key]] <- data.frame(
      condition = character(0), beta_level = numeric(0), x_value = numeric(0),
      rmse = numeric(0), bias = numeric(0), variance = numeric(0),
      x_level = character(0), bin_id = integer(0),
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
          for (beta_sel in beta_levels_select) {
            out <- extract_parameter_metrics_by_bins(
              file_path = file_path,
              beta_level_select = beta_sel,
              x_col = x_col,
              x_level = x_level,
              bin_spec = bin_spec,
              parameter = parameter
            )
            if (!is.null(out) && nrow(out) > 0) {
              out$condition <- condition
              out$beta_level <- as.numeric(beta_sel)
              panel_data_raw[[panel_key]] <- rbind(
                panel_data_raw[[panel_key]],
                out[, c("condition", "beta_level", "x_value", "rmse", "bias", "variance", "x_level", "bin_id")]
              )
              all_rmse <- c(all_rmse, out$rmse)
              all_bias <- c(all_bias, out$bias)
              all_var <- c(all_var, out$variance)
              all_x <- c(all_x, out$x_value)
            }
            progress_step <- progress_step + 1
            setTxtProgressBar(progress_bar, progress_step)
          }
        }
      }
    }
  }

  close(progress_bar)

  if (length(all_rmse) == 0 || length(all_x) == 0) {
    stop("No valid data points were found for the selected T and beta levels (check nested paths and RData).")
  }

  x_pad <- (max(all_x) - min(all_x)) * 0.05
  x_range <- c(min(all_x) - x_pad, max(all_x) + x_pad)
  x_ticks <- sort(unique(c(bin_spec$breaks[[x_levels[1]]], bin_spec$breaks[[x_levels[2]]])))

  panel_data <- list()
  for (fixed_level in fixed_levels) {
    dd <- panel_data_raw[[fixed_level]]
    if (nrow(dd) == 0) {
      panel_data[[fixed_level]] <- dd
    } else {
      agg <- aggregate(list(rmse = dd$rmse, bias = dd$bias, variance = dd$variance),
                       by = list(condition = dd$condition, beta_level = dd$beta_level,
                                 x_level = dd$x_level, bin_id = dd$bin_id, x_value = dd$x_value),
                       FUN = mean, na.rm = TRUE)
      panel_data[[fixed_level]] <- agg[order(agg$x_value),]
    }
  }
  panel_data <- lapply(panel_data, normalize_nested_panel_df)

  if (is.null(y_range_rmse)) {
    y_pad <- max((max(all_rmse, na.rm = TRUE) - min(all_rmse, na.rm = TRUE)) * 0.08, 0.01)
    y_range_rmse <- c(max(0, min(all_rmse, na.rm = TRUE) - y_pad), max(all_rmse, na.rm = TRUE) + y_pad)
  }
  if (is.null(y_range_bias)) {
    br <- max(all_bias, na.rm = TRUE) - min(all_bias, na.rm = TRUE)
    pad <- max(br * 0.1, 0.01)
    y_range_bias <- c(min(all_bias, na.rm = TRUE) - pad, max(all_bias, na.rm = TRUE) + pad)
  }
  if (is.null(y_range_variance)) {
    y_range_variance <- c(0, max(all_var, na.rm = TRUE) * 1.1)
  }

  y_ticks_rmse <- pretty(y_range_rmse, n = 6)
  y_ticks_bias <- pretty(y_range_bias, n = 6)
  y_ticks_var <- pretty(y_range_variance, n = 6)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  suffix <- if (is.null(custom_title_label)) "" else paste0("_", custom_title_label)
  beta_tag <- paste0("beta", paste(beta_levels_select, collapse = "_"))

  metric_specs <- list(
    list(name = "rmse", label = "RMSE", col = "rmse", y_range = y_range_rmse, y_ticks = y_ticks_rmse, zero = FALSE),
    list(name = "bias", label = "Bias", col = "bias", y_range = y_range_bias, y_ticks = y_ticks_bias, zero = TRUE),
    list(name = "variance", label = "Variance", col = "variance", y_range = y_range_variance, y_ticks = y_ticks_var, zero = FALSE)
  )

  n_beta <- length(beta_levels_select)
  pdf_height <- max(6.6, 2.45 * n_beta + 1.85)

  for (ms in metric_specs) {
    output_filename <- paste0(ms$label, "_by_", x_param, "_", parameter, "_T", t_level, ".pdf")
    output_path <- file.path(output_dir, output_filename)

    pdf(output_path, width = 10.6, height = pdf_height)
    par(mfrow = c(n_beta, 2), oma = c(3.5, 3.15, 1.95, 2.35),      
        mar = c(1.55, 2, 0.35, 0.45), mgp = c(1.5, 0.95, 0), cex = 1.0)

    for (row_idx in seq_len(n_beta)) {
      beta_row <- beta_levels_select[row_idx]
      for (col_idx in seq_along(fixed_levels)) {
        panel_key <- fixed_levels[col_idx]
        show_x_axis <- (row_idx == n_beta)
        show_y_axis <- (col_idx == 1)
        show_legend <- (row_idx == n_beta && col_idx == 2)
        plot_rmse_nested_single_panel(
          panel_df = panel_data[[panel_key]],
          conditions = conditions,
          condition_labels = condition_labels,
          beta_level = as.numeric(beta_row),
          x_range = x_range,
          y_range = ms$y_range,
          x_ticks = x_ticks,
          y_ticks = ms$y_ticks,
          metric_col = ms$col,
          show_x_axis = show_x_axis,
          show_y_axis = show_y_axis,
          show_legend = show_legend,
          add_zero_line = ms$zero,
          legend_position = "bottomleft"
        )
      }
    }

    mtext(panel_titles[1], side = 3, line = -0.3, at = 0.25, outer = TRUE, cex = 2.2, font = 2)
    mtext(panel_titles[2], side = 3, line = -0.3, at = 0.75, outer = TRUE, cex = 2.2, font = 2)
    ylab_outer <- plot_RMSEgrid_outer_ylab(ms$name, parameter)
    mtext(ylab_outer, side = 2, line = 0.7, outer = TRUE, cex = 2.4, font = 1)
    mtext(x_axis_label, side = 1, line = 2.4, outer = TRUE, cex = 2.2, font = 1)

    for (r in seq_len(n_beta)) {
      at_r <- 1 - (r - 0.5) / n_beta
      mtext(
        as.expression(substitute(beta == x, list(x = beta_levels_select[r]))),
        side = 4, line = 1.3, at = at_r, outer = TRUE, las = 3, cex = 2.4, font = 2)
    }
    dev.off()

    cat(ms$label, "nested-by-parameter plot saved to:", output_path, "\n")
  }
}

# One panel with all simulation beta levels overlaid (line type by beta, color by condition).
plot_rmse_nested_multibeta_panel <- function(panel_df, conditions, condition_labels, beta_levels_select,
                                             x_range, y_range, x_ticks, y_ticks, metric_col,
                                             show_x_axis, show_y_axis, show_legend,
                                             add_zero_line = FALSE, legend_position = "bottomleft") {
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  if (add_zero_line) {
    abline(h = 0, lty = 2, col = "gray50", lwd = 1.2)
  }
  if (show_x_axis) {
    axis(1, at = x_ticks, cex.axis = 1.15)
  }
  if (show_y_axis) {
    axis(2, at = y_ticks, las = 1, cex.axis = 1.15)
  }

  condition_colors <- c("#d3540b", "#160f0fea", "#47D647", "#E982FF")
  widths <- c(4, 4, 3, 3)
  beta_lty <- seq_along(beta_levels_select)
  pch_beta <- c(19, 17, 15, 18, 4, 8)

  tol_beta <- 1e-6
  if (nrow(panel_df) > 0) {
    panel_df$beta_level <- as.numeric(as.character(panel_df$beta_level))
    panel_df$condition <- as.character(panel_df$condition)
  }
  for (i in seq_along(conditions)) {
    cond <- as.character(conditions[i])
    for (j in seq_along(beta_levels_select)) {
      bj <- as.numeric(beta_levels_select[j])
      dd <- panel_df[abs(panel_df$beta_level - bj) < tol_beta & panel_df$condition == cond, ]
      if (nrow(dd) == 0) {
        next
      }
      for (lev in unique(dd$x_level)) {
        dd_lev <- dd[dd$x_level == lev, ]
        ord <- order(dd_lev$x_value)
        lines(
          dd_lev$x_value[ord], dd_lev[[metric_col]][ord],
          col = condition_colors[i], lwd = widths[i], lty = beta_lty[j]
        )
        points(
          dd_lev$x_value[ord], dd_lev[[metric_col]][ord],
          col = condition_colors[i],
          pch = pch_beta[((j - 1) %% length(pch_beta)) + 1],
          cex = 1.25
        )
      }
    }
  }

  if (show_legend) {
    n_b <- length(beta_levels_select)
    n_c <- length(conditions)
    leg_lab <- character(n_b * n_c)
    leg_col <- character(n_b * n_c)
    leg_lty <- integer(n_b * n_c)
    leg_pch <- integer(n_b * n_c)
    k <- 0
    for (j in seq_len(n_b)) {
      for (i in seq_len(n_c)) {
        k <- k + 1
        leg_lab[k] <- paste0(condition_labels[i], ", beta = ", beta_levels_select[j])
        leg_col[k] <- condition_colors[i]
        leg_lty[k] <- beta_lty[j]
        leg_pch[k] <- pch_beta[((j - 1) %% length(pch_beta)) + 1]
      }
    }
    legend(
      legend_position,
      legend = leg_lab,
      col = leg_col,
      lty = leg_lty,
      pch = leg_pch,
      lwd = 2,
      bty = "n",
      cex = 0.58,
      ncol = 1
    )
  }
}

# Nested layout, one row x two columns: pool all simulation beta curves in each panel.
# split_by: drift_mean -> columns low vs high drift; bound_mean -> columns low vs high bound.
# x_axis_param: drift_mean or bound_mean (binned on x, same bin rules as other nested plots).
#
# Curves for different simulation beta strata need not coincide: each point is the metric
# conditional on rows with that true betaweight. Small vertical shifts are easier to see when
# overlaid than when spread across separate rows with the same y scale (plot_RMSE_nested_by_param).
plot_RMSE_nested_onerow_by_split <- function(main_dir, output_dir, t_level_select = 40,
                                             beta_levels_select = c(0, 0.2, 0.4),
                                             split_by = "bound_mean",
                                             x_axis_param = "drift_mean",
                                             parameter = "betaweight",
                                             custom_title_label = NULL,
                                             y_range_rmse = NULL, y_range_bias = NULL, y_range_variance = NULL) {
  valid_params <- c("betaweight", "drift_mean", "bound_mean", "nondt_mean")
  if (!parameter %in% valid_params) {
    stop(paste("Invalid parameter. Choose one of:", paste(valid_params, collapse = ", ")))
  }
  if (!split_by %in% c("drift_mean", "bound_mean")) {
    stop("split_by must be \"drift_mean\" or \"bound_mean\" (defines the two columns).")
  }
  if (!x_axis_param %in% c("drift_mean", "bound_mean")) {
    stop("x_axis_param must be \"drift_mean\" or \"bound_mean\" (what is binned on the x-axis).")
  }

  condition_info <- read_conditions(main_dir = main_dir)
  if (condition_info$layout != "nested") {
    stop("plot_RMSE_nested_onerow_by_split() requires a nested simulation layout.")
  }
  combinations <- condition_info$parameter_cells
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
  size_info <- read_size(parent_dirs = nested_paths)
  p_level <- size_info$p_levels[1]
  t_levels <- size_info$t_levels
  t_level <- if (t_level_select %in% t_levels) t_level_select else t_levels[1]

  if (identical(split_by, "drift_mean")) {
    col_split_levels <- c("lowDrift", "highDrift")
    panel_title_exprs <- list(
      expression(paste("Low ", mu[nu])),
      expression(paste("High ", mu[nu]))
    )
  } else {
    col_split_levels <- c("lowBound", "highBound")
    panel_title_exprs <- list(
      expression(paste("Low ", mu[alpha])),
      expression(paste("High ", mu[alpha]))
    )
  }

  if (identical(x_axis_param, "drift_mean")) {
    x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
    x_levels_axis <- c("lowDrift", "highDrift")
  } else {
    x_axis_label <- expression(paste("True population mean boundary separation (", mu[alpha], ")"))
    x_levels_axis <- c("lowBound", "highBound")
  }
  x_axis_short <- if (identical(x_axis_param, "drift_mean")) "drift" else "bound"
  bin_spec_x <- get_param_bin_spec_nested(x_axis_short)

  infer_x_level_bins <- function(combo) {
    if (identical(x_axis_param, "drift_mean")) {
      if (grepl("lowDrift", combo)) {
        return("lowDrift")
      }
      if (grepl("highDrift", combo)) {
        return("highDrift")
      }
    } else {
      if (grepl("lowBound", combo)) {
        return("lowBound")
      }
      if (grepl("highBound", combo)) {
        return("highBound")
      }
    }
    NA_character_
  }

  n_files <- 0
  for (sk in col_split_levels) {
    matched_combos <- combinations[grepl(sk, combinations)]
    for (combo in matched_combos) {
      for (condition in conditions) {
        condition_path <- file.path(main_dir, combo, condition)
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        n_files <- n_files + length(list.files(condition_path, pattern = pattern, full.names = TRUE))
      }
    }
  }
  progress_max <- max(1, n_files * length(beta_levels_select))
  progress_bar <- txtProgressBar(min = 0, max = progress_max, style = 3)
  progress_step <- 0

  init_panel_df <- function() {
    data.frame(
      condition = character(0), beta_level = numeric(0), x_value = numeric(0),
      rmse = numeric(0), bias = numeric(0), variance = numeric(0),
      x_level = character(0), bin_id = integer(0),
      stringsAsFactors = FALSE
    )
  }
  panel_raw <- list()
  for (sk in col_split_levels) {
    panel_raw[[sk]] <- init_panel_df()
  }

  all_rmse <- numeric(0)
  all_bias <- numeric(0)
  all_var <- numeric(0)
  all_x <- numeric(0)

  for (sk in col_split_levels) {
    matched_combos <- combinations[grepl(sk, combinations)]
    for (combo in matched_combos) {
      x_level_bins <- infer_x_level_bins(combo)
      if (is.na(x_level_bins)) {
        next
      }
      for (condition in conditions) {
        condition_path <- file.path(main_dir, combo, condition)
        if (!dir.exists(condition_path)) {
          next
        }
        pattern <- paste0("sim_P", p_level, "T", t_level, "_.*\\.RData$")
        file_list <- list.files(condition_path, pattern = pattern, full.names = TRUE)
        if (length(file_list) == 0) {
          next
        }
        for (file_path in file_list) {
          for (beta_sel in beta_levels_select) {
            out <- extract_parameter_metrics_by_bins(
              file_path = file_path,
              beta_level_select = beta_sel,
              x_col = x_axis_param,
              x_level = x_level_bins,
              bin_spec = bin_spec_x,
              parameter = parameter
            )
            if (!is.null(out) && nrow(out) > 0) {
              out$condition <- condition
              out$beta_level <- as.numeric(beta_sel)
              panel_raw[[sk]] <- rbind(
                panel_raw[[sk]],
                out[, c("condition", "beta_level", "x_value", "rmse", "bias", "variance", "x_level", "bin_id")]
              )
              all_rmse <- c(all_rmse, out$rmse)
              all_bias <- c(all_bias, out$bias)
              all_var <- c(all_var, out$variance)
              all_x <- c(all_x, out$x_value)
            }
            progress_step <- progress_step + 1
            setTxtProgressBar(progress_bar, progress_step)
          }
        }
      }
    }
  }

  close(progress_bar)

  if (length(all_x) == 0 || length(all_rmse) == 0) {
    stop("No valid data for plot_RMSE_nested_onerow_by_split() (paths, T, or columns).")
  }

  agg_panel <- function(dd) {
    if (nrow(dd) == 0) {
      return(dd)
    }
    aggregate(
      list(rmse = dd$rmse, bias = dd$bias, variance = dd$variance),
      by = list(
        condition = dd$condition,
        beta_level = dd$beta_level,
        x_level = dd$x_level,
        bin_id = dd$bin_id,
        x_value = dd$x_value
      ),
      FUN = mean,
      na.rm = TRUE
    )
  }
  panel_data <- list()
  for (sk in col_split_levels) {
    dd <- agg_panel(panel_raw[[sk]])
    panel_data[[sk]] <- if (nrow(dd) > 0) dd[order(dd$x_value), ] else dd
  }
  panel_data <- lapply(panel_data, normalize_nested_panel_df)

  x_pad <- (max(all_x) - min(all_x)) * 0.05
  x_range <- c(min(all_x) - x_pad, max(all_x) + x_pad)
  x_ticks <- sort(unique(c(
    bin_spec_x$breaks[[x_levels_axis[1]]],
    bin_spec_x$breaks[[x_levels_axis[2]]]
  )))

  if (is.null(y_range_rmse)) {
    y_pad <- max((max(all_rmse, na.rm = TRUE) - min(all_rmse, na.rm = TRUE)) * 0.08, 0.01)
    y_range_rmse <- c(max(0, min(all_rmse, na.rm = TRUE) - y_pad), max(all_rmse, na.rm = TRUE) + y_pad)
  }
  if (is.null(y_range_bias)) {
    br <- max(all_bias, na.rm = TRUE) - min(all_bias, na.rm = TRUE)
    pad <- max(br * 0.1, 0.01)
    y_range_bias <- c(min(all_bias, na.rm = TRUE) - pad, max(all_bias, na.rm = TRUE) + pad)
  }
  if (is.null(y_range_variance)) {
    y_range_variance <- c(0, max(all_var, na.rm = TRUE) * 1.1)
  }

  y_ticks_rmse <- pretty(y_range_rmse, n = 6)
  y_ticks_bias <- pretty(y_range_bias, n = 6)
  y_ticks_var <- pretty(y_range_variance, n = 6)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  suffix <- if (is.null(custom_title_label)) "" else paste0("_", custom_title_label)
  beta_tag <- paste0("beta", paste(beta_levels_select, collapse = "_"))
  split_tag <- if (identical(split_by, "drift_mean")) "driftSplit" else "boundSplit"
  x_tag <- if (identical(x_axis_param, "drift_mean")) "driftX" else "boundX"

  metric_specs <- list(
    list(name = "rmse", label = "RMSE", col = "rmse", y_range = y_range_rmse, y_ticks = y_ticks_rmse, zero = FALSE),
    list(name = "bias", label = "Bias", col = "bias", y_range = y_range_bias, y_ticks = y_ticks_bias, zero = TRUE),
    list(name = "variance", label = "Variance", col = "variance", y_range = y_range_variance, y_ticks = y_ticks_var, zero = FALSE)
  )

  pdf_height <- max(4.8, 3.6)
  for (ms in metric_specs) {
    output_filename <- paste0(
      ms$label, "_", split_tag, "_", x_tag, "_", parameter, "_T", t_level, "_", beta_tag, suffix, ".pdf"
    )
    output_path <- file.path(output_dir, output_filename)

    pdf(output_path, width = 10.8, height = pdf_height)
    par(
      mfrow = c(1, 2),
      oma = c(3.2, 3.35, 2.0, 0.6),
      mar = c(1.55, 1.95, 0.45, 0.45),
      mgp = c(1.5, 0.95, 0),
      cex = 1.0
    )

    for (col_idx in seq_along(col_split_levels)) {
      sk <- col_split_levels[col_idx]
      show_x_axis <- TRUE
      show_y_axis <- (col_idx == 1)
      show_legend <- (col_idx == 2)
      plot_rmse_nested_multibeta_panel(
        panel_df = panel_data[[sk]],
        conditions = conditions,
        condition_labels = condition_labels,
        beta_levels_select = beta_levels_select,
        x_range = x_range,
        y_range = ms$y_range,
        x_ticks = x_ticks,
        y_ticks = ms$y_ticks,
        metric_col = ms$col,
        show_x_axis = show_x_axis,
        show_y_axis = show_y_axis,
        show_legend = show_legend,
        add_zero_line = ms$zero,
        legend_position = "bottomleft"
      )
    }

    at_col <- (seq_len(2) - 0.5) / 2
    mtext(panel_title_exprs[[1]], side = 3, line = 0.05, at = at_col[1], outer = TRUE, cex = 1.95, font = 2)
    mtext(panel_title_exprs[[2]], side = 3, line = 0.05, at = at_col[2], outer = TRUE, cex = 1.95, font = 2)
    ylab_outer <- plot_RMSEgrid_outer_ylab(ms$name, parameter)
    mtext(ylab_outer, side = 2, line = 1.05, outer = TRUE, cex = 2.0, font = 1)
    mtext(x_axis_label, side = 1, line = 2.25, outer = TRUE, cex = 2.0, font = 1)

    dev.off()
    cat(ms$label, "nested one-row by split saved to:", output_path, "\n")
  }

  invisible(NULL)
}

# Nested cell study: 3 rows (simulation beta), 4 columns = drift estimand (low / high bound) then
# betaweight (low / high bound). X-axis is binned true drift (mu_nu) in every panel.
# Bin edges (this figure only): five bins per drift tier, lowDrift [0,1] and highDrift [2,3].
# Writes three PDFs (RMSE, Bias, Variance). Assign rbind directly to panel_raw_*[[fixed_level]]
# so list updates are not lost (tmp <- list; tmp[[i]] <- ... can fail to update the original).
# annotate_bias_highbound_segment: if TRUE, Bias PDF only, bracket + label on high bound + drift_mean
# panels (column 2) for each beta in annotate_bias_segment_beta (default c(0.4, 0.2): rows beta 0.4 and 0.2).
# annotate_bias_segment_bin_pairs: named list by beta string; each value is c(bin_i, bin_j) indexing
# highDrift bin centers (1 = leftmost on the high segment). Default list("0.2" = c(1, 2)) for adjacent
# bins 0.2 apart on mu_nu; betas not named use c(1, 3).
# Values from annotate_bias_segment_condition (default EZ_contaminated).
# print_bias_EZ_highbound_console: if TRUE, print a short console summary:
# EZ x Contaminated only, high population bound (figure columns 2 and 4), x_level = highDrift, all bin
# centers and all simulation beta levels. For bin index > 1, append | Diff = bias minus bias at bin 1.
# Set FALSE to suppress.
# print_bias_console_condition: condition folder id to filter (default EZ_contaminated).
# Return value (invisible list): bias_highDrift$drift_mean / $betaweight / $combined; plus metadata;
# bias_EZ_highBound_highDrift = same slices as printed (for programmatic use).
plot_RMSE_meanDrift_beta_sideBySide <- function(main_dir, output_dir, t_level_select = 40,
                                                beta_levels_select = c(0, 0.2, 0.4),
                                                custom_title_label = NULL,
                                                y_range_rmse_drift = NULL, y_range_rmse_beta = NULL,
                                                y_range_bias_drift = NULL, y_range_bias_beta = NULL,
                                                y_range_var_drift = NULL, y_range_var_beta = NULL,
                                                annotate_bias_highbound_segment = FALSE,
                                                annotate_bias_segment_beta = c(0.4, 0.2),
                                                annotate_bias_segment_condition = "EZ_contaminated",
                                                annotate_bias_segment_bin_pairs = list("0.2" = c(1L, 2L)),
                                                print_bias_EZ_highbound_console = TRUE,
                                                print_bias_console_condition = "EZ_contaminated") {
  condition_info <- read_conditions(main_dir = main_dir)
  if (condition_info$layout != "nested") {
    stop("plot_bias_nested_drift_betaweight_by_bound() requires a nested simulation layout.")
  }
  combinations <- condition_info$parameter_cells
  conditions <- condition_info$conditions
  condition_labels <- condition_info$condition_labels
  nested_paths <- as.vector(outer(file.path(main_dir, combinations), conditions, file.path))
  size_info <- read_size(parent_dirs = nested_paths)
  p_level <- size_info$p_levels[1]
  t_levels <- size_info$t_levels
  t_level <- if (t_level_select %in% t_levels) t_level_select else t_levels[1]

  x_levels <- c("lowDrift", "highDrift")
  fixed_levels <- c("lowBound", "highBound")
  x_col <- "drift_mean"
  x_axis_label <- expression(paste("True population intercept drift (", mu[nu], ")"))
  # Five bins: low drift tier 0–0.2, …, 0.8–1; high drift tier 2–2.2, …, 2.8–3 (not get_param_bin_spec_nested()).
  low_breaks_biasfig <- seq(0, 1, length.out = 6)
  high_breaks_biasfig <- seq(2, 3, length.out = 6)
  bin_spec <- list(
    breaks = list(lowDrift = low_breaks_biasfig, highDrift = high_breaks_biasfig),
    centers = list(
      lowDrift = (low_breaks_biasfig[-1] + low_breaks_biasfig[-length(low_breaks_biasfig)]) / 2,
      highDrift = (high_breaks_biasfig[-1] + high_breaks_biasfig[-length(high_breaks_biasfig)]) / 2
    )
  )

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
  n_beta <- length(beta_levels_select)
  progress_max <- max(1, n_files * n_beta * 2)
  progress_bar <- txtProgressBar(min = 0, max = progress_max, style = 3)
  progress_step <- 0

  init_panel_df <- function() {
    data.frame(
      condition = character(0), beta_level = numeric(0), x_value = numeric(0),
      rmse = numeric(0), bias = numeric(0), variance = numeric(0),
      x_level = character(0), bin_id = integer(0),
      stringsAsFactors = FALSE
    )
  }
  panel_raw_drift <- list(lowBound = init_panel_df(), highBound = init_panel_df())
  panel_raw_beta <- list(lowBound = init_panel_df(), highBound = init_panel_df())
  all_x <- numeric(0)
  all_rmse_drift <- numeric(0)
  all_rmse_beta <- numeric(0)
  all_bias_drift <- numeric(0)
  all_bias_beta <- numeric(0)
  all_var_drift <- numeric(0)
  all_var_beta <- numeric(0)

  for (fixed_level in fixed_levels) {
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
          for (beta_sel in beta_levels_select) {
            for (param in c("drift_mean", "betaweight")) {
              out <- extract_parameter_metrics_by_bins(
                file_path = file_path,
                beta_level_select = beta_sel,
                x_col = x_col,
                x_level = x_level,
                bin_spec = bin_spec,
                parameter = param
              )
              if (!is.null(out) && nrow(out) > 0) {
                out$condition <- condition
                out$beta_level <- as.numeric(beta_sel)
                row_block <- out[, c("condition", "beta_level", "x_value", "rmse", "bias", "variance", "x_level", "bin_id")]
                if (identical(param, "drift_mean")) {
                  panel_raw_drift[[fixed_level]] <- rbind(panel_raw_drift[[fixed_level]], row_block)
                  all_x <- c(all_x, out$x_value)
                  all_rmse_drift <- c(all_rmse_drift, out$rmse)
                  all_bias_drift <- c(all_bias_drift, out$bias)
                  all_var_drift <- c(all_var_drift, out$variance)
                } else {
                  panel_raw_beta[[fixed_level]] <- rbind(panel_raw_beta[[fixed_level]], row_block)
                  all_x <- c(all_x, out$x_value)
                  all_rmse_beta <- c(all_rmse_beta, out$rmse)
                  all_bias_beta <- c(all_bias_beta, out$bias)
                  all_var_beta <- c(all_var_beta, out$variance)
                }
              }
              progress_step <- progress_step + 1
              setTxtProgressBar(progress_bar, progress_step)
            }
          }
        }
      }
    }
  }

  close(progress_bar)

  if (length(all_x) == 0 || length(all_rmse_drift) == 0) {
    stop("No valid drift_mean data (check nested paths, RData, and column names).")
  }
  if (length(all_rmse_beta) == 0) {
    stop("No valid betaweight data (check nested paths, RData, and column names).")
  }

  agg_one <- function(dd) {
    if (nrow(dd) == 0) {
      return(dd)
    }
    aggregate(
      list(rmse = dd$rmse, bias = dd$bias, variance = dd$variance),
      by = list(
        condition = dd$condition,
        beta_level = dd$beta_level,
        x_level = dd$x_level,
        bin_id = dd$bin_id,
        x_value = dd$x_value
      ),
      FUN = mean,
      na.rm = TRUE
    )
  }
  panel_drift <- lapply(fixed_levels, function(fl) {
    dd <- agg_one(panel_raw_drift[[fl]])
    if (nrow(dd) > 0) dd[order(dd$x_value), ] else dd
  })
  names(panel_drift) <- fixed_levels
  panel_beta <- lapply(fixed_levels, function(fl) {
    dd <- agg_one(panel_raw_beta[[fl]])
    if (nrow(dd) > 0) dd[order(dd$x_value), ] else dd
  })
  names(panel_beta) <- fixed_levels

  # Full bias summary for highDrift bins only (true drift in [2,3]): drift_mean vs betaweight estimands,
  # all bounds and simulation beta levels. Returned invisibly at the end.
  bias_highDrift_drift_mean <- {
    pieces <- list()
    for (fl in fixed_levels) {
      dd <- panel_drift[[fl]]
      dd <- dd[dd$x_level == "highDrift", ]
      if (nrow(dd) > 0) {
        dd$estimand <- "drift_mean"
        dd$bound <- fl
        pieces[[fl]] <- dd
      }
    }
    if (length(pieces) == 0L) {
      data.frame(
        condition = character(0), beta_level = numeric(0), x_value = numeric(0),
        rmse = numeric(0), bias = numeric(0), variance = numeric(0),
        x_level = character(0), bin_id = integer(0), estimand = character(0),
        bound = character(0), stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, pieces)
    }
  }
  bias_highDrift_betaweight <- {
    pieces <- list()
    for (fl in fixed_levels) {
      dd <- panel_beta[[fl]]
      dd <- dd[dd$x_level == "highDrift", ]
      if (nrow(dd) > 0) {
        dd$estimand <- "betaweight"
        dd$bound <- fl
        pieces[[fl]] <- dd
      }
    }
    if (length(pieces) == 0L) {
      data.frame(
        condition = character(0), beta_level = numeric(0), x_value = numeric(0),
        rmse = numeric(0), bias = numeric(0), variance = numeric(0),
        x_level = character(0), bin_id = integer(0), estimand = character(0),
        bound = character(0), stringsAsFactors = FALSE
      )
    } else {
      do.call(rbind, pieces)
    }
  }
  add_condition_labels <- function(df) {
    if (nrow(df) == 0) {
      df$condition_label <- character(0)
      return(df)
    }
    ib <- match(df$condition, conditions)
    df$condition_label <- ifelse(!is.na(ib), condition_labels[ib], as.character(df$condition))
    df
  }
  bias_highDrift_drift_mean <- add_condition_labels(bias_highDrift_drift_mean)
  bias_highDrift_betaweight <- add_condition_labels(bias_highDrift_betaweight)
  reorder_bias_cols <- function(df) {
    if (nrow(df) == 0) {
      return(df)
    }
    cn <- c(
      "estimand", "bound", "condition", "condition_label", "beta_level",
      "x_level", "bin_id", "x_value", "bias", "rmse", "variance"
    )
    cn <- cn[cn %in% names(df)]
    df <- df[, cn]
    df[order(df$estimand, df$bound, df$beta_level, df$condition, df$bin_id), ]
  }
  bias_highDrift_drift_mean <- reorder_bias_cols(bias_highDrift_drift_mean)
  bias_highDrift_betaweight <- reorder_bias_cols(bias_highDrift_betaweight)
  bias_highDrift_combined <- if (nrow(bias_highDrift_drift_mean) == 0L && nrow(bias_highDrift_betaweight) == 0L) {
    bias_highDrift_drift_mean
  } else if (nrow(bias_highDrift_drift_mean) == 0L) {
    bias_highDrift_betaweight
  } else if (nrow(bias_highDrift_betaweight) == 0L) {
    bias_highDrift_drift_mean
  } else {
    rbind(bias_highDrift_drift_mean, bias_highDrift_betaweight)
  }

  tol_b <- 1e-8
  cond_console <- print_bias_console_condition
  ib_lab <- match(cond_console, conditions)
  cond_label_console <- if (!is.na(ib_lab)) condition_labels[ib_lab] else cond_console

  bias_EZ_highBound_highDrift_drift <- bias_highDrift_drift_mean[
    as.character(bias_highDrift_drift_mean$condition) == cond_console &
      bias_highDrift_drift_mean$bound == "highBound",
  ]
  bias_EZ_highBound_highDrift_beta <- bias_highDrift_betaweight[
    as.character(bias_highDrift_betaweight$condition) == cond_console &
      bias_highDrift_betaweight$bound == "highBound",
  ]

  if (isTRUE(print_bias_EZ_highbound_console)) {
    cat(
      "\n",
      "========================================================================\n",
      "  Bias on the high drift segment (x_level = highDrift, true drift on [2, 3])\n",
      "  Condition: ", cond_label_console, " (", cond_console, ")\n",
      "  Figure column 2 = drift_mean estimand | column 4 = betaweight (both high population bound)\n",
      "  Bin centers on x: ",
      paste(sprintf("%.3f", bin_spec$centers$highDrift), collapse = ", "),
      "\n",
      "  Simulation betaweight (rows): ",
      paste(beta_levels_select, collapse = ", "),
      "\n",
      "========================================================================\n",
      sep = ""
    )
    betas_ord <- sort(beta_levels_select)
    for (bj in betas_ord) {
      cat("\n  --- Simulation beta = ", bj, " ---\n", sep = "")
      cat("  Column 2 (drift_mean bias, high bound)\n")
      dd <- bias_EZ_highBound_highDrift_drift[abs(bias_EZ_highBound_highDrift_drift$beta_level - bj) < tol_b, ]
      if (nrow(dd) == 0L) {
        cat("    (no rows)\n")
      } else {
        dd <- dd[order(dd$bin_id), ]
        ref_bias_d <- {
          w1 <- dd$bias[dd$bin_id == 1L]
          if (length(w1) > 0L) w1[1L] else NA_real_
        }
        for (ri in seq_len(nrow(dd))) {
          tr <- dd$x_value[ri]
          b <- dd$bias[ri]
          est <- tr + b
          diff_s <- ""
          if (!is.na(ref_bias_d) && dd$bin_id[ri] > 1L) {
            diff_s <- sprintf(" | Diff from bin 1 = %s", sprintf("%.6f", b - ref_bias_d))
          }
          cat(
            " | True drift bin center = ", sprintf("%.3f", tr),
            " | Bias = ", sprintf("%.6f", b),
            " | Estimated = ", sprintf("%.6f", est),
            diff_s,
            "\n",
            sep = ""
          )
        }
      }
      cat("  Column 4 (betaweight bias, high bound)\n")
      db <- bias_EZ_highBound_highDrift_beta[abs(bias_EZ_highBound_highDrift_beta$beta_level - bj) < tol_b, ]
      if (nrow(db) == 0L) {
        cat("    (no rows)\n")
      } else {
        db <- db[order(db$bin_id), ]
        ref_bias_b <- {
          w1 <- db$bias[db$bin_id == 1L]
          if (length(w1) > 0L) w1[1L] else NA_real_
        }
        for (ri in seq_len(nrow(db))) {
          tr_beta <- bj
          b <- db$bias[ri]
          est <- tr_beta + b
          diff_s <- ""
          if (!is.na(ref_bias_b) && db$bin_id[ri] > 1L) {
            diff_s <- sprintf(" | Diff from bin 1 = %s", sprintf("%.6f", b - ref_bias_b))
          }
          cat(
            " | True drift bin center = ", sprintf("%.3f", db$x_value[ri]),
            " | Bias = ", sprintf("%.6f", b),
            " | Estimated = ", sprintf("%.6f", est),
            diff_s,
            "\n",
            sep = ""
          )
        }
      }
    }
    cat("\n")
    flush.console()
  }

  x_pad <- (max(all_x) - min(all_x)) * 0.05
  x_range <- c(min(all_x) - x_pad, max(all_x) + x_pad)
  x_ticks <- sort(unique(c(bin_spec$breaks[[x_levels[1]]], bin_spec$breaks[[x_levels[2]]])))

  default_y_rmse_drift <- {
    pad_d <- max(diff(range(all_rmse_drift, na.rm = TRUE)) * 0.08, 0.01)
    c(max(0, min(all_rmse_drift, na.rm = TRUE) - pad_d), max(all_rmse_drift, na.rm = TRUE) + pad_d)
  }
  default_y_rmse_beta <- {
    pad_b <- max(diff(range(all_rmse_beta, na.rm = TRUE)) * 0.08, 0.01)
    c(max(0, min(all_rmse_beta, na.rm = TRUE) - pad_b), max(all_rmse_beta, na.rm = TRUE) + pad_b)
  }
  default_y_bias_drift <- {
    pad_d <- max(diff(range(all_bias_drift, na.rm = TRUE)) * 0.1, 0.01)
    c(min(all_bias_drift, na.rm = TRUE) - pad_d, max(all_bias_drift, na.rm = TRUE) + pad_d)
  }
  default_y_bias_beta <- {
    pad_b <- max(diff(range(all_bias_beta, na.rm = TRUE)) * 0.1, 0.01)
    c(min(all_bias_beta, na.rm = TRUE) - pad_b, max(all_bias_beta, na.rm = TRUE) + pad_b)
  }
  default_y_var_drift <- c(0, max(all_var_drift, na.rm = TRUE) * 1.1)
  default_y_var_beta <- c(0, max(all_var_beta, na.rm = TRUE) * 1.1)

  y_range_rmse_drift <- if (is.null(y_range_rmse_drift)) default_y_rmse_drift else y_range_rmse_drift
  y_range_rmse_beta <- if (is.null(y_range_rmse_beta)) default_y_rmse_beta else y_range_rmse_beta
  y_range_bias_drift <- if (is.null(y_range_bias_drift)) default_y_bias_drift else y_range_bias_drift
  y_range_bias_beta <- if (is.null(y_range_bias_beta)) default_y_bias_beta else y_range_bias_beta
  y_range_var_drift <- if (is.null(y_range_var_drift)) default_y_var_drift else y_range_var_drift
  y_range_var_beta <- if (is.null(y_range_var_beta)) default_y_var_beta else y_range_var_beta

  y_ticks_rmse_drift <- pretty(y_range_rmse_drift, n = 6)
  y_ticks_rmse_beta <- pretty(y_range_rmse_beta, n = 6)
  y_ticks_bias_drift <- pretty(y_range_bias_drift, n = 6)
  y_ticks_bias_beta <- pretty(y_range_bias_beta, n = 6)
  y_ticks_var_drift <- pretty(y_range_var_drift, n = 6)
  y_ticks_var_beta <- pretty(y_range_var_beta, n = 6)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  suffix <- if (is.null(custom_title_label)) "" else paste0("_", custom_title_label)
  beta_tag <- paste0("beta", paste(beta_levels_select, collapse = "_"))

  metric_specs <- list(
    list(
      name = "rmse", label = "RMSE", col = "rmse",
      y_drift = y_range_rmse_drift, y_beta = y_range_rmse_beta,
      ticks_drift = y_ticks_rmse_drift, ticks_beta = y_ticks_rmse_beta,
      zero = FALSE
    ),
    list(
      name = "bias", label = "Bias", col = "bias",
      y_drift = y_range_bias_drift, y_beta = y_range_bias_beta,
      ticks_drift = y_ticks_bias_drift, ticks_beta = y_ticks_bias_beta,
      zero = TRUE
    ),
    list(
      name = "variance", label = "Variance", col = "variance",
      y_drift = y_range_var_drift, y_beta = y_range_var_beta,
      ticks_drift = y_ticks_var_drift, ticks_beta = y_ticks_var_beta,
      zero = FALSE
    )
  )

  pdf_height <- max(6.6, 2.45 * n_beta + 1.95)

  for (ms in metric_specs) {
    output_filename <- paste0(
      ms$label, "_drift_and_betaweight_by_drift_bins_T", t_level, ".pdf"
    )
    output_path <- file.path(output_dir, output_filename)

    pdf(output_path, width = 16.2, height = pdf_height)
    par(
      mfrow = c(n_beta, 4),
      oma = c(3.5, 5.3, 5, 2.3),
      mgp = c(1.5, 0.95, 0),
      cex = 1.0
    )

    col_bound <- c("lowBound", "highBound", "lowBound", "highBound")
    col_is_drift <- c(TRUE, TRUE, FALSE, FALSE)
    mar_base <- c(1.55, 0.5, 0.35, 0.2)

    row_annot_targets <- if (isTRUE(annotate_bias_highbound_segment) && identical(ms$name, "bias")) {
      betas_req <- unique(as.numeric(annotate_bias_segment_beta))
      rows <- integer(0)
      for (bj in betas_req) {
        w <- which(abs(beta_levels_select - bj) < 1e-8)
        if (length(w) == 1L) {
          rows <- c(rows, as.integer(w))
        } else {
          warning(
            "annotate_bias_highbound_segment: annotate_bias_segment_beta value ",
            bj,
            " must match exactly one entry in beta_levels_select; that target skipped."
          )
        }
      }
      unique(sort(rows))
    } else {
      integer(0)
    }

    for (row_idx in seq_len(n_beta)) {
      beta_row <- beta_levels_select[row_idx]
      for (col_idx in seq_len(4)) {
        mar_this <- mar_base
        if (col_idx == 2) {
          mar_this[4] <- 1.4
        }
        if (col_idx == 3) {
          mar_this[2] <- 2.5
        }
        par(mar = mar_this)

        bound_key <- col_bound[col_idx]
        use_drift <- col_is_drift[col_idx]
        panel_df <- if (use_drift) panel_drift[[bound_key]] else panel_beta[[bound_key]]
        y_range <- if (use_drift) ms$y_drift else ms$y_beta
        y_ticks <- if (use_drift) ms$ticks_drift else ms$ticks_beta
        show_x_axis <- (row_idx == n_beta)
        show_y_axis <- (col_idx %in% c(1, 3))
        show_legend <- (row_idx == 1L && col_idx == 4L)
        plot_rmse_nested_single_panel(panel_df = panel_df, conditions = conditions,
                                      condition_labels = condition_labels,
                                      beta_level = as.numeric(beta_row),
                                      x_range = x_range, y_range = y_range,
                                      x_ticks = x_ticks, y_ticks = y_ticks,
                                      metric_col = ms$col, legend_cex = 1.2,
                                      show_x_axis = show_x_axis,
                                      show_y_axis = show_y_axis,
                                      show_legend = show_legend,
                                      add_zero_line = ms$zero,
                                      legend_position = "bottomleft")

        # Bias PDF only: mark two highDrift bins for one condition and this row's beta (pair from annotate_bias_segment_bin_pairs).
        if (isTRUE(annotate_bias_highbound_segment) && identical(ms$name, "bias") &&
            row_idx %in% row_annot_targets &&
            col_idx == 2L && use_drift && bound_key == "highBound") {
            tol_b <- 1e-6
            tol_x <- 1e-5
            dd_annot <- panel_drift[["highBound"]]
            dd_annot <- dd_annot[
              abs(dd_annot$beta_level - as.numeric(beta_row)) < tol_b &
                dd_annot$x_level == "highDrift" &
                as.character(dd_annot$condition) == annotate_bias_segment_condition,
            ]
            hi_centers <- bin_spec$centers$highDrift
            bin_pair_default <- c(1L, 3L)
            bin_pair <- bin_pair_default
            if (length(annotate_bias_segment_bin_pairs) > 0L) {
              for (nm in names(annotate_bias_segment_bin_pairs)) {
                if (abs(as.numeric(nm) - beta_row) < 1e-8) {
                  bin_pair <- as.integer(annotate_bias_segment_bin_pairs[[nm]])
                  break
                }
              }
            }
            if (length(bin_pair) != 2L || any(bin_pair < 1L) || any(bin_pair > length(hi_centers))) {
              bin_pair <- bin_pair_default
            }
            if (length(hi_centers) >= 1L && max(bin_pair) <= length(hi_centers)) {
              xc_mark <- hi_centers[bin_pair]
              y_mark <- vapply(xc_mark, function(xc) {
                sub <- dd_annot[abs(dd_annot$x_value - xc) < tol_x, ]
                if (nrow(sub) == 0) {
                  return(NA_real_)
                }
                mean(sub$bias, na.rm = TRUE)
              }, numeric(1))
              if (!any(is.na(y_mark))) {
                bias_diff <- y_mark[2] - y_mark[1]

                x_span <- diff(x_range)
                x_join <- max(
                  x_range[1] + 0.04 * x_span,
                  min(1.95, min(xc_mark) - 0.05 * x_span)
                )
                y_low <- min(y_mark[1], y_mark[2])
                y_high <- max(y_mark[1], y_mark[2])

                marker_color <- gray(0.2)
                # Two short horizontal segments from the selected points joined by one vertical segment.
                segments(x_join, y_mark[1], xc_mark[1], y_mark[1], col = marker_color, lwd = 1.05, lty = 1)
                segments(x_join, y_mark[2], xc_mark[2], y_mark[2], col = marker_color, lwd = 1.05, lty = 1)
                segments(x_join, y_low, x_join, y_high, col = marker_color, lwd = 1.05, lty = 1)

                text(
                  x_join - 0.012 * x_span,
                  (y_mark[1] + y_mark[2]) / 2,
                  labels = sprintf("Bias difference = %.3f", bias_diff),
                  adj = c(1, 0.5),
                  cex = 0.82,
                  col = "gray35"
                )
              }
            }
        }
      }
    }

    # Columns 1 and 2 outer label: "Estimation [metric] for drift_mean"
    out_top_line1 <- 2.2
    out_top_cex1 <- 2.7

    label_cols12 <- plot_RMSEgrid_outer_ylab(ms$name, "drift_mean")
    mtext(label_cols12, side = 3, line = out_top_line1, at = 0.25, outer = TRUE, cex = out_top_cex1, font = 2)
    
    # Columns 3 and 4 outer label: "Estimation [metric] for betaweight"
    label_cols34 <- plot_RMSEgrid_outer_ylab(ms$name, "betaweight")    
    mtext(label_cols34, side = 3, line = out_top_line1, at = 0.75, outer = TRUE, cex = out_top_cex1, font = 2)

    # Column labels: "Low" or "High" population bound
    col_at <- (seq_len(4) - 0.5) / 4
    for (j in seq_len(4)) {
      lab_col <- if (j %% 2 == 1L) {
        expression(paste("Low ", mu[alpha]))
      } else {
        expression(paste("High ", mu[alpha]))
      }
      mtext(lab_col, side = 3, line = -0.35, at = col_at[j], outer = TRUE, cex = 2, font = 2)
    }
    
    # X-axis label: "True population [x-axis parameter]"
    mtext(x_axis_label, side = 1, line = 2.5, outer = TRUE, cex = 2.8, font = 1)

    # Y-axis label for this combined figure (both drift_mean and betaweight columns).
    ylab_outer <- switch(
      ms$name,
      "rmse" = "RMSE",
      "bias" = "Estimation bias",
      "variance" = "Estimation variance",
      ms$label
    )
    mtext(ylab_outer, side = 2, line = 3.2, at = 0.5, outer = TRUE, cex = 2.8, font = 1)

    for (r in seq_len(n_beta)) {
      at_r <- 1 - (r - 0.5) / n_beta
      beta_text <- bquote("True " * beta == .(beta_levels_select[r]))
      mtext(beta_text, side = 4, line = 1.3, at = at_r, outer = TRUE, las = 3, cex = 2.2, font = 2)
    }

    dev.off()
    cat(ms$label, "drift + betaweight nested plot saved to:", output_path, "\n")
  }

  invisible(list(
    bias_highDrift = list(
      drift_mean = bias_highDrift_drift_mean,
      betaweight = bias_highDrift_betaweight,
      combined = bias_highDrift_combined
    ),
    bias_EZ_highBound_highDrift = list(
      drift_mean = bias_EZ_highBound_highDrift_drift,
      betaweight = bias_EZ_highBound_highDrift_beta
    ),
    beta_levels_select = beta_levels_select,
    highDrift_bin_centers = bin_spec$centers$highDrift,
    highDrift_bin_breaks = bin_spec$breaks$highDrift,
    t_level = t_level
  ))
}

# Bin definitions aligned with plot_AUCgrid.R::get_param_bin_spec().
get_param_bin_spec_nested <- function(x_param) {
  if (x_param %in% c("drift", "drift_mean")) {
    low_breaks <- seq(0, 1, length.out = 5)
    high_breaks <- seq(2, 3, length.out = 5)
    return(list(
      breaks = list(lowDrift = low_breaks, highDrift = high_breaks),
      centers = list(
        lowDrift = (low_breaks[-1] + low_breaks[-length(low_breaks)]) / 2,
        highDrift = (high_breaks[-1] + high_breaks[-length(high_breaks)]) / 2
      )
    ))
  }
  if (!x_param %in% c("bound", "bound_mean")) {
    stop("Invalid x_param in get_param_bin_spec_nested(). Use drift or bound.")
  }
  low_breaks <- seq(2, 2.5, length.out = 3)
  high_breaks <- seq(3.5, 4, length.out = 3)
  list(
    breaks = list(lowBound = low_breaks, highBound = high_breaks),
    centers = list(
      lowBound = (low_breaks[-1] + low_breaks[-length(low_breaks)]) / 2,
      highBound = (high_breaks[-1] + high_breaks[-length(high_breaks)]) / 2
    )
  )
}


# Bin-wise RMSE, bias, and variance for `parameter` estimates; rows still stratify by simulation betaweight.
extract_parameter_metrics_by_bins <- function(file_path, beta_level_select, x_col, x_level, bin_spec, parameter) {
  e <- new.env(parent = emptyenv())
  load(file_path, envir = e)
  if (!exists("simStudy_Beta", envir = e, inherits = FALSE)) {
    return(NULL)
  }

  sim_beta <- get("simStudy_Beta", envir = e)
  if (!parameter %in% colnames(sim_beta$true) || !parameter %in% colnames(sim_beta$estimates)) {
    return(NULL)
  }

  true_betas <- as.numeric(sim_beta$true[, "betaweight"])
  x_vals <- as.numeric(sim_beta$true[, x_col])
  true_param <- as.numeric(sim_beta$true[, parameter])
  est_param <- as.numeric(sim_beta$estimates[, parameter])

  n <- min(length(true_betas), length(x_vals), length(true_param), length(est_param))
  if (n == 0) {
    return(NULL)
  }
  true_betas <- true_betas[seq_len(n)]
  x_vals <- x_vals[seq_len(n)]
  true_param <- true_param[seq_len(n)]
  est_param <- est_param[seq_len(n)]

  breaks <- bin_spec$breaks[[x_level]]
  centers <- bin_spec$centers[[x_level]]
  if (is.null(breaks) || is.null(centers)) {
    return(NULL)
  }

  out <- data.frame(
    x_value = numeric(0), rmse = numeric(0), bias = numeric(0), variance = numeric(0),
    x_level = character(0), bin_id = integer(0),
    stringsAsFactors = FALSE
  )

  for (b in seq_len(length(breaks) - 1)) {
    in_bin <- if (b < (length(breaks) - 1)) {
      x_vals >= breaks[b] & x_vals < breaks[b + 1]
    } else {
      x_vals >= breaks[b] & x_vals <= breaks[b + 1]
    }

    idx <- in_bin & (abs(true_betas - beta_level_select) < 1e-8)
    if (!any(idx)) {
      next
    }

    true_s <- true_param[idx]
    est_s <- est_param[idx]
    valid <- !is.na(true_s) & !is.na(est_s)
    if (!any(valid)) {
      next
    }

    err <- est_s[valid] - true_s[valid]
    rmse_val <- sqrt(mean(err^2))
    bias_val <- mean(err)
    variance_val <- if (length(err) >= 2) stats::var(err) else NA_real_

    out <- rbind(out, data.frame(
      x_value = centers[b],
      rmse = rmse_val,
      bias = bias_val,
      variance = variance_val,
      x_level = x_level,
      bin_id = b,
      stringsAsFactors = FALSE
    ))
  }

  out
}

# One facet of the nested RMSE grid: single simulation beta, four conditions (color + lty for print).
plot_rmse_nested_single_panel <- function(panel_df, conditions, condition_labels, beta_level,
                                          x_range, y_range, x_ticks, y_ticks, metric_col,
                                          show_x_axis, show_y_axis, show_legend, legend_cex = 1.4,
                                          add_zero_line = FALSE, legend_position = "bottomleft",
                                          inner_ylab = NULL, draw_inner_ylab = FALSE) {
  if (nrow(panel_df) > 0) {
    panel_df$beta_level <- as.numeric(as.character(panel_df$beta_level))
    panel_df$condition <- as.character(panel_df$condition)
  }
  tol_beta <- 1e-6
  dd_all <- panel_df[abs(panel_df$beta_level - as.numeric(beta_level)) < tol_beta, ]
  plot(NA, NA, xlim = x_range, ylim = y_range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o")
  if (add_zero_line) {
    abline(h = 0, lty = 2, col = "gray50", lwd = 1.2)
  }
  if (show_x_axis) {
    axis(1, at = x_ticks, cex.axis = 1.4)
  }
  if (show_y_axis) {
    axis(2, at = y_ticks, las = 1, cex.axis = 1.4)
  }
  if (draw_inner_ylab && !is.null(inner_ylab)) {
    mtext(inner_ylab, side = 2, line = 2.45, cex = 0.95, las = 0)
  }

  # Match plotCell_RMSEmetric_by_beta: distinct lty / lwd per condition; pch 18 (diamonds)
  # needs larger cex than filled disc / square so symbols read at similar visual weight.
  condition_colors <- c("#d3540b", "#160f0fea", "#47D647", "#E982FF")
  widths <- c(3, 3, 3, 4)
  styles <- c(2, 1, 4, 3)
  point_pch <- c(19, 17, 15, 18)
  point_cex <- c(1.7, 1.9, 1.9, 2.4)

  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    dd <- dd_all[dd_all$condition == cond, ]
    if (nrow(dd) == 0) {
      next
    }

    for (lev in unique(dd$x_level)) {
      dd_lev <- dd[dd$x_level == lev, ]
      ord <- order(dd_lev$x_value)
      lines(
        dd_lev$x_value[ord], dd_lev[[metric_col]][ord],
        col = condition_colors[i], lwd = widths[i], lty = styles[i]
      )
      points(
        dd_lev$x_value[ord], dd_lev[[metric_col]][ord],
        col = condition_colors[i], pch = point_pch[i], cex = point_cex[i]
      )
    }
  }

  if (show_legend) {
    n_leg <- length(condition_labels)
    lines_order <- if (n_leg == 4L) c(3L, 4L, 1L, 2L) else seq_len(n_leg)
    legend_point_cex <- c(1.6, 2, 2, 2.6)[seq_len(n_leg)]
    legend(
      legend_position,
      legend = condition_labels[lines_order],
      col = condition_colors[lines_order],
      lty = styles[lines_order],
      lwd = widths[lines_order],
      pch = point_pch[lines_order],
      pt.cex = legend_point_cex[lines_order],
      seg.len = 2.2,
      bty = "n",
      cex = legend_cex,
      xpd = FALSE
    )
  }
}
