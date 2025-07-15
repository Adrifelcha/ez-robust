# Function to plot ROC curves for each beta level against Beta = 0
plot_cellROC <- function(resultsFile = NULL, epsilon = 0.05, levelsM = 1000, 
                         show_legend = TRUE, output_dir = NULL, bty = "n",
                         show_x_axis = TRUE, show_y_axis = TRUE) {
  
  # Compute ROC data
  roc_data <- get_cellROCs(resultsFile = resultsFile, epsilon = epsilon, levelsM = levelsM)
  
  # Plotting settings
  if (!is.null(output_dir)) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Generate PDF filename
    filename <- basename(resultsFile)
    pdfFile <- ifelse(grepl("-Outliers\\.RData$", filename), 
                      gsub("-Outliers\\.RData$", "Outlier_ROC.pdf", filename),
                      gsub("-Clean\\.RData$", "Clean_ROC.pdf", filename))
    pdfFile <- file.path(output_dir, pdfFile)
    pdf(pdfFile)
    
    # Set margins
    par(mar = c(4, 4, 2, 1))
    par(oma = c(0, 0, 0, 0.5))
  }
  
  # Set up the plotting area
  plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), 
       ann = FALSE, axes = FALSE, type = "n")
  box(bty = bty)
  
  # Define colors for different beta levels
  beta_colors <- c("#555555", "#2ea02d", "#de8520", "#d62728", "#9467bd")
  
  # Plot ROC curves for each beta level
  for (i in seq_along(roc_data$beta_levels)) {
    beta_level <- roc_data$beta_levels[i]
    tpr <- roc_data$tpr_list[[as.character(beta_level)]]
    
    # Plot the ROC curve (FPR on x-axis, TPR on y-axis)
    lines(roc_data$fpr, tpr, col = beta_colors[i+1], lwd = 3, lty = 1)
  }
  
  # Add diagonal reference line (random classifier)
  abline(0, 1, col = "gray", lty = 2, lwd = 1)
  
  # Add axes
  if(show_x_axis){
  axis(1, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 0.8)
  }
  if(show_y_axis){
    axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), cex.axis = 0.8, las = 2)
  } 
  # Add legend
  if (show_legend) {
    # Create legend labels with proper beta symbol
    legend_labels <- sapply(roc_data$beta_levels, function(x) 
      as.expression(bquote(beta == .(format(x, digits = 1)))))
    
    legend("bottomright", 
           legend = legend_labels,
           col = beta_colors[2:(length(roc_data$beta_levels) + 1)], 
           lwd = 3, lty = 1, cex = 1.3, bty = "n")
  }
  
  # Add title
  if (!is.null(output_dir)) {
     # Add axis labels
    mtext(expression(paste("False Positive Rate (", beta, " = 0)")), 1, line = 2.5, cex = 1.1, font = 2)
    mtext("True Positive Rate", 2, line = 2.5, cex = 1.1, font = 2)
  
    title_part <- gsub("\\.RData$", "", basename(resultsFile))
    mtext(title_part, 3, line = 0, adj = 1, cex = 1.2, col = "darkgray")
  }
  
  if (!is.null(output_dir)) {
    dev.off()
    cat("ROC plot saved to:", pdfFile, "\n")
  }
}