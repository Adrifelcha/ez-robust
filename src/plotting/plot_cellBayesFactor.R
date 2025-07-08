# Function to plot Bayes Factors as jittered dots across beta levels
plot_BF_jitter <- function(resultsFile = NULL, show_x_axis = TRUE, 
show_y_axis = TRUE, output_dir = NULL) {

  # Load results file
  load(resultsFile)

  # Extract true betas and chains
  true_betas <- as.vector(unlist(simStudy_Beta$true[,"betaweight"]))
  chains <- simStudy_Beta$beta_chains

  # Compute Bayes Factors
  bf_results <- compute_BF_ROPE(true_betas, chains)
  log_bf <- bf_results$bayes_factors$log_bayes_factor

  # Get unique beta levels and sort them
  beta_levels <- sort(unique(bf_results$bayes_factors$true_beta))

  # Plotting settings
  spacing <- 2.5  # Space between beta levels
  xlim <- c((spacing/2), length(beta_levels) * spacing + (spacing/2))
  y_range <- range(log_bf[!is.infinite(log_bf)])
  n_inf <- sum(is.infinite(log_bf))

  # Extend y-range slightly for better visualization
  y_range <- y_range + c(-0.05, 0.05) * diff(y_range)
  
  if (!is.null(output_dir)) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Generate PDF filename
    filename <- basename(resultsFile)
    pdfFile <- ifelse(grepl("-Outliers\\.RData$", filename), 
                      gsub("-Outliers\\.RData$", "Outlier_BF.pdf", filename),
                      gsub("-Clean\\.RData$", "Clean_BF.pdf", filename))
    pdfFile <- file.path(output_dir, pdfFile)
    pdf(pdfFile)
    
    # Set margins
    par(mar = c(4, 3.5, 2, 0.5))
    par(oma = c(0, 0, 0, 0.5))
  }
  
  # Set up the plotting area
  plot(NA, NA, xlim = xlim, ylim = c(y_range[1], y_range[2] + 0.5), ann = F, bty = "n", axes = F)
  
  # Add background polygons for negative and positive regions
  # Negative region (dark gray)
  polygon(c(xlim[1], xlim[2], xlim[2], xlim[1]), 
          c(y_range[1], y_range[1], 0, 0), 
          col = "darkgray", border = NA)
  
  # Positive region (light gray)
  polygon(c(xlim[1], xlim[2], xlim[2], xlim[1]), 
          c(0, 0, y_range[2] + 0.5, y_range[2] + 0.5), 
          col = "lightgray", border = NA)
  
  # Add axes
  if (show_x_axis) {
    axis(1, at = seq(1, length(beta_levels)) * spacing, labels = beta_levels, line = -1)
  }
  if (show_y_axis) {
    y_axis <- seq(y_range[1], y_range[2], length.out = 7)
    # Draw axis without labels first
    axis(2, at = c(y_axis, y_range[2] + 0.5), labels = FALSE, las = 2, line = -1)
    
    # Add labels manually with different colors
    text(par("usr")[1] - 0.2, y_axis, labels = sprintf("%.1f", y_axis), 
         las = 2, col = "black", xpd = TRUE)
    text(par("usr")[1] - 0.2, y_range[2] + 0.5, labels = "Inf", 
         las = 2, col = "red", xpd = TRUE)
  }
  
  if (!is.null(output_dir)) {
    mtext("True fixed effect", 1, line = 1.75, f = 2, cex = 1.2)
    mtext("Log Bayes Factor", 2, line = 1.75, f = 2, cex = 1.2)
    
    title_part <- gsub("\\.RData$", "", basename(resultsFile))
    mtext(title_part, 3, line = 0, adj = 1, cex = 1.2, col = "darkgray")
    text(xlim[1] + 1.2, y_range[2] - 0.3, expression(BF["01"]), cex = 3, col = "gray20")
  }
  
  # Add jittered dots for each beta level
  for (i in seq_along(beta_levels)) {
    current_level <- beta_levels[i]
    level_data <- bf_results$bayes_factors[bf_results$bayes_factors$true_beta == current_level, ]
    
    # Position for this beta level
    x_pos <- i * spacing
    
    # Add jittered points
    n_points <- nrow(level_data)
    jitter_width <- 0.8  # Width of jitter
    
    for (j in 1:n_points) {
      # Add horizontal jitter
      jitter_x <- runif(1, -jitter_width/2, jitter_width/2)
      
      # Add the point
      if (!is.infinite(level_data$log_bayes_factor[j])) {
      points(x_pos + jitter_x, level_data$log_bayes_factor[j], 
             col = i, pch = 16, cex = 0.6)
      }else{
        points(x_pos + jitter_x, y_range[2] + 0.5, 
             col = "red", pch = 16, cex = 0.6)
      }
    }
    
    # Add horizontal line at log(BF) = 0 (BF = 1, no evidence either way)
    segments(x_pos - 0.4, 0, x_pos + 0.4, 0, 
             col = "red", lwd = 2, lty = 2)
  }
  
  # Add reference line at log(BF) = 0
  abline(h = 0, col = "red", lwd = 2, lty = 2)
  
  if (!is.null(output_dir)) {
    dev.off()
    cat("Bayes Factor jitter plot saved to:", pdfFile, "\n")
  }
}

# Example usage:
resultsFile <- here("output", "simStudy_results", "EZ_clean", "sim_P20T20_ez-Clean.RData")
output_dir <- here("output", "figures_BayesFactors", "EZ_clean")
plot_BF_jitter(resultsFile = resultsFile, output_dir = output_dir)

