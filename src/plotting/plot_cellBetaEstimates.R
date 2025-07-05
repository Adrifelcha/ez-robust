plot_cellBetaEstimates <- function(resultsFile, show_frequency_bars = FALSE,
                                   show_x_axis = TRUE, show_y_axis = TRUE, y_range = NA,
                                   single_plot = FALSE) {
  # Extract the filename part from the full path
  filename <- basename(resultsFile)
  # Remove the .RData extension
  title_part <- gsub("\\.RData$", "", filename)

  # Identify true beta levels
  load(resultsFile)
  betas <- as.vector(unlist(simStudy_Beta$true[,"betaweight"]))
  beta_levels <- sort(unique(betas))  # Sort in ascending order
  
  # Get betaweight estimates (mean posteriors)
  estimates <- unlist(simStudy_Beta$estimates[,"betaweight"])
  
  # Plotting settings
  spacing <- 2.5  # Space between distributions
  last_x_pos <- length(beta_levels) * spacing
  xlim <- c((spacing/2), last_x_pos + (spacing/2))
  if(is.na(y_range)){
    y_range <- c(-0.4, 1)
  }

  if(single_plot){
    # Set margins (bottom, left, top, right)
    par(mar = c(4, 4, 2, 0.5))  # Inner margins
    par(oma = c(2, 1, 0, 0.5))  # Outer margins
  }

  # Set up the plotting area to be square with no margin lines
  plot(NA, NA, xlim = xlim, ylim = y_range, ann = F, bty = "n", axes = F) 

  if(show_x_axis){
      axis(1, at = seq(1, length(beta_levels)) * spacing, labels = beta_levels, line = 1)
  }
  if(show_y_axis){  
      y_axis <- seq(y_range[1], y_range[2], length.out = 7)
      axis(2, at = y_axis, labels = sprintf("%.2f", y_axis), las=2, line = 0)  
  }
  if(single_plot){
      mtext("True fixed effect",1, line=3.5, f=2, cex=1.2)
      mtext("Recovered value (Mean posterior)",2, line=3.5, f=2, cex=1.2)
      text(xlim[1] + 0.75, y_range[2] - 0.1, expression(beta), cex = 4, col = "darkgray")
      mtext(title_part, 3, line=0, adj=1, cex=1.2, col="darkgray")
  }

  # Create violin plots with density curves and jittered points for each beta level
  for (i in seq_along(beta_levels)) {
    current_level <- beta_levels[i]
    level_estimates <- estimates[betas == current_level]
    
    # Calculate density
    dens <- density(level_estimates, adjust = 1.5)
    
    # Scale the density to preserve actual variability
    # Use a fixed scale factor instead of normalizing to max
    scale_factor <- 0.3  # Adjust this to control overall width
    scaled_dens <- dens$y * scale_factor
    
    # Position for this beta level
    x_pos <- i * spacing
    
    # Add frequency bars if requested - draw them AFTER the violin plot
    if (show_frequency_bars) {
      # Create histogram data for frequency bars
      hist_data <- hist(level_estimates, plot = FALSE, breaks = 20)  # Reduced breaks for larger bars
      
      # Scale the frequency bars - make them much more visible
      max_freq <- max(hist_data$counts)
      scaled_freq <- hist_data$counts / max_freq * 1  # Increased from 0.3 to 0.5
      
      # Draw frequency bars as filled rectangles
      for (j in seq_along(hist_data$mids)) {
        # Calculate bin boundaries
        bin_left <- hist_data$breaks[j]
        bin_right <- hist_data$breaks[j + 1]
        bar_width <- scaled_freq[j]
        
        # Draw filled rectangle for each bin - make them more prominent
        rect(x_pos - bar_width, bin_left, 
             x_pos + bar_width, bin_right, 
             col = adjustcolor("black", alpha = 0.7),  # Increased alpha from 0.4 to 0.7
             border = "black", lwd = 1)  # Increased border width
      }
    }


    # Draw the violin plot outline (density curve) - preserve actual width
    polygon(c(x_pos - scaled_dens, rev(x_pos + scaled_dens)), 
            c(dens$x, rev(dens$x)), 
            col = adjustcolor(i, alpha = 0.3), 
            border = i, lwd = 2)
    
    # Add jittered points within the density boundaries
    n_points <- length(level_estimates)
    for (j in 1:n_points) {
      # Find the density value at this estimate
      dens_idx <- which.min(abs(dens$x - level_estimates[j]))
      max_jitter <- scaled_dens[dens_idx]
      
      # Jitter horizontally within the density boundary
      jitter_x <- runif(1, -max_jitter, max_jitter)
      
      # Add the point
      points(x_pos + jitter_x, level_estimates[j], 
             col = i, pch = 16, cex = 0.6)
    }
    
    # Add individual horizontal line for this beta level (not full width)
    line_width <- 0.8  # Width of the line segment
    segments(x_pos - line_width/2, current_level, 
             x_pos + line_width/2, current_level, 
             col = i, lwd = 3, lty = 2)
  }
}

# Example usage (commented out):
# rm(list = ls())
 library(here)
 pdf(here("output", "simStudy_results", "EZ_clean", "sim_P160T80_ezClean.pdf"))
 plot_cellBetaEstimates(here("output", "simStudy_results", "EZ_clean", "sim_P160T80_ez-Clean.RData"), 
                        show_frequency_bars = TRUE, show_x_axis = TRUE, single_plot = TRUE)
 dev.off()