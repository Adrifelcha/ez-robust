# Auxiliary function to generate a sequence of 6 colors from dark to light,
# centered around the input solid_color.
makeColors <- function(solid_color){    
    nColors <- 6
    
    col_dark <- colorRampPalette(c("black", solid_color))(3)[2]
    col_light <- colorRampPalette(c(solid_color, "white"))(3)[2]
    main_palette <- colorRampPalette(c(col_dark, solid_color, col_light))

    # --- Generate the nColors ---
    color_sequence <- main_palette(nColors)

    colors <- matrix(color_sequence, ncol = 2, byrow = TRUE)

    # Return the sequence of 6 hex color codes
    return(colors)
}

# --- Example Usage ---
# my_color_ramp <- makeColors("darkgreen")
# print(my_color_ramp)
# # Expected output: A vector of 6 hex codes, starting dark green,
# # passing through darkgreen around the 3rd/4th position, ending light green.
# # Example:
# # [1] "#003200" "#004B00" "#006400" "#559F55" "#AADFAA" "#FFFFFF"
# # (Exact hex codes depend on interpolation details, but this is the idea)
#
# # You can visualize the ramp:
# # scales::show_col(my_color_ramp) # Requires the 'scales' package


# Internal function to plot the posterior distribution of a parameter
plot_posterior <- function(density_data, parameter_name, color = c("blue", "darkblue"), max_density=NA, true_value=NA){
                # Calculate density once
                dens <- density(density_data)

                if(is.na(max_density)){    max_density <- max(dens$y)        }

                # Calculate plot limits, ensuring true_value is included if present
                minX <- min(dens$x)
                maxX <- max(dens$x)
                if(!is.na(true_value)){    minX <- min(minX, true_value); maxX <- max(maxX, true_value)    }
                # Add a small buffer to limits
                x_buffer <- (maxX - minX) * 0.06
                xlim <- c(minX - x_buffer,maxX + x_buffer)

                # Determine main title based on parameter name
                if(parameter_name == "bound_mean"){
                    main_title <- "Mean boundary separation"
                    parameter <- expression(paste(alpha))
                } else if(parameter_name == "nondt_mean"){
                    main_title <- "Mean non-decision time"
                    parameter <- expression(paste(tau[0]))
                } else if(parameter_name == "drift_mean"){
                    main_title <- "Mean drift rate"
                    parameter <- expression(paste(delta))
                } else {
                    main_title <- parameter_name # Default title if name not matched
                }

                # Set up the plot area WITHOUT drawing the density line yet
                plot(dens, main="", ylim=c(0, max_density), xlim=xlim,
                     ann=FALSE, axes=FALSE, type='n') # type='n' prevents plotting

                # Draw the polygon with TRANSPARENT fill and SOLID border                
                polygon(dens, col=color[2], border = NA)
                lines(dens, col=color[1], lwd=4)
                #polygon(dens, col="green", border="purple")

                # Add axes
                x_label <- round(seq(minX,maxX, length.out=7),1)
                axis(1, x_label, x_label, line=-0.9) # X-axis                

                # Add title and axis labels (optional, customize as needed)
                mtext(main_title, side=3, cex=0.8, line=0, f=2)
                mtext(parameter, side=1, cex=1.4, line=3)

                # Add true value line if provided
                if(!is.na(true_value)){
                    abline(v = true_value, col = "black", lty = 2, lwd = 1.5) # Make line slightly thicker
                }
}


# Main function to plot the posterior distributions of the three parameters
# color_scheme: a vector of two colors, the first is the solid color, the second is the transparent color
# true_value: a list of true values for the parameters
plot_posteriorDistributions <- function(samples, color_scheme=c("orange"), true_value=NA){
                bound_mean1 <- samples$BUGSoutput$sims.list$bound_mean
                nondt_mean1 <- samples$BUGSoutput$sims.list$nondt_mean
                drift_mean1 <- samples$BUGSoutput$sims.list$drift_mean
                
                colors <- makeColors(color_scheme)
                # Find the maximum density value across all three parameters
                max_density <- max(
                max(density(bound_mean1)$y),
                max(density(nondt_mean1)$y),
                max(density(drift_mean1)$y)
                )

                par(mfrow=c(1,3), mar=c(4.5,0,3,0))                
                plot_posterior(density_data = bound_mean1, parameter_name = "bound_mean", 
                               color = colors[1,], max_density = max_density, true_value = true_value$bound_mean)
                plot_posterior(density_data = nondt_mean1, parameter_name = "nondt_mean", 
                               color = colors[2,], max_density = max_density, true_value = true_value$nondt_mean)
                plot_posterior(density_data = drift_mean1, parameter_name = "drift_mean",
                               color = colors[3,], max_density = max_density, true_value = true_value$drift_mean)

                par(mfrow=c(1,1))
}


#plot_posteriorDistributions(samples, color_scheme = c("darkorange"), true_value = trueVals)
