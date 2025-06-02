# This R script contains four functions:
# sample_hierarchical_parameters(): Define  hierarchical parameters used to sample individual parameters
# sample_drift(): Samples individual drift rates according to the simulation design
# get_simulation_parameters(): Obtain true parameter values used to generate DDM data
#################################################################################
# Check if truncnorm package is installed, install if needed
if (!requireNamespace("truncnorm", quietly = TRUE)) {
  install.packages("truncnorm")
}
# Load the truncnorm package
library(truncnorm)

################################################################################
# Function 1: Sample hierarchical parameters
################################################################################
# This function samples hierarchical parameter values used to later sample
# individual parameters
# Inputs:
# - true_means: A list containing uniform ranges or point values for the true hierarchical means
# - true_sdevs: A list containing uniform ranges or point values for the true hierarchical standard deviations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_hierarchical_parameters <- function(true_means, true_sdevs) {
      # Function requires a list specifying the uniform distributions for each parameter
      if(is.null(true_means)){
            cat("No ranges for the true hierarchical means found. Please fill in the true_means list.")
            return(list(
                  "bound_mean" = c(NA, NA), "nondt_mean" = c(NA, NA), "drift_mean" = c(NA, NA)))
      } else {
        ##########################################################
        # Sample HIERARCHICAL MEANS
        ##########################################################
        hierarchical_parameters <- list(
          "bound_mean" = ifelse(length(true_means$bound_mean) == 2,
                                runif(1, true_means$bound_mean[1], true_means$bound_mean[2]),
                                true_means$bound_mean),
          "drift_mean" = ifelse(length(true_means$drift_mean) == 2,
                                runif(1, true_means$drift_mean[1], true_means$drift_mean[2]),
                                true_means$drift_mean),
          "nondt_mean" = ifelse(length(true_means$nondt_mean) == 2,
                            runif(1, true_means$nondt_mean[1], true_means$nondt_mean[2]),
                            true_means$nondt_mean),
        ##########################################################
        # Sample HIERARCHICAL STANDARD DEVIATIONS
        ##########################################################
          "bound_sdev" = ifelse(is.null(true_sdevs$bound_sdev),
                                bound_mean / 5,
                                ifelse(length(true_sdevs$bound_sdev) == 2,
                                      runif(1, true_sdevs$bound_sdev[1], true_sdevs$bound_sdev[2]),
                                      true_sdevs$bound_sdev)),
          "nondt_sdev" = ifelse(is.null(true_sdevs$nondt_sdev),
                                nondt_mean / 5,
                                ifelse(length(true_sdevs$nondt_sdev) == 2,
                                    runif(1, true_sdevs$nondt_sdev[1], true_sdevs$nondt_sdev[2]),
                                    true_sdevs$nondt_sdev)),
          "drift_sdev" = ifelse(is.null(true_sdevs$drift_sdev),
                                drift_mean / 5,
                                ifelse(length(true_sdevs$drift_sdev) == 2,
                                    runif(1, true_sdevs$drift_sdev[1], true_sdevs$drift_sdev[2]),
                                    true_sdevs$drift_sdev)))
        return(hierarchical_parameters)
      }
}

# Test the function
#true_means <- list("bound_mean" = c(0, 1), "drift_mean" = c(0, 1), "nondt_mean" = c(0, 1))
#true_sdevs <- list("bound_sdev" = c(0, 1), "drift_sdev" = c(0, 1), "nondt_sdev" = c(0, 1))
#hierarchical_parameters <- sample_hierarchical_parameters(true_means = true_means, true_sdevs = true_sdevs)
#print(hierarchical_parameters)

################################################################################
# Function 2: Sample individual drift rates
################################################################################
# This function samples individual drift rates according to the simulation design
# Inputs:
# - nPart: Number of participants to simulate
# - drift_mean: Mean value for the drift rate distribution
# - drift_sdev: Standard deviation for the drift rate distribution
# - betaweight: Optional fixed value for the regression coefficient
# - X: Vector of predictors for regression models
# - withinSubject: Whether to use a within-subject design (default: FALSE)
# - modelType: Type of model ("hierarchical" or "ttest")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_drift <- function(nPart, drift_mean, drift_sdev, betaweight, X, 
                         withinSubject = FALSE, modelType = NULL) {
      if(withinSubject) {   drift <- rnorm(nPart*2, drift_mean + (betaweight*X), drift_sdev)
      }else{  
              if(modelType == "ttest") {
                drift <- rnorm(nPart, drift_mean + (betaweight*X), drift_sdev)
              } else {
                drift <- rnorm(nPart, drift_mean, drift_sdev)
              }
      }
return(drift)
}


################################################################################
# Function 3: Main function - Obtain true parameter values for the simulation study
################################################################################
# Inputs:
# - true_means: A list containing the uniform ranges or point values for the true hierarchical means
# - true_sdevs: A list containing the uniform ranges or point values for the true hierarchical standard deviations
# - Show: Whether to display the sampled parameters in the console
# - nPart: Number of participants to simulate
# - modelType: Type of model ("hierarchical" or "ttest")
# - X: Vector of predictors for regression models
# - fixedBeta: Optional fixed value for the regression coefficient
# - withinSubject: Whether to use a within-subject design (default: FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_simulation_parameters <- function(true_means = NULL, true_sdevs = NULL, Show = TRUE,
                              nPart, modelType = NULL, X = NULL, fixedBeta = NA, withinSubject = FALSE) {
      # Defensive programming: Default modelType depends on whether design is within-subject
      if(is.null(modelType)){
            if(withinSubject){
                  modelType <- "ttest"
            } else {
                  modelType <- "hierarchical"
            }
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      # Sample hierarchical parameters (means and standard deviations)
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      hierarchical_parameters <- sample_hierarchical_parameters(true_means = true_means,
                                                                true_sdevs = true_sdevs)
            drift_mean <- hierarchical_parameters$drift_mean
            nondt_mean <- hierarchical_parameters$nondt_mean
            bound_mean <- hierarchical_parameters$bound_mean
            drift_sdev <- hierarchical_parameters$drift_sdev
            nondt_sdev <- hierarchical_parameters$nondt_sdev
            bound_sdev <- hierarchical_parameters$bound_sdev
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      # Sample non-hierarchical parameters
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      betaweight <- ifelse(is.na(fixedBeta) && modelType == "ttest", runif(1, -1, 1), fixedBeta)
      bound <- rtruncnorm(n = nPart, a = 0, mean = bound_mean, sd = bound_sdev)
      nondt <- rtruncnorm(n = nPart, a = 0, mean = nondt_mean, sd = nondt_sdev)       
      drift <- sample_drift(nPart = nPart, drift_mean = drift_mean, drift_sdev = drift_sdev, 
                            betaweight = betaweight, X = X, withinSubject = withinSubject, modelType = modelType)  
      # Create parameter set list with all generated values
      parameter_set <- list(
        "bound_mean" = bound_mean, "drift_mean" = drift_mean, "nondt_mean" = nondt_mean, 
        "bound_sdev" = bound_sdev, "drift_sdev" = drift_sdev, "nondt_sdev" = nondt_sdev,
        "bound" = bound, "drift" = drift, "nondt" = nondt)
        if(modelType == "ttest"||!is.na(betaweight)) {
          parameter_set <- c(parameter_set, list("betaweight" = betaweight))
        }

      if(Show) {      show_parameters(parameter_set)          }
return(parameter_set)
}