# This R script contains SEVEN functions to automatically write the JAGS model file
# format_truncation(): A custom function to format truncation limits for JAGS model
# default_truncation_list(): A function to provide default truncation limits for JAGS model
# We start writing the model components:
# model_hierarchical_means(): We write the hierarchical means first
# model_hierarchical_sdevs(): We write the hierarchical standard deviations next
# model_individual_parameters(): We specify the individual parameters
# model_EZequations(): We write the EZ equations last
# We build the .BUGS model file altogether:
# JAGS_writeModel(): Writes the .BUGS model file
#################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Function to handle truncation limits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
format_truncation <- function(limits) {
    if(limits[1] == "" && limits[2] == "") {
      return("")  # No truncation
    } else if(limits[1] == "") {
      return(paste0("T(,", limits[2], ")"))
    } else if(limits[2] == "") {
      return(paste0("T(", limits[1], ",)"))
    } else {
      return(paste0("T(", limits[1], ",", limits[2], ")"))
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Function to provide default truncation limits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
default_truncation_list <- function(){
   return(list("bound_mean" = c(0.1, ""),
               "nondt_mean" = c(0.05, ""),
               "drift_mean" = c(-3, 3),
               "bound" = c(0.1, ""),
               "nondt" = c(0.05, ""),
               "drift" = c(-3, 3),
               "betaweight" = c(-3, 3)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Function to write the hierarchical means
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
model_hierarchical_means <- function(priors, t = NULL){
  # Define prior distributions for hierarchical means
  # Each parameter has a normal prior with specified mean and precision (1/variance)
  # T() indicates truncation to keep parameters in reasonable ranges
  priors.bound_m  <- paste("          bound_mean ~ dnorm(", priors$bound_mean_mean,",pow(",priors$bound_mean_sdev,",-2))", 
                        format_truncation(t$bound_mean), sep="")
  priors.nondt_m  <- paste("          nondt_mean ~ dnorm(", priors$nondt_mean_mean,",pow(",priors$nondt_mean_sdev,",-2))", 
                        format_truncation(t$nondt_mean), sep="")
  priors.drift_m  <- paste("          drift_mean ~ dnorm(", priors$drift_mean_mean,",pow(",priors$drift_mean_sdev,",-2))", 
                        format_truncation(t$drift_mean), sep="")
  return(c(priors.bound_m, priors.nondt_m, priors.drift_m))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Function to write the hierarchical standard deviations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
model_hierarchical_sdevs <- function(priors, t = NULL){
  # Define prior distributions for hierarchical standard deviations
  # Check if we're using uniform or inverse gamma priors
  if(all(c("bound_sdev_lower", "bound_sdev_upper") %in% names(priors))){
    # Uniform priors
    priors.bound_sd <- paste("          bound_sdev ~ dunif(", priors$bound_sdev_lower,",",priors$bound_sdev_upper,")", sep="")
    priors.nondt_sd <- paste("          nondt_sdev ~ dunif(", priors$nondt_sdev_lower,",",priors$nondt_sdev_upper,")", sep="")
    priors.drift_sd <- paste("          drift_sdev ~ dunif(", priors$drift_sdev_lower,",",priors$drift_sdev_upper,")", sep="")
  } else if(all(c("bound_sdev_shape", "bound_sdev_scale") %in% names(priors))){
    # Inverse gamma priors
    priors.bound_sd <- paste("          bound_sdev ~ dgamma(", priors$bound_sdev_shape,",",priors$bound_sdev_scale,")", 
                        format_truncation(t$bound_sdev), sep="")
    priors.nondt_sd <- paste("          nondt_sdev ~ dgamma(", priors$nondt_sdev_shape,",",priors$nondt_sdev_scale,")", 
                        format_truncation(t$nondt_sdev), sep="")
    priors.drift_sd <- paste("          drift_sdev ~ dgamma(", priors$drift_sdev_shape,",",priors$drift_sdev_scale,")", 
                        format_truncation(t$drift_sdev), sep="")
  } else {
    cat("Unknown prior distribution for the hierarchical standard deviations")
  }
  return(c(priors.bound_sd, priors.nondt_sd, priors.drift_sd))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Function to write the individual parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
model_individual_parameters <- function(t = NULL, modelType = "hierarchical", withinSubject = FALSE){
  if(withinSubject==TRUE){
      return(paste0("
  
                  for(p in 1:nParticipants) {
                      bound[p] ~ dnorm(bound_mean, pow(bound_sdev, -2))", 
                        format_truncation(t$bound),
                  "\n                      nondt[p] ~ dnorm(nondt_mean, pow(nondt_sdev, -2))", 
                        format_truncation(t$nondt),
                  "\n                      for(j in 1:2){
                          drift[p,j] ~ dnorm(drift_mean+betaweight*(j-1), pow(drift_sdev, -2))", 
                        format_truncation(t$drift),
                  "\n                      }
                  }"))
  } else {
    if(modelType == "hierarchical"){
      return(paste0(
                "\n                # Sampling model",
                "\n                for (p in 1:nParticipants){",
                "\n                    bound[p] ~ dnorm(bound_mean, pow(bound_sdev, -2))", 
                        format_truncation(t$bound),
                "\n                    nondt[p] ~ dnorm(nondt_mean, pow(nondt_sdev, -2))", 
                        format_truncation(t$nondt),
                "\n                    drift[p] ~ dnorm(drift_mean, pow(drift_sdev, -2))", 
                        format_truncation(t$drift)))
    } else {
         return(paste0(
              "\n                  # Sampling model",
              "\n                  for (p in 1:nParticipants){",
              "\n                      drift[p] ~ dnorm(drift_mean + betaweight*X[p], pow(drift_sdev, -2))", 
                    format_truncation(t$drift),
              "\n                      bound[p] ~ dnorm(bound_mean, pow(bound_sdev, -2))", 
                    format_truncation(t$bound),
              "\n                      nondt[p] ~ dnorm(nondt_mean, pow(nondt_sdev, -2))", 
                    format_truncation(t$nondt)))
    }    
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Function to write the EZ equations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
model_EZequations <- function(withinSubject = FALSE){
  if(withinSubject==TRUE){
    content.end <- "
              # Forward equations from EZ Diffusion
              for (k in 1:length(meanRT)) {
                  ey[k]  = exp(-bound[P[k]] * drift[P[k],(X[k])+1])
                  Pc[k]  = 1 / (1 + ey[k])
                  PRT[k] = 2 * pow(drift[P[k],(X[k]+1)], 3) / bound[P[k]] * pow(ey[k] + 1, 2) / (2 * -bound[P[k]] * drift[P[k],(X[k]+1)] * ey[k] - ey[k]*ey[k] + 1)
                  MDT[k] = (bound[P[k]] / (2 * drift[P[k],(X[k]+1)])) * (1 - ey[k]) / (1 + ey[k])
                  MRT[k] = MDT[k] + nondt[P[k]]

                  # Loss functions using MRT, PRT, and Pc
                  correct[k] ~ dbin(Pc[k], nTrialsPerCondition)
                  meanRT[k]  ~ dnorm(MRT[k], PRT[k] * nTrialsPerCondition)
                  varRT[k]   ~ dnorm(1/PRT[k], 0.5*(nTrialsPerCondition-1) * PRT[k] * PRT[k])
              }
      }"
  } else {
     content.end <- "
                  # Forward equations from EZ Diffusion
                  ey[p]  = exp(-bound[p] * drift[p])
                  Pc[p]  = 1 / (1 + ey[p])
                  PRT[p] = 2 * pow(drift[p], 3) / bound[p] * pow(ey[p] + 1, 2) / (2 * -bound[p] * drift[p] * ey[p] - ey[p]*ey[p] + 1)
                  MDT[p] = (bound[p] / (2 * drift[p])) * (1 - ey[p]) / (1 + ey[p])
                  MRT[p] = MDT[p] + nondt[p]

                  # Loss functions using MRT, PRT, and Pc
                  correct[p] ~ dbin(Pc[p], nTrialsPerPerson)
                  meanRT[p]  ~ dnorm(MRT[p], PRT[p] * nTrialsPerPerson)
                  varRT[p]   ~ dnorm(1/PRT[p], 0.5*(nTrialsPerPerson-1) * PRT[p] * PRT[p])
              }
      }"
  }
}

######################################################################################
# Main function to write the .BUGS model file
######################################################################################
JAGS_writeModel <- function(priors, modelType, withinSubject = FALSE, modelFile=NA, custom_truncation_list = NULL){

  # If no custom truncation list is provided, use these default values
  if(is.null(custom_truncation_list)){
        custom_truncation_list <- default_truncation_list()
  }  

  # If no model file is provided, use a generic model file name
  if(is.na(modelFile)){
    modelFile <- here("output", "BUGS-models", "JAGS_model.txt")
  }

  opening <- "model{"
  hierarchical_means <- model_hierarchical_means(priors, t = custom_truncation_list)
  hierarchical_sdevs <- model_hierarchical_sdevs(priors, t = custom_truncation_list)
  individual_parameters <- model_individual_parameters(t = custom_truncation_list, modelType, withinSubject)
  EZequations <- model_EZequations(withinSubject)

  start <- c(opening, hierarchical_means, hierarchical_sdevs)
  if((modelType != "hierarchical")||(withinSubject==TRUE)){
    priors.beta     <- paste("          betaweight ~ dnorm(", priors$betaweight_mean,",pow(",priors$betaweight_sdev,",-2))", 
                            format_truncation(custom_truncation_list$betaweight), sep="")
    start <- c(start, priors.beta)
  }
    
  # Write the complete model to file
  final_file <- file(modelFile)
  writeLines(c(start, individual_parameters, EZequations), final_file)
  close(final_file)
}
