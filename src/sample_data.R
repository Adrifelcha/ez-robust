#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function generates a complete dataset for simulation studies using the DDM model
# Ex: It creates data for multiple participants with individual parameter values.
#
# Inputs:
# - nPart: Number of participants to simulate
# - nTrials: Total number of trials per participant (used when no conditions)
# - parameter_set: A list containing parameter values for each participant:
#   * bound: Decision threshold (a) for each participant
#   * drift: Drift rate (v) for each participant/condition
#   * nondt: Non-decision time (t) for each participant
# - nTrialsPerCondition: Number of trials per condition (used when conditions exist)
#
# Returns:
# - A matrix with columns for participant ID, reaction time, accuracy
#   (and condition, if applicable)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Sample data using simulation settings and true parameter values sampled
sample_data <- function(nPart, nTrials = NA, parameter_set, 
                        nTrialsPerCondition = NA, prevent_zero_accuracy = TRUE,
                        outlier_probability = 0){
        # Calculate total number of observations and initialize data matrix
        nObs <- nPart*nTrials
        data <- matrix(NA, ncol=3, nrow=nObs)
        
        # Assign participant IDs to each row
        data[,1] <- rep(1:nPart, each=nTrials)
        
        # Generate data for each participant
        for(i in 1:nPart){
          # Identify rows for current participant
          this.sub <- which(data[,1]==i)
          
          # Generate dataset for this participant using their specific parameters
          # First generate the dataset once
          temp <- sample_dataset(a = parameter_set$bound[i], 
                                 v = parameter_set$drift[i], 
                                 t = parameter_set$nondt[i], 
                                 n = nTrials,
                                 outlier_probability = outlier_probability)
          accuracy <- temp$accuracy
          
          # If prevent_zero_accuracy is TRUE and we got all zeros, keep trying
          while(sum(accuracy)==0 && prevent_zero_accuracy){
            temp <- sample_dataset(a = parameter_set$bound[i], 
                                   v = parameter_set$drift[i], 
                                   t = parameter_set$nondt[i], 
                                   n = nTrials,
                                   outlier_probability = outlier_probability)
            accuracy <- temp$accuracy
          }
          
          # Store reaction times and accuracy in the data matrix
          data[this.sub,3] <- accuracy
          data[this.sub,2] <- temp$RT
        }
        
        # Convert to matrix and add column names
        data <- as.matrix(data)
        colnames(data) <- c("sub", "rt", "accuracy")
  
  return(data)
}