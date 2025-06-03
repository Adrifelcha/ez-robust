# This R script contains four functions:
# identify_design(): Determines whether we are running a within- or between-subject design
# sample_data(): Ensures trial data is generated according to the prevent_zero_accuracy flag
# get_simulation_data(): Generates a complete dataset according to the simulation study cell design
# sample_summaryStats(): Uses the EZ-DDM equations to generate EZ summary statistics
#################################################################################


################################################################################
# Function 1: Identify design
################################################################################
# This function identifies the design (between- vs within-subject) 
# Inputs:
# - nPart: Number of participants to simulate
# - nTrials: Total number of trials per participant (used when no conditions)
# - nTrialsPerCondition: Number of trials per condition (used when conditions exist)
# - parameter_set: A list containing parameter values for each participant
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
identify_design <- function(nPart, nTrials, nTrialsPerCondition, parameter_set){
      withinSubject <- !is.na(nTrialsPerCondition)
      if(withinSubject){
          nObs <- nPart*nTrialsPerCondition*2
          nCols <- 4
          colNames <- c("sub","cond","rt", "accuracy")
          n_subIndex <- nTrialsPerCondition*2
          col_accuracy <- 4
          col_rt <- 3
          N <- nTrialsPerCondition
          adjusted_parameter_set <- list(
                  bound = rep(parameter_set$bound, each=2),
                  drift = parameter_set$drift,
                  nondt = rep(parameter_set$nondt, each=2))
      } else {
          nObs <- nPart*nTrials
          nCols <- 3
          colNames <- c("sub","rt", "accuracy")
          n_subIndex <- nTrials
          col_accuracy <- 3
          col_rt <- 2
          N <- nTrials
          adjusted_parameter_set <- parameter_set
      }

      data <- matrix(NA, ncol=nCols, nrow=nObs)  
      data[,1] <- rep(1:nPart, each=n_subIndex)
      cell_index <- data[,1]

      if(withinSubject){
          data[,2] <- rep(rep(c(1,0), each=nTrialsPerCondition), nPart)
          cell_index <- rep(1:(nPart*2), each=nTrialsPerCondition)
      }

  return(list(data = data, nObs = nObs,
              nCols = nCols, N = N,
              n_subIndex = n_subIndex,
              col_accuracy = col_accuracy,
              col_rt = col_rt,
              colNames = colNames,
              adjusted_parameter_set = adjusted_parameter_set,
              cell_index = cell_index))
}

################################################################################
# Function 2: Sample data
################################################################################
# This function generates a complete dataset for simulation studies using the DDM model
# Inputs:
# - nPart: Number of participants to simulate
# - nTrials: Total number of trials per participant (used when no conditions)
# - nTrialsPerCondition: Number of trials per condition (used when conditions exist)
# - parameter_set: A list containing parameter values for each participant
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_data <- function(nPart, N, params, contamination_probability = 0, prevent_zero_accuracy = FALSE){
      # Generate dataset for this participant using their specific parameters
      # First generate the dataset once
      temp <- get_DDM_data(a = params$bound, v = params$drift, t = params$nondt, n = N,
                           contamination_probability = contamination_probability)
      accuracy <- temp$accuracy

      # If prevent_zero_accuracy is TRUE and we got all zeros, keep trying
      while(sum(accuracy)==0 && prevent_zero_accuracy){
        temp <- get_DDM_data(a = params$bound, v = params$drift, t = params$nondt, n = N,
                            contamination_probability = contamination_probability)
        accuracy <- temp$accuracy
      }

return(temp)
}


################################################################################
# Function 3: Main function to generate a dataset
################################################################################
# This function combines the above functions to generate the appropriate dataset
# Inputs:
# - nPart: Number of participants to simulate
# - nTrials: Total number of trials per participant (used when no conditions)
# - nTrialsPerCondition: Number of trials per condition (used when conditions exist)
# - parameter_set: A list containing parameter values for each participant
# - contamination_probability: probability of contamination
# - prevent_zero_accuracy: A flag (TRUE/FALSE) indicating if zero accuracy should be prevented
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_simulation_data <- function(nPart, nTrials = NA, parameter_set, 
                                nTrialsPerCondition = NA, contamination_probability = 0, prevent_zero_accuracy = FALSE){
  
  design <- identify_design(nPart, nTrials, nTrialsPerCondition, parameter_set)
            data <- design$data
            nObs <- design$nObs
            nCols <- design$nCols
            n_subIndex <- design$n_subIndex
            col_accuracy <- design$col_accuracy
            col_rt <- design$col_rt
            N <- design$N
            colNames <- design$colNames
            adjusted_parameter_set <- design$adjusted_parameter_set
            cell_index <- design$cell_index
            n_par_sets <- length(adjusted_parameter_set$bound)

  for(i in 1:n_par_sets){
        # Identify rows for current participant
        this.cell <- which(cell_index==i)

        params <- list(bound = adjusted_parameter_set$bound[i],
                       drift = adjusted_parameter_set$drift[i],
                       nondt = adjusted_parameter_set$nondt[i])
        # Generate dataset for this participant using their specific parameters
        # First generate the dataset once
        temp <- sample_data(nPart, N, params, contamination_probability, prevent_zero_accuracy)
        data[this.cell,col_accuracy] <- temp$accuracy
        data[this.cell,col_rt] <- temp$RT
  }
  # Convert to matrix and add column names
  data <- as.matrix(data)
  colnames(data) <- colNames
  
return(data)
}


################################################################################
# Function 4: Generate summary statistics
################################################################################
# This function generates summary statistics using the EZ-DDM equations
# Inputs:
# - parameter_set: A list containing parameter values for each participant
# - n_trials: Number of trials per participant used for sampling
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_summaryStats <- function(parameter_set, n_trials){
  # Extract parameters from the input list
  v <- parameter_set$drift    # Drift rate for each participant
  a <- parameter_set$bound    # Decision threshold for each participant
  t <- parameter_set$nondt    # Non-decision time for each participant
  N <- n_trials            # Number of trials per participant
  nP <- length(v)          # Number of participants
  
  # Calculate intermediate value used in EZ-diffusion equations
  y <- exp(-a * v)
  
  # Forward EZ equations to predict theoretical values
  PredAccuracyRate <- 1 / (1 + y)  # Equation 1: Predicted accuracy rate
  PredMean <- t + ((a / (2 * v)) * ((1 - y) / (1 + y)))  # Equation 2: Predicted mean RT
  PredVariance <- (a / (2 * v^3)) * (( 1 - 2 * a * v * y - exp(-a*2*v)) / ((y + 1)^2))  # Equation 3: Predicted RT variance
  
  # Sample observed summary statistics from theoretical predictions
  # Sample total number of accurate responses using binomial distribution
  ObservedAccuracyTotal <- rbinom(nP, size = N, prob = PredAccuracyRate)
  
  # Sample observed mean RT using normal distribution
  # Standard error of the mean = sqrt(variance/N)
  ObservedMean <- rnorm(nP, mean = PredMean, sd = sqrt(PredVariance / N))
  
  # Sample observed RT variance using normal distribution
  # Standard error of variance = sqrt(2*variance^2/(N-1))
  ObservedVariance <- rnorm(nP, mean = PredVariance, sd = sqrt((2 * (PredVariance^2)) / (N - 1)))
  
  # Return a data frame with the sampled summary statistics
  return(data.frame("A" = ObservedAccuracyTotal, 
                    "Mrt" = ObservedMean, 
                    "Vrt" = ObservedVariance))
}