#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function generates summary statistics for simulated DDM data without generating
# the full trial-by-trial dataset. It uses the sampling distributions implied by the EZ-DDM equations
# to sample the total number of correct responses, mean RT, and RT variance.
#
# Inputs:
# - indiv_pars: A list containing parameter values for each participant:
#   * bound: Decision threshold (a) for each participant
#   * drift: Drift rate (v) for each participant
#   * nondt: Non-decision time (t) for each participant
# - n_trials: Number of trials per participant used for sampling
#
# Returns:
# - A data frame with columns for:
#   * A: Total number of accurate responses (out of n_trials)
#   * Mrt: Mean reaction time
#   * Vrt: Variance of reaction time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_summaryStats <- function(indiv_pars, n_trials){
  # Extract parameters from the input list
  v <- indiv_pars$drift    # Drift rate for each participant
  a <- indiv_pars$bound    # Decision threshold for each participant
  t <- indiv_pars$nondt    # Non-decision time for each participant
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