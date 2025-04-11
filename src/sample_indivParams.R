#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function generates individual participant parameters for DDM simulations
# by sampling from truncated normal distributions.
#
# Inputs:
# - n_participants: Number of participants to simulate
# - dmean: Mean value for the drift rate (v) distribution
# - bmean: Mean value for the decision boundary/threshold (a) distribution
# - nmean: Mean value for the non-decision time (t) distribution
#
# Returns:
# - A data frame with columns for:
#   * bound: Decision threshold (a) for each participant
#   * drift: Drift rate (v) for each participant
#   * nondt: Non-decision time (t) for each participant
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

# Check if truncnorm package is installed, install if needed
if (!requireNamespace("truncnorm", quietly = TRUE)) {
  install.packages("truncnorm")
}
# Load the truncnorm package
library(truncnorm)

sample_indivParams <- function(n_participants, dmean, bmean, nmean){
  # Sample drift rates from truncated normal distribution
  # No lower bound specified, allowing for negative drift rates
  # SD = 1 allows for reasonable individual differences
  d <- rtruncnorm(n = n_participants, mean = dmean, sd = 0.75)
  
  # Sample decision boundaries from truncated normal distribution
  # Lower bound of 0 ensures boundaries are positive
  # SD = 0.5 provides moderate individual differences
  b <- rtruncnorm(n = n_participants, a = 0, mean = bmean, sd = 0.75)
  
  # Sample non-decision times from truncated normal distribution
  # Lower bound of 0 ensures non-decision times are positive
  # SD = 0.15 provides smaller individual differences for this parameter
  n <- rtruncnorm(n = n_participants, a = 0, mean = nmean, sd = 0.05)
  
  # Return a data frame with the sampled parameters for each participant
  return(data.frame("bound" = b, "drift" = d, "nondt" = n))
}