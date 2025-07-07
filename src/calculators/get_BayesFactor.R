# Script to compute Bayes Factors against Beta = 0
###################################################

# Function to compute Bayes Factors using ROPE approach
# (Same as in our paper)
compute_BF_ROPE <- function(true_betas, chains, epsilon = 0.1) {
  
  # Get unique beta levels
  beta_levels <- sort(unique(true_betas))
  
  # Define region of practical equivalence (ROPE) for Bayes factor calculation
  prior_constant <- pnorm(epsilon) - pnorm(-epsilon)  # Prior probability mass in ROPE
  
  # Initialize results as data frame
  results <- data.frame(
    true_beta = numeric(length(chains)),
    post_mass_rope = numeric(length(chains)),
    log_bayes_factor = numeric(length(chains)),
    bayes_factor = numeric(length(chains))
  )
  
  # Compute Bayes Factor for each chain
  for (i in 1:length(chains)) {
    current_chain <- chains[[i]]
    current_true_beta <- true_betas[i]
    
    # Calculate posterior probability mass in ROPE
    post_mass <- mean(current_chain > -epsilon & current_chain < epsilon)
    
    # Calculate log Bayes factor (log of prior/posterior odds)
    if (post_mass == 0) {
      log_bf <- Inf  # Handle cases where posterior mass is zero
    } else {
      log_bf <- log(prior_constant / post_mass)
    }
    
    # Store results in data frame
    results$true_beta[i] <- current_true_beta
    results$post_mass_rope[i] <- post_mass
    results$log_bayes_factor[i] <- log_bf
    results$bayes_factor[i] <- exp(log_bf)
  }
  
  return(list(
    bayes_factors = results,
    beta_levels = beta_levels,
    true_betas = true_betas,
    epsilon = epsilon,
    prior_constant = prior_constant
  ))
}
