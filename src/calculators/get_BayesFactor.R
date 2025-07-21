# Script to compute Bayes Factors against Beta = 0
###################################################
# Function to compute Bayes Factors using ROPE approach
# (Same as in our paper)
compute_BF_ROPE <- function(true_betas, chains, epsilon = 0.05) {
  
  # First, try to clean the chains if needed
  cleaned_result <- diagnose_and_clean_chains(chains, true_betas)
  
  if (is.list(cleaned_result) && "chains" %in% names(cleaned_result)) {
    # Cleaned result contains both chains and true_betas
    chains <- cleaned_result$chains
    true_betas <- cleaned_result$true_betas
  } else {
    # Cleaned result is just the chains
    chains <- cleaned_result
  }
  
  # Get unique beta levels
  beta_levels <- sort(unique(true_betas))
  
  # Define region of practical equivalence (ROPE) for Bayes factor calculation
  prior_constant <- pnorm(epsilon) - pnorm(-epsilon)  # Prior probability mass in ROPE
  
  # Initialize results as data frame
  results <- data.frame(
    true_beta = numeric(length(chains)),
    mean_posterior = numeric(length(chains)),
    post_mass_rope = numeric(length(chains)),
    log_bayes_factor = numeric(length(chains)),
    bayes_factor = numeric(length(chains))
  )
  
  # Compute Bayes Factor for each chain
  for (i in 1:length(chains)) {
    current_chain <- unlist(chains[[i]])

    if(length(current_chain) == 1){
        current_chain <- chains[i,]
    }

    current_chain <- as.vector(current_chain)
    current_chain <- current_chain[!is.na(current_chain)]

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
    results$mean_posterior[i] <- mean(current_chain)
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
