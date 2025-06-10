# This R script contains four functions:
# run_random_walk(): Emulates a single-trial Wiener diffusion random walk
# sample_trials(): Simulates multiple trials (with contaminant data)
# add_contaminant(): Generates contaminant data (RT outliers or accuracy guesses)
# get_DDM_data(): Combines the above functions to generate a dataset
#################################################################################


################################################################################
# Function 1: Simulate single Random Walk trial
################################################################################
# This function generates a SINGLE TRIAL observation from the DDM model
# Inputs:
# - a, v: boundary and drift rate
# - dt: time step size
# - max_steps: maximum number of steps
# This function emulates the random walk process implied by the DDM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
run_random_walk <- function(a, v, dt, max_steps){
  # Initialize the evidence accumulator at 0 (starting point)
  x <- 0
  
  # Generate  random  step size noise in advance 
  # (Stochastic component of evidence accumulation - Noise around the drift rate)
  random_dev <- rnorm(max_steps)  
  
  # Scale the random noise and drift rate
  # (This ensures the discrete process approximates a continuous Wiener process)
  noise <- random_dev * sqrt(dt)
  drift <- v * dt
  
  # Start the random walk!
  for(i in 2:max_steps){
    # Calculate the evidence 'sampled' in this step
    this_step = drift + noise[i]
    
    # Update the evidence accumulated so far
    x = x + this_step
    
    # Check if a decision boundary has been reached
    # The boundaries are at +a/2 and -a/2
    if(abs(x)>=(a/2)){  
      break  # Stop the simulation if a boundary is reached
    }
  }
  
  # Create the output:
  # - RT: Reaction time (scaled step count)
  #   Note: We subtract 1 to avoid counting the starting point as a step
  # - C: Final evidence accumulated (x)
  #   Note: This value gets discretized in function sample_dataset()
  output <- list("RT" = (i-1)*dt, "C"  = x)
  
  return(output)
}

# Test the function
# a <- 1; v <- 1; dt <- 0.001; max_steps <- 10000
# run_random_walk(a, v, dt, max_steps)


################################################################################
# Function 2: Sample multiple trials from the DDM model
################################################################################
# This function runs *n* random walk trials by calling run_random_walk() n times
# Inputs:
# - a, v, t: boundary, drift rate, and non-decision time
# - n: number of trials
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
sample_trials <- function(a,v,t,n){
  # Set random walk emulation parameters
  dt = 0.001  # Time step size (in seconds)
  max_steps = 10 / dt  # Maximum number of steps (~10 seconds)
  
  # Initialize vectors to store results
  rt = rep(NA,n)        # Reaction times
  accuracy = rep(NA,n)  # Accuracy (1 = correct, 0 = incorrect)
  
  # Generate n trials
  for(i in 1:n){
      # Generate a single DDM trial
      X <- run_random_walk(a, v, dt, max_steps)
      
      # Store the reaction time (without non-decision time for now)
      rt[i] <- X$RT 
      
      # Determine accuracy based on the final evidence value
      if(X$C > 0){  
          # Hitting the upper boundary is considered a "correct" response
          accuracy[i] <- 1
      } else {      
          # Hitting the lower boundary is considered an "incorrect" response
          accuracy[i] <- 0  
      }
  }
  
  # Create output data frame with reaction times and accuracy
  output <- data.frame("RT" = rt + t, "accuracy" = accuracy)
 
return(output)
}  

# Test the function
#a <- 1; v <- 1; t <- 0; n <- 1000
#sample_trials(a, v, t, n)


################################################################################
# Function 3: Add contaminant data to the dataset
################################################################################
# This function adds contaminant data to the dataset
# Inputs:
# - a, t: boundary and non-decision time
# - ddm_data: data generated from the DDM model
# - contamination_probability: probability of contamination
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
add_contaminant <- function(a,t, ddm_data, contamination_probability){
      # Step 1: Decide which trials will be contaminated
      n <- length(ddm_data$RT)
      # We fix the number of contaminant trials based on the contamination probability
      contaminant_count <- n*contamination_probability
      # The first few trials will be contaminated
      replace_trial <- 1:contaminant_count      
      # We determine what type of contaminant each trial is
      contaminant_type <- rbinom(contaminant_count, 1, 0.5) # 0 = participant guessed, 1= RT outlier
      
      outlier_contaminants <- sum(contaminant_type)
      guess_contaminants <- contaminant_count - outlier_contaminants

      # Step 3: Replace RTs with outliers
      if(outlier_contaminants>0){
        RT_replace <- ddm_data$RT[replace_trial][as.logical(contaminant_type)]         
        ddm_data$RT[replace_trial][as.logical(contaminant_type)]  <- RT_replace + runif(outlier_contaminants, 2, 3)
      }

      # Step 4: Replace accuracy with guess data
      if(guess_contaminants>0){
        tmp_x <- sample_trials(a,v=0,t,guess_contaminants)
        ddm_data$accuracy[replace_trial][!as.logical(contaminant_type)]  <- tmp_x$accuracy
        ddm_data$RT[replace_trial][!as.logical(contaminant_type)]  <- tmp_x$RT + t
      }

      return(ddm_data)
}


################################################################################
# Function 4: Main function to generate a dataset
################################################################################
# This function combines the above functions to generate the appropriate dataset
# Inputs:
# - a, v, t: boundary, drift rate, and non-decision time
# - n: number of trials
# - contamination_probability: probability of contamination
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_DDM_data <- function(a, v, t, n, contamination_probability = 0, separate_datasets = FALSE){
  clean_data <- sample_trials(a, v, t, n)
  if(contamination_probability > 0){
      output <- add_contaminant(a, t, ddm_data = clean_data, contamination_probability)
  }else{
      output <- clean_data
  }

  if(separate_datasets){
     return(list("clean_data" = clean_data, "contaminated_data" = output))
  } else {
     return(output)
  }
}