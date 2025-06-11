###################################################################################
# MODULAR FUNCTION CODING AHEAD:
#
# HDDM_runJAGS is a higher-level function that manages the MCMC sampling process
# using JAGS to estimate hierarchical drift diffusion model parameters.
###################################################################################
# This function handles the complete JAGS workflow for HDDM parameter estimation:
# 1. Prepares the data for JAGS
# 2. Runs the MCMC sampling process
# 3. Extracts and processes the posterior samples
# 4. Calculates summary statistics and diagnostic metrics
#
# Inputs:
# - summaryData: Data frame containing EZ-DDM summary statistics per cell design
# - nTrials: Number of trials per participant
# - X: Predictor vector for models with a regression structure
# - jagsData: List of data to be passed to JAGS
# - jagsParameters: Vector of parameter names to be monitored
# - jagsInits: List of initial values for MCMC chains
# - n.chains: Number of MCMC chains to run (default: 4)
# - n.burnin: Number of burn-in iterations (default: 250)
# - n.iter: Total number of iterations including burn-in (default: 2000)
# - n.thin: Thinning interval for samples (default: 1)
# - modelFile: Path to the JAGS model file (default: "./EZHBDDM.bug")
# - Show: Whether to display diagnostic plots (default: TRUE)
# - track_allParameters: Whether to track all or only hierarchical parameters
#
# Returns a list containing:
#   * estimates: Posterior means for all monitored parameters
#   * sd: Posterior standard deviations
#   * credInterval: 95% credible intervals
#   * rhats: Convergence diagnostics
#   * clock: Computation time in seconds
#   * n.iter: Number of iterations used
###################################################################################

###################################################################################
simStudy_runJAGS <- function(summaryData, nTrials, X, jagsData, jagsParameters, jagsInits, 
                         n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile=NA, Show = TRUE,
                         track_allParameters = track_allParameters, separate_datasets = FALSE,
                         include_EZ_Robust = FALSE){  
      # If modelFile is not provided, use the default model file
      if(is.na(modelFile)){
        modelFile <- here("output", "BUGS-models", "EZHBDDM.bug")
      }

      # Step 1: Extract key summary statistics from the data
      # These are the sufficient statistics for the EZ-DDM
      if(separate_datasets){
          EZ_stats_contaminated <- getEZ_stats(summaryData$contaminated_summary, nTrials)
          EZ_stats_clean <- getEZ_stats(summaryData$clean_summary, nTrials)
          
      }else{
          EZ_stats <- getEZ_stats(summaryData, nTrials)
      }

      model  <- "EZ"
      if(include_EZ_Robust){
         model <- c(model, "EZRobust")
      }

      results <- list()
      for(m in model){
        if(separate_datasets){
          results[[m]] <- list("contaminated" = HDDM_runJAGS(EZ_stats_contaminated, nTrials, X, jagsData[[m]], jagsParameters, jagsInits,
                                                        n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile=modelFile[m], Show = TRUE,
                                                        track_allParameters = track_allParameters),
                               "clean" = HDDM_runJAGS(EZ_stats_clean, nTrials, X, jagsData[[m]], jagsParameters, jagsInits,
                                                      n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile=modelFile[m], Show = TRUE,
                                                      track_allParameters = track_allParameters))
        }else{
          results[[m]] <- HDDM_runJAGS(EZ_stats, nTrials, X, jagsData, jagsParameters, jagsInits,
                                       n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile=modelFile[m], Show = TRUE,
                                       track_allParameters = track_allParameters)
        }
      }
return(results)
}




HDDM_runJAGS <- function(EZ_stats, nTrials, X, jagsData, jagsParameters, jagsInits, 
                         n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile=NA, Show = TRUE,
                         track_allParameters = track_allParameters){
  
  # If modelFile is not provided, use the default model file
  if(is.na(modelFile)){
    modelFile <- here("output", "BUGS-models", "EZHBDDM.bug")
  }

  # Step 1: Extract key summary statistics from the data
  # These are the sufficient statistics for the EZ-DDM
  correct <- EZ_stats$correct
  varRT   <- EZ_stats$varRT
  meanRT  <- EZ_stats$meanRT
  medianRT <- EZ_stats$medianRT
  iqrVarRT <- EZ_stats$iqrVarRT
  nTrialsPerCondition <- EZ_stats$nTrialsPerCondition
  nParticipants    <- EZ_stats$nParticipants
  
  # Step 2: Run JAGS to sample from the posterior distribution
  # Record start time for performance tracking
  tic <- clock::date_now(zone="UTC")
  
  # Call JAGS with the specified model and parameters
  # suppressMessages reduces console output during sampling
  suppressMessages(samples <- jags(data=jagsData, 
                                   parameters.to.save=jagsParameters, 
                                   model=modelFile, 
                                   n.chains=n.chains, 
                                   n.iter=n.iter, 
                                   n.burnin=n.burnin, 
                                   n.thin=n.thin, 
                                   DIC=T,           # Calculate Deviance Information Criterion
                                   inits=jagsInits))
  
  # Record end time and calculate duration
  toc <- clock::date_now(zone="UTC")
  clock <- as.numeric(toc-tic, units="secs")  # Computation time in seconds
  
  # Step 3: Extract MCMC samples and calculate convergence diagnostics
  object <- samples$BUGSoutput$sims.array  # 3D array: [iterations, chains, parameters]
  rhats <- apply(object, 3, getRhat)       # Calculate R-hat for each parameter
  
  # Step 4: Generate diagnostic plots if requested
  if(Show){  
    JAGS_plotChain(samples = samples)   
  }
  
  # Step 5: Process posterior samples for each monitored parameter
  estimates <- list()      # Will store posterior means
  error <- list()          # Will store posterior standard deviations
  credInterval <- list()   # Will store 95% credible intervals
  
  # Loop through each parameter to calculate summary statistics
  for(i in 1:length(jagsParameters)){
      # Extract samples for the current parameter
      posteriorParameters <- JAGS_extractSamples(jagsParameters[i], samples)
      
      # Handle differently depending on parameter dimensionality
      if(length(dim(posteriorParameters))==3){
         # For vector parameters (e.g., participant-level parameters)
         meanPost    <- apply(posteriorParameters, 3, mean)         # Mean for each element
         sdPost      <- apply(posteriorParameters, 3, sd)           # SD for each element
         percentiles <- apply(posteriorParameters, 3, quantile, probs=c(0.025, 0.975))  # 95% CI
      } else {
         # For scalar parameters (e.g., group-level parameters)
         meanPost    <- mean(posteriorParameters)                   # Overall mean
         sdPost      <- sd(posteriorParameters)                     # Overall SD
         percentiles <- quantile(posteriorParameters, probs = c(0.025, 0.975))  # 95% CI
      }
      
      # Store results in respective lists
      estimates    <- c(estimates, list(meanPost))
      error        <- c(error, list(sdPost))
      credInterval <- c(credInterval, list(percentiles))
  }
  
  # Name the list elements according to parameter names for easy access
  names(estimates)    <- jagsParameters
  names(error)        <- jagsParameters
  names(credInterval) <- jagsParameters
  
  # Step 6: Return comprehensive results object
  return(list("estimates" = estimates,       # Posterior means
              "sd" = error,                  # Posterior standard deviations
              "credInterval" = credInterval, # 95% credible intervals
              "rhats" = rhats,               # Convergence diagnostics
              "clock" = clock,               # Computation time
              "n.iter" = n.iter))            # Number of iterations
}

# Helper function
getEZ_stats <- function(sumData, nTrials){
   return(list("correct" = sumData[,"sum_correct"],  # Number of correct responses per participant
               "varRT" = sumData[,"varRT"],          # Variance of response times
               "meanRT" = sumData[,"meanRT"],        # Mean response time
               "medianRT" = sumData[,"medianRT"],    # Median response time
               "iqrVarRT" = sumData[,"iqrVarRT"],    # Interquartile range of response times
               "nTrialsPerCondition" = nTrials,         # Number of trials per participant
               "nParticipants" = nrow(sumData),
               "P" = sumData[,"sub"]))     # Number of participants
}

