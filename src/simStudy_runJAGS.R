###################################################################################
# MODULAR FUNCTION CODING AHEAD:
#
# This is a higher-level function that manages the MCMC sampling process
# using JAGS to estimate hierarchical drift diffusion model parameters.
###################################################################################
# Contains three functions:
# 1. simStudy_runJAGS: Main function that orchestrates the MCMC sampling process
# 2. runJAGS: Function that runs the MCMC sampling process locally (for a gvien dataset and modelFile)
# 3. getEZ_stats: Function that extracts the sufficient statistics from the data
###################################################################################

###################################################################################
# Main function that orchestrates the MCMC sampling process
###################################################################################
simStudy_runJAGS <- function(summaryData, modelType, nTrials, X, jagsParameters, jagsInits, this.seed,
                             n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, modelFile=NA, Show = TRUE,
                             track_allParameters = track_allParameters, redo_if_bad_rhat=FALSE, rhat_cutoff= 1.05, 
                             separate_datasets = FALSE, include_EZ_Robust = FALSE, withinSubject = FALSE, contamination_probability=0){  
      # If modelFile is not provided, use the default model file
      if(length(modelFile) == 1 && is.na(modelFile)){
        modelFile <- here("output", "BUGS-models", "EZHBDDM.bug")
      }

      # Step 1: Extract key summary statistics from the data
      # These are the sufficient statistics for the EZ-DDM
      if(separate_datasets){
          EZ_stats_contaminated <- getEZ_stats(sumData = summaryData$contaminated_summary, 
                                               nTrials = nTrials, 
                                               withinSubject = withinSubject)
          EZ_stats_clean <- getEZ_stats(sumData = summaryData$clean_summary, 
                                        nTrials = nTrials, 
                                        withinSubject = withinSubject)
          
      }else{
          EZ_stats <- getEZ_stats(sumData = summaryData, nTrials = nTrials, withinSubject = withinSubject)
      }

      model  <- "EZ"
      if(include_EZ_Robust){
         model <- c(model, "EZRobust")
      }

      results <- list()
      for(m in model){
        if(separate_datasets){
          results[[m]] <- list("contaminated" = runJAGS(EZ_stats = EZ_stats_contaminated, 
                                                        EZmodel = m,
                                                        modelType = modelType,
                                                        withinSubject = withinSubject, this.seed = this.seed,
                                                        jagsParameters = jagsParameters, 
                                                        jagsInits = jagsInits,
                                                        n.chains = n.chains, n.burnin = n.burnin, 
                                                        n.iter = n.iter, n.thin = n.thin, 
                                                        modelFile=modelFile[m], Show = TRUE, 
                                                        track_allParameters = track_allParameters,
                                                        redo_if_bad_rhat=redo_if_bad_rhat, rhat_cutoff= rhat_cutoff),
                               "clean" = runJAGS(EZ_stats = EZ_stats_clean, 
                                                        EZmodel = m,
                                                        modelType = modelType,
                                                        withinSubject = withinSubject, this.seed = this.seed,
                                                        jagsParameters = jagsParameters, 
                                                        jagsInits = jagsInits,
                                                        n.chains = n.chains, n.burnin = n.burnin, 
                                                        n.iter = n.iter, n.thin = n.thin, 
                                                        modelFile=modelFile[m], Show = TRUE, 
                                                        track_allParameters = track_allParameters,
                                                        redo_if_bad_rhat=redo_if_bad_rhat, rhat_cutoff= rhat_cutoff))
        }else{
          if(contamination_probability==0){
            data_type <- "clean"
          }else{
            data_type <- "contaminated"
          }
          results[[m]] <- runJAGS(EZ_stats = EZ_stats,
                                  EZmodel = m,
                                  modelType = modelType,
                                  withinSubject = withinSubject, this.seed = this.seed,
                                  jagsParameters = jagsParameters, jagsInits = jagsInits,
                                  n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter, 
                                  n.thin = n.thin, modelFile = modelFile[m], Show = TRUE,
                                  track_allParameters = track_allParameters,
                                  redo_if_bad_rhat=redo_if_bad_rhat, rhat_cutoff= rhat_cutoff)
          results[[m]][["data_type"]] <- data_type
        }
      }
return(results)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Local function that runs the MCMC sampling process
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
runJAGS <- function(EZ_stats, EZmodel, modelType, withinSubject, this.seed,
                    jagsParameters, jagsInits, n.chains=4, n.burnin=250, n.iter=2000, n.thin=1, 
                    modelFile=NA, Show = TRUE, track_allParameters = track_allParameters,
                    redo_if_bad_rhat=FALSE, rhat_cutoff= 1.05){
  
  # If modelFile is not provided, use the default model file
  if(is.na(modelFile)){
    modelFile <- here("output", "BUGS-models", "EZHBDDM.bug")
  }

  # Step 1: Combine base data with summary statistics required by the model
  jagsData <- JAGS_passData(EZ_stats, modelType=modelType, EZmodel=EZmodel, withinSubject=withinSubject)

  # Step 2: Run JAGS to sample from the posterior distribution  
  # Flag to control R-hat checking loop
  rhat_not_verified <- TRUE
  this.seed <- seed
  nIter <- n.iter
  nBurnin <- n.burnin
  nThin <- n.thin  
  bad_rhat_count <- 0
  while(rhat_not_verified){
                    # Set up start time
                    tic <- clock::date_now(zone="UTC")
                    # Call JAGS (suppressMessages reduces console output during sampling
                    suppressMessages(samples <- jags(data=jagsData, 
                                                    parameters.to.save=jagsParameters, 
                                                    model=modelFile, 
                                                    n.chains=n.chains, 
                                                    n.iter=nIter, 
                                                    n.burnin=nBurnin, 
                                                    n.thin=nThin, 
                                                    DIC=T, inits=jagsInits))                    
                    # Record end time
                    toc <- clock::date_now(zone="UTC")
                    # Compute computation time
                    clock <- as.numeric(toc-tic, units="secs")  # Computation time in seconds
                    
                    # Extract MCMC samples and calculate convergence diagnostics
                    object <- samples$BUGSoutput$sims.array  # 3D array: [iterations, chains, parameters]
                    rhats <- apply(object, 3, getRhat)       # Calculate R-hat for each parameter

                    # Check if R-hats are too high
                    if(redo_if_bad_rhat){
                            if(any(rhats > rhat_cutoff)){
                                  bad_rhat_count <- bad_rhat_count + 1
                                  rhat_not_verified <- TRUE
                                  this.seed <- this.seed + 10000
                                  nIter <- round(nIter * 1.5,0)
                                  nBurnin <- round(nBurnin * 1.5,0)
                                  if(bad_rhat_count > 2){      nThin <- round(nThin * 2,0)      }
                                  set.seed(this.seed)
                            }else{    rhat_not_verified <- FALSE   }
                    }else{            rhat_not_verified <- FALSE   }
  }
  # Step 3: Generate diagnostic plots if requested
  if(Show){  
    JAGS_plotChain(samples = samples)
  }
  
  # Step 4: Process posterior samples for each monitored parameter
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
  
  # Step 5: Return comprehensive results object
  return(list("summStats" = EZ_stats,
              "beta_chains" = JAGS_extractSamples("betaweight", samples),
              "estimates" = estimates,       # Posterior means
              "sd" = error,                  # Posterior standard deviations
              "credInterval" = credInterval, # 95% credible intervals
              "rhats" = rhats,               # Convergence diagnostics
              "clock" = clock,               # Computation time
              "bad_rhat_count" = bad_rhat_count,
              "nIter" = nIter,
              "nBurnin" = nBurnin,
              "nThin" = nThin,
              "this.seed" = this.seed,
              "jags_data" = list("nTrialsPerCondition" = EZ_stats$nTrialsPerCondition,
                                 "nParticipants" = EZ_stats$nParticipants)))            # Number of iterations
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper function that extracts the sufficient statistics from the data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getEZ_stats <- function(sumData, nTrials, withinSubject = FALSE){
   return(list("correct" = sumData[,"sum_correct"],  # Number of correct responses per participant
               "varRT" = sumData[,"varRT"],          # Variance of response times
               "meanRT" = sumData[,"meanRT"],        # Mean response time
               "medianRT" = sumData[,"medianRT"],    # Median response time
               "iqrVarRT" = sumData[,"iqrVarRT"],    # Interquartile range of response times
               "nTrialsPerCondition" = nTrials,         # Number of trials per participant
               "nParticipants" = length(unique(sumData[,"sub"])),
               "P" = sumData[,"sub"],     # Number of participants
               "X" = if(withinSubject) sumData[,"cond"] else NULL))
}

