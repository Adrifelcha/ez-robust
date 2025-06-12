simStudy_runCell <- function(p, t, nTPC, d, X, b, settings, Show, prevent_zero_accuracy, this.seed){
            set.seed(this.seed)
            # Generate dataset with known parameters
            design <- simStudy_setup(nPart = p, nTrials = t, nTrialsPerCondition = nTPC, 
                                     true_sdevs = settings$true_sdevs, true_means = settings$true_means, 
                                     modelType = d, X = X, Show = Show, prevent_zero_accuracy = prevent_zero_accuracy,
                                     fixedBeta = b, withinSubject = settings$withinSubject,
                                     contamination_probability = settings$contaminant_prob,
                                     separate_datasets = settings$separate_datasets)
                                
            # Attempt to run JAGS with error handling
            start_time <- Sys.time()
            runJags <- NULL # Default to NULL
            # Run JAGS
            z <- try(runJags <- simStudy_runJAGS(
                                    summaryData = design$sumData, 
                                    nTrials = t, 
                                    X = X, 
                                    jagsData = settings$jagsData[[d]], 
                                    jagsParameters = settings$jagsParameters[,d], 
                                    jagsInits = settings$jagsInits[[as.character(p)]], 
                                    n.chains = settings$n.chains, 
                                    n.burnin = nBurnin, 
                                    n.iter = nIter, 
                                    n.thin = nThin, 
                                    modelFile = settings$modelFile[,d], 
                                    Show = Show, 
                                    track_allParameters = FALSE,
                                    separate_datasets = settings$separate_datasets,
                                    include_EZ_Robust = settings$include_EZ_Robust,
                                    withinSubject = settings$withinSubject))
            end_time <- Sys.time()
    TimeTaken <- difftime(end_time, start_time, units = "secs")
    JAGS_error <- inherits(z, "try-error")

    # Return a list containing the design, runJags, start_time, and end_time
    return(list("design" = design, "runJags" = runJags, "TimeTaken" = TimeTaken, "JAGS_error" = JAGS_error))
}