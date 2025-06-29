###############################################################################################################
# This function runs a single cell of the simulation study
#
# Inputs:
# - p: number of participants
# - t: number of trials
# - nTPC: number of trials per condition
# - d: model type
# - X: design matrix
# - b: fixed effects (if available)
# - settings: settings for the simulation study
# - Show: whether to show the output
# - prevent_zero_accuracy: whether to prevent zero accuracy
# - this.seed: the seed for the random number generator
###############################################################################################################
simStudy_runCell <- function(p, t, nTPC, d, X, b = NA, settings, Show, prevent_zero_accuracy, redo_if_bad_rhat=FALSE, 
                             rhat_cutoff= 1.05, this.seed, nIter, nBurnin, nThin, nChains){
            set.seed(this.seed)
            # Generate dataset with known parameters
            design <- simStudy_setup(nPart = p, nTrials = t, nTrialsPerCondition = nTPC, 
                                     true_sdevs = settings$true_sdevs, true_means = settings$true_means, 
                                     modelType = d, X = X, Show = Show, prevent_zero_accuracy = prevent_zero_accuracy,
                                     fixedBeta = b, withinSubject = settings$withinSubject,
                                     contamination_probability = settings$contaminant_prob,
                                     separate_datasets = settings$separate_datasets)
                                
            # Attempt to run JAGS with error handling            
            runJags <- NULL # Default to NULL
            # Run JAGS
            z <- try(runJags <- simStudy_runJAGS(
                                    summaryData = design$sumData, 
                                    modelType = d,
                                    nTrials = t, 
                                    X = X,                                     
                                    jagsParameters = settings$jagsParameters[,d], 
                                    jagsInits = settings$jagsInits[[as.character(p)]], 
                                    this.seed = this.seed,
                                    n.chains = nChains, 
                                    n.burnin = nBurnin, 
                                    n.iter = nIter, 
                                    n.thin = nThin, 
                                    modelFile = settings$modelFile[,d], 
                                    Show = Show, 
                                    track_allParameters = FALSE,
                                    redo_if_bad_rhat=redo_if_bad_rhat, 
                                    rhat_cutoff=rhat_cutoff,
                                    separate_datasets = settings$separate_datasets,
                                    include_EZ_Robust = settings$include_EZ_Robust,
                                    withinSubject = settings$withinSubject,
                                    contamination_probability = settings$contaminant_prob))    
    JAGS_error <- inherits(z, "try-error")

    # Return a list containing the design, runJags, start_time, and end_time
    return(list("design" = design, "runJags" = runJags, "JAGS_error" = JAGS_error))
}