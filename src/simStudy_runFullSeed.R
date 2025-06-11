###################################################################################
# MODULAR FUNCTION CODING AHEAD:
#
# HDDM_runFullSeed is a higher-level function intended for running
#         SIMULATION STUDIES with PARALLEL COMPUTING.
# It generates a single parameter recovery check across all simulation cells.
###################################################################################
# Inputs:
# - seed: Master random seed to use for all simulations
# - settings: List containing all simulation settings and design factors
# - forceRun: Whether to force new simulations even if results exist
# - redo_if_bad_rhat: Whether to repeat simulations with poor convergence
# - rhat_cutoff: Threshold for acceptable R-hat values (default: 1.05)
#
# Returns:
# - A list containing:
#   * hierarchical: Results from hierarchical model simulations
#   * betaEffect: Results from models with predictor effects (metaregression, t-test)
#   * reps: Count of repeated simulations due to errors or poor convergence
###################################################################################

simStudy_runFullSeed <- function(seed, settings, forceRun, prevent_zero_accuracy=TRUE,
                                 redo_if_bad_rhat=FALSE, rhat_cutoff=NA, Show=FALSE,
                                 include_EZ_Robust=FALSE){
  # Start timing the entire simulation study
  grand_tic <- clock::date_now(zone="UTC")  
  # Set master random seed for reproducibility
  set.seed(seed)
  
  # Load required libraries with suppressed messages
  suppressMessages(library(R2jags))
  suppressMessages(library(rstan))
  
  # Determine output file name to store seed-specific results
  fileName <- paste(settings$output.folder, "seed-", seed, ".RData", sep="")
  
  # Check if this seed has already been run
  if(file.exists(fileName) & !forceRun){
     stop("Seed already run")
  }
  
  # Set default R-hat cutoff if not specified
  if(redo_if_bad_rhat & is.na(rhat_cutoff)){
    rhat_cutoff <- 1.05
  }

  # Defensive coding: Simulation study may contain only one model type
  if(is.null(settings$design_levels)&&length(settings$modelType)==1){
    settings$design_levels <- settings$modelType
  }

  if(is.null(settings$separate_datasets)){ settings$separate_datasets <- FALSE }
  if(is.null(settings$contaminant_prob)){ settings$contaminant_prob <- 0 }

  if(settings$contaminant_prob==0){
        settings$separate_datasets <- FALSE
  }

  EZ <- ifelse(include_EZ_Robust, 2, 1)
  
  # Create a marker file to indicate simulation has started
  write('Seed has been initiated', paste(settings$output.folder, "seed-", seed, "_start.txt", sep=""))
  
  # Initialize storage for results - separate lists for hierarchical and regression structured models
  out_H <- list()     # Will store results from hierarchical models
  out_Beta <- list()  # Will store results from models with random Beta value
  out_NoEffect <- list()  # Will store results from models with no effect (Beta fixed to 0)
  out_Effect <- list()  # Will store results from models with effect (Beta fixed to anything other than 0)
  nTPC <- NULL        # Number of trials per condition - Default to NULL for between-subjects designs
  cell <- 1           # Counter for current cell
  redo_JAGS <- 0      # Counter for JAGS errors
  redo_Rhat <- 0      # Counter for R-hat convergence issues
 
   # Loop through all design factors (hierarchical, metaregression, t-test)
  for(d in settings$design_levels){
      
      # Loop through all participant count levels
      for(p in settings$participant_levels){
          # Create participant-specific covariates/predictors
          X <- get_X_covariate(p, modelType = d)
          
          # Loop through all trial count levels
          for(t in settings$trial_levels){

                if(d == "hierarchical" || is.null(settings$beta_level)){
                    beta_levels <- NA
                }else{
                    beta_levels <- settings$beta_level
                }

                if(!is.null(settings$withinSubject)){
                    if(settings$withinSubject){       nTPC <- t         }
                }
                
                for(b in beta_levels){
                        # Flag to control R-hat checking loop
                        rhat_not_verified <- TRUE
                        
                        # Initialize default seed to master seed
                        this.seed <- seed
                        
                        # Display progress information
                        cat("Running cell", cell, "of", settings$nCells, "\n")
                        
                        nIter <- settings$n.iter
                        nBurnin <- settings$n.burnin
                        nThin <- settings$n.thin

                        # Keep generating and analyzing datasets until R-hat criteria are met
                        while(rhat_not_verified){
                                # Set seed for this attempt
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
                                z <- try(runJags <- simStudy_runJAGS(
                                                        summaryData = design$sumData, 
                                                        nTrials = t, 
                                                        X = X, 
                                                        jagsData = settings$jagsData[[d]], 
                                                        jagsParameters = settings$jagsParameters[[d]], 
                                                        jagsInits = settings$jagsInits[[as.character(p)]], 
                                                        n.chains = settings$n.chains, 
                                                        n.burnin = nBurnin, 
                                                        n.iter = nIter, 
                                                        n.thin = nThin, 
                                                        modelFile = settings$modelFile[,d], 
                                                        Show = Show, 
                                                        track_allParameters = FALSE,
                                                        separate_datasets = settings$separate_datasets,
                                                        include_EZ_Robust = settings$include_EZ_Robust))
                                end_time <- Sys.time()
                                if(Show){
                                    cat("Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds\n")
                                }
                                # If JAGS error occurs, retry with different seed
                                if(inherits(z, "try-error")){ 
                                    cat("Repeating cell", cell, "of", settings$nCells, "due to a JAGS error \n")
                                    this.seed <- this.seed + 10000  # Change seed by 10,000
                                    set.seed(this.seed)
                                    
                                    # Generate new dataset and try again
                                    design <- HDDM_setup(priors = settings$priors[[d]], nPart = p, nTrials = t, 
                                                modelType = d, X = X[,d], criterion = c, 
                                                fromPrior = settings$fromPrior, Show = FALSE, 
                                                prevent_zero_accuracy = prevent_zero_accuracy,
                                                generative_uniforms = settings$generative_uniforms)
                                    z <- try(runJags <- HDDM_runJAGS(
                                        summaryData = design$sumData, 
                                        nTrials = t, 
                                        X = X[,d], 
                                        jagsData = settings$jagsData[[d]], 
                                        jagsParameters = settings$jagsParameters[[d]], 
                                        jagsInits = settings$jagsInits[[as.character(p)]], 
                                        n.chains = settings$n.chains, 
                                        n.burnin = nBurnin, 
                                        n.iter = nIter, 
                                        n.thin = nThin, 
                                        modelFile = settings$modelFile[d,c], 
                                        Show = Show, 
                                        track_allParameters = FALSE))

                                    # Increment error counter and break if too many errors
                                    redo_JAGS <- redo_JAGS + 1
                                    if(redo_JAGS > 5){ 
                                        break  # Give up after 5 attempts
                                    }
                                }

                                # Check if R-hat values indicate good convergence
                                count_bad_rhats <- sum(runJags$rhats[settings$jagsParameters[[d]]] > rhat_cutoff, na.rm = TRUE)
                                
                                # Exit loop if R-hat check is disabled or all R-hats are good
                                if((!redo_if_bad_rhat) | (count_bad_rhats == 0)){ 
                                    rhat_not_verified <- FALSE
                                } else { 
                                    # Otherwise, try again with different seed
                                    cat("Repeating cell", cell, "of", settings$nCells, "due to bad Rhats \n")
                                    this.seed <- this.seed + 10000
                                    redo_Rhat <- redo_Rhat + 1
                                }       
                                if(redo_Rhat>0){
                                            nIter <- nIter * 2
                                            nBurnin <- nBurnin * 2
                                            nThin <- nThin * 2
                                }                        
                        } # Close while() loop for R-hat verification
                }

                if(d != "hierarchical"){
                # Store results for this regression structured design cell
                out_Beta <- rbind(out_Beta, list(seed = this.seed,           # Seed used for this cell (could change if R-hats were bad)
                                                 p = p,                      # Number of participants
                                                 t = t,                      # Number of trials
                                                 d = d,                      # Design type
                                                 c = c,                      # Criterion parameter
                                                 rhats = runJags$rhats,      # Convergence diagnostics
                                                 true.values = design$parameter_set,      # True parameter values
                                                 mean.estimates = runJags$estimates,      # Posterior means
                                                 std.estimates = runJags$estimates,       # Posterior SDs
                                                 elapsed.time = runJags$clock             # Computation time
                                                ))            
                } else {
                    
                    # Store results for this hierarchical design cell
                    out_H <- rbind(out_H, list(
                        seed = this.seed,           # Seed used for this cell
                        p = p,                      # Number of participants
                        t = t,                      # Number of trials
                        d = d,                      # Design type
                        rhats = runJags$rhats,      # Convergence diagnostics
                        true.values = design$parameter_set,      # True parameter values
                        mean.estimates = runJags$estimates,      # Posterior means
                        std.estimates = runJags$estimates,       # Posterior SDs
                        elapsed.time = runJags$clock             # Computation time
                    ))                    
                }
                                # Increment cell counter
                cell <- 1 + cell

          }
      }
  }
  
  # Create a marker file to indicate simulation has completed
  grand_toc <- clock::date_now(zone="UTC")
  total_time <- difftime(grand_toc, grand_tic, units="mins")
  write(paste("Running this seed took ", total_time, "minutes.\n"), 
        paste(settings$output.folder, "seed-", seed, "_end.txt", sep=""))
  
  # Create and save output object
  # Start by storing the number of JAGS errors and R-hat issues
  output <- list("reps" = data.frame("bad_JAGS" = redo_JAGS,          # Count of JAGS errors
                                     "bad_Rhat" = redo_Rhat),         # Count of R-hat issues
                "settings" = settings)
  # Add the results from the hierarchical models (if any)
  if("hierarchical" %in% settings$design_levels){
    output <- c(output, list("hierarchical" = out_H))
  }
  # Add the results from the regression structured models (if any)
  if("ttest" %in% settings$design_levels | "metaregression" %in% settings$design_levels){
    output <- c(output, list("betaEffect" = out_Beta))
  }
  
  # Save results to file
  save(output, file=fileName)  
  
  # Return results
  return(output)
}

#x <- HDDM_runFullSeed(seed = 1, settings = settings, 
#                      forceRun = TRUE, redo_if_bad_rhat = TRUE, rhat_cutoff = 1.05)