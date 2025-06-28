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
  suppressMessages(library(truncnorm))

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
          X <- get_X_covariate(p, modelType = d, withinSubject = settings$withinSubject)
          
          # Loop through all trial count levels
          for(t in settings$trial_levels){

                if(d == "hierarchical" || is.null(settings$beta_levels)){
                    beta_levels <- NA
                }else{
                    beta_levels <- settings$beta_levels
                }

                if(!is.null(settings$withinSubject)){
                    if(settings$withinSubject){       nTPC <- t         }
                }
                
                for(b in beta_levels){
                        write('Cell is running', paste(settings$log.folder, "seed-", seed, "_cell-", cell, "_p", p, "t", t, "b", b,".txt", sep=""))
                        # Flag to control R-hat checking loop
                        rhat_not_verified <- TRUE
                        
                        # Initialize default seed to master seed
                        this.seed <- seed
                        
                        # Display progress information
                        cat("Running cell", cell, "of", settings$nCells, "\n")
                        
                        # Run the cell
                        start_time <- Sys.time()
                        runCell <- simStudy_runCell(p = p, t = t, nTPC = nTPC, d = d, X = X, b = b, 
                                                    settings = settings, Show = Show, this.seed = this.seed,
                                                    prevent_zero_accuracy = prevent_zero_accuracy,
                                                    redo_if_bad_rhat=redo_if_bad_rhat, rhat_cutoff=rhat_cutoff,
                                                    nIter = nIter, nBurnin = nBurnin, nThin = nThin, nChains = nChains)
                                                        
                        # If JAGS error occurs, retry with different seed
                        while(runCell$JAGS_error){ 
                            cat("Repeating cell", cell, "of", settings$nCells, "due to a JAGS error \n")
                            this.seed <- this.seed + 10000  # Change seed by 10,000
                            
                            runCell <- simStudy_runCell(p = p, t = t, nTPC = nTPC, d = d, X = X, b = b, 
                                                        settings = settings, Show = Show, this.seed = this.seed,
                                                        prevent_zero_accuracy = prevent_zero_accuracy,
                                                        redo_if_bad_rhat=redo_if_bad_rhat, rhat_cutoff=rhat_cutoff,
                                                        nIter = nIter, nBurnin = nBurnin, nThin = nThin, nChains = nChains)

                            # Increment error counter and break if too many errors
                            redo_JAGS <- redo_JAGS + 1
                            if(redo_JAGS > 5){ 
                                break  # Give up after 5 attempts
                            }
                        }
                        end_time <- Sys.time()
                        if(Show){       
                                    cat("Time taken: ", difftime(end_time, start_time, units = "secs"), " seconds\n")        
                                }

                        runJAGS_output <- load_JAGS_cellResults(runJags = runCell$runJags,
                                                                parameter_set = runCell$design$parameter_set)
                        results_cell <- runJAGS_output$output

                        # Store results for this design cell by stacking them
                        if(d != "hierarchical"){ # Designs with betaweight parameter
                                if(is.na(b)){
                                    # Beta is random accross datasets
                                    out_Beta <- rbind(out_Beta, results_cell)            
                                } else {
                                    # Beta is fixed accross datasets
                                    if(b == 0){  # No effect
                                        out_NoEffect <- rbind(out_NoEffect, results_cell)
                                    } else {     # Fixed effect
                                        out_Effect <- rbind(out_Effect, results_cell)
                                    }
                                }
                        } else { # Hierarchical designs
                            out_H <- rbind(out_H, results_cell)
                        }
                        # Increment cell counter
                        cell <- 1 + cell
                        file.remove(paste(settings$log.folder, "seed-", seed, "_cell-", cell-1, "_p", p, "t", t, "b", b,".txt", sep=""))
                }
          }
      }
  }
  
  # Create a marker file to indicate simulation has completed
  grand_toc <- clock::date_now(zone="UTC")
  total_time <- difftime(grand_toc, grand_tic, units="mins")
  write(paste("Running this seed took ", total_time, "minutes.\n"), 
        paste(settings$log.folder, "seed-", seed, "_end.txt", sep=""))
  
  # Create and save output object
  # Start by storing the number of JAGS errors and R-hat issues
  output <- list("reps" = data.frame("bad_JAGS" = redo_JAGS,          # Count of JAGS errors
                                     "bad_Rhat" = redo_Rhat,
                                     "total_time" = total_time),         # Count of R-hat issues
                "settings" = settings)

    if(length(out_Effect) > 0){
        output$fixedEffect <- out_Effect
    }
    if(length(out_NoEffect) > 0){
        output$noEffect <- out_NoEffect
    }
    if(length(out_Beta) > 0){
        output$betaEffect <- out_Beta
    }
    if(length(out_H) > 0){
        output$hierarchical <- out_H
    }

  # Save results to file
  save(output, file=fileName)
  
  # Return results
  return(output)
}

# Helper function to combine results from multiple seeds
# Simply using rbind() doesn't work because the results are stored in a list
combine_results <- function(...) {
    # Get all results
    results <- list(...)
    
    # Initialize combined output
    combined <- list(
        reps = do.call(rbind, lapply(results, function(x) x$reps)),
        settings = results[[1]]$settings  # settings should be the same for all
    )
    
    # Combine noEffect if it exists in any result
    noEffect_list <- lapply(results, function(x) if(!is.null(x$noEffect)) x$noEffect else NULL)
    if(any(!sapply(noEffect_list, is.null))) {
        combined$noEffect <- do.call(rbind, noEffect_list[!sapply(noEffect_list, is.null)])
    }
    
    # Combine fixedEffect if it exists in any result
    fixedEffect_list <- lapply(results, function(x) if(!is.null(x$fixedEffect)) x$fixedEffect else NULL)
    if(any(!sapply(fixedEffect_list, is.null))) {
        combined$fixedEffect <- do.call(rbind, fixedEffect_list[!sapply(fixedEffect_list, is.null)])
    }
    
    return(combined)
}
