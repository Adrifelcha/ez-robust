jags_localResults <- function(runJags_i, EZ_variation, data_variation, true_parameter_set){
    localResults <- list(seed = runJags_i$this.seed,           # Seed used for this cell (could change if R-hats were bad)
                         p = runJags_i$jags_data$nParticipants,                      # Number of participants
                         t = runJags_i$jags_data$nTrialsPerCondition,                # Number of trials
                         ez_type = EZ_variation,                                           # Design type                         
                         data_type = data_variation,
                         rhats = runJags_i$rhats,      # Convergence diagnostics
                         true.values = true_parameter_set,      # True parameter values
                         mean.estimates = runJags_i$estimates,         # Posterior means
                         std.estimates = runJags_i$estimates,          # Posterior SDs
                         elapsed.time = runJags_i$clock,               # Computation time
                         nIter = runJags_i$nIter,                      # Final used number of iterations
                         nBurnin = runJags_i$nBurnin,                 # Final used number of burn-in iterations
                         nThin = runJags_i$nThin,                     # Final used number of thinning iterations
                         bad_rhat_count = runJags_i$bad_rhat_count    # Number of times R-hats were bad                         
                        )
                        
    return(localResults)
}

load_JAGS_cellResults <- function(runJags, parameter_set){
    
    EZ_variations <- names(runJags)
    data_variations <- names(runJags$EZ)
    
    output <- list()
    for(EZ_var in EZ_variations){
        if(length(runJags[[EZ_var]]) > 1){            
            for(data_var in data_variations){                
                # Create name for this combination
                name <- paste0(EZ_var, "_", data_var)
                # Store results in named list element
                output[[name]] <- jags_localResults(runJags_i = runJags[[EZ_var]][[data_var]], 
                                                  EZ_variation = EZ_var, 
                                                  data_variation = data_var,
                                                  true_parameter_set = parameter_set)
            }            
        }else{
            name <- paste0(EZ_var, "_", runJags[[EZ_var]]$data_type)
            output[[name]] <- jags_localResults(runJags_i = runJags[[EZ_var]], 
                                              EZ_variation = EZ_var, 
                                              data_variation = runJags[[EZ_var]]$data_type,
                                              true_parameter_set = parameter_set)
        }     
    }

    return(output)
}