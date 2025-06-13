jags_localResults <- function(runJags_i, EZ_variation, data_variation, true_parameter_set, this.seed){
    localResults <- list(seed = this.seed,           # Seed used for this cell (could change if R-hats were bad)
                         p = runJags_i$jags_data$nParticipants,                      # Number of participants
                         t = runJags_i$jags_data$nTrialsPerCondition,                # Number of trials
                         ez_type = EZ_variation,                                           # Design type                         
                         data_type = data_variation,
                         rhats = runJags_i$rhats,      # Convergence diagnostics
                         true.values = true_parameter_set,      # True parameter values
                         mean.estimates = runJags_i$estimates,      # Posterior means
                         std.estimates = runJags_i$estimates,       # Posterior SDs
                         elapsed.time = runJags_i$clock             # Computation time
                        )
    return(localResults)
}

load_JAGS_cellResults <- function(runJags, this.seed, parameter_set){
    
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
                                                  this.seed = this.seed,
                                                  true_parameter_set = parameter_set)
            }            
        }else{
            name <- paste0(EZ_var, "_", runJags[[EZ_var]]$data_type)
            output[[name]] <- jags_localResults(runJags_i = runJags[[EZ_var]], 
                                              EZ_variation = EZ_var, 
                                              data_variation = runJags[[EZ_var]]$data_type,
                                              this.seed = this.seed, 
                                              true_parameter_set = parameter_set)
        }     
    }

    rhats <- c()
    for(name in names(output)){
        rhats <- c(rhats, output[[name]]$rhats)
    }    

    final_output <- list(rhats = rhats, output = output)
    return(final_output)
}