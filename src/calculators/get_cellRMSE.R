# Custom Function to compute RMSE for parameter estimation
# Inputs:
# - resultsFile: the path to the results file
#   Ex: "output/RData_simStudy_results/EZ_clean/sim_P20T20_ez-Clean.RData"
# - parameter: the parameter to compute RMSE for
#   Ex: "bound_mean", "drift_mean", "nondt_mean", "betaweight"
# Output:
# - A list with the following elements:
#   - parameter: the parameter to compute RMSE for
#   - beta_levels: the beta levels used in the simulation study
#   - rmse_by_beta: a vector of RMSEs for each beta level
###################################################################

get_cellRMSE <- function(resultsFile = NULL, parameter = "bound_mean") {
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load results file
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    load(resultsFile)
    
    # Extract true values and estimates for the specified parameter
    true_values <- as.vector(unlist(simStudy_Beta$true[, parameter]))
    estimates <- as.vector(unlist(simStudy_Beta$estimates[, parameter]))
    
    # Extract true beta values for grouping
    true_betas <- as.vector(unlist(simStudy_Beta$true[, "betaweight"]))
    # Get unique beta levels and sort them
    beta_levels <- sort(unique(true_betas))
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate RMSE for each beta level
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Start empty vector to store RMSEs for each beta level
    rmse_by_beta <- rep(NA, length(beta_levels))
    names(rmse_by_beta) <- as.character(beta_levels)
    
    # Fill the RMSE vector with the RMSE for each beta level
    for (i in seq_along(beta_levels)) {        
            # Extract the current beta level
            beta_val <- beta_levels[i]

            # Find indices for the current beta level
            beta_indices <- which(true_betas == beta_val)
            
            # Isolate the true and estimated values for the current beta level
            true_subset <- true_values[beta_indices]
            est_subset <- estimates[beta_indices]
            
            # Calculate RMSE for this beta level, ignoring NA values
            # Check available true and estimated values
            valid_indices <- !is.na(true_subset) & !is.na(est_subset)
            if (sum(valid_indices) > 0) {
                # If there are at least one valid true and estimated value, calculate RMSE
                squared_errors <- (true_subset[valid_indices] - est_subset[valid_indices])^2
                rmse_by_beta[i] <- sqrt(mean(squared_errors))
            } else {
                # If there are no valid true and estimated values, set RMSE to NA
                rmse_by_beta[i] <- NA
            }
    }
    
    # Return the RMSE data
    return(list(parameter = parameter,
                beta_levels = beta_levels,
                rmse_by_beta = rmse_by_beta))
}

#resultsFile <- "output/RData_simStudy_results/EZ_clean/sim_P20T20_ez-Clean.RData"
#parameter <- "bound_mean"
#rmse_data <- get_cellRMSE(resultsFile = resultsFile, parameter = parameter)
#print(rmse_data)