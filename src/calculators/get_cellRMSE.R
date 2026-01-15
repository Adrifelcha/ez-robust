# Custom Function to compute RMSE, bias, and variance for parameter estimation
#   MSE = Bias² + Variance
#   RMSE = sqrt(MSE)
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
#   - mse_by_beta: a vector of MSEs for each beta level
#   - bias_by_beta: a vector of bias values for each beta level (mean(estimates - true))
#   - variance_by_beta: a vector of variance values for each beta level var(estimates - true)
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
    # Calculate RMSE, bias, and variance for each beta level
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Start empty vectors to store metrics for each beta level
    mse_by_beta <- rep(NA, length(beta_levels))
    rmse_by_beta <- rep(NA, length(beta_levels))    
    bias_by_beta <- rep(NA, length(beta_levels))
    variance_by_beta <- rep(NA, length(beta_levels))    
    names(mse_by_beta) <- as.character(beta_levels)
    names(rmse_by_beta) <- as.character(beta_levels)
    names(bias_by_beta) <- as.character(beta_levels)
    names(variance_by_beta) <- as.character(beta_levels)
    
    # Fill the vectors with metrics for each beta level
    for (i in seq_along(beta_levels)) {        
            # Extract the current beta level
            beta_val <- beta_levels[i]

            # Find indices for the current beta level
            beta_indices <- which(true_betas == beta_val)
            
            # Isolate the true and estimated values for the current beta level
            true_subset <- true_values[beta_indices]
            est_subset <- estimates[beta_indices]
            
            # Check available true and estimated values
            valid_indices <- !is.na(true_subset) & !is.na(est_subset)
            if (sum(valid_indices) > 0) {
                # Extract valid values
                true_valid <- true_subset[valid_indices]
                est_valid <- est_subset[valid_indices]
                
                error <- est_valid - true_valid
                
                # Calculate bias: mean(estimates - true)
                bias_by_beta[i] <- mean(error)
                
                # Calculate variance: var(estimates - true)                                
                variance_by_beta[i] <- var(error)
                
                # Calculate RMSE: sqrt(mean((estimates - true)^2))
                # Note: MSE = Bias² + Variance, so RMSE = sqrt(MSE)
                squared_errors <- error^2                
                mean_squared_errors <- mean(squared_errors)
                mse_by_beta[i] <- mean_squared_errors
                rmse_by_beta[i] <- sqrt(mean_squared_errors)
            } else {
                # If there are no valid true and estimated values, set all to NA
                mse_by_beta[i] <- NA
                rmse_by_beta[i] <- NA
                bias_by_beta[i] <- NA
                variance_by_beta[i] <- NA
            }
    }
    
    # Return the data
    return(list(parameter = parameter,
                beta_levels = beta_levels,
                rmse_by_beta = rmse_by_beta,
                mse_by_beta = mse_by_beta,
                bias_by_beta = bias_by_beta,
                variance_by_beta = variance_by_beta))
}

#resultsFile <- "output/RData_simStudy_results/EZ_clean/sim_P20T20_ez-Clean.RData"
#parameter <- "bound_mean"
#rmse_data <- get_cellRMSE(resultsFile = resultsFile, parameter = parameter)
#print(rmse_data)