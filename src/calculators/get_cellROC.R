# Custom Function to compute ROC curves for each beta level 
# against Beta = 0 for a given simulation study design cell
###################################################################
get_cellROCs <- function(resultsFile = NULL, epsilon = 0.05, levelsM = 1000) {
    
    # Load results file
    load(resultsFile)
    # Extract true betas and chains
    true_betas <- as.vector(unlist(simStudy_Beta$true[,"betaweight"]))
    chains <- simStudy_Beta$beta_chains
    
    # Compute Bayes Factors
    bf_results <- compute_BF_ROPE(true_betas, chains, epsilon = epsilon)
    # Get unique beta levels and sort them
    beta_levels <- sort(unique(bf_results$bayes_factors$true_beta))
    
    # Define range of log Bayes factor thresholds for ROC analysis
    m <- seq(-10, 10, length.out = levelsM)  # Range of log Bayes factor thresholds
    m[1] <- -Inf  # Set first threshold to negative infinity
    m[levelsM] <- Inf  # Set last threshold to positive infinity
    
    # Get null effect data for false positive rate calculation
    null_data <- bf_results$bayes_factors[bf_results$bayes_factors$true_beta == 0, ]
    null_log_bf <- null_data$log_bayes_factor
    
    # Calculate false positive rate (null effect)
    fpr <- numeric(levelsM)
    for(j in 1:levelsM){
        fpr[j] <- mean(null_log_bf > m[j], na.rm = TRUE)
    }
    
    # Initialize a list to store TPR for each beta level
    tpr_list <- list()
    
    # Calculate true positive rate for each beta level
    for(i in seq_along(beta_levels)){
        if(beta_levels[i] == 0){   next    }
        
        beta_data <- bf_results$bayes_factors[bf_results$bayes_factors$true_beta == beta_levels[i], ]
        log_bf <- beta_data$log_bayes_factor
        
        tpr <- numeric(levelsM)
        for(j in 1:levelsM){
            tpr[j] <- mean(log_bf > m[j], na.rm = TRUE)
        }
        
        tpr_list[[as.character(beta_levels[i])]] <- tpr
    }
    
    # Return the ROC data
    return(list(beta_levels = beta_levels[beta_levels != 0],
                thresholds = m,
                fpr = fpr,
                tpr_list = tpr_list,
                epsilon = epsilon))
} 