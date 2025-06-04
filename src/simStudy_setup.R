###################################################################################
# MACRO-FUNCTION: simStudy_setup
#
# This higher-level function orchestrates all data generation processes by calling
# the sub-functions that generate the parameter set, the raw data, and the summary data.
###################################################################################
# Inputs:
# - nPart: Number of participants to simulate
# - nTrials: Number of trials per participant
# - nTrialsPerCondition: Number of trials per condition (only used for hierarchical models)
# - true_sdevs: Standard deviation of the true parameter values
# - true_means: Mean of the true parameter values
# - modelType: Type of model ("hierarchical", "metaregression", "ttest", etc.)
# - X: Design matrix for models with predictors
# - Show: Whether to display the sampled parameters
# - prevent_zero_accuracy: Whether to prevent zero accuracy (TRUE) or allow it (FALSE)
# - fixedBeta: Fixed value of the beta parameter (only used for metaregression models)
# - withinSubject: Whether to simulate within-subject data (TRUE) or between-subject data (FALSE)
###################################################################################
simStudy_setup <- function(nPart, nTrials, nTrialsPerCondition=NULL, true_sdevs, true_means, modelType=NA, 
                           X=NA, Show=TRUE, prevent_zero_accuracy=TRUE, fixedBeta=NA, withinSubject=FALSE,
                           contamination_probability = 0, separate_datasets = FALSE){

    if(Show){
        local_settings = list("nPart" = nPart, "nTrials" = nTrials, "modelType" = modelType, "nTrialsPerCondition" = nTrialsPerCondition)
        show_design(local_settings, withinSubjectDesign = withinSubject)
    }
    # Step 1: Obtain parameter values to be used as ground truth in simulation studies    
    parameter_set <- get_simulation_parameters(true_means = true_means, true_sdevs = true_sdevs, Show = Show,
                                               nPart = nPart, modelType = modelType, X = X, 
                                               fixedBeta = fixedBeta, withinSubject = withinSubject)

    # Step 2: Generate hierarchical DDM data from the parameter set, using simulation settings
    rawData = get_simulation_data(nPart = nPart, nTrials = nTrials, parameter_set = parameter_set,
                                  nTrialsPerCondition = nTrialsPerCondition,
                                  contamination_probability = contamination_probability,
                                  separate_datasets = separate_datasets,
                                  prevent_zero_accuracy = prevent_zero_accuracy)

    # Step 3: Calculate EZ summary statistics from the raw data
    summData = get_summaryStats(data = rawData) 
    
    # Return all components needed for subsequent analysis
    return(list("parameter_set" = parameter_set,
                "rawData" = rawData, 
                "sumData" = summData))
}