#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function specifies which data objects need to be passed to the JAGS model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Inputs:
# - modelType: Type of model being used (e.g., "hierarchical", or "ttest")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
JAGS_passData <- function(modelType=NA, EZRobust=FALSE, withinSubject=FALSE){
  # Define the basic data variables needed for all model types
  passData <- list("nParticipants", "correct")

  if(EZRobust){
    passData <- c(passData, "medianRT", "iqrVarRT")
  }else{
    passData <- c(passData, "meanRT", "varRT")
  }

  if(withinSubject){
    passData <- c("nTrialsPerCondition", "P", passData)
  }else{
    passData <- c("nTrialsPerPerson", passData)    
  }
  
  # For models that include predictors (i.e., metaregression, t-test),
  # we also need to pass the predictor vector X
  if(!(is.na(modelType)|modelType=="hierarchical")){
          # Add the predictor variable X to the list of data to pass
          passData <- c(passData, "X")
  }
  
  # Return the list of data to pass to JAGS
  return(passData)
}