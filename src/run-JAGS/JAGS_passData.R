#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function specifies which data objects need to be passed to the JAGS model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Inputs:
# - modelType: Type of model being used (e.g., "hierarchical", or "ttest")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
JAGS_passData <- function(EZ_stats, modelType=NA, EZmodel=FALSE, withinSubject=FALSE){
  # Start with the basic data
  passData <- list(
    nParticipants = EZ_stats$nParticipants,
    correct = EZ_stats$correct
  )

  # Add RT summary statistics depending on model type
  if(EZmodel == "EZRobust"){
    passData$medianRT <- EZ_stats$medianRT
    passData$iqrVarRT <- EZ_stats$iqrVarRT
  } else {
    passData$meanRT <- EZ_stats$meanRT
    passData$varRT <- EZ_stats$varRT
  }

  # Add trial/participant indexing depending on design
  if(withinSubject){
    passData$nTrialsPerCondition <- EZ_stats$nTrialsPerCondition
    passData$P <- EZ_stats$P
  } else {
    passData$nTrialsPerPerson <- EZ_stats$nTrialsPerPerson
  }

  # Add predictor variable X if appropriate
  if(!(is.na(modelType) || modelType == "hierarchical")){
    passData$X <- EZ_stats$X
  }

  return(passData)
}