#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function specifies which data objects need to be passed to the JAGS model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# Inputs:
# - modelType: Type of model being used (e.g., "hierarchical", or "ttest")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
JAGS_passData <- function(EZ_stats, modelType=NA, EZmodel=FALSE, withinSubject=FALSE){
  # Define the basic data variables needed for all model types
  passData <- list("nParticipants" = EZ_stats$nParticipants, 
                   "correct" = EZ_stats$correct)

  if(EZmodel == "EZRobust"){
    passData <- c(passData, "medianRT" = EZ_stats$medianRT, 
                            "iqrVarRT" = EZ_stats$iqrVarRT)
  }else{
    passData <- c(passData, "meanRT" = EZ_stats$meanRT, 
                            "varRT" = EZ_stats$varRT)
  }

  if(withinSubject){
    passData <- c("nTrialsPerCondition" = EZ_stats$nTrialsPerCondition, 
                  "P" = EZ_stats$P, passData)
  }else{
    passData <- c("nTrialsPerPerson" = EZ_stats$nTrialsPerPerson, passData)    
  }
  
  # For models that include predictors (i.e., metaregression, t-test),
  # we also need to pass the predictor vector X
  if(!(is.na(modelType)|modelType=="hierarchical")){
          # Add the predictor variable X to the list of data to pass
          passData <- c(passData, "X" = EZ_stats$X)
  }
  
  # Return the list of data to pass to JAGS
  return(passData)
}