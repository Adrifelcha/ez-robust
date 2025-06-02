# This function creates an X object with participant-specific covariates
# according to the model type
# Inputs:
# - p: number of participants
# - modelType: type of model (hierarchical, ttest, metaregression)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_X_covariate <- function(p, modelType){
    if(modelType == "hierarchical"){
        X <- cbind(rep(NA, p))
        colnames(X) <- c("hierarchical")
    }
    if(modelType == "ttest"){
        X <- cbind((0:(p-1))%%2)
        colnames(X) <- c("ttest")
    }
    if(modelType == "metaregression"){
        X <- cbind((0:(p-1))/p)
        colnames(X) <- c("metaregression")
    }
    return(X)
}
