# This function creates an X object with participant-specific covariates
# according to the model type
# Inputs:
# - p: number of participants
# - modelType: type of model (hierarchical, ttest, metaregression)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
get_X_covariate <- function(p, modelType, withinSubject = FALSE){
    if(modelType == "hierarchical"){
        X <- cbind(rep(NA, p))
        colnames(X) <- c("hierarchical")
    }
    if(modelType == "ttest"){
        X <- if(withinSubject) {
            cbind(rep(c(1,0),p))  # Within-subject design
        } else {
            cbind((0:(p-1))%%2)   # Between-subject design
        }
        colnames(X) <- c("ttest")
    }
    if(modelType == "metaregression"){
        X <- cbind((0:(p-1))/p)
        colnames(X) <- c("metaregression")
    }
    return(X)
}
