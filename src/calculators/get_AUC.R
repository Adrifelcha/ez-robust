# Custom Function to compute the Area Under the Curve (AUC)
# from a vector of false positive rates (fpr) and a vector 
# of true positive rates (tpr)
###################################################################
compute_AUC <- function(fpr, tpr){

  # Sort by false positive rate
  sorted_fpr <- order(fpr)
  fpr <- fpr[sorted_fpr]
  tpr <- tpr[sorted_fpr]

  # Calculate AUC using the trapezoidal rule
  auc <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)])) / 2

return(auc)
}