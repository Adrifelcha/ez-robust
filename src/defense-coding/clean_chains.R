
# Function to diagnose and clean problematic beta_chains
diagnose_and_clean_chains <- function(chains, true_betas, show_message = TRUE) {
  if(show_message){
      cat("Diagnosing chains structure...\n")
      cat("Class of chains:", class(chains), "\n")
      cat("Dimensions:", if(is.matrix(chains)) paste(dim(chains)) else paste(length(chains)), "\n")
  }
  # Handle the case where the chains object has a corrupted structure
  if (is.matrix(chains) && nrow(chains) != 1) {
    if(show_message){
        cat("WARNING: Corrupted structure detected!\n")
        cat("Matrix dimensions:", dim(chains), "\n")
    }    
    # Check if we can reconstruct the list structure
    n_seeds <- length(true_betas)
    validated_chains <- list()
    for(i in 1:n_seeds){
      validated_chains[[i]] <- matrix(chains[i, ], ncol = 3)
    }
    if(show_message){
       cat("Successfully converted matrix to list structure\n")
      }
  }else{
    if(show_message){    
       cat("No corrupted structure detected\n")    
       }
    validated_chains <- chains
  }
return(validated_chains)
}