#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function generates initial values for the individual drift rates
# to pass to JAGS.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
JAGS_inits <- function(n.chains, nParticipants, custom_sd=0.3, withinSubject=FALSE){
  # Create a list to hold initial values for each chain
  myinits <- rep(list(list()), n.chains)
  
  if(withinSubject==FALSE){
      # For each chain, generate random initial values
      for(i in 1:n.chains){
        # Initialize the `drift` parameter MCMC chain for each participant
        myinits[[i]] <- list(drift = rnorm(nParticipants, 0, custom_sd))
      }
  } else {
      for(i in 1:n.chains){
        myinits[[i]] <- list(drift = matrix(rep(rnorm(nParticipants,0,custom_sd),2),ncol=2, byrow = FALSE))
      }
  }
  return(myinits)
}