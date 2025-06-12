#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function calculates the R-hat convergence diagnostic
#
# R-hat measures the convergence of MCMC chains by comparing the variance between
# chains to the variance within chains. Values close to 1.0 indicate good convergence.
# Generally, we consider values below 1.05 to be acceptable.
#
# Inputs:
# - posterior_chains: A matrix where columns store MCMC chains and rows store iterations
# - n.chains: Number of chains (only needed if posterior_chains has been reshaped into a vector)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
JAGS_Rhat <- function(runJags, withinSubject = FALSE, EZRobust = FALSE, separate_datasets = FALSE) {
    if(EZRobust){
        if(withinSubject){
            rhats <- runJags$BUGSoutput$Rhat[c("bound_mean", "drift_mean", "nondt_mean", "betaweight")]
        }else{
            rhats <- runJags$BUGSoutput$Rhat[c("bound_mean", "drift_mean", "nondt_mean", "betaweight")]
        }
    }
}

# Reference: https://medium.com/@nialloulton/understanding-the-r-hat-statistic-d83b3b5ca162

runJags$EZ$clean
runJags$EZ$contaminated
runJags$EZRobust$clean
runJags$EZRobust$contaminated
