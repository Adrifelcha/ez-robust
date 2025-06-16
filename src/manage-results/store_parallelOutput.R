

get_settings <- function(x){
    # Get the settings from the first result
    isolate_settings <- x[,"settings"]
    sample_results <- isolate_settings[[1]]
    return(sample_results)
}



store_parallelOutput <- function(output, settings, saveTo = "./", simStudyName = "genUniform"){
  #################################################################################
  # Identify relevant properties of the simulation study
  #################################################################################   
   nCells               <- settings$nCells
   output.folder        <- settings$output.folder
   allP   <- settings$participant_levels
   allT   <- settings$trial_levels
   allD   <- settings$design_levels
   allC   <- settings$criterion_levels
   nDatasets  <- nrow(output)

   
   # Is a simple Hierarchial study included?
   Hierarchical <- "hierarchical" %in% settings$design_levels
   # Is a Beta-effect design included?
   BetaEffect <- "ttest" %in% settings$design_levels | "metaregression" %in% settings$design_levels
   if(Hierarchical){   outH <- output[,"hierarchical"]    }
   if(BetaEffect){     outB <- output[,"betaEffect"]      }
   Ttest <- "ttest" %in% settings$design_levels
   Meta <- "metaregression" %in% settings$design_levels
   
   ################################################################################
   # Empty lists to story P-specific arrays with T-specific pages an nDataset rows
   ################################################################################
   if(Hierarchical){      trueVals_Hier <- list()                 
                          meanPosts_Hier <- list()
                          sdevPosts_Hier <- list()
                          rhats_Hier <- list()                    }
   if(Ttest){             trueVals_Ttst_nondt <- list()           
                          trueVals_Ttst_drift <- list()
                          trueVals_Ttst_bound <- list()
                          meanPosts_Ttst_nondt <- list()           
                          meanPosts_Ttst_drift <- list()
                          meanPosts_Ttst_bound <- list()
                          sdevPosts_Ttst_nondt <- list()
                          sdevPosts_Ttst_drift <- list()
                          sdevPosts_Ttst_bound <- list()
                          rhats_Ttst_nondt <- list()
                          rhats_Ttst_drift <- list()
                          rhats_Ttst_bound <- list()               }
   if(Meta){              trueVals_Meta_nondt <- list()
                          trueVals_Meta_drift <- list()    
                          trueVals_Meta_bound <- list()
                          meanPosts_Meta_nondt <- list()
                          meanPosts_Meta_drift <- list()    
                          meanPosts_Meta_bound <- list()
                          sdevPosts_Meta_nondt <- list()
                          sdevPosts_Meta_drift <- list()
                          sdevPosts_Meta_bound <- list()
                          rhats_Meta_nondt <- list()
                          rhats_Meta_drift <- list()   
                          rhats_Meta_bound <- list()               }
   # Store running times (in seconds)
   clock_base <- array(NA, dim=c(nDatasets,length(allT),length(allP)), 
                       dimnames = list(paste("seed", 1:nDatasets), paste("T", allT, sep=""), paste("P", allP,sep="")))
   if(Hierarchical){     clock_Hier <- clock_base                   }
   if(Ttest){ clock_Ttst <- list("bound" = clock_base, "drift" = clock_base, "nondt" = clock_base)}
   if(Meta){  clock_Meta <- list("bound" = clock_base, "drift" = clock_base, "nondt" = clock_base)}
   
   ###############################################################################################################################
   # Fill the empty storage objects
   ###############################################################################################################################
   i <- 1
   for(p in allP){       
       # No. parameters depend on No. of participants       
       if(Hierarchical){      nParamsH <- length(unlist(settings$jagsParameters["hierarchical"]))     }
       if(BetaEffect){                
            if(Ttest){ nParamsB <- length(unlist(settings$jagsParameters["ttest"]))     
            }else{     nParamsB <- length(unlist(settings$jagsParameters["metaregression"]))     }       
       }       
       
       # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Fill the empty lists created above with empty arrays to store output
       #################################################################################################
       if(Hierarchical){
           # Identify parameter names
           Hnames_true      <- sub('.*true.values.','',names(unlist(outH[[1]][which(outH[[1]][,"p"]==p)[1],"true.values"])))
           Hnames_estimates <- sub(".*\\.", "", sub('.*mean.estimates.','',names(unlist(outH[[1]][which(outH[[1]][,"p"]==p)[1],"mean.estimates"]))))
           Hnames_rhats     <- sub(".*\\.", "", names(unlist(outH[[1]][which(outH[[1]][,"p"]==p)[1],"rhats"])))
           # Create empty arrays to store True values
           emptyObj_Htrue <- array(NA, dim=c(nDatasets,nParamsH,length(allT)), 
                                   dimnames = list(paste("seed", 1:nDatasets), Hnames_true, paste("T",allT,sep="")))
           # Create empty arrays for estimates and errors
           emptyObj_H <- array(NA, dim=c(nDatasets,nParamsH,length(allT)), 
                               dimnames = list(paste("seed", 1:nDatasets), Hnames_estimates, paste("T",allT,sep="")))
           # Create empty arrays for Rhats
           emptyObj_R_H <- array(NA, dim=c(nDatasets,nParamsH+1,length(allT)), 
                                 dimnames = list(paste("seed", 1:nDatasets), Hnames_rhats, paste("T",allT,sep="")))
           # Add these arrays (specific number of columns) to each list
           trueVals_Hier <- c(trueVals_Hier, list(emptyObj_Htrue))
           meanPosts_Hier <- c(meanPosts_Hier, list(emptyObj_H))
           sdevPosts_Hier <- c(sdevPosts_Hier, list(emptyObj_H))
           rhats_Hier <- c(rhats_Hier, list(emptyObj_R_H))
       }
       if(BetaEffect){
           # Identify parameter names
           Bnames_true      <- sub('.*true.values.','',names(unlist(outB[[1]][which(outB[[1]][,"p"]==p)[1],"true.values"])))
           Bnames_estimates <- sub(".*\\.", "", sub('.*mean.estimates.','',names(unlist(outB[[1]][which(outB[[1]][,"p"]==p)[1],"mean.estimates"]))))
           Bnames_rhats     <- sub(".*\\.", "", names(unlist(outB[[1]][which(outB[[1]][,"p"]==p)[1],"rhats"])))
           # Create empty arrays to store True values
           emptyObj_Btrue <- array(NA, dim=c(nDatasets,length(Bnames_true),length(allT)), 
                               dimnames = list(paste("seed", 1:nDatasets), Bnames_true, paste("T",allT,sep="")))
           # Create empty arrays for estimates and errors
           emptyObj_B <- array(NA, dim=c(nDatasets,nParamsB,length(allT)), 
                               dimnames = list(paste("seed", 1:nDatasets), Bnames_estimates, paste("T",allT,sep="")))
           # Create empty arrays for Rhats
           emptyObj_R_B <- array(NA, dim=c(nDatasets,nParamsB+1,length(allT)), 
                                 dimnames = list(paste("seed", 1:nDatasets), Bnames_rhats, paste("T",allT,sep="")))
           # Add these arrays (specific number of columns) to each list
           if(Ttest){
                      trueVals_Ttst_nondt <- c(trueVals_Ttst_nondt, list(emptyObj_Btrue))     # True values
                      trueVals_Ttst_drift <- c(trueVals_Ttst_drift, list(emptyObj_Btrue))
                      trueVals_Ttst_bound <- c(trueVals_Ttst_bound, list(emptyObj_Btrue))
                      meanPosts_Ttst_nondt <- c(meanPosts_Ttst_nondt, list(emptyObj_B))       # Mean Posteriors
                      meanPosts_Ttst_drift <- c(meanPosts_Ttst_drift, list(emptyObj_B))
                      meanPosts_Ttst_bound <- c(meanPosts_Ttst_bound, list(emptyObj_B))
                      sdevPosts_Ttst_nondt <- c(sdevPosts_Ttst_nondt, list(emptyObj_B))       # Posterior variance
                      sdevPosts_Ttst_drift <- c(sdevPosts_Ttst_drift, list(emptyObj_B))
                      sdevPosts_Ttst_bound <- c(sdevPosts_Ttst_bound, list(emptyObj_B))
                      rhats_Ttst_nondt <- c(rhats_Ttst_nondt, list(emptyObj_R_B))             # Rhats
                      rhats_Ttst_drift <- c(rhats_Ttst_drift, list(emptyObj_R_B))
                      rhats_Ttst_bound <- c(rhats_Ttst_bound, list(emptyObj_R_B))
           }
           if(Meta){
                       trueVals_Meta_nondt <- c(trueVals_Meta_nondt, list(emptyObj_Btrue))    # True values
                       trueVals_Meta_drift <- c(trueVals_Meta_drift, list(emptyObj_Btrue))
                       trueVals_Meta_bound <- c(trueVals_Meta_bound, list(emptyObj_Btrue))
                       meanPosts_Meta_nondt <- c(meanPosts_Meta_nondt, list(emptyObj_B))      # Mean Posteriors
                       meanPosts_Meta_drift <- c(meanPosts_Meta_drift, list(emptyObj_B))
                       meanPosts_Meta_bound <- c(meanPosts_Meta_bound, list(emptyObj_B))
                       sdevPosts_Meta_nondt <- c(sdevPosts_Meta_nondt, list(emptyObj_B))      # Posterior variance
                       sdevPosts_Meta_drift <- c(sdevPosts_Meta_drift, list(emptyObj_B))
                       sdevPosts_Meta_bound <- c(sdevPosts_Meta_bound, list(emptyObj_B))
                       rhats_Meta_nondt <- c(rhats_Meta_nondt, list(emptyObj_R_B))            # Rhats
                       rhats_Meta_drift <- c(rhats_Meta_drift, list(emptyObj_R_B))
                       rhats_Meta_bound <- c(rhats_Meta_bound, list(emptyObj_R_B))             
           }
       }
       
       # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Fill the arrays corresponding to each level 'p'
       #################################################################################################
       for(k in 1:nDatasets){
           j <- 1
           for(t in allT){
               if(Hierarchical){
               thisH <- which(outH[[k]][,"p"]==p&outH[[k]][,"t"]==t)
                        trueVals_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"true.values"])
                        meanPosts_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"mean.estimates"])
                        sdevPosts_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"std.estimates"])
                        rhats_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"rhats"])
                        clock_Hier[k,j,i] <- as.numeric(outH[[k]][thisH,"elapsed.time"])
               }
               if(Meta){
                       # Meta-regression - Criterion: nondt
                       thisMn <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="nondt"&outB[[k]][,"d"]=="metaregression")
                       if(length(thisMn)>0){
                         trueVals_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"true.values"])
                         meanPosts_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"mean.estimates"])
                         sdevPosts_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"std.estimates"])
                         rhats_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"rhats"])
                         clock_Meta[["nondt"]][k,j,i] <- as.numeric(outB[[k]][thisMn,"elapsed.time"])
                       }
                       # Meta-regression - Criterion: drift
                       thisMd <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="drift"&outB[[k]][,"d"]=="metaregression")
                       if(length(thisMd)>0){
                         trueVals_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"true.values"])
                         meanPosts_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"mean.estimates"])
                         sdevPosts_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"std.estimates"])
                         rhats_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"rhats"])
                         clock_Meta[["drift"]][k,j,i] <- as.numeric(outB[[k]][thisMd,"elapsed.time"])
                       }
                       # Meta-regression - Criterion: bound
                       thisMb <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="bound"&outB[[k]][,"d"]=="metaregression")
                       if(length(thisMb)>0){
                         trueVals_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"true.values"])
                         meanPosts_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"mean.estimates"])
                         sdevPosts_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"std.estimates"])
                         rhats_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"rhats"])
                         clock_Meta[["bound"]][k,j,i] <- as.numeric(outB[[k]][thisMb,"elapsed.time"])
                       }
               }
               if(Ttest){
                       # ttest - Criterion: nondt
                       thisTn <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="nondt"&outB[[k]][,"d"]=="ttest")
                       if(length(thisTn)>0){
                         trueVals_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"true.values"])
                         meanPosts_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"mean.estimates"])
                         sdevPosts_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"std.estimates"])
                         rhats_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"rhats"])
                         clock_Ttst[["nondt"]][k,j,i] <- as.numeric(outB[[k]][thisTn,"elapsed.time"])
                       }
                       # ttest - Criterion: drift
                       thisTd <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="drift"&outB[[k]][,"d"]=="ttest")
                       if(length(thisTd)>0){
                         trueVals_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"true.values"])
                         meanPosts_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"mean.estimates"])
                         sdevPosts_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"std.estimates"])
                         rhats_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"rhats"])
                         clock_Ttst[["drift"]][k,j,i] <- as.numeric(outB[[k]][thisTd,"elapsed.time"])
                       }
                       # ttest - Criterion: bound
                       thisTb <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="bound"&outB[[k]][,"d"]=="ttest")
                       if(length(thisTb)>0){
                         trueVals_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"true.values"])
                         meanPosts_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"mean.estimates"])
                         sdevPosts_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"std.estimates"])
                         rhats_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"rhats"])
                         clock_Ttst[["bound"]][k,j,i] <- as.numeric(outB[[k]][thisTb,"elapsed.time"])
                       }
                }
           j <- j+1 
           }
       }
       i <- i + 1
   }

   ################################################################################
   # Name each array by the number-of-participants included
   ################################################################################
   if(Hierarchical){
                      names(trueVals_Hier) <- paste("P",allP,sep="") 
                      names(meanPosts_Hier) <- paste("P",allP,sep="")
                      names(sdevPosts_Hier) <- paste("P",allP,sep="")
                      names(rhats_Hier) <- paste("P",allP,sep="")
   }
   if(Ttest){
                       names(trueVals_Ttst_nondt) <- paste("P",allP,sep="")     # True values
                       names(trueVals_Ttst_drift) <- paste("P",allP,sep="")
                       names(trueVals_Ttst_bound) <- paste("P",allP,sep="")
                       names(meanPosts_Ttst_nondt) <- paste("P",allP,sep="")    # Mean Posteriors
                       names(meanPosts_Ttst_drift) <- paste("P",allP,sep="")
                       names(meanPosts_Ttst_bound) <- paste("P",allP,sep="")
                       names(sdevPosts_Ttst_nondt) <- paste("P",allP,sep="")    # Posterior variance
                       names(sdevPosts_Ttst_drift) <- paste("P",allP,sep="")
                       names(sdevPosts_Ttst_bound) <- paste("P",allP,sep="")
                       names(rhats_Ttst_nondt) <- paste("P",allP,sep="")        # Rhats
                       names(rhats_Ttst_drift) <- paste("P",allP,sep="")
                       names(rhats_Ttst_bound) <- paste("P",allP,sep="")
   }
   if(Meta){
                       names(trueVals_Meta_nondt) <- paste("P",allP,sep="")     # True values
                       names(trueVals_Meta_drift) <- paste("P",allP,sep="")
                       names(trueVals_Meta_bound) <- paste("P",allP,sep="")
                       names(meanPosts_Meta_nondt) <- paste("P",allP,sep="")    # Mean Posteriors
                       names(meanPosts_Meta_drift) <- paste("P",allP,sep="")   
                       names(meanPosts_Meta_bound) <- paste("P",allP,sep="")
                       names(sdevPosts_Meta_nondt) <- paste("P",allP,sep="")    # Posterior variance
                       names(sdevPosts_Meta_drift) <- paste("P",allP,sep="")
                       names(sdevPosts_Meta_bound) <- paste("P",allP,sep="")
                       names(rhats_Meta_nondt) <- paste("P",allP,sep="")        # Rhats
                       names(rhats_Meta_drift) <- paste("P",allP,sep="")
                       names(rhats_Meta_bound) <- paste("P",allP,sep="")
   }
   
   ################################################################################
   # FINAL: Store the output per design-parameter by writing an .RData object
   ################################################################################
   if(Hierarchical){
                     simStudy_Hierarchical <- list("true" = trueVals_Hier,
                                                   "recovered" = meanPosts_Hier,
                                                   "estimates_sdev" = sdevPosts_Hier,
                                                   "rhats" = rhats_Hier)
                     outputFile = paste(saveTo,"/simStudy_Hierarchical.RData",sep="")
                     save(simStudy_Hierarchical, file=outputFile)
   }
   if(Meta){
            if("nondt" %in% settings$criterion_levels){
                     simStudy_Meta_nondt <- list("true" = trueVals_Meta_nondt, 
                                                 "recovered" = meanPosts_Meta_nondt,
                                                 "estimates_sdev" = sdevPosts_Meta_nondt,
                                                 "rhats" = rhats_Meta_nondt)
                     save(simStudy_Meta_nondt, file=paste(saveTo,"/sim_",simStudyName,"_Meta_nondt.RData",sep=""))          }
            if("drift" %in% settings$criterion_levels){
                     simStudy_Meta_drift <- list("true" = trueVals_Meta_drift,
                                                 "recovered" = meanPosts_Meta_drift,
                                                 "estimates_sdev" = sdevPosts_Meta_drift,
                                                 "rhats" = rhats_Meta_drift)
                     save(simStudy_Meta_drift, file=paste(saveTo,"/sim_",simStudyName,"_Meta_drift.RData",sep=""))           }
            if("bound" %in% settings$criterion_levels){
                     simStudy_Meta_bound <- list("true" = trueVals_Meta_bound,
                                                 "recovered" = meanPosts_Meta_bound,
                                                 "estimates_sdev" = sdevPosts_Meta_bound,
                                                 "rhats" = rhats_Meta_bound)
                     save(simStudy_Meta_bound, file=paste(saveTo,"/sim_",simStudyName,"_Meta_bound.RData",sep=""))           }
   }
   if(Ttest){
                     simStudy_Ttst_nondt <- list("true" = trueVals_Ttst_nondt,
                                                 "recovered" = meanPosts_Ttst_nondt,
                                                 "estimates_sdev" = sdevPosts_Ttst_nondt,
                                                 "rhats" = rhats_Ttst_nondt)
                     save(simStudy_Ttst_nondt, file=paste(saveTo,"/sim_",simStudyName,"_Ttest_nondt.RData", sep=""))
                     
                     simStudy_Ttest_drift <- list("true" = trueVals_Ttst_drift,
                                                  "recovered" = meanPosts_Ttst_drift,
                                                  "estimates_sdev" = sdevPosts_Ttst_drift,
                                                  "rhats" = rhats_Ttst_drift)
                     save(simStudy_Ttest_drift, file=paste(saveTo,"/sim_",simStudyName,"_Ttest_drift.RData",sep=""))
                     
                     simStudy_Ttest_bound <- list("true" = trueVals_Ttst_bound, 
                                                  "recovered" = meanPosts_Ttst_bound,
                                                  "estimates_sdev" = sdevPosts_Ttst_bound,
                                                  "rhats" = rhats_Ttst_bound)
                     save(simStudy_Ttest_bound, file=paste(saveTo,"/sim_",simStudyName,"_Ttest_bound.RData",sep=""))
   }
} 



#######################################################
# S T O R E   P A R A L L E L   O U T P U T
#################################################################################
# This function processes and stores the results from a parallel simulation study
# 
# Inputs:
#   "output": Combined results from parallel processing (from foreach %dopar%)
#   "settings": List containing simulation parameters
#   "saveTo": Directory path where to save the processed results
#################################################################################

store_BetaParallelOutput <- function(output, settings, saveTo = NA){  
  # Extract simulation settings from the settings list
  nDatasets      <- nrow(output)
  nParams        <- length(settings$jagsParameters)  
  beta_levels    <- settings$beta_levels    # The effect sizes tested
  nchain         <- settings$n.chains       # Number of MCMC chains used
  n.iter         <- settings$n.iter         # Number of MCMC iterations
  n.burnin       <- settings$n.burnin       # Number of MCMC burn-in iterations
  n.thin         <- settings$n.thin         # Number of MCMC thinning
  nIter <- (n.iter - n.burnin) / n.thin
  
  # Create file names for different beta levels (B, B1, B2, etc.)
  B.files <- paste("B", c("", 1:(length(settings$beta_levels)-1)), sep="")
    
  
  # Process each beta level separately
  i <- 1
  for(b in beta_levels){
      # Select appropriate data based on whether this is a null effect (beta=0) or not
      if(b == 0){ 
          out = output[,"noEffect"]  
      } else {  
          out = output[,"fixedEffect"]  
      }
      
      # Initialize arrays to store timing information
      clock <- rep(NA, nDatasets)   # JAGS runtime
      clock2 <- rep(NA, nDatasets)  # Total runtime
      
      # Initialize arrays to store re-run information
      rerun_jags <- rep(NA, nDatasets)  # Count of JAGS failures
      rerun_rhat <- rep(NA, nDatasets)  # Count of Rhat failures
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Extract parameter names from the first dataset to use as array dimensions
      #################################################################################################
      # Find the first result with this beta level
      first_result_idx <- which(out[[1]][,"beta"] == b)[1]
      
      # Extract parameter names from different components of the results
      names_true      <- sub('.*true.values.', '', names(unlist(out[[1]][first_result_idx, "true.values"])))
      names_estimates <- sub(".*\\.", "", sub('.*mean.estimates.', '', names(unlist(out[[1]][first_result_idx, "mean.estimates"]))))
      names_rhats     <- sub(".*\\.", "", names(unlist(out[[1]][first_result_idx, "rhats"])))
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Create empty arrays to store aggregated results
      #################################################################################################
      # Array for beta parameter chains across all datasets
      betaChains <- array(NA, dim=c(nIter, nchain, nDatasets))
      
      # Array for true parameter values
      trueVals <- array(NA, dim=c(nDatasets, nParams), 
                       dimnames = list(paste("seed", 1:nDatasets), names_true))
      
      # Array for posterior mean estimates
      meanPosts <- array(NA, dim=c(nDatasets, nParams), 
                        dimnames = list(paste("seed", 1:nDatasets), names_estimates))
      
      # Array for posterior standard deviations
      sdevPosts <- array(NA, dim=c(nDatasets, nParams), 
                        dimnames = list(paste("seed", 1:nDatasets), names_estimates))
      
      # Array for MCMC convergence diagnostics (Rhat values)
      rhats <- array(NA, dim=c(nDatasets, nParams+1), 
                    dimnames = list(paste("seed", 1:nDatasets), names_rhats))
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Extract data from each dataset and fill the arrays
      #################################################################################################
      for(k in 1:nDatasets){
          # Find which row in this result corresponds to the current beta level
          thisH <- as.numeric(which(out[[k]][,"beta"] == b))
          
          # Extract values and store in the appropriate arrays
          trueVals[k,] <- unlist(out[[k]][thisH, "true.values"])
          meanPosts[k,] <- unlist(out[[k]][thisH, "mean.estimates"])
          sdevPosts[k,] <- unlist(out[[k]][thisH, "std.estimates"])
          rhats[k,] <- unlist(out[[k]][thisH, "rhats"])
          betaChains[,,k] <- as.matrix(out[[k]][thisH, "beta_chains"]$beta_chains)
          
          # Store timing information
          clock[k] <- as.numeric(out[[k]][thisH, "jags.time"])
          clock2[k] <- as.numeric(out[[k]][thisH, "total.time"])
          
          # Extract re-run information if available
          if("reps" %in% names(output[k,])) {
            rerun_jags[k] <- output[k,"reps"]$bad_JAGS
            rerun_rhat[k] <- output[k,"reps"]$bad_Rhat
          }
      }
   
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Create a list with all processed results for this beta level
      #################################################################################################
      simStudy_Beta <- list(
          "true" = trueVals,                # True parameter values
          "estimates" = meanPosts,          # Posterior mean estimates
          "variance" = sdevPosts^2,         # Posterior variance (square of SD)
          "rhats" = rhats,                  # MCMC convergence diagnostics
          "jagsTime" = clock,               # JAGS runtime
          "totalTime" = clock2,             # Total runtime
          "settings" = settings,            # Simulation settings
          "beta_chain" = betaChains,        # Beta parameter chains
          "reruns" = data.frame(            # New: add re-run information
              jags = rerun_jags,
              rhat = rerun_rhat
          )
      )
      
      # Create the output file name with participant and trial counts
      # Normalize the path to remove any double slashes
      if(is.na(saveTo)){
        saveTo <- here("output", "RData-results")
        dir.create(saveTo, recursive = TRUE, showWarnings = FALSE)
      }

      saveTo <- gsub("/$", "", saveTo)  # Remove trailing slash if present
      outputFile <- file.path(saveTo, paste0("simHypTesting_", B.files[i], ".RData"))
      
      # Save the processed results to a file
      save(simStudy_Beta, file=outputFile)
      
      # Move to the next beta level
      i <- i + 1
   }
} 





store_parallelOutput <- function(output, settings, saveTo = "./", simStudyName = "genUniform"){
  #################################################################################
  # Identify relevant properties of the simulation study
  #################################################################################   
   nCells               <- settings$nCells
   output.folder        <- settings$output.folder
   allP   <- settings$participant_levels
   allT   <- settings$trial_levels
   allD   <- settings$design_levels
   allC   <- settings$criterion_levels
   nDatasets  <- nrow(output)

   
   # Is a simple Hierarchial study included?
   Hierarchical <- "hierarchical" %in% settings$design_levels
   # Is a Beta-effect design included?
   BetaEffect <- "ttest" %in% settings$design_levels | "metaregression" %in% settings$design_levels
   if(Hierarchical){   outH <- output[,"hierarchical"]    }
   if(BetaEffect){     outB <- output[,"betaEffect"]      }
   Ttest <- "ttest" %in% settings$design_levels
   Meta <- "metaregression" %in% settings$design_levels
   
   ################################################################################
   # Empty lists to story P-specific arrays with T-specific pages an nDataset rows
   ################################################################################
   if(Hierarchical){      trueVals_Hier <- list()                 
                          meanPosts_Hier <- list()
                          sdevPosts_Hier <- list()
                          rhats_Hier <- list()                    }
   if(Ttest){             trueVals_Ttst_nondt <- list()           
                          trueVals_Ttst_drift <- list()
                          trueVals_Ttst_bound <- list()
                          meanPosts_Ttst_nondt <- list()           
                          meanPosts_Ttst_drift <- list()
                          meanPosts_Ttst_bound <- list()
                          sdevPosts_Ttst_nondt <- list()
                          sdevPosts_Ttst_drift <- list()
                          sdevPosts_Ttst_bound <- list()
                          rhats_Ttst_nondt <- list()
                          rhats_Ttst_drift <- list()
                          rhats_Ttst_bound <- list()               }
   if(Meta){              trueVals_Meta_nondt <- list()
                          trueVals_Meta_drift <- list()    
                          trueVals_Meta_bound <- list()
                          meanPosts_Meta_nondt <- list()
                          meanPosts_Meta_drift <- list()    
                          meanPosts_Meta_bound <- list()
                          sdevPosts_Meta_nondt <- list()
                          sdevPosts_Meta_drift <- list()
                          sdevPosts_Meta_bound <- list()
                          rhats_Meta_nondt <- list()
                          rhats_Meta_drift <- list()   
                          rhats_Meta_bound <- list()               }
   # Store running times (in seconds)
   clock_base <- array(NA, dim=c(nDatasets,length(allT),length(allP)), 
                       dimnames = list(paste("seed", 1:nDatasets), paste("T", allT, sep=""), paste("P", allP,sep="")))
   if(Hierarchical){     clock_Hier <- clock_base                   }
   if(Ttest){ clock_Ttst <- list("bound" = clock_base, "drift" = clock_base, "nondt" = clock_base)}
   if(Meta){  clock_Meta <- list("bound" = clock_base, "drift" = clock_base, "nondt" = clock_base)}
   
   ###############################################################################################################################
   # Fill the empty storage objects
   ###############################################################################################################################
   i <- 1
   for(p in allP){       
       # No. parameters depend on No. of participants       
       if(Hierarchical){      nParamsH <- length(unlist(settings$jagsParameters["hierarchical"]))     }
       if(BetaEffect){                
            if(Ttest){ nParamsB <- length(unlist(settings$jagsParameters["ttest"]))     
            }else{     nParamsB <- length(unlist(settings$jagsParameters["metaregression"]))     }       
       }       
       
       # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Fill the empty lists created above with empty arrays to store output
       #################################################################################################
       if(Hierarchical){
           # Identify parameter names
           Hnames_true      <- sub('.*true.values.','',names(unlist(outH[[1]][which(outH[[1]][,"p"]==p)[1],"true.values"])))
           Hnames_estimates <- sub(".*\\.", "", sub('.*mean.estimates.','',names(unlist(outH[[1]][which(outH[[1]][,"p"]==p)[1],"mean.estimates"]))))
           Hnames_rhats     <- sub(".*\\.", "", names(unlist(outH[[1]][which(outH[[1]][,"p"]==p)[1],"rhats"])))
           # Create empty arrays to store True values
           emptyObj_Htrue <- array(NA, dim=c(nDatasets,nParamsH,length(allT)), 
                                   dimnames = list(paste("seed", 1:nDatasets), Hnames_true, paste("T",allT,sep="")))
           # Create empty arrays for estimates and errors
           emptyObj_H <- array(NA, dim=c(nDatasets,nParamsH,length(allT)), 
                               dimnames = list(paste("seed", 1:nDatasets), Hnames_estimates, paste("T",allT,sep="")))
           # Create empty arrays for Rhats
           emptyObj_R_H <- array(NA, dim=c(nDatasets,nParamsH+1,length(allT)), 
                                 dimnames = list(paste("seed", 1:nDatasets), Hnames_rhats, paste("T",allT,sep="")))
           # Add these arrays (specific number of columns) to each list
           trueVals_Hier <- c(trueVals_Hier, list(emptyObj_Htrue))
           meanPosts_Hier <- c(meanPosts_Hier, list(emptyObj_H))
           sdevPosts_Hier <- c(sdevPosts_Hier, list(emptyObj_H))
           rhats_Hier <- c(rhats_Hier, list(emptyObj_R_H))
       }
       if(BetaEffect){
           # Identify parameter names
           Bnames_true      <- sub('.*true.values.','',names(unlist(outB[[1]][which(outB[[1]][,"p"]==p)[1],"true.values"])))
           Bnames_estimates <- sub(".*\\.", "", sub('.*mean.estimates.','',names(unlist(outB[[1]][which(outB[[1]][,"p"]==p)[1],"mean.estimates"]))))
           Bnames_rhats     <- sub(".*\\.", "", names(unlist(outB[[1]][which(outB[[1]][,"p"]==p)[1],"rhats"])))
           # Create empty arrays to store True values
           emptyObj_Btrue <- array(NA, dim=c(nDatasets,length(Bnames_true),length(allT)), 
                               dimnames = list(paste("seed", 1:nDatasets), Bnames_true, paste("T",allT,sep="")))
           # Create empty arrays for estimates and errors
           emptyObj_B <- array(NA, dim=c(nDatasets,nParamsB,length(allT)), 
                               dimnames = list(paste("seed", 1:nDatasets), Bnames_estimates, paste("T",allT,sep="")))
           # Create empty arrays for Rhats
           emptyObj_R_B <- array(NA, dim=c(nDatasets,nParamsB+1,length(allT)), 
                                 dimnames = list(paste("seed", 1:nDatasets), Bnames_rhats, paste("T",allT,sep="")))
           # Add these arrays (specific number of columns) to each list
           if(Ttest){
                      trueVals_Ttst_nondt <- c(trueVals_Ttst_nondt, list(emptyObj_Btrue))     # True values
                      trueVals_Ttst_drift <- c(trueVals_Ttst_drift, list(emptyObj_Btrue))
                      trueVals_Ttst_bound <- c(trueVals_Ttst_bound, list(emptyObj_Btrue))
                      meanPosts_Ttst_nondt <- c(meanPosts_Ttst_nondt, list(emptyObj_B))       # Mean Posteriors
                      meanPosts_Ttst_drift <- c(meanPosts_Ttst_drift, list(emptyObj_B))
                      meanPosts_Ttst_bound <- c(meanPosts_Ttst_bound, list(emptyObj_B))
                      sdevPosts_Ttst_nondt <- c(sdevPosts_Ttst_nondt, list(emptyObj_B))       # Posterior variance
                      sdevPosts_Ttst_drift <- c(sdevPosts_Ttst_drift, list(emptyObj_B))
                      sdevPosts_Ttst_bound <- c(sdevPosts_Ttst_bound, list(emptyObj_B))
                      rhats_Ttst_nondt <- c(rhats_Ttst_nondt, list(emptyObj_R_B))             # Rhats
                      rhats_Ttst_drift <- c(rhats_Ttst_drift, list(emptyObj_R_B))
                      rhats_Ttst_bound <- c(rhats_Ttst_bound, list(emptyObj_R_B))
           }
           if(Meta){
                       trueVals_Meta_nondt <- c(trueVals_Meta_nondt, list(emptyObj_Btrue))    # True values
                       trueVals_Meta_drift <- c(trueVals_Meta_drift, list(emptyObj_Btrue))
                       trueVals_Meta_bound <- c(trueVals_Meta_bound, list(emptyObj_Btrue))
                       meanPosts_Meta_nondt <- c(meanPosts_Meta_nondt, list(emptyObj_B))      # Mean Posteriors
                       meanPosts_Meta_drift <- c(meanPosts_Meta_drift, list(emptyObj_B))
                       meanPosts_Meta_bound <- c(meanPosts_Meta_bound, list(emptyObj_B))
                       sdevPosts_Meta_nondt <- c(sdevPosts_Meta_nondt, list(emptyObj_B))      # Posterior variance
                       sdevPosts_Meta_drift <- c(sdevPosts_Meta_drift, list(emptyObj_B))
                       sdevPosts_Meta_bound <- c(sdevPosts_Meta_bound, list(emptyObj_B))
                       rhats_Meta_nondt <- c(rhats_Meta_nondt, list(emptyObj_R_B))            # Rhats
                       rhats_Meta_drift <- c(rhats_Meta_drift, list(emptyObj_R_B))
                       rhats_Meta_bound <- c(rhats_Meta_bound, list(emptyObj_R_B))             
           }
       }
       
       # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       # Fill the arrays corresponding to each level 'p'
       #################################################################################################
       for(k in 1:nDatasets){
           j <- 1
           for(t in allT){
               if(Hierarchical){
               thisH <- which(outH[[k]][,"p"]==p&outH[[k]][,"t"]==t)
                        trueVals_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"true.values"])
                        meanPosts_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"mean.estimates"])
                        sdevPosts_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"std.estimates"])
                        rhats_Hier[[i]][k,,j] <- unlist(outH[[k]][thisH,"rhats"])
                        clock_Hier[k,j,i] <- as.numeric(outH[[k]][thisH,"elapsed.time"])
               }
               if(Meta){
                       # Meta-regression - Criterion: nondt
                       thisMn <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="nondt"&outB[[k]][,"d"]=="metaregression")
                       if(length(thisMn)>0){
                         trueVals_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"true.values"])
                         meanPosts_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"mean.estimates"])
                         sdevPosts_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"std.estimates"])
                         rhats_Meta_nondt[[i]][k,,j] <- unlist(outB[[k]][thisMn,"rhats"])
                         clock_Meta[["nondt"]][k,j,i] <- as.numeric(outB[[k]][thisMn,"elapsed.time"])
                       }
                       # Meta-regression - Criterion: drift
                       thisMd <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="drift"&outB[[k]][,"d"]=="metaregression")
                       if(length(thisMd)>0){
                         trueVals_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"true.values"])
                         meanPosts_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"mean.estimates"])
                         sdevPosts_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"std.estimates"])
                         rhats_Meta_drift[[i]][k,,j] <- unlist(outB[[k]][thisMd,"rhats"])
                         clock_Meta[["drift"]][k,j,i] <- as.numeric(outB[[k]][thisMd,"elapsed.time"])
                       }
                       # Meta-regression - Criterion: bound
                       thisMb <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="bound"&outB[[k]][,"d"]=="metaregression")
                       if(length(thisMb)>0){
                         trueVals_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"true.values"])
                         meanPosts_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"mean.estimates"])
                         sdevPosts_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"std.estimates"])
                         rhats_Meta_bound[[i]][k,,j] <- unlist(outB[[k]][thisMb,"rhats"])
                         clock_Meta[["bound"]][k,j,i] <- as.numeric(outB[[k]][thisMb,"elapsed.time"])
                       }
               }
               if(Ttest){
                       # ttest - Criterion: nondt
                       thisTn <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="nondt"&outB[[k]][,"d"]=="ttest")
                       if(length(thisTn)>0){
                         trueVals_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"true.values"])
                         meanPosts_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"mean.estimates"])
                         sdevPosts_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"std.estimates"])
                         rhats_Ttst_nondt[[i]][k,,j] <- unlist(outB[[k]][thisTn,"rhats"])
                         clock_Ttst[["nondt"]][k,j,i] <- as.numeric(outB[[k]][thisTn,"elapsed.time"])
                       }
                       # ttest - Criterion: drift
                       thisTd <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="drift"&outB[[k]][,"d"]=="ttest")
                       if(length(thisTd)>0){
                         trueVals_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"true.values"])
                         meanPosts_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"mean.estimates"])
                         sdevPosts_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"std.estimates"])
                         rhats_Ttst_drift[[i]][k,,j] <- unlist(outB[[k]][thisTd,"rhats"])
                         clock_Ttst[["drift"]][k,j,i] <- as.numeric(outB[[k]][thisTd,"elapsed.time"])
                       }
                       # ttest - Criterion: bound
                       thisTb <- which(outB[[k]][,"p"]==p&outB[[k]][,"t"]==t&outB[[k]][,"c"]=="bound"&outB[[k]][,"d"]=="ttest")
                       if(length(thisTb)>0){
                         trueVals_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"true.values"])
                         meanPosts_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"mean.estimates"])
                         sdevPosts_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"std.estimates"])
                         rhats_Ttst_bound[[i]][k,,j] <- unlist(outB[[k]][thisTb,"rhats"])
                         clock_Ttst[["bound"]][k,j,i] <- as.numeric(outB[[k]][thisTb,"elapsed.time"])
                       }
                }
           j <- j+1 
           }
       }
       i <- i + 1
   }

   ################################################################################
   # Name each array by the number-of-participants included
   ################################################################################
   if(Hierarchical){
                      names(trueVals_Hier) <- paste("P",allP,sep="") 
                      names(meanPosts_Hier) <- paste("P",allP,sep="")
                      names(sdevPosts_Hier) <- paste("P",allP,sep="")
                      names(rhats_Hier) <- paste("P",allP,sep="")
   }
   if(Ttest){
                       names(trueVals_Ttst_nondt) <- paste("P",allP,sep="")     # True values
                       names(trueVals_Ttst_drift) <- paste("P",allP,sep="")
                       names(trueVals_Ttst_bound) <- paste("P",allP,sep="")
                       names(meanPosts_Ttst_nondt) <- paste("P",allP,sep="")    # Mean Posteriors
                       names(meanPosts_Ttst_drift) <- paste("P",allP,sep="")
                       names(meanPosts_Ttst_bound) <- paste("P",allP,sep="")
                       names(sdevPosts_Ttst_nondt) <- paste("P",allP,sep="")    # Posterior variance
                       names(sdevPosts_Ttst_drift) <- paste("P",allP,sep="")
                       names(sdevPosts_Ttst_bound) <- paste("P",allP,sep="")
                       names(rhats_Ttst_nondt) <- paste("P",allP,sep="")        # Rhats
                       names(rhats_Ttst_drift) <- paste("P",allP,sep="")
                       names(rhats_Ttst_bound) <- paste("P",allP,sep="")
   }
   if(Meta){
                       names(trueVals_Meta_nondt) <- paste("P",allP,sep="")     # True values
                       names(trueVals_Meta_drift) <- paste("P",allP,sep="")
                       names(trueVals_Meta_bound) <- paste("P",allP,sep="")
                       names(meanPosts_Meta_nondt) <- paste("P",allP,sep="")    # Mean Posteriors
                       names(meanPosts_Meta_drift) <- paste("P",allP,sep="")   
                       names(meanPosts_Meta_bound) <- paste("P",allP,sep="")
                       names(sdevPosts_Meta_nondt) <- paste("P",allP,sep="")    # Posterior variance
                       names(sdevPosts_Meta_drift) <- paste("P",allP,sep="")
                       names(sdevPosts_Meta_bound) <- paste("P",allP,sep="")
                       names(rhats_Meta_nondt) <- paste("P",allP,sep="")        # Rhats
                       names(rhats_Meta_drift) <- paste("P",allP,sep="")
                       names(rhats_Meta_bound) <- paste("P",allP,sep="")
   }
   
   ################################################################################
   # FINAL: Store the output per design-parameter by writing an .RData object
   ################################################################################
   if(Hierarchical){
                     simStudy_Hierarchical <- list("true" = trueVals_Hier,
                                                   "recovered" = meanPosts_Hier,
                                                   "estimates_sdev" = sdevPosts_Hier,
                                                   "rhats" = rhats_Hier)
                     outputFile = paste(saveTo,"/simStudy_Hierarchical.RData",sep="")
                     save(simStudy_Hierarchical, file=outputFile)
   }
   if(Meta){
            if("nondt" %in% settings$criterion_levels){
                     simStudy_Meta_nondt <- list("true" = trueVals_Meta_nondt, 
                                                 "recovered" = meanPosts_Meta_nondt,
                                                 "estimates_sdev" = sdevPosts_Meta_nondt,
                                                 "rhats" = rhats_Meta_nondt)
                     save(simStudy_Meta_nondt, file=paste(saveTo,"/sim_",simStudyName,"_Meta_nondt.RData",sep=""))          }
            if("drift" %in% settings$criterion_levels){
                     simStudy_Meta_drift <- list("true" = trueVals_Meta_drift,
                                                 "recovered" = meanPosts_Meta_drift,
                                                 "estimates_sdev" = sdevPosts_Meta_drift,
                                                 "rhats" = rhats_Meta_drift)
                     save(simStudy_Meta_drift, file=paste(saveTo,"/sim_",simStudyName,"_Meta_drift.RData",sep=""))           }
            if("bound" %in% settings$criterion_levels){
                     simStudy_Meta_bound <- list("true" = trueVals_Meta_bound,
                                                 "recovered" = meanPosts_Meta_bound,
                                                 "estimates_sdev" = sdevPosts_Meta_bound,
                                                 "rhats" = rhats_Meta_bound)
                     save(simStudy_Meta_bound, file=paste(saveTo,"/sim_",simStudyName,"_Meta_bound.RData",sep=""))           }
   }
   if(Ttest){
                     simStudy_Ttst_nondt <- list("true" = trueVals_Ttst_nondt,
                                                 "recovered" = meanPosts_Ttst_nondt,
                                                 "estimates_sdev" = sdevPosts_Ttst_nondt,
                                                 "rhats" = rhats_Ttst_nondt)
                     save(simStudy_Ttst_nondt, file=paste(saveTo,"/sim_",simStudyName,"_Ttest_nondt.RData", sep=""))
                     
                     simStudy_Ttest_drift <- list("true" = trueVals_Ttst_drift,
                                                  "recovered" = meanPosts_Ttst_drift,
                                                  "estimates_sdev" = sdevPosts_Ttst_drift,
                                                  "rhats" = rhats_Ttst_drift)
                     save(simStudy_Ttest_drift, file=paste(saveTo,"/sim_",simStudyName,"_Ttest_drift.RData",sep=""))
                     
                     simStudy_Ttest_bound <- list("true" = trueVals_Ttst_bound, 
                                                  "recovered" = meanPosts_Ttst_bound,
                                                  "estimates_sdev" = sdevPosts_Ttst_bound,
                                                  "rhats" = rhats_Ttst_bound)
                     save(simStudy_Ttest_bound, file=paste(saveTo,"/sim_",simStudyName,"_Ttest_bound.RData",sep=""))
   }
} 