##############################################
# Print parameter values to console
##############################################
show_parameters <- function(parameter_set){  
      cat("===== EZBHDDM True Parameters: ==============\n")
      
      # Print drift rate hierarchical parameters
      cat("Drift Mean:   ", parameter_set$drift_mean,"\n")
      cat("Drift SD:     ", parameter_set$drift_sdev,"\n")
      # Print decision boundary hierarchical parameters
      cat("Bound Mean:   ", parameter_set$bound_mean,"\n")
      cat("Bound SD:     ", parameter_set$bound_sdev,"\n")
      # Print non-decision time hierarchical parameters
      cat("Non-decision Time Mean:", parameter_set$nondt_mean,"\n")
      cat("Non-decision Time SD:  ", parameter_set$nondt_sdev,"\n")
      # Print betaweight parameter if it exists
      if(!is.null(parameter_set$betaweight)){
        cat("Betaweight:   ", parameter_set$betaweight,"\n")
      }
      cat("Individual drift range: ", min(parameter_set$drift), " to ", max(parameter_set$drift), "\n")
      cat("Individual non-decision time range: ", min(parameter_set$nondt), " to ", max(parameter_set$nondt), "\n")
      cat("Negative non-decision time parameters: ", sum(parameter_set$nondt < 0), "\n")
      cat("Individual boundary range: ", min(parameter_set$bound), " to ", max(parameter_set$bound), "\n")
      cat("Negative boundary parameters: ", sum(parameter_set$bound < 0), "\n")
      cat("=============================================\n")
}