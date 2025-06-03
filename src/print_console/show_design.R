#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function displays the simulation design parameters to the console
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
show_design = function(local_settings, withinSubjectDesign = FALSE){
      
      # Print header
      cat("========== EZBHDDM Design Parameters: ==========\n")

      # Print model information based on design type
      # Case 1: Within-subject t-test on drift rate
      if(withinSubjectDesign){
            cat("Model type:  Within-subject t-test on drift rate\n")
      } else {           
      # Case 2: Between-subjects designs
            # modelType is either "hierarchical", "metaregression", or "ttest"
            cat("Model type:              ", local_settings$modelType,"\n")
            cat("Testing for an effect on the drift rate\n")
      }
            
      # Display number of participants
      cat("Number of Participants:  ", local_settings$nPart,"\n")
      
      # Display trial information
      if(withinSubjectDesign){
            # Within-subject design
            cat("Trials Per Condition:    ", local_settings$nTrialsPerCondition,"\n")
            cat("Number of Conditions:    ", 2,"\n")
      } else {
            # Between-subjects design
            cat("Trials Per Person:       ", local_settings$nTrials,"\n")
      }
}
