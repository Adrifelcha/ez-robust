#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function generates standardized filenames for simulation output files
# Note: This function ISN'T made to work for within-subject designs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
name_simStudyCell <- function(nTrials, nParticipants, nDatasets, modelType, beta_level = NA, output.folder = NA) {
    # Set default output folder using here() if not specified
    if(is.na(output.folder)) {
        output.folder <- here::here("output", "simStudy_results")
    }
    
    # Ensure output folder exists
    if(!dir.exists(output.folder)){
        dir.create(output.folder, recursive = TRUE)
    }
    
    beta <- gsub("\\.", "", as.character(beta_level))

    # Create a compact base filename
    filename <- paste("sim_P", nParticipants, "T", nTrials, "D", nDatasets, "_Beta-", beta, modelType, ".RData", sep="")
            
    # Combine folder and filename, ensuring proper path separator
    if(substr(output.folder, nchar(output.folder), nchar(output.folder)) != "/" && 
       substr(output.folder, nchar(output.folder), nchar(output.folder)) != "\\") {
        output.folder <- paste0(output.folder, "/")
    }
    outputFile <- paste0(output.folder, filename)
    
    return(outputFile)
}
