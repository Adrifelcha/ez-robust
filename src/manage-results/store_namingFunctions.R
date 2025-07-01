#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function generates standardized filenames for simulation output files
# Note: This function ISN'T made to work for within-subject designs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
name_simStudyCell <- function(nTrials, nParticipants, nDatasets, data_type, ez_model, output.folder = NA) {
    # Set default output folder using here() if not specified
    if(is.na(output.folder)) {
        output.folder <- here::here("output", "simStudy_results")
    }
    
    # Ensure output folder exists
    if(!dir.exists(output.folder)){
        dir.create(output.folder, recursive = TRUE)
    }
    
    # Create a compact base filename
    filename <- paste("sim_P", nParticipants, "T", nTrials, "_", ez_model, "-", data_type,".RData", sep="")
            
    # Combine folder and filename, ensuring proper path separator
    if(substr(output.folder, nchar(output.folder), nchar(output.folder)) != "/" && 
       substr(output.folder, nchar(output.folder), nchar(output.folder)) != "\\") {
        output.folder <- paste0(output.folder, "/")
    }
    outputFile <- paste0(output.folder, filename)
    
    return(outputFile)
}


name_subCell_ID <- function(subCell_name){

    # Split the colname string at the underscore
    split_result <- strsplit(subCell_name, "_")

    if(split_result[[1]][2] == "clean"){
        data_type <- "Clean"
    }else{
        if(split_result[[1]][2] == "contaminated"){
            data_type <- "Outliers"
        }else{
            data_type <- split_result[[1]][2]
        }
    }

    if(split_result[[1]][1] == "EZ"){
        EZ_model <- "ez"
    }else{
        if(split_result[[1]][1] == "EZRobust"){
            EZ_model <- "ezRob"
        }else{
            EZ_model <- split_result[[1]][1]
        }
    }

    return(list("EZ_model" = EZ_model, "data_type" = data_type))
}