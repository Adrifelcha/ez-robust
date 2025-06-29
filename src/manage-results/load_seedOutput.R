#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
# This function loads the seed-specific .RData files stored in the /samples directory
# It combines them into a single list object, with a structure identical to the one
# produced by running the simulation in a single `foreach` loop with our custom combine function.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
load_seedOutput <- function(directory = NA, object_name = "output") {  
    # seed-specific results are stored in .RData files
    pattern <- ".RData$"
    # List all files matching the pattern in the directory
    files <- list.files(directory, pattern = pattern, full.names = TRUE)
    
    # First check if there are any files to load
    if(length(files) == 0) {
      stop("No '", pattern, "' files found in directory: ", directory)    
    }
    
    # Store all results in a temporary list
    all_seed_results <- list()
    valid_count <- 0
    for(i in seq_along(files)) {
          # Load the file into a temporary environment to avoid namespace conflicts
          archive <- files[i]
          cat("Loading file:", basename(archive), "\n")        
          temp_env <- new.env()
          load(archive, envir = temp_env)
          
          # Check if the loaded file contains the expected object
          if(!exists(object_name, envir = temp_env)) {
            warning("RData file ", archive, " does not contain a '", object_name, "' object")
            next  # Skip this file
          }
          
          # Get the object from the environment and add to list
          valid_count <- valid_count + 1
          all_seed_results[[valid_count]] <- get(object_name, envir = temp_env)
    }
  
    # Get the number of valid seeds
    n_seeds <- length(all_seed_results)

    # Find the first valid result to use as a reference for structure
    detect <- NULL
    for (res in all_seed_results) {
      # A valid result will have more than "reps" and "settings"
      if (!is.null(res) && length(res) > 2) {
        detect <- res
        break
      }
    }

    # If no comprehensive result was found, fall back to the first loaded result.
    if(is.null(detect)) {
      warning("Could not find any seed results with data to combine. Returning only settings and reps if available.")
      detect <- all_seed_results[[1]]
    }

    # Extract settings and find the names of the result components
    settings <- detect$settings
    available_results <- setdiff(names(detect), c("reps", "settings"))
    
    # Initialize the final combined output list
    final_output <- list(settings = settings)

    # Combine 'reps' data from all seeds
    all_reps <- lapply(all_seed_results, function(x) if(!is.null(x$reps)) x$reps else NULL)
    all_reps <- all_reps[!sapply(all_reps, is.null)] # Filter out NULLs
    if(length(all_reps) > 0) {
      final_output$reps <- do.call(rbind, all_reps)
    }

    # Combine each available result component from all seeds
    for(res_name in available_results) {    
      component_list <- lapply(all_seed_results, function(x) if(!is.null(x[[res_name]])) x[[res_name]] else NULL)
      component_list <- component_list[!sapply(component_list, is.null)] # Filter out NULLs    
      if(length(component_list) > 0) {
        final_output[[res_name]] <- do.call(rbind, component_list)
      }
    }

    cat("Successfully combined", n_seeds, "seed-specific output files into a single list object.\n")
    
    return(final_output)
}
##########################