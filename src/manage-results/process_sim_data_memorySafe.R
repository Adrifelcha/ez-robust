process_sim_data_by_cell <- function(seed_dir, output_dir) {
  # List all seed-specific .RData files
  all_seed_files <- list.files(seed_dir, pattern = "\\.RData$", full.names = TRUE)
  n_seeds <- length(all_seed_files)
  cat("Found", n_seeds, "seed-result files to process.\n")
  
  # --- Get metadata from the first seed file ---
  cat("Extracting metadata from the first seed-result file...\n")
  temp_env <- new.env()
  load(all_seed_files[1], envir = temp_env)
  sample_output <- get("output", envir = temp_env)
  
  # Extract settings and design space
  settings <- sample_output$settings
  participant_levels <- settings$participant_levels
  trial_levels <- settings$trial_levels
  param_names <- c(settings$jagsParameters)
  
  # Identify effects (i.e., noEffect vs fixedEffect)
  available_effects <- setdiff(names(sample_output), c("reps", "settings"))
  # Identify conditions (i.e., Model by Data type; e.g., "EZ_clean")
  available_conditions <- colnames(sample_output[[available_effects[1]]])
  
  cat("Found", length(available_conditions), "simulation conditions:\n", 
      paste(1:length(available_conditions),") ", available_conditions, collapse=",\n ", sep=""), "\n")
  rm(temp_env, sample_output) # Clean up
  
  
  # Loop over each model/data condition (e.g., "EZ_clean")
  for(condition in available_conditions){
      cat(paste0("\n--- Processing Condition: ", condition, " ---\n"))

      target_folder <- here(output_dir, condition)
      dir.create(target_folder, recursive = TRUE, showWarnings = FALSE)
    
      for(p in participant_levels){         # For every participant level
          for(t in trial_levels){           # And for every trial level
              cat(paste0("  Processing cell: P = ", p, ", T = ", t, "\n"))
              # Initialize progress bar to be printed to the console
              progress_bar <- txtProgressBar(min = 0, max = n_seeds,
                                            style = 3, width = 50, char = "=")
              
              # Empty objects
              rhats_matrix <- c()
              true_values_matrix <- c()
              mean_estimates_matrix <- c()
              std_estimates_matrix <- c()
              credInterval_matrix <- c()                                          
              beta_chains_matrix <- c()
              summary_matrix <- c()
              seed_names <- c()

              # Loop through every single seed file to find data for this cell
              for(i in seq_along(all_seed_files)){
                  seed_file <- all_seed_files[i]
                  temp_env_seed <- new.env()
                  
                  load_successful <- FALSE
                  tryCatch({
                    load(seed_file, envir = temp_env_seed)
                    load_successful <- TRUE
                  }, error = function(e) {
                    cat(paste("\nCould not load file:", seed_file, "- skipping.\n"))
                  })
                  
                  if(load_successful){
                    seed_output <- get("output", envir = temp_env_seed)
                  
                    # We check the dataframes for all effects
                    effects_in_seed <- setdiff(names(seed_output), c("reps", "settings"))
                    for(effect in effects_in_seed){

                        effect_df <- seed_output[[effect]]
                        # Verify that we have results for this condition and effect
                        
                        if(condition %in% colnames(effect_df)) {
                          these_results <- effect_df[,condition]
                          
                          if (!is.null(these_results) && length(these_results) > 0) {
                            # The actual list of all PxT results is the first (and only) element.
                            matching_results <- Filter(function(x) x$p == p && x$t == t, these_results)

                            # Extract rhats
                            these_rhats <- t(sapply(matching_results, function(x) x$rhats))              
                            rhats_matrix <- rbind(rhats_matrix, these_rhats[,param_names])
                            # Extract true values
                            these_truevals <- t(sapply(matching_results, function(x) x$true.values))              
                            true_values_matrix <- rbind(true_values_matrix, these_truevals[,param_names])
                            # Extract mean estimates
                            these_estimates <- t(sapply(matching_results, function(x) x$mean.estimates))              
                            mean_estimates_matrix <- rbind(mean_estimates_matrix, these_estimates[,param_names])
                            # Extract std estimates
                            these_std <- t(sapply(matching_results, function(x) x$std.estimates))              
                            std_estimates_matrix <- rbind(std_estimates_matrix, these_std[,param_names])
                            # Extract credible intervals
                            these_CI <- t(sapply(matching_results, function(x) x$credInterval))
                            credInterval_matrix <- rbind(credInterval_matrix, these_CI[,param_names])
                            # Extract beta chains
                            getbetas <- t(sapply(matching_results, function(x) x$beta_chains))
                            if(nrow(getbetas) == 1){
                               beta_chains_matrix <- cbind(beta_chains_matrix, getbetas)
                            }else{
                               list_of_chains <- split(getbetas, 1:nrow(getbetas))
                              for(fila in 1:nrow(getbetas)){
                                list_of_chains[[fila]] <- matrix(list_of_chains[[fila]], ncol = seed_output$settings$n.chains)
                              }                            
                              beta_chains_row <- matrix(list_of_chains, nrow = 1)                                  
                            }

                            current_summary <- data.frame("seed" = sapply(matching_results, function(x) x$seed),
                                                        "jagsTime" = c(t(sapply(matching_results, function(x) x$jagsTime))),
                                                        "nIter" = c(t(sapply(matching_results, function(x) x$nIter))),
                                                        "nBurnin" = c(t(sapply(matching_results, function(x) x$nBurnin))),
                                                        "nThin" = c(t(sapply(matching_results, function(x) x$nThin))),
                                                        "bad_rhats" = c(t(sapply(matching_results, function(x) x$bad_rhat_count))))
                            summary_matrix <- rbind(summary_matrix, current_summary)
                          }
                        }
                    }
                    rm(temp_env_seed, seed_output) # Clean up memory
                  } else {
                    rm(temp_env_seed)
                  }

                  # Update the progress bar
                  setTxtProgressBar(progress_bar, i)
              }
              
              # Close the progress bar
              close(progress_bar)
                            
              cat("    -> Found", nrow(rhats_matrix), "matching results across seeds. Collating and saving...\n")
                            
              simStudy_Beta <- list("true" = true_values_matrix,
                                    "estimates" = mean_estimates_matrix,
                                    "error.sd" = std_estimates_matrix,
                                    "rhats" = rhats_matrix,
                                    "credInterval" = credInterval_matrix,
                                    "beta_chains" = beta_chains_matrix,                                    
                                    "summary" = summary_matrix)
            
              subCell_ID <- name_subCell_ID(subCell_name = condition)
              outputFile <- name_simStudyCell(nTrials = t, nParticipants = p, nDatasets = n_seeds, 
                                              data_type = subCell_ID$data_type, 
                                              ez_model = subCell_ID$EZ_model, output.folder = target_folder)
            
              save(simStudy_Beta, file = outputFile)
              
          } # end trial loop
    } # end participant loop
  } # end condition loop
  
  cat("\nProcessing complete. All cell-specific files have been generated in:", output_dir, "\n")
}