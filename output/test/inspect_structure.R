library(here)
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

# Define directories
seed_dir <- here("demos", "simulation-study", "samplesSummStats")

# Load the first seed file
all_seed_files <- list.files(seed_dir, pattern = "\\.RData$", full.names = TRUE)[1]
cat("Loading from:", all_seed_files, "\n\n")
temp_env <- new.env()
load(all_seed_files, envir = temp_env)
seed_output <- get("output", envir = temp_env)

# Get first effect and first condition
effects_in_seed <- setdiff(names(seed_output), c("reps", "settings"))
effect <- effects_in_seed[1]
effect_df <- seed_output[[effect]]

conditions <- colnames(effect_df)
condition <- conditions[1]

cat("=== Inspecting effect_df structure ===\n")
cat("Effect:", effect, "\n")
cat("Condition:", condition, "\n")
cat("Class of effect_df:", class(effect_df), "\n")
cat("Dimensions of effect_df:", dim(effect_df), "\n\n")

# Get the results for this condition
these_results <- effect_df[,condition]

cat("=== Inspecting these_results ===\n")
cat("Class of these_results:", class(these_results), "\n")
cat("Length of these_results:", length(these_results), "\n")
cat("Is list?", is.list(these_results), "\n\n")

if (is.list(these_results) && length(these_results) > 0) {
  # Check first element
  first_result <- these_results[[1]]
  
  cat("=== Inspecting first result element ===\n")
  cat("Names:", names(first_result), "\n\n")
  
  # Check p and t
  cat("=== Inspecting p and t ===\n")
  cat("x$p:\n")
  print(first_result$p)
  cat("Class:", class(first_result$p), "\n")
  cat("Length:", length(first_result$p), "\n")
  cat("Type:", typeof(first_result$p), "\n")
  cat("Is atomic?", is.atomic(first_result$p), "\n\n")
  
  cat("x$t:\n")
  print(first_result$t)
  cat("Class:", class(first_result$t), "\n")
  cat("Length:", length(first_result$t), "\n")
  cat("Type:", typeof(first_result$t), "\n")
  cat("Is atomic?", is.atomic(first_result$t), "\n\n")
  
  # Check summStats structure
  cat("=== Inspecting summStats structure ===\n")
  cat("x$summStats:\n")
  print(str(first_result$summStats))
  cat("\n")
  
  # Check true.values structure
  cat("=== Inspecting true.values structure ===\n")
  cat("x$true.values:\n")
  print(str(first_result$true.values))
  cat("\n")
  
  # Check mean.estimates structure
  cat("=== Inspecting mean.estimates structure ===\n")
  cat("x$mean.estimates:\n")
  print(str(first_result$mean.estimates))
  cat("\n")
  
  # Try to understand what happens when we try to compare
  cat("=== Testing comparison ===\n")
  p <- 20  # Example value
  t <- 20  # Example value
  
  cat("Comparing first_result$p == p:\n")
  tryCatch({
    result <- first_result$p == p
    cat("Result:", result, "\n")
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
  })
  
  cat("\nComparing first_result$t == t:\n")
  tryCatch({
    result <- first_result$t == t
    cat("Result:", result, "\n")
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
  })
  
  # Test Filter function with these_results
  cat("\n=== Testing Filter function ===\n")
  cat("Testing with p=20, t=20:\n")
  tryCatch({
    filtered <- Filter(function(x) x$p == p && x$t == t, these_results)
    cat("Filtered results count:", length(filtered), "\n")
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
    cat("Error class:", class(e), "\n")
  })
  
  # Check if all results have summStats
  cat("\n=== Checking summStats across all results ===\n")
  summStats_count <- sum(sapply(these_results, function(x) !is.null(x$summStats)))
  cat("Results with non-NULL summStats:", summStats_count, "out of", length(these_results), "\n")
  
  # Find one with summStats if it exists
  if (summStats_count > 0) {
    result_with_stats <- these_results[[which(sapply(these_results, function(x) !is.null(x$summStats)))[1]]]
    cat("\nFound result with summStats. Structure:\n")
    print(str(result_with_stats$summStats))
    
    # Show how it compares to true.values
    cat("\n=== Comparing summStats structure to true.values ===\n")
    cat("summStats is a list:", is.list(result_with_stats$summStats), "\n")
    cat("true.values is a list:", is.list(result_with_stats$true.values), "\n")
    cat("summStats names:", names(result_with_stats$summStats), "\n")
    cat("true.values names:", names(result_with_stats$true.values), "\n")
  }
}
