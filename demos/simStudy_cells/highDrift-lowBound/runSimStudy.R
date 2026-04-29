##########################################################
# LOAD FUNCTIONS/PACKAGES
##########################################################
######## Load required R packages for parallel processing
# Set working directory to repo root (this script is in demos/simulation-study/)
library(here)
library(foreach)
library(doParallel)

# Call the function within the src directory
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

##########################################################
# START OUTPUT DIRECTORIES IF NEEDED
##########################################################
# A directory to store the simulation results
output_dir <- here("demos", "simStudy_cells", "highDrift-lowBound", "samples", "samples-T40-fullP")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# A directory to store the simulation logs
log_dir <- here("demos", "simStudy_cells", "highDrift-lowBound", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

##########################################################
# LOAD DEFAULT SIMULATION SETTINGS
##########################################################
source(here("demos", "simStudy_setup", "simulation_settings.R"))

##########################################################
# SET CUSTOM SIMULATION SETTINGS
##########################################################
settings$participant_levels <- 160
#settings$trial_levels <- c(40,160)
settings$trial_levels <- 40
settings$true_means$drift_mean <- c(2, 3)
settings$true_means$bound_mean <- c(2, 2.5)

settings$jagsParameters <- rbind("bound_mean", "drift_mean", "nondt_mean", "betaweight",
                                 "bound", "drift", "nondt")

################################################################
# Run simulation study
################################################################
cores       <-  detectCores()
my.cluster  <-  makeCluster(cores[1]-4)

registerDoParallel(cl = my.cluster)
# Then modify your foreach call to use this combine function
resultado <- foreach(seed = 1:1000, 
                     .errorhandling = "pass",
                     .combine = combine_results
) %dopar% {
  Z <- simStudy_runFullSeed(seed = seed,
                            settings = settings,
                            forceRun = forceRun,
                            include_EZ_Robust = settings$include_EZ_Robust,
                            redo_if_bad_rhat = TRUE,
                            rhat_cutoff = 1.05,
                            prevent_zero_accuracy = FALSE,
                            Show = TRUE)
}
stopCluster(cl = my.cluster)