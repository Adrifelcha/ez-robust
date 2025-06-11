# Do you want to force re-run the simulations?
forceRun = TRUE


##########################################################
# LOAD FUNCTIONS/PACKAGES
##########################################################
######## Load required R packages for parallel processing
library(here)
library(foreach)
library(doParallel)

# Call the function within the src directory
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()

##########################################################
# SIMULATION SETTINGS
##########################################################
# Create output directory if it doesn't exist
output_dir <- here("demos", "simStudy_EZ-Clean", "samples")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Fixed simulation design variables
settings <- list("output.folder" = file.path(output_dir, "/"),
                 "participant_levels" = c(20,40,80,160,320),
                 "trial_levels" = c(20,40,80,160,320),
                 "beta_levels" = c(0, 0.1, 0.2, 0.4, 0.8),
                 "separate_datasets" = TRUE,
                 "contaminant_prob" = 0.05,
                 "include_EZ_Robust" = TRUE,
                 "nDatasets" = 1000,
                 "modelType" = "ttest",
                 "withinSubject" = TRUE,
                 "n.chains" = 4,
                 "n.burnin" = 500,
                 "n.iter" = 4000,
                 "n.thin" = 1,
                 "jagsParameters" = rbind("bound_mean", "drift_mean", "nondt_mean", "betaweight"),
                 "modelFile" = rbind(here("output", "BUGS-models", "EZHBDDM-withinSubject.bug"),
                                     here("output", "BUGS-models", "EZHBDDM-withinSubject-Robust.bug")))

# Implied number of cells
settings <- c(settings,
              list("nCells" = prod(length(settings$participant_levels),length(settings$trial_levels), length(settings$beta_levels)),
                   "jagsData" = list(list("EZ" = JAGS_passData(settings$modelType, withinSubject = settings$withinSubject),
                                          "EZRobust" = JAGS_passData(settings$modelType, EZRobust = settings$include_EZ_Robust, withinSubject = settings$withinSubject)))))

# Prepare JAGS objects
jagsInits <- list()
for(i in settings$participant_levels){
    jagsInits <- c(jagsInits, list(JAGS_inits(n.chains = settings$n.chains, nParticipants = i, 
                   custom_sd = 0.25, withinSubject = settings$withinSubject)))
}

# Add JAGS objects to settings
settings <- c(settings, list("jagsInits" = jagsInits))
# Change names so they can be called more easily
colnames(settings$jagsParameters) <- settings$modelType
names(settings$jagsData) <- settings$modelType
names(settings$jagsInits) <- settings$participant_levels
rownames(settings$modelFile) <- c("EZ", "EZRobust")
colnames(settings$modelFile) <- settings$modelType


##########################################################
# Write JAGS models
##########################################################
# Prepare specific prior distribution parameters
custom_priors_list <- list(
                      "bound_mean_mean" = 2.5,    "bound_mean_sdev" = 2.00,
                      "drift_mean_mean" = 0.00,    "drift_mean_sdev" = 3.00,
                      "nondt_mean_mean" = 0.5,    "nondt_mean_sdev" = 0.20,
                      "bound_sdev_lower" = 0.01,   "bound_sdev_upper" = 2.00,
                      "drift_sdev_lower" = 0.01,   "drift_sdev_upper" = 2.00,
                      "nondt_sdev_lower" = 0.01,   "nondt_sdev_upper" = 0.50,
                      "betaweight_mean" = 0,       "betaweight_sdev" = 1)

# Add JAGS objects to settings
settings <- c(settings, list("priors" = list(JAGS_priors(Show=FALSE, "ttest", custom_prior_list = custom_priors_list))))
names(settings$priors) <- settings$modelType

# Define custom truncation list
custom_truncation_list <- list(
        "bound_mean" = c(0.1, ""), "nondt_mean" = c(0.05, ""), "drift_mean" = c("", ""),
        "bound_sdev" = c(0.01, ""), "nondt_sdev" = c(0.01, ""), "drift_sdev" = c(0.01, ""),
        "drift" = c("", ""), "bound" = c(0.0001, ""), "nondt" = c(0.0001, ""), "betaweight" = c("-3", "3"))

for(EZR in c(FALSE, TRUE)){
        modelFile <- settings$modelFile[EZR+1]
        JAGS_writeModel(priors = settings$priors[[settings$modelType]], modelType = settings$modelType, 
                        withinSubject = settings$withinSubject, EZRobust = EZR, modelFile = modelFile,
                        custom_truncation_list = custom_truncation_list)    
}

settings <- c(settings, list("true_sdevs" = list("bound_sdev" = 0.5, "nondt_sdev" = 0.1, "drift_sdev" = 0.75),
                             "true_means" = list("bound_mean" = c(2, 4), "nondt_mean" = c(0.2, 0.4), "drift_mean" = c(-3, 3))))

################################################################
# Define simulation functions
################################################################
cores       <-  detectCores()
my.cluster  <-  makeCluster(cores[1]-4)

registerDoParallel(cl = my.cluster)
resultado <- foreach(seed = 201:210, 
                  .errorhandling = "pass",
                  .combine = 'rbind'
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



resultado <- load_seedOutput(directory = here("demos", "simulation-studies", "generative_uniforms", "samples"),
                             object_name = "output")
settings$nDatasets <- nrow(resultado$results)
# Store the results to repo-root/output/RData-results
store_parallelOutput(output = resultado$results, saveTo = here("output", "RData-results"),
                     settings = settings, simStudyName = "genUniform")

source(here("src", "plotting", "plot_simStudyOutput.R"))
makeSimStudyPlot(here("output", "RData-results", "sim_genUniform_Meta_drift.RData"), plotType = 2)
