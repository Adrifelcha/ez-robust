
rm(list = ls())
library(here)

# Call the function within the src directory
source(here("src", "load_allFunctions.R"))
load_allCustomFunctions()
directory <- here("demos", "simulation-study", "samples")
object_name = "output"

x <- load_seedOutput(directory, object_name)


store_simStudyResults(output = x, saveTo = here("output", "simStudy_results"))


rm(list = ls())
library(here)
load(here("demos", "simulation-study", "samples", "seed-20.RData"))
ls()
names(output)
