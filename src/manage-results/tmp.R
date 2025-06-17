library(here)
directory <- here("demos", "simulation-study", "samples")

x <- load_seedOutput(directory)

length(x)

names(x)

x$reps

length(x$settings)

rm(list = ls())
library(here)
load(here("demos", "simulation-study", "samples", "seed-1.RData"))
ls()
names(output)


rm(list = ls())
library(here)
load(here("demos", "simulation-study", "samples", "seed-2.RData"))
ls()
names(output)


