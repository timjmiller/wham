library(TMB)
dyn.load(dynlib("wham_v1"))
#Southern New England - Mid-Atlantic yellowtail flounder 1973-2016
load("wham_v1.RData")
source("other_wham_functions.r")

#full state-space model, abundance is the state vector, B-H S-R model estimated.
#age comp likelihoods are logit normal
mod = fit.wham.fn(input)
#gdbsource("wham_v1.r", TRUE)
save(mod, file = "mod.RData")
#rho.fn(mod)
#x = retro.res.fn(mod)
