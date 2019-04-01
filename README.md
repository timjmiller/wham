### WHAM: a state-space age-structured assessment model

This repository has R and TMB code for a general model with many similarities to those in [Miller et al. 2016](https://doi.org/10.1139/cjfas-2015-0339), [Miller and Hyun 2018](https://doi.org/10.1139/cjfas-2017-0035), and [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124).
The model also has many similarities to ASAP [Legault and Restreppo]. Indeed, there are R functions provided to generate a WHAM input file (with some default settings) from an ASAP3 dat file, and many of the plotting functions for input data, results, and diagnostics are borrowed and modified from [ASAP](https://github.com/cmlegault/ASAPplots).
The WHAM model can be configured to estimate a range of assessment models from a traditional statistical catch-at-age model with recruitments as (possibly penalized) fixed effects, SCAA with recruitments as random effects, or abundance at all ages treated as random effects.

### Selectivity
* A number of selectivity "blocks" is specified and a pointer to each fleet and index each year is specified.
* Selectivity models for each block can be composed of age-specific parameters, logistic (increasing or decreasing), or double-logistic. 

### Observation model components:
## Aggregate relative abundance indices:
* These can be in numbers or biomass, determined by an indicator.
* A normal distribution is assumed for the log-indices. 
* Input standard errors for log-indices are year-specific, but index-specific scalars may be estimated.
* There can be any number of indices.
* There are indicators for whether to use the index in a given year.
* The time of year the index occurs is specified for each year.
* The selectivity "block" to use for an index is specified for each year.
  
## Index age composition:
* These can also be in numbers or biomass, determined by an indicator.
* Several distribution assumptions available: multinomial, Dirichlet, Dirichlet-multinomial, logistic normal, two types of zero-one inflated logistic normal.
* Effective sample sizes for multinomial are input and year-specific, but variance parameters may be estimated for other distributions.
* The order is the same as the relative abundance indices.
* There are indicators for whether to use the age composition in a given year.
  
## Aggregate fleet catch:
* This is in biomass.
* A normal distribution is assumed for the log-catch. 
* Input standard errors for log-catch are year-specific, but fleet-specific scalars may be estimated.
* There can be any number of fleets.
* The selectivity "block" to use for an fleet is specified for each year.
  
## Fleet age composition:
* These can only be in numbers.
* Same distribution assumptions as indices available: multinomial, Dirichlet, Dirichlet-multinomial, logistic normal, two types of zero-one inflated logistic normal.
* Effective sample sizes for multinomial are input and year-specific, but variance parameters may be estimated for other distributions.
* The order is the same as the fleet catches.
* There are indicators for whether to use the age composition in a given year.

### Inputs:

## Weight-at-age:
* Like ASAP weight-at-age are input and pointers are specified for them for fleet catches, relative abundance indices (if necessary), Spawning Stock Biomass, etc.
* Future versions may include options for internal estimation of weight at age, like [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124).

### Biological Reference Points:
* The user specifies a percentage of unfished spawning biomass per recruit to make estimates of fishing mortality and SSB reference points.
* Beverton-Holt and Ricker spawner-recruit models may be assumed. Traditional "alpha" and "beta" parameters or steepness and $$R_0$$ may be estimated.
* Annual conditional reference points are estimated in case weight, maturity, natural mortality, or selectivity at age change over time. 
* If a spawner-recruit model is assumed, MSY-based reference points are also estimated. Both of these classes of reference points employ a Newton method internally to determine the fishing mortality reference points, thereby propogating uncertainty of inputs as in [Miller et al. 2016,](https://doi.org/10.1139/cjfas-2015-0339) and [Miller et al. 2018,](https://doi.org/10.1139/cjfas-2017-0124)

### Parameters:

## Natural Mortality:
* In general, this can be year- and age-specific. However, by default it is assumed known at the initial parameter values using the map argument in TMB.
* It can also be treated as a random walk of random effects, or an allometric function of weight at age.
  
