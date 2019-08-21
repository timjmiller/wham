## WHAM: a state-space age-structured assessment model

The Woods Hole Assessment Model (WHAM) is a state-space age-structured stock assessment model that can include environmental effects on population processes. WHAM can be configured to estimate a range of assessment models from a traditional statistical catch-at-age (SCAA) model with recruitments as (possibly penalized) fixed effects, SCAA with recruitments as random effects, or abundance at all ages treated as random effects. WHAM is a generalization of the R and TMB code from [Miller et al. (2016)](https://doi.org/10.1139/cjfas-2015-0339), [Miller and Hyun 2018](https://doi.org/10.1139/cjfas-2017-0035), and [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124). WHAM has many similarities to ASAP ([code](https://www.nefsc.noaa.gov/nft/ASAP.html), [Legault and Restrepo 1999](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2007/WGMHSA/Annex%203%20-%20ICCAT%20Working%20Document.pdf)), including the input data file structure. Many of the plotting functions for input data, results, and diagnostics are modified from ASAP code written by Chris Legault and Liz Brooks ([ASAPplots](https://github.com/cmlegault/ASAPplots)).

### Installation and basic use

Install WHAM:
```
# if you want vignettes
# not sure why, but these look better if you *build using R Studio*
devtools::install_github("timjmiller/wham", dependencies=TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))

# if you do NOT want vignettes
devtools::install_github("timjmiller/wham", dependencies=TRUE)

# if you're having problems with dependencies not installing:
devtools::install_github("kaskr/adcomp/TMB", dependencies=TRUE)
to.install <- c("plotrix","ellipse","Hmisc","gplots","fields","RColorBrewer","colorspace","mnormt","Deriv","tidyr","dplyr","ggplot2"))
new.packages <- to.install[!(to.install %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

View installed vignettes:
```
browseVignettes("wham")
```

For a clean, runnable `.R` script, look at `ex1_SNEMA_yellowtail_flounder.R` in the `example_scripts` folder of the `wham` package install:
```
library(wham)
wham.dir <- find.package("wham")
file.path(wham.dir, "example_scripts")
```

You can run this entire example script with:
```
write.dir <- "choose/where/to/save/output" # otherwise will be saved in working directory
source(file.path(wham.dir, "example_scripts", "ex1_SNEMA_yellowtail_flounder.R"))
```

### Observation model components:
#### Aggregate relative abundance indices:
* These can be in numbers or biomass, determined by an indicator.
* A normal distribution is assumed for the log-indices. 
* Input standard errors for log-indices are year-specific, but index-specific scalars may be estimated.
* There can be any number of indices.
* There are indicators for whether to use the index in a given year.
* The time of year the index occurs is specified for each year.
* The selectivity "block" to use for an index is specified for each year.
  
#### Index age composition:
* These can also be in numbers or biomass, determined by an indicator.
* Several distribution assumptions available: multinomial, Dirichlet, Dirichlet-multinomial, logistic normal, two types of zero-one inflated logistic normal.
* Effective sample sizes for multinomial are input and year-specific, but variance parameters may be estimated for other distributions.
* The order is the same as the relative abundance indices.
* There are indicators for whether to use the age composition in a given year.
  
#### Aggregate fleet catch:
* This is in biomass.
* A normal distribution is assumed for the log-catch. 
* Input standard errors for log-catch are year-specific, but fleet-specific scalars may be estimated.
* There can be any number of fleets.
* The selectivity "block" to use for an fleet is specified for each year.
  
#### Fleet age composition:
* These can only be in numbers.
* Same distribution assumptions as indices available: multinomial, Dirichlet, Dirichlet-multinomial, logistic normal, two types of zero-one inflated logistic normal.
* Effective sample sizes for multinomial are input and year-specific, but variance parameters may be estimated for other distributions.
* The order is the same as the fleet catches.
* There are indicators for whether to use the age composition in a given year.

##### Hidden process model components:
* Abundance at age:
  * As mentioned above, the model can be configured to have no hidden processes (random effects), 
  * or log-recruitment can be treated as a normal random variable,
  * or all log-abundances at age can be treated as normal random variables. Abundances at age are currently conditionally independent of each other.
* Natural mortality:
  * can be assumed known at age, or estimated as fixed effects parameters.
  * if an allometric function of size is estimated, the (log of the) exponential parameter can be estimated as a normal random effect with a prior based on Lorenzen's work.
  * or annual log-natural mortality parameters can be treated as a normal random walk.

### Other Inputs:

#### Weight-at-age:
* Like ASAP weight-at-age are input and pointers are specified for them for fleet catches, relative abundance indices (if necessary), Spawning Stock Biomass, etc.
* Future versions may include options for internal estimation of weight at age, like [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124).

## Maturity-at-age:
* One matrix of annual maturity at age is input and used to estimate annual spawning stock biomass.


### Parameter options:

* Natural Mortality:
  * In general, this can be year- and age-specific. However, by default it is assumed known at the initial parameter values using the map argument in TMB.
  * It can also be treated as a random walk of random effects, or an allometric function of weight at age.
* Variance parameters for natural mortality:
  * If a random walk in annual natural mortality parameters, then corresponding variance parameters can be estimated.
* Initial numbers at age:
  * Estimate between 1 and n_ages age-specfic parameters
  * Estimate an initial recruitment and an "initial F" and rest of numbers at age are based on equilibrium calculations.
* Recruitment:
  * Can be estimated as annual fixed effects, a random walk, or random effects around a mean (potentially related to spawning stock biomass)
* Abundance at age:
  * Numbers at age after recruitment can be deterministically determined by numbers at previous age in previous year (statistical catch-at-age)
  * Transitions in abundance at age can be stochastic
* Variance parameters for abundance at age
  * If transitions in abundance at age are stochastic, estimate between between 1 and n_ages corresponding variance parameters. 
* Catchability is estimated on a logit scale with the user specifying upper and lower bounds.
* Selectivity parameters: 
* A set of selectivity "blocks" is specified and a pointer to each fleet and index each year is specified.
* Selectivity models for each block can be composed of age-specific parameters, logistic (increasing or decreasing), or double-logistic. 
* Parameters estimated on a logit scale with the user specifying upper and lower bounds.

* Variance parameters for aggregate catch and indices
  * Fleet- or index-specific parameters may be estimated to scale corresponding input standard errors.
* Variance parameters for age composition observations
  * Depends on models assumed, but parameters for each age composition set (i.e., fleet or index) can be estimated.
    
### Ouput Estimates:

* Annual abundance-at-age
* Annual spawning stock biomass
* Selectivity at age for each selectivity "block"
* Natural mortality at age
* Catchability for each relative abundance index, and catchability at age
* Annual fully-selected fishing mortalty for each fleet
* Annual fishing mortality at age for each fleet, and total F at age
* Predicted catch for each fleet and corresponding age composition
* Predicted indices and corresponding age composition
* Estimates of unfished spawning biomass per recruit each year based on corresponding weight-, maturity-, and natural mortality-at-age.
* Biological Reference Points:
  * The user specifies a percentage of unfished spawning biomass per recruit to make estimates of fishing mortality and SSB reference points.
  * Beverton-Holt and Ricker spawner-recruit models may be assumed. Traditional "alpha" and "beta" parameters or steepness and $$R_0$$ may be estimated.
  * Annual conditional reference points are estimated in case weight, maturity, natural mortality, or selectivity at age change over time. 
  * If a spawner-recruit model is assumed, MSY-based reference points are also estimated. 
  * Both of these classes of reference points employ a Newton method internally to determine the fishing mortality reference points, thereby propogating uncertainty of inputs as in [Miller et al. 2016](https://doi.org/10.1139/cjfas-2015-0339) and [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124).
