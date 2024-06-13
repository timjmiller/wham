# WHAM: a state-space age-structured assessment model

The Woods Hole Assessment Model (WHAM) is a general state-space age-structured stock assessment framework designed to include environmental effects on population processes. The state-space framework is attractive because it can estimate observation and process error, as well as naturally propagate random effect parameters in stock projections. WHAM can be configured to estimate a range of assessment models (see [Ex 1](https://timjmiller.github.io/wham/articles/ex1_basics.html) and [Ex 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html)):

- statistical catch-at-age (SCAA) model with recruitments as fixed effects, 
- SCAA with recruitments as random effects
- "full state-space model", abundance at all ages are random effects
- multi-stock and/or multi-region abundance at age as random effects

WHAM advances fisheries assessment because it can estimate constrained random deviations, i.e. random effects, on parameters such as:

- recruitment / numbers-at-age ([Ex 2](https://timjmiller.github.io/wham/articles/ex2_CPI_recruitment.html) and [Ex 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html)),
- selectivity ([Ex 4](https://timjmiller.github.io/wham/articles/ex4_selectivity.html)),
- natural mortality ([Ex 5](https://timjmiller.github.io/wham/articles/ex5_GSI_M.html)), 
- index catchability ([Ex 11](https://timjmiller.github.io/wham/articles/ex11_catchability.html)), 
- movement between regions, and
- environmental effects on the above (e.g., [Ex 2](https://timjmiller.github.io/wham/articles/ex2_CPI_recruitment.html) and [Ex 5](https://timjmiller.github.io/wham/articles/ex5_GSI_M.html))

A nice property of treating population and environmental processes as random effects is that their uncertainty is naturally propagated in projections/forecasts ([Ex 3](https://timjmiller.github.io/wham/articles/ex3_projections.html)).

Overview of WHAM presentation (Jan 8 2021): https://www.youtube.com/watch?v=o8vJvbIaOdE

## Background

WHAM generalizes and extends R and TMB code from [Miller et al. (2016)](https://doi.org/10.1139/cjfas-2015-0339), [Miller and Hyun 2018](https://doi.org/10.1139/cjfas-2017-0035), and [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124). WHAM has many similarities to ASAP ([code](https://www.nefsc.noaa.gov/nft/ASAP.html), [Legault and Restrepo 1998](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2007/WGMHSA/Annex%203%20-%20ICCAT%20Working%20Document.pdf)), including the input data file structure. Many of the plotting functions for input data, results, and diagnostics are modified from ASAP code written by Chris Legault and Liz Brooks ([ASAPplots](https://github.com/cmlegault/ASAPplots)).

A paper describing WHAM has been published, which includes the model equations, simulation tests, and demos of random effects options for numbers-at-age, *M*, selectivity, and environment-recruitment: [https://doi.org/10.1016/j.fishres.2021.105967](https://doi.org/10.1016/j.fishres.2021.105967).

[Stock et al. (2021)](https://doi.org/10.1016/j.fishres.2021.105873) describes the 2D (year x age) AR(1) correlation structure that can be used on numbers-at-age, *M*, and selectivity in WHAM (as in [Berg and Nielsen 2016](https://doi.org/10.1093/icesjms/fsw046), [Cadigan 2016](https://doi.org/10.1139/cjfas-2015-0047), and [Xu et al. 2019](https://doi.org/10.1139/cjfas-2017-0446)).

As mentioned above, WHAM has also been extended (after release 1.0.9) to allow multiple stocks and/or multiple regions to be modeled with movement between regions. Seasonal changes in fleet fishing effort and movement are also possible. The single stock version of WHAM remains available as the "single_wham" branch, but it will not be developed further.

WHAM is written in R and TMB, and would not be possible without these superb open-source tools. For more information, see:

- R: https://www.r-project.org/
- TMB: https://kaskr.github.io/adcomp/_book/Introduction.html

## Installation

For the latest stable, tested release:

```
devtools::install_github("timjmiller/wham", dependencies=TRUE)
```

For the development version with recent bug fixes and features (potentially untested):

```
devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
```

### ON WINDOWS

If you get an error about cc1plus.exe running out of memory during installation, try installing only 64bit:
```
devtools::install_github("timjmiller/wham", dependencies=TRUE, INSTALL_opts=c("--no-multiarch"))
```
or for the devel branch:
```
devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "devel", INSTALL_opts=c("--no-multiarch"))
```

Also consider using the "pak" package for installation:
```
pak::pkg_install("timjmiller/wham")
```
or for the devel branch:
```
pak::pkg_install("timjmiller/wham@devel")
```
Using "pak" seems to avoid many installation hurdles.

If you're having problems with dependencies not installing. It is probably because some are being used in one or more R sessions. After closing all R sessions and restarting R without any packages first check make sure no packages are loaded (even by e.g. .Rprofile):
```
ls() ## no variables
search() ## no packages other than Base
```
Then:
```
to.install <- c("plotrix","ellipse","Hmisc","gplots","fields","RColorBrewer","colorspace","mnormt","Deriv","tidyr","dplyr","ggplot2","viridis", "abind", "rmarkdown", "pander", "kableExtra")
new.packages <- to.install[!(to.install %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```
Then use one of the calls above to install wham.

If you want pdfs of parameter tables that are generted by plot_wham_output you will need a tex installation. If you do not use RStudio, use the tinytex package:
```
install.packages("tinytex")
tinytex::install_tinytex()
```
and add the path to pandoc in your .Rprofile so Rmarkdown can find your pandoc
```
Sys.setenv(RSTUDIO_PANDOC = "path/to/your/pandoc") 
```

## Tutorial

We suggest walking through the vignettes to familiarize yourself with WHAM: https://timjmiller.github.io/wham/articles.

Clean, runnable `.R` scripts for most vignettes are also available in the `example_scripts` folder of the `wham` package install:
```
library(wham)
wham.dir <- find.package("wham")
file.path(wham.dir, "example_scripts")
```

You can then run the entire first example script with:
```
setwd("choose/where/to/save/output")
source(file.path(wham.dir, "example_scripts", "ex1_basics.R"))
```

You can run ALL examples with (takes 1 hour):
```
library(wham)
wham.dir <- find.package("wham")
source(file.path(wham.dir, "example_scripts", "run_all_examples.R"))
```

## Short-course materials

A short course was given in Woods Hole in June 2024 on using the WHAM package. Slides and corresponding R scripts are available in this [repository](https://github.com/timjmiller/wham_course_WH_2024)


## Installing vignettes

Installation from GitHub does not include the vignettes by default because they can be accessed online anytime at https://timjmiller.github.io/wham/articles. If you want to build the vignettes locally, they look best if you *build using R Studio*:
```
devtools::install_github("timjmiller/wham", dependencies=TRUE, build_vignettes = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
```

View installed vignettes:
```
browseVignettes("wham")
```

## References

> Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and Bell, B. M. 2016. TMB: Automatic differentiation and Laplace approximation. Journal of Statistical Software 70(5): 1–21. doi: [10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05).

> Legault, C. M. and Restrepo, V. R. 1998. A flexible forward age-structured assessment program. ICCAT. Col. Vol. Sci. Pap. 49: 246-253. [link](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2007/WGMHSA/Annex%203%20-%20ICCAT%20Working%20Document.pdf)

> Miller, T. J., Hare, J. A., and Alade, L. A. 2016. A state-space approach to incorporating environmental effects on recruitment in an age-structured assessment model with an application to Southern New England yellowtail flounder. Canadian Journal of Fisheries and Aquatic Sciences 73(8): 1261-1270. doi: [10.1139/cjfas-2015-0339](https://doi.org/10.1139/cjfas-2015-0339)

> Miller, T. J. and Hyun, S-Y. 2018. Evaluating evidence for alternative natural mortality and process error assumptions using a state-space, age-structured assessment model. Canadian Journal of Fisheries and Aquatic Sciences 75(5): 691-703. doi: [10.1139/cjfas-2017-0035](https://doi.org/10.1139/cjfas-2017-0035)

> Miller, T. J., O’Brien, L., and Fratantoni, P. S. 2018. Temporal and environmental variation in growth and maturity and effects on management reference points of Georges Bank Atlantic cod. Canadian Journal of Fisheries and Aquatic Sciences 75(12): 2159-2171. doi: [10.1139/cjfas-2017-0124](https://doi.org/10.1139/cjfas-2017-0124)

> R Core Team. 2020. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

> Stock, B. C., and Miller, T. J. 2021. The Woods Hole Assessment Model (WHAM): A general state-space assessment framework that incorporates time- and age-varying processes via random effects and links to environmental covariates. Fisheries Research, 240: 105967. doi: [10.1016/j.fishres.2021.105967](https://doi.org/10.1016/j.fishres.2021.105967)

> Stock, B. C., Xu, H., Miller, T. J., Thorson, J. T., and Nye, J. A. 2021. Implementing two-dimensional autocorrelation in either survival or natural mortality improves a state-space assessment model for Southern New England-Mid Atlantic yellowtail flounder. Fisheries Research, 237: 105873. doi: [10.1016/j.fishres.2021.105873](https://doi.org/10.1016/j.fishres.2021.105873)

## NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.

****************************

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [NOAA Fisheries](https://www.fisheries.noaa.gov/)

