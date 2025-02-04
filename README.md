# WHAM: a state-space age-structured assessment model

The Woods Hole Assessment Model (WHAM) is a general age-structured stock assessment framework that can be configured to estimate assessment models that range in complexity from statistical catch-at-age (SCAA) model with annual recruitments as fixed effects, to state-space, multi-stock, multi-region, age-structured models where many parameters can be treated as time- and age-varying process errors and/or allowing effects of environmental covariates. 

The default configuration of WHAM will estimate a traditional SCAA model with no random effects (see [Example 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html), [Example 8](https://timjmiller.github.io/wham/articles/ex8_compare.html), and [Example 13](https://timjmiller.github.io/wham/articles/ex13_no_ASAP.html))

WHAM can estimate process errors as random effects on:

- recruitment ([Example 1](https://timjmiller.github.io/wham/articles/ex1_basics.html), [Example 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html), and [Example 8](https://timjmiller.github.io/wham/articles/ex8_compare.html))
- transitions in numbers-at-age over time ([Example 1](https://timjmiller.github.io/wham/articles/ex1_basics.html), [Example 6](https://timjmiller.github.io/wham/articles/ex6_NAA.html), [Example 8](https://timjmiller.github.io/wham/articles/ex8_compare.html), and [Example 13](https://timjmiller.github.io/wham/articles/ex13_no_ASAP.html),
- selectivity ([Example 4](https://timjmiller.github.io/wham/articles/ex4_selectivity.html)),
- natural mortality ([Example 5](https://timjmiller.github.io/wham/articles/ex5_GSI_M.html)), 
- index catchability ([Example 11](https://timjmiller.github.io/wham/articles/ex11_catchability.html)), and
- movement of stocks between regions ([Example 12](https://timjmiller.github.io/wham/articles/ex12_multistock.html))

WHAM can also include explicit environmental effects on all parameters above (except selectivity) and uses a state-space treatment of the covariates so that observation error of the covariates is considered (e.g., [Example 2](https://timjmiller.github.io/wham/articles/ex2_CPI_recruitment.html) and [Example 5](https://timjmiller.github.io/wham/articles/ex5_GSI_M.html))

A nice property of treating population and environmental processes as random effects is that their uncertainty is naturally propagated in projections/forecasts ([Example 3](https://timjmiller.github.io/wham/articles/ex3_projections.html)).

A presentation providing an overview of an early version of WHAM (Jan 8 2021): https://www.youtube.com/watch?v=o8vJvbIaOdE

## Background

WHAM generalizes and extends R and TMB code from [Miller et al. (2016)](https://doi.org/10.1139/cjfas-2015-0339), [Miller and Hyun 2018](https://doi.org/10.1139/cjfas-2017-0035), and [Miller et al. 2018](https://doi.org/10.1139/cjfas-2017-0124). WHAM has many similarities to ASAP ([code](https://www.nefsc.noaa.gov/nft/ASAP.html), [Legault and Restrepo 1998](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2007/WGMHSA/Annex%203%20-%20ICCAT%20Working%20Document.pdf)), including the input data file structure. Many of the plotting functions for input data, results, and diagnostics are modified from code written by Chris Legault and Liz Brooks ([ASAPplots](https://github.com/cmlegault/ASAPplots)).

A paper describing WHAM has been published, which includes the model equations, simulation tests, and demos of random effects options for numbers-at-age, *M*, selectivity, and environment-recruitment: [https://doi.org/10.1016/j.fishres.2021.105967](https://doi.org/10.1016/j.fishres.2021.105967). Stay tuned for a paper describing generalizations to multiple stocks and regions and movement.

[Stock et al. (2021)](https://doi.org/10.1016/j.fishres.2021.105873) describes the 2D (year x age) AR(1) correlation structure that can be used on numbers-at-age, *M*, and selectivity in WHAM (as in [Berg and Nielsen 2016](https://doi.org/10.1093/icesjms/fsw046), [Cadigan 2016](https://doi.org/10.1139/cjfas-2015-0047), and [Xu et al. 2019](https://doi.org/10.1139/cjfas-2017-0446)).

As mentioned above, WHAM has also been extended (from release 2.0.0) to allow multiple stocks and/or multiple regions to be modeled with movement between regions. Seasonal changes in fleet fishing effort and movement are also possible. The version of WHAM that only allows a single stock and region remains available as the "single_wham" branch, but it will not be developed further.

WHAM is an R package with computation and estimation based on the TMB package [Kristensen et al. 2016](https://doi.org/10.18637/jss.v070.i05), and would not be possible without these superb open-source tools. For more information, see:

- R: https://www.r-project.org/
- TMB: https://kaskr.github.io/adcomp/_book/Introduction.html

## Installation

For the latest stable, tested release, with any bug fixes:
```
devtools::install_github("timjmiller/wham", dependencies=TRUE)
```

For the development version with recent extensions that passes all existing tests, but with potentially untested new features:
```
devtools::install_github("timjmiller/wham@devel", dependencies=TRUE)
```

For the single stock version:
```
devtools::install_github("timjmiller/wham@single_wham")
```

You can even install a specific release or commit. E.g.,
```
devtools::install_github("timjmiller/wham@v1.0.9")
#or equivalently the associated commit
devtools::install_github("timjmiller/wham@40cc14b")
```

### ON WINDOWS

The least frustrating installation of WHAM is via the "pak" package:
```
pak::pkg_install("timjmiller/wham")
```
or for the devel branch:
```
pak::pkg_install("timjmiller/wham@devel")
```
or for single stock wham:
```
pak::pkg_install("timjmiller/wham@single_wham")
```
Using "pak" seems to avoid many installation hurdles.

You may still try to install with devtools, but if you get an error about cc1plus.exe running out of memory during installation, try installing only 64bit:
```
devtools::install_github("timjmiller/wham", dependencies=TRUE, INSTALL_opts=c("--no-multiarch"))
```
or for the devel branch:
```
devtools::install_github("timjmiller/wham", dependencies=TRUE, ref = "devel", INSTALL_opts=c("--no-multiarch"))
```

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

NOTE: If you specify `pak::pkg_install` to install wham to library directory that is not in your `.libPaths()`. It is best to add the directory because `pak` installs all dependencies to the same directory and may not be found for example when creating model output with `prepare_wham_input`. E.g.,
```
pak::pkg_install("timjmiller/wham@single_wham", lib = "single_wham")
library(wham, lib.loc = "single_wham")
.libPaths("single_wham")
```

## Releases and Versions

WHAM is now used for management for several stocks and to reduce confusion during the assessment process the latest commit pushed to the master and devel branch will get a unique version number from release v2.0.0 forward.

All releases of the package have versions denoted as (v)a.b.0. (e.g., the first release is v1.0.0). Releases that have lesser changes will increment the second number (e.g., v2.1.0) and releases with large changes will increment the first number (e.g. v3.0.0). When there are commits for bug fixes to the latest release in the master branch of the repo they will increment the third number (e.g., v2.0.1). 

The devel branch will always be up-to-date with the master branch and versioned initially as (v)a.b.c.9000. For example, when the released version is changed from v2.0.0 to v2.0.1, devel branch would also go from v2.0.0.9xxx to v2.0.1.9000. The latest pushed commits to devel branch between fixes to the released version will also have a unique version. For example, the next pushed change to devel branch would be v2.0.1.9001.

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

You can run ALL examples with (this is not quick):
```
library(wham)
wham.dir <- find.package("wham")
source(file.path(wham.dir, "example_scripts", "run_all_examples.R"))
```

## Short-course materials

A short course was given in Woods Hole in June 2024 on using the WHAM package. Slides and corresponding R scripts are available in this [repository](https://github.com/timjmiller/wham_course_WH_2024)

A workshop was given at Memorial University in September 2024 on using the WHAM package, that expanded on the previous short course and included making inputs without ASAP dat files and fitting multi-stock models. Slides and corresponding R scripts are available in this [repository](https://github.com/timjmiller/wham_workshop_MUN_2024).

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

