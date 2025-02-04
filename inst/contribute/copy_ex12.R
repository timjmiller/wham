# after running the example scripts, copy output to correct locations:
#   plot files used by vignettes (vignettes/exY_plots/)
#   data files used by vignettes (data/)
#   test results (inst/extdata/exY_test_results.rds)

# This file uses the 'here' package to work across computers, as long as you 
#   *open from the top-level wham folder* (where you git cloned wham into).
# If you do not have 'here' installed:
# install.packages("here")
library(here)
pkg.dir <- here()
pkgbuild::compile_dll(debug = FALSE); 
pkgload::load_all(compile=FALSE)

# assumes you have a .gitignore'd folder 'sandbox' in the wham folder
# assumes you ran the ex today. if not, need to change
main.dir <- here("sandbox", "pkg_example_results", "ex12")
write.dir <- tempdir(check=TRUE) #will be passed to the ex*.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")
#to.dir <- here("vignettes", "ex12_plots")

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex12_multistock.R"))
setwd(pkg.dir)

vign_12_info <- list()

vign_12_info$move <- list(
	nofit_seasonal_Ps_terminal = nofit_move$rep$seasonal_Ps_terminal, 	
	nofit_annual_Ps = nofit_move$rep$annual_Ps,
	obj = fit_move$opt$obj,
	mu = fit_move$rep$mu
)
vign_12_info$move_age <- list(
	n_unique_mu_re = length(unique(input_move_age$map$mu_re)),
	nofit_mu = nofit_move_age$rep$mu,
	obj = fit_move_age$opt$obj,
	mu = fit_move_age$rep$mu
)
vign_12_info$move_year <- list(
	n_unique_mu_re = length(unique(input_move_year$map$mu_re)),
	nofit_mu = nofit_move_year$rep$mu,
	obj = fit_move_year$opt$obj,
	mu = fit_move_year$rep$mu,
	mu_repars = fit_move_year$parList$mu_repars
)
vign_12_info$move_prior <- list(
	mean_vals = move_prior$mean_vals,
	n_mu_re = length(unique(input_move_prior$map$mu_re)), #+1 for NAs
	n_mu_prior_re = length(unique(input_move_prior$map$mu_prior_re)), #+1 for NAs
	mu = fit_move_prior$rep$mu
)
save(vign_12_info, file=file.path(pkg.dir, "data","vign_12_info.RData"))
