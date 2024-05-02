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
main.dir <- here("sandbox", "pkg_example_results", "ex11")
write.dir <- tempdir(check=TRUE) #will be passed to the ex*.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")
to.dir <- here("vignettes", "ex11_plots")

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex11_catchability.R"))
setwd(pkg.dir)

file.copy(from = file.path(write.dir, "prior_posterior_q.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "prior_posterior_q.png"), to = to.dir, overwrite = TRUE)

for(i in 1:4){
	file.copy(from = file.path(write.dir, paste0("q_time_series_",i,".png")), to = main.dir, overwrite = TRUE)
	file.copy(from = file.path(write.dir, paste0("q_time_series_",i,".png")), to = to.dir, overwrite = TRUE)
}
