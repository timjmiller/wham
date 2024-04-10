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
main.dir <- here("sandbox", "pkg_example_results", "ex09")
write.dir <- tempdir(check=TRUE) #will be passed to the ex9_retro_pred.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

# copy to vignette plots folder
to.dir <- file.path(pkg.dir, "vignettes", "ex9_plots")
file.copy(from=file.path(write.dir,"retro_pred_R.png"), to=file.path(to.dir,"retro_pred_R.png"), overwrite = TRUE)

