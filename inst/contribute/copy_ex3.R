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
#pkgload::load_all(compile=FALSE)

# assumes you have a .gitignore'd folder 'sandbox' in the wham folder
# assumes you ran ex1 today. if not, need to change
main.dir <- here("sandbox", "pkg_example_results", "ex03")
write.dir <- tempdir(check=TRUE) #will be passed to the ex3_projections.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex3_projections.R"))
setwd(pkg.dir)

mod <- readRDS(file.path(write.dir,"m5.rds"))
mod_proj <- readRDS(file.path(write.dir,"m5_proj.rds"))

ex3_test_results <- list(pars=NULL, nll=NULL)
ex3_test_results$pars <- mod$opt$par
ex3_test_results$nll <- as.numeric(mod$opt$obj)
saveRDS(ex3_test_results, file=here("tests", "testthat", "fixtures", "ex03_test_results.RDS"))

file.copy(from=file.path(write.dir,"mod.rds"), to=main.dir, overwrite = T)
file.copy(from=file.path(write.dir,"mod_proj.rds"), to=main.dir, overwrite = T)
for(i in 1:11) file.copy(from = file.path(write.dir, paste0("proj_",i)), to = main.dir, recursive = TRUE, overwrite = TRUE, copy.date = TRUE)

vign3_nll_orig <- mod$opt$obj
vign3_nll_proj <- nll_proj

to.dir <-here("vignettes", "ex03_files")
dir.create(to.dir)

save(vign3_nll_proj, vign3_nll_orig, file=file.path(to.dir, "ex03_nlls.RData"))

# copy plots from sandbox to ex3_plots for vignette
from.files <- c(file.path(write.dir,paste0("proj_",1:5), "plots_png", "results", "Ecov_1_CPI.png"),
  file.path(write.dir,paste0("proj_",1:11), "plots_png", "results", "SSB_F_trend.png"),
  file.path(write.dir,paste0("proj_",1:11), "plots_png", "ref_points", "Kobe_status.png"),
  file.path(write.dir,paste0("proj_",1:11), "plots_png", "ref_points", "Kobe_msy_status.png"))
to.files <- c(
  file.path(to.dir, paste0("proj_",1:5, "_Ecov_1_CPI.png")),
  file.path(to.dir, paste0("proj_",1:11, "_SSB_F_trend.png")),
  file.path(to.dir, paste0("proj_",1:11, "_Kobe_status.png")),
  file.path(to.dir, paste0("proj_",1:11, "_Kobe_msy_status.png")))

file.copy(from=from.files, to=to.files, overwrite = T)
