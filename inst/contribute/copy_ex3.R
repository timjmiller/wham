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
# assumes you ran ex1 today. if not, need to change
main.dir <- here("sandbox", "pkg_example_results", "ex03")
write.dir <- tempdir(check=TRUE) #will be passed to the ex3_projections.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex3_projections.R"))


mod <- readRDS(file.path(write.dir,"m5.rds"))
mod_proj <- readRDS(file.path(write.dir,"m5_proj.rds"))

ex3_test_results <- list(pars=NULL, nll=NULL)
ex3_test_results$pars <- mod$opt$par
ex3_test_results$nll <- as.numeric(mod$opt$obj)
saveRDS(ex3_test_results, file=file.path(pkg.dir, "inst","extdata","ex3_test_results.rds"))

file.copy(from=file.path(write.dir,"m5.rds"), to=main.dir, overwrite = T)
file.copy(from=file.path(write.dir,"m5_proj.rds"), to=main.dir, overwrite = T)
for(i in 1:11) file.copy(from = file.path(write.dir, paste0("proj_",i)), to = main.dir, recursive = TRUE, overwrite = TRUE, copy.date = TRUE)

mod <- readRDS(file.path(main.dir,"m5.rds"))
mod_proj <- readRDS(file.path(main.dir,"m5_proj.rds"))

#this is now provided when running project_wham
# for(i in 1:length(mod_proj)){
#   mod_proj[[i]]$marg_nll <- mod_proj[[i]]$fn()
# }
#mod_proj <- saveRDS(mod_proj, file.path(main.dir,"m5_proj.rds")) 

nll_proj <-  sapply(mod_proj, function(x) x$marg_nll)

vign3_nll_orig <- mod$opt$obj
vign3_nll_proj <- nll_proj
save(vign3_nll_proj, vign3_nll_orig, file=file.path(pkg.dir, "data", "vign3_nlls.RData"))

# copy plots from sandbox to ex3_plots for vignette
from.files <- c(file.path(write.dir,paste0("proj_",1:5), "plots_png", "results", "Ecov_1_CPI.png"),
  file.path(write.dir,paste0("proj_",1:11), "plots_png", "results", "SSB_F_trend.png"),
  file.path(write.dir,paste0("proj_",1:11), "plots_png", "ref_points", "Kobe_status.png"),
  file.path(write.dir,paste0("proj_",1:11), "plots_png", "ref_points", "Kobe_msy_status.png"))
to.files <- c(
  file.path(pkg.dir, "vignettes", "ex3_plots", paste0("proj_",1:5, "_Ecov_1_CPI.png")),
  file.path(pkg.dir, "vignettes", "ex3_plots", paste0("proj_",1:11, "_SSB_F_trend.png")),
  file.path(pkg.dir, "vignettes", "ex3_plots", paste0("proj_",1:11, "_Kobe_status.png")),
  file.path(pkg.dir, "vignettes", "ex3_plots", paste0("proj_",1:11, "_Kobe_msy_status.png")))

file.copy(from=from.files, to=to.files, overwrite = T)
