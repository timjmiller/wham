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
main.dir <- here("sandbox", "pkg_example_results", "ex02")
write.dir <- tempdir(check=TRUE) #will be passed to the ex2_CPI_recruitment.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex2_CPI_recruitment.R"))

# save results table
save(df.mods, file=here("data","vign2_res.RData"))
save(df.mods, file=here("sandbox", "pkg_example_results","ex02", "vign2_res.RData"))
save(mods, file = here("sandbox", "pkg_example_results","ex02","ex2_models.RData"))
# write.csv(df.mods, file="vign2_res.csv",row.names=F,quote=F)
for(m in 1:length(mods)){
  plot_wham_output(mod=mods[[m]], dir.main=file.path(write.dir,paste0("m",m)))
  file.copy(from = file.path(write.dir, paste0("m",m)), to = main.dir, recursive = TRUE, overwrite = TRUE, copy.date = TRUE)
}

# pkgload::load_all(pkg.dir, compile=FALSE)
# wham:::plot.MSY.annual(mods[[7]])
# wham:::kobe.plot(mods[[5]], status.years=NULL, static = FALSE, msy = TRUE)

#already made when sourcing ex2_CPI_recruitment.R
# load(file.path(write.dir,"vign2_res.RData")) # get 'df.mods'
# mod.list <- file.path(write.dir, paste0(df.mods$Model,".rds"))
# mod.list <- file.path(write.dir, paste0("m",1:dim(df.mods)[1],".rds"))
# mods <- lapply(mod.list, readRDS)
# vign2_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))

save(vign2_conv, file=here("data", "vign2_conv.RData"))
vign2_res <- df.mods
save(vign2_res, file=here("data", "vign2_res.RData"))

# made in test_ex02 script
# ex2_test_results <- list()
# ex2_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
# ex2_test_results$nll <- sapply(mods, function(x) x$opt$obj)
# saveRDS(ex2_test_results, file=file.path(pkg.dir, "inst", "extdata", "ex2_test_results.rds"))

# copy plots from sandbox to ex2_plots for vignette
from.files <- here("sandbox", "pkg_example_results","ex02", paste0("m",1:7), "plots_png", "ref_points", "Kobe_status.png")
to.files <- here("vignettes", "ex2_plots", paste0("Kobe_status_m", 1:7, ".png"))
file.copy(from=from.files, to=to.files, overwrite = T)

from.files <- here("sandbox", "pkg_example_results","ex02", paste0("m",3:7), "plots_png", "ref_points", "Kobe_msy_status.png")
to.files <- here("vignettes", "ex2_plots", paste0("Kobe_msy_status_m", 3:7, ".png"))
file.copy(from=from.files, to=to.files, overwrite = T)

from.files <- here("sandbox", "pkg_example_results","ex02", paste0("m",c(1,6)), "plots_png", "retro", "stock_1_region_1_NAA_age_1_retro_relative.png")
to.files <- here("vignettes", "ex2_plots", paste0("NAA_age1_retro_relative_m", c(1,6), ".png"))
file.copy(from=from.files, to=to.files, overwrite = T)

from.files <- here("sandbox", "pkg_example_results","ex02", paste0("m",4:7), "plots_png", "results", "SSB_Rec_stock_1_fit.png")
to.files <- here("vignettes", "ex2_plots", paste0("SSB_Rec_fit_m", 4:7, ".png"))
file.copy(from=from.files, to=to.files, overwrite = T)

from.files <- here("sandbox", "pkg_example_results","ex02", paste0("m",4:5), "plots_png", "diagnostics", "OSA_resid_ecov_4panel_CPI.png")
to.files <- here("vignettes", "ex2_plots", paste0("OSA_resid_ecov_4panel_CPI_m", 4:5, ".png"))
file.copy(from=from.files, to=to.files, overwrite = T)

