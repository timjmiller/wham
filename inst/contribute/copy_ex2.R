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
devtools::load_all(path = pkg.dir)

# assumes you have a .gitignore'd folder 'sandbox' in the wham folder
# assumes you ran ex1 today. if not, need to change
main.dir <- file.path(pkg.dir, "sandbox", paste0("runall-",format(Sys.Date(), "%Y%m%d")))
write.dir <- file.path(main.dir,"ex2")
load(file.path(write.dir,"vign2_res.RData")) # get 'df.mods'

mod.list <- file.path(write.dir, paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)
vign2_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))

save(vign2_conv, file=file.path(pkg.dir, "data", "vign2_conv.RData"))
vign2_res <- df.mods
save(vign2_res, file=file.path(pkg.dir, "data", "vign2_res.RData"))

ex2_test_results <- list()
ex2_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
ex2_test_results$nll <- sapply(mods, function(x) x$opt$obj)
saveRDS(ex2_test_results, file=file.path(pkg.dir, "inst", "extdata", "ex2_test_results.rds"))

# copy plots from sandbox to ex2_plots for vignette
to.files <- list.files(file.path(pkg.dir, "vignettes", "ex2_plots"), full.names = T)
from.files <- c(file.path(write.dir,"m1","plots_png","ref_points","Kobe_status.png"),
                file.path(write.dir,"/m2/plots_png/ref_points/Kobe_status.png"),
                file.path(write.dir,"/m3/plots_png/ref_points/Kobe_status.png"),
                file.path(write.dir,"/m4/plots_png/ref_points/Kobe_status.png"),
                file.path(write.dir,"/m5/plots_png/ref_points/Kobe_status.png"),
                file.path(write.dir,"/m1/plots_png/retro/NAA_age1_retro_relative.png"),
                file.path(write.dir,"/m5/plots_png/retro/NAA_age1_retro_relative.png"),
                file.path(write.dir,"/m4/plots_png/diagnostics/OSAresid_ecov_4panel_1.png"),
                file.path(write.dir,"/m5/plots_png/diagnostics/OSAresid_ecov_4panel_1.png"),
                file.path(write.dir,"/m3/plots_png/results/SSB_Rec_fit.png"),
                file.path(write.dir,"/m4/plots_png/results/SSB_Rec_fit.png"))
file.copy(from=from.files, to=to.files, overwrite = T)
