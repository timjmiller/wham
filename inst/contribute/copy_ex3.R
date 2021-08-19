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
write.dir <- file.path(main.dir,"ex3")

mod <- readRDS(file.path(write.dir,"m5.rds"))
mod_proj <- readRDS(file.path(write.dir,"m5_proj.rds"))
nll_proj <-  sapply(mod_proj, function(x) x$fn(x$opt$obj)

vign3_nll_orig <- mod$opt$obj
vign3_nll_proj <- nll_proj
save(vign3_nll_orig, file=file.path(pkg.dir, "data", "vign3_nll_orig.RData"))
save(vign3_nll_proj, file=file.path(pkg.dir, "data", "vign3_nll_proj.RData"))

# copy plots from sandbox to ex3_plots for vignette
to.files <- list.files(file.path(pkg.dir, "vignettes", "ex3_plots"), full.names = T)
from.files <- c(paste0(write.dir,"/proj_1/plots_png/results/Ecov_1.png"),
                paste0(write.dir,"/proj_2/plots_png/results/Ecov_1.png"),
                paste0(write.dir,"/proj_3/plots_png/results/Ecov_1.png"),
                paste0(write.dir,"/proj_4/plots_png/results/Ecov_1.png"),
                paste0(write.dir,"/proj_5/plots_png/results/Ecov_1.png"),
                paste0(write.dir,"/proj_1/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_2/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_3/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_4/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_5/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_6/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_7/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_8/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_9/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_10/plots_png/ref_points/Kobe_status.png"),
                paste0(write.dir,"/proj_1/plots_png/results/SSB_F_trend.png"),
                paste0(write.dir,"/proj_6/plots_png/results/SSB_F_trend.png"),
                paste0(write.dir,"/proj_7/plots_png/results/SSB_F_trend.png"),
                paste0(write.dir,"/proj_8/plots_png/results/SSB_F_trend.png"),
                paste0(write.dir,"/proj_9/plots_png/results/SSB_F_trend.png"),
                paste0(write.dir,"/proj_10/plots_png/results/SSB_F_trend.png"))
file.copy(from=from.files, to=to.files, overwrite = T)

