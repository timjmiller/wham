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
main.dir <- here("sandbox", "pkg_example_results", "ex08")
write.dir <- tempdir(check=TRUE) #will be passed to the ex5_M_GSI.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex8_compare_asap.R"))
setwd(pkg.dir)

saveRDS(mods, file = file.path(main.dir, "mods.rds"))

# copy to vignette plots folder
to.dir <- file.path(pkg.dir, "vignettes", "ex8_plots")
from.dir <- file.path(write.dir,"compare_png")
from.list <- list.files(from.dir,full.names=T)
to.list <- file.path(to.dir, list.files(from.dir))
file.copy(from=from.list, to=to.dir, overwrite = TRUE)
file.copy(from=file.path(write.dir,"model_comparison.csv"), to=to.dir, overwrite = TRUE)
file.copy(from=file.path(write.dir,"res.rds"), to=to.dir, overwrite = TRUE)

x <- list.files(write.dir)
special.plots <- grep(x, pattern = ".png", fixed = T, value = T)
file.copy(from=file.path(write.dir,special.plots), to=to.dir, overwrite = TRUE)

compare_wham_models(mods, fdir=main.dir, do.table=F, plot.opts=list(which=10, kobe.yr=2010, kobe.prob=F))
file.copy(from=file.path(main.dir,"compare_png","compare_rel_status_kobe.png"), to = file.path(to.dir, "which10_2010.png"), overwrite = TRUE)

res <- readRDS(file.path(write.dir,"res.rds"))
vign8_res <- round(res$tab,2)
save(vign8_res, file=file.path(pkg.dir, "data","vign8_res.RData"))
