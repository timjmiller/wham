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

res <- compare_wham_models(mods, fdir=main.dir, plot.opts = list(browse = FALSE))
compare.png.files <- grep(list.files(file.path(main.dir, "compare_png")), pattern = ".png", fixed = T, value = T)

# copy to vignette plots folder
to.dir <- file.path(pkg.dir, "vignettes", "ex08_files")
dir.create(to.dir)
file.copy(from=file.path(main.dir,"compare_png",compare.png.files), to=to.dir, overwrite = TRUE)

x <- list.files(write.dir)
special.plots <- grep(x, pattern = ".png", fixed = T, value = T)
file.copy(from=file.path(write.dir,special.plots), to=to.dir, overwrite = TRUE)

res <- readRDS(file.path(write.dir,"res.rds"))
res_table <- round(res$tab,2)
saveRDS(res_table, file=file.path(to.dir, "res_table.RDS"))

compare_wham_models(mods, fdir=write.dir, do.table=F, plot.opts=list(which=10, kobe.yr=2010, kobe.prob=F, browse = FALSE))
file.copy(from=file.path(main.dir,"compare_png","compare_rel_status_kobe.png"), to = file.path(to.dir, "which10_2010.png"), overwrite = TRUE)

