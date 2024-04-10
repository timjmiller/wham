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
main.dir <- here("sandbox", "pkg_example_results", "ex10")
write.dir <- tempdir(check=TRUE) #will be passed to the ex9_retro_pred.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")
to.dir <- here("vignettes", "ex10_plots")

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex10_simulation.R"))
setwd(pkg.dir)

#Section 2

file.copy(from = file.path(write.dir, "self_fit_1_ssb.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "self_fit_1_ssb.png"), to = to.dir, overwrite = TRUE)

file.copy(from = file.path(write.dir, "ssb_rel_bias.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "ssb_rel_bias.png"), to = to.dir, overwrite = TRUE)


#Section 3

file.copy(from = file.path(write.dir, "sim_fit_3a.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "sim_fit_3a.png"), to = to.dir, overwrite = TRUE)

file.copy(from = file.path(write.dir, "sim_fit_3b.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "sim_fit_3b.png"), to = to.dir, overwrite = TRUE)

file.copy(from = file.path(write.dir, "sim_fit_3c.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "sim_fit_3c.png"), to = to.dir, overwrite = TRUE)

file.copy(from = file.path(write.dir, "sim_fit_3d_1.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "sim_fit_3d_1.png"), to = to.dir, overwrite = TRUE)

file.copy(from = file.path(write.dir, "sim_fit_3d_2.png"), to = main.dir, overwrite = TRUE)
file.copy(from = file.path(write.dir, "sim_fit_3d_2.png"), to = to.dir, overwrite = TRUE)

save(input_3, sim_3a, sim_3b, file = file.path(main.dir, "vign_10_3.RData"))
save(input_3, sim_3a, sim_3b, file = here("data", "vign_10_3.RData"))

#Section 4
save(stock_om_4, file = file.path(main.dir, "vign_10_4_stock_om.RData"))
save(stock_om_4, file = here("data", "vign_10_4_stock_om.RData"))

