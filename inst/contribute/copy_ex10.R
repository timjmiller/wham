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
main.dir <- here("sandbox", "pkg_example_results", "ex10")
write.dir <- tempdir(check=TRUE) #will be passed to the ex9_retro_pred.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")
to.dir <- here("vignettes", "ex10_files")
dir.create(to.dir)

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex10_simulation.R"))
setwd(pkg.dir)

#save RData file
mod_1$sdrep = summary(mod_1$sdrep, "fixed")
mod_1_simulate <- paste(capture.output(mod_1$simulate), collapse = "\n")
save.objects <- c("mod_1", "mod_1_simulate", "vign_10_3_input", "vign_10_3a_sim", "vign_10_3b_sim", 
	"vign_10_4_om_input", "vign_10_4a_om_input", "vign_10_4b_om_input", "vign_10_4c_om_sel_repars", "vign_10_4d_om_q_repars", "vign_10_4e_om_Ecov_process_pars", "vign_10_4e_om_Ecov_beta_R")
save(list = save.objects, file = file.path(main.dir, "ex_10.RData"))
save(list = save.objects, file = file.path(to.dir, "ex_10.RData"))
#cat(mod_1_simulate, sep="\n")

files <- c("self_fit_1_ssb.png", "ssb_rel_bias.png", "sim_fit_3a.png", "sim_fit_3b.png", "sim_fit_3c.png", "sim_fit_3d_1.png", "sim_fit_3d_2.png",
	"sim_plot_4d.png", "sim_plot_4e_1.png", "sim_plot_4e_2.png", "sim_fit_4f_1.png", "sim_fit_4f_2.png", "looped_rep_5_1.png", "looped_rep_5_2.png")
from.files <- file.path(write.dir, files)
file.copy(from = from.files, to = main.dir, overwrite = TRUE)
file.copy(from = from.files, to = to.dir, overwrite = TRUE)

