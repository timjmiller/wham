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
write.dir <- file.path(main.dir,"ex1")
load(file.path(write.dir,"ex1_models.RData")) # get 'mods'
res <- compare_wham_models(mods, table.opts=list(fname="ex1_table", sort=TRUE))

# sapply(mods, function(x) x$runtime)
# sum(sapply(mods, function(x) x$runtime))
#   m1   m2   m3   m4 
# 0.64 1.10 0.13 0.21
# total = 2.08

ex1_test_results <- list(nll=sapply(mods, function(x) x$opt$obj),
						 par=lapply(mods, function(x) as.numeric(x$opt$par)))
saveRDS(ex1_test_results, file=file.path(pkg.dir, "inst", "extdata", "ex1_test_results.rds"))

# data files used by vignettes must be .RData with one object saved (eg 'vign1_m1_conv'), and that object's name must match the filename ('vign1_m1_conv.RData')
# for example, the commented-out line below will not work (can't just save 'res' as 'vign1_res.RData')
#   save(res, file=file.path(pkg.dir, "data", "vign1_res.RData"))
vign1_res <- res
save(vign1_res, file=file.path(pkg.dir, "data", "vign1_res.RData"))

vign1_m1_conv <- capture.output(check_convergence(mods[[1]]))
vign1_m2_conv <- capture.output(check_convergence(mods[[2]]))
vign1_m3_conv <- capture.output(check_convergence(mods[[3]]))
vign1_m4_conv <- capture.output(check_convergence(mods[[4]]))
save(vign1_m1_conv, file=file.path(pkg.dir, "data", "vign1_m1_conv.RData"))
save(vign1_m2_conv, file=file.path(pkg.dir, "data", "vign1_m2_conv.RData"))
save(vign1_m3_conv, file=file.path(pkg.dir, "data", "vign1_m3_conv.RData"))
save(vign1_m4_conv, file=file.path(pkg.dir, "data", "vign1_m4_conv.RData"))

# copy plots from sandbox to ex1_plots for vignette
# these are hard-coded linux file paths... 
from.dir <- file.path(write.dir,"plots_png")
to.dir <- file.path(pkg.dir, "vignettes", "ex1_plots")
plot.files <- c("diagnostics/Catch_age_comp_index1.png",
                "diagnostics/Catch_age_comp_resids_index1.png",
                "misc/catch_at_age_consistency_obs_fleet1.png",
                "input_data/catch_by_fleet.png",
                "misc/catch_curves_fleet1_obs.png",
                "ref_points/FSPR_relative.png",
                "input_data/index.png",
                "diagnostics/Index_4panel_2.png",
                "ref_points/Kobe_status.png",
                "diagnostics/likelihood.png",
                "diagnostics/NAA_4panel_1.png",
                "diagnostics/NAA_4panel_5.png",
                "diagnostics/OSAresid_catch_4panel_fleet1.png",
                "diagnostics/OSAresid_catch_4panel_index1.png",
                "diagnostics/OSAresid_catch_4panel_index2.png",
                "results/Selectivity_fleet1.png",
                "ref_points/SPR_targets_ave_plot.png",
                "results/SSB_at_age_proportion.png",
                "results/SSB_F_trend.png",
                "results/SSB_Rec_loglog.png",
                "retro/SSB_retro.png",
                "retro/SSB_retro_relative.png",
                "input_data/weight_at_age_fleet1.png")
from.files <- file.path(from.dir, plot.files)
to.files <- paste0(to.dir, sapply(strsplit(plot.files,"/"), function(x) x[2]))
file.copy(from=from.files, to=to.files, overwrite = T)
