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
#main.dir <- file.path(pkg.dir, "sandbox", paste0("runall-",format(Sys.Date(), "%Y%m%d")))
main.dir <- here("sandbox", "pkg_example_results", "ex01")
write.dir <- tempdir(check=TRUE) #will be passed to the ex1_basics.R script for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex1_basics.R"))

mods$m4_bad <- m4 #mods[[5]]
mods$m4_proj <- m4_proj #mods[[6]]
save(mods, file = here("sandbox", "pkg_example_results","ex01","ex1_models.RData"))
#load(here("sandbox", "pkg_example_results","ex01","ex1_models.RData")) # get 'mods'
res <- compare_wham_models(mods[1:4], table.opts=list(fname="ex1_table", sort=TRUE), fdir = here("sandbox", "pkg_example_results","ex01"))
setwd(main.dir)

# sapply(mods, function(x) x$runtime)
# sum(sapply(mods, function(x) x$runtime))
#   m1   m2   m3   m4 
# 0.64 1.10 0.13 0.21
# total = 2.08

# made in test_ex01 script
# ex1_test_results <- list(nll=sapply(mods, function(x) x$opt$obj),
# 						 par=lapply(mods, function(x) as.numeric(x$opt$par)))
# saveRDS(ex1_test_results, file=file.path(pkg.dir, "inst", "extdata", "ex1_test_results.rds"))

# data files used by vignettes must be .RData with one object saved (eg 'vign1_m1_conv'), and that object's name must match the filename ('vign1_m1_conv.RData')
# for example, the commented-out line below will not work (can't just save 'res' as 'vign1_res.RData')
#   save(res, file=file.path(pkg.dir, "data", "vign1_res.RData"))
vign1_res <- res
save(vign1_res, file=here("data", "vign1_res.RData"))

vign1_m1_conv <- capture.output(check_convergence(mods[[1]]))
vign1_m2_conv <- capture.output(check_convergence(mods[[2]]))
vign1_m3_conv <- capture.output(check_convergence(mods[[3]]))
vign1_m4_conv <- capture.output(check_convergence(mods[[4]]))
vign1_m4_bad_conv <- capture.output(check_convergence(mods[[5]]))
save(vign1_m1_conv,vign1_m2_conv,vign1_m3_conv,vign1_m4_conv,vign1_m4_bad_conv, file=here("data", "vign1_conv.RData"))
# save(vign1_m2_conv, file=here("data", "vign1_m2_conv.RData"))
# save(vign1_m3_conv, file=here("data", "vign1_m3_conv.RData"))
# save(vign1_m4_conv, file=here("data", "vign1_m4_conv.RData"))
# save(vign1_m4_bad_conv, file=here("data", "vign1_m4_bad_conv.RData"))

# copy plots from sandbox to ex1_plots for vignette
# these are hard-coded linux file paths... 
from.dir <- here("sandbox", "pkg_example_results", "ex01","plots_png") #file.path(write.dir,"plots_png")
to.dir <- here("vignettes","ex1_plots") #file.path(pkg.dir, "vignettes", "ex1_plots")
plot.files <- c("diagnostics/Catch_age_comp_index_1_region_1.png",
                "diagnostics/Catch_age_comp_resids_index_1.png",
                "diagnostics/Index_4panel_index_2_region_1.png",
                "diagnostics/NAA_4panel_stock_1_region_1_age_1.png",
                "diagnostics/NAA_4panel_stock_1_region_1_age_5.png",
                "diagnostics/likelihood.png",
                "diagnostics/OSA_resid_catch_4panel_fleet_1.png",
                "diagnostics/OSA_resid_catch_4panel_index_1.png",
                "diagnostics/OSA_resid_catch_4panel_index_2.png",
                "diagnostics/OSA_resid_paa_6panel_fleet_1.png",
                "diagnostics/OSA_resid_paa_6panel_index_1.png",
                "diagnostics/OSA_resid_paa_6panel_index_2.png",
                "misc/catch_at_age_consistency_obs_fleet_1_region_1.png",
                "misc/catch_curves_fleet_1_region_1_obs.png",
                "ref_points/FSPR_absolute.png",
                "ref_points/FSPR_relative.png",
                "ref_points/Kobe_status.png",
                #"ref_points/SPR_targets_ave_plot.png",
                "results/Selectivity_fleet_1_region_1.png",
                "results/SSB_at_age_proportion_stock_1.png",
                "results/SSB_F_trend.png",
                "results/SSB_Rec_loglog_stock_1.png",
                "retro/stock_1_SSB_retro.png",
                "retro/stock_1_SSB_retro_relative.png",
                "input_data/catch_by_fleet.png",
                "input_data/index.png",
                "input_data/weight_at_age_fleet_1_fleet.png")
from.files <- file.path(from.dir, plot.files)
to.files <- paste0(to.dir, "/", sapply(strsplit(plot.files,"/"), function(x) x[2]))
file.copy(from=from.files, to=to.files, overwrite = T)

file.copy(from = here("sandbox", "pkg_example_results", "ex01", "wham_html_diagnostics.png"), 
  to = here("vignettes", "wham_html_diagnostics.png"), overwrite = T)
file.copy(from = here("sandbox", "pkg_example_results", "ex01", "wham_html_tables.png"), 
  to = here("vignettes", "wham_html_tables.png"), overwrite = T)
