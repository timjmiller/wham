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
#pkgload::load_all() #this is done in the ex1_basics.R script and running multiple times causes multiple dlls to get loaded and confuse make_osa_residuals()

# assumes you have a .gitignore'd folder 'sandbox' in the wham folder
# assumes you ran ex1 today. if not, need to change
main.dir <- here("sandbox", "pkg_example_results", "ex01")
write.dir <- tempdir(check=TRUE) #will be passed to the ex1_basics.R script for writing all results
path_to_scripts <- here("inst","example_scripts")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex1_basics.R"))
setwd(pkg.dir)

mods[["m4_fixed"]] <- m4_fixed #mods[[5]]

mods[["m4_bad"]] <- m4 #mods[[5]]
mods[["m4_proj"]] <- m4_proj #mods[[6]]
saveRDS(mods, file = file.path(main.dir, "ex1_models.RDS"))
#mods <- readRDS(file.path(main.dir, "ex1_models.RDS")) # get 'mods'
res <- compare_wham_models(mods[c("m1","m2","m3","m4")], table.opts=list(fname="ex1_table", sort=TRUE), fdir = main.dir)
plot_wham_output(mod=m4_proj, dir.main = main.dir)

#setwd(main.dir)

# sapply(mods, function(x) x$runtime)
# sum(sapply(mods, function(x) x$runtime))
#   m1   m2   m3   m4 
# 0.64 1.10 0.13 0.21
# total = 2.08

# test results are made in test_ex01_SNEMAYT.R because bias correction is still used for it. 

to.dir <-here("vignettes", "ex01_files")
dir.create(to.dir)

vign1_res <- res
#save(vign1_res, file=here("data", "vign1_res.RData"))
vign1_res <- res[c("tab","best")]

vign1_m1_conv <- capture.output(check_convergence(mods[["m1"]]), type = "message")
vign1_m2_conv <- capture.output(check_convergence(mods[["m2"]]), type = "message")
vign1_m3_conv <- capture.output(check_convergence(mods[["m3"]]), type = "message")
vign1_m4_conv <- capture.output(check_convergence(mods[["m4"]]), type = "message")
vign1_m4_fixed <- capture.output(check_convergence(mods[["m4_fixed"]]), type = "message")
vign1_m4_bad_conv <- capture.output(check_convergence(mods[["m4_bad"]]), type = "message")
save(vign1_res, vign1_m1_conv,vign1_m2_conv,vign1_m3_conv,vign1_m4_conv,vign1_m4_bad_conv,vign1_m4_fixed, file=file.path(to.dir, "ex01_res.RData"))

# copy plots from sandbox to ex01_files for vignette
# these are hard-coded linux file paths... 
from.dir <- file.path(write.dir,"plots_png") #here("sandbox", "pkg_example_results", "ex01","plots_png") #file.path(write.dir,"plots_png")
#to.dir <- here("vignettes","ex01_files") #file.path(pkg.dir, "vignettes", "ex01_files")
plot.files <- c("diagnostics/Catch_age_comp_index_1_region_1_a.png",
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
                "results/Selectivity_fleet_1_region_1.png",
                "results/SSB_at_age_proportion_stock_1.png",
                "results/SSB_F_trend.png",
                "results/SSB_Rec_loglog_stock_1.png",
                "retro/stock_1_SSB_retro.png",
                "retro/stock_1_SSB_retro_relative.png",
                "input_data/weight_at_age_fleet_1_fleet.png",
                "input_data/scaled_indices.png",
                "input_data/weight_at_age_fleet_1_fleet.png")
from.files <- file.path(from.dir, plot.files)
to.files <- paste0(to.dir, "/", sapply(strsplit(plot.files,"/"), function(x) x[2]))
  # img <- magick::image_read(from.files[1])
  # img_small <- magick::image_scale(img, "600") # Resizes width to 600px
  # magick::image_write(img_small, path = file.path(to.dir, "temp.jpg"), format = "jpg", quality = 75)

# sapply(1:length(plot.files), function(x){
#   img <- magick::image_read(from.files[x])
#   img_small <- magick::image_scale(img, "600") # Resizes width to 600px
#   magick::image_write(img_small, path = to.files[x], format = "png")

#   })

to.files <- paste0(to.dir, "/", sapply(strsplit(plot.files,"/"), function(x) x[2]))
file.copy(from=from.files, to=to.files, overwrite = TRUE)

#use snipping tool to screen grab the html.
file.copy(from = file.path(main.dir, "wham_html_diagnostics.png"), 
  to = file.path(to.dir, "wham_html_diagnostics.png"), overwrite = TRUE)
file.copy(from =file.path(main.dir, "wham_html_tables.png"), 
  to = file.path(to.dir, "wham_html_tables.png"), overwrite = TRUE)

#save input for use in help file examples
saveRDS(input3, here::here("inst","extdata", "ytf_rec1.RDS"))

setwd(here())