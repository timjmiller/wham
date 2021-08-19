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
# assumes you ran the ex today. if not, need to change
main.dir <- file.path(pkg.dir, "sandbox", paste0("runall-",format(Sys.Date(), "%Y%m%d")))
write.dir <- file.path(main.dir,"ex6")

df.mods <- read.csv(file.path(write.dir,"ex6_table.csv"))
mod.list <- file.path(write.dir, paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# copy to vignette plots folder
to.dir <- file.path(pkg.dir, "vignettes", "ex6_plots")
for(m in c(1,7,13)){
  from.dir <- paste0(write.dir,"/m",m,"/plots_png")
  file.copy(from=file.path(from.dir,"ref_points","Kobe_status.png"),
            to=paste0(to.dir,"/m",m,"_Kobe_status.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"results","SSB_F_trend.png"),
            to=paste0(to.dir,"/m",m,"_SSB_F_trend.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"results","Ecov_1.png"),
            to=paste0(to.dir,"/m",m,"_Ecov_1.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"retro","Fbar_retro_relative.png"),
            to=paste0(to.dir,"/m",m,"_Fbar_retro_relative.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"retro","NAA_retro_relative.png"),
            to=paste0(to.dir,"/m",m,"_NAA_retro_relative.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"retro","SSB_retro_relative.png"),
            to=paste0(to.dir,"/m",m,"_SSB_retro_relative.png"), overwrite = TRUE)
}

vign6_res <- df.mods
save(vign6_res, file=file.path(pkg.dir, "data","vign6_res.RData"))

ex6_test_results <- list(pars=NULL, nll=NULL)
ex6_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
ex6_test_results$nll <- sapply(mods, function(x) x$opt$objective)
saveRDS(ex6_test_results, file=file.path(pkg.dir, "inst","extdata","ex6_test_results.rds"))

file.copy(from=file.path(write.dir,"NAA_devs.png"), to=file.path(pkg.dir, "vignettes","ex6_plots","NAA_devs.png"), overwrite=T)

