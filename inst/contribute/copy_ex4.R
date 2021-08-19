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
write.dir <- file.path(main.dir,"ex4")

mod.list <- file.path(write.dir,paste0("m",1:9,".rds"))
mods <- lapply(mod.list, readRDS)
vign4_conv <- lapply(mods, function(x) capture.output(check_convergence(x)))
df.mods <- read.csv(file.path(write.dir,"ex4_table.csv"))

save(vign4_conv, file=file.path(pkg.dir, "data","vign4_conv.RData"))
vign4_res <- df.mods
save(vign4_res, file=file.path(pkg.dir, "data","vign4_res.RData"))

selAA <- lapply(mods, function(x) x$report()$selAA[[1]])
vign4_selAA <- selAA
save(vign4_selAA, file=file.path(pkg.dir, "data","vign4_selAA.RData"))

# copy selAA.png to vignette plots folder
file.copy(from=file.path(write.dir,"selAA.png"), to=file.path(pkg.dir, "vignettes","ex4_plots","selAA.png"), overwrite=T)

# save results for testing
ex4_test_results <- list(pars=NULL, nll=NULL)
ex4_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
ex4_test_results$nll <- sapply(mods, function(x) x$opt$objective)
saveRDS(ex4_test_results, file=file.path(pkg.dir, "inst","extdata","ex4_test_results.rds"))

