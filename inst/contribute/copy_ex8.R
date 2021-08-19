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

# copy to vignette plots folder
to.dir <- file.path(pkg.dir, "vignettes", "ex8_plots")
from.dir <- file.path(write.dir,"compare_png")
from.list <- list.files(from.dir,full.names=T)
to.list <- file.path(to.dir, list.files(from.dir))
file.copy(from=from.list, to=to.list, overwrite = TRUE)
file.copy(from=file.path(write.dir,"model_comparison.csv"), to=file.path(to.dir,"model_comparison.csv"), overwrite = TRUE)
file.copy(from=file.path(write.dir,"res.rds"), to=file.path(to.dir,"res.rds"), overwrite = TRUE)

res <- readRDS(file.path(write.dir,"res.rds"))
vign8_res <- round(res$tab,2)
save(vign8_res, file=file.path(pkg.dir, "data","vign8_res.RData"))
