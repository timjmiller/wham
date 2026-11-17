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
main.dir <- here("sandbox", "pkg_example_results", "ex06")
write.dir <- tempdir(check=TRUE) #will be passed to the ex5_M_GSI.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex6_NAA.R"))
setwd(pkg.dir)

to.dir <-here("vignettes", "ex06_files")
dir.create(to.dir)

# df.mods <- read.csv(file.path(write.dir,"ex6_table.csv"))
# mod.list <- file.path(write.dir, paste0(df.mods$Model,".rds"))
# mods <- lapply(mod.list, readRDS)

is_conv <- !not_conv
saveRDS(mods, file = file.path(main.dir, "mods.rds"))
# mods <- readRDS(file = file.path(main.dir, "mods.rds"))
for(m in which(is_conv)) file.copy(from = file.path(write.dir, paste0("m",m)), to = main.dir, recursive = TRUE, overwrite = TRUE, copy.date = TRUE)
file.copy(from = file.path(write.dir, "ex6_table.csv"), to = main.dir)
file.copy(from = file.path(write.dir, "NAA_devs.png"), to = main.dir, overwrite = T)
file.copy(from=file.path(write.dir,"NAA_devs.png"), to=to.dir, overwrite=T)

# copy to vignette plots folder
to.dir.m <- file.path(to.dir,paste0("m", c(1,7,13)))
from.dir <- file.path(main.dir, paste0("m", c(1,7,13)),"plots_png")

# for(m in c(1,7,13)){
#   from.dir <- paste0(write.dir,"/m",m,"/plots_png")
  file.copy(from=file.path(from.dir,"ref_points","Kobe_status.png"),
            to=paste0(to.dir.m,"_Kobe_status.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"results","SSB_F_trend.png"),
            to=paste0(to.dir.m,"_SSB_F_trend.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"results","Ecov_1_GSI.png"),
            to=paste0(to.dir.m,"_Ecov_1.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"retro","region_1_Fbar_retro_relative.png"),
            to=paste0(to.dir.m,"_Fbar_retro_relative.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"retro","stock_1_region_1_NAA_retro_relative.png"),
            to=paste0(to.dir.m,"_NAA_retro_relative.png"), overwrite = TRUE)
  file.copy(from=file.path(from.dir,"retro","stock_1_SSB_retro_relative.png"),
            to=paste0(to.dir.m,"_SSB_retro_relative.png"), overwrite = TRUE)
# }

res_table <- df.mods
thebest = c(11,13)
res_table[1,"dAIC"] = '---'
res_table[1,"AIC"] = '---'
res_table <- res_table |>
  dplyr::rename("Converged"="conv", "PD Hessian"="pdHess", "Runtime (min)"="runtime", "$\\Delta AIC$"="dAIC",
         "$\\rho_{R}$"="rho_R", "$\\rho_{SSB}$"="rho_SSB", "$\\rho_{\\overline{F}}$"="rho_Fbar") |>
  knitr::kable(escape = F) |>
  kableExtra::kable_styling(bootstrap_options = c("condensed","responsive")) |>
  kableExtra::row_spec(thebest, bold=TRUE)

saveRDS(res_table, file=file.path(to.dir,"res_table.RDS"))
#res_table <- readRDS(file=file.path(to.dir, "res_table.RDS"))
#res_table

# save results for testing
make_test_res <- FALSE
if(make_test_res){
  ex6_test_results <- list()
  ex6_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
  ex6_test_results$nll <- sapply(mods, function(x) x$opt$objective)
  ex6_test_results$is_conv <- vign6_res$pdHess
  saveRDS(ex6_test_results, file=file.path(pkg.dir, "tests","testthat","fixtures","ex06_test_results.RDS"))
}
