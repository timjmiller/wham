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
main.dir <- here("sandbox", "pkg_example_results", "ex05")
write.dir <- tempdir(check=TRUE) #will be passed to the ex5_M_GSI.R script to set working directory and for writing all results
path_to_scripts <- system.file("example_scripts", package="wham")

#do bias-correction results?
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

#should write everything to write.dir and create any R objects in the session
source(file.path(path_to_scripts, "ex5_M_GSI.R"))
setwd(pkg.dir)
to.dir <- here("vignettes", "ex05_files")
dir.create(to.dir)

#df.mods <- read.csv(file.path(write.dir,"ex5_table.csv"))
#mod.list <- file.path(write.dir, paste0(df.mods$Model,".rds"))
#mods <- lapply(mod.list, readRDS)

res_table <- df.mods
posdef <- which(res_table$pdHess == TRUE)
thebest <- c(12)
res_table$Ecov_how [grep("linear", res_table$Ecov_how)] <- "linear"
res_table$Ecov_how [grep("poly-2", res_table$Ecov_how)] <- "poly-2"
res_table <- res_table |>
  dplyr::select(!(c(age_specific,mean_model))) |>
  dplyr::rename("Mean model"="M_model", "RE model" = "M_re", "GSI model"="Ecov_process", "GSI link"="Ecov_how",
         "Converged"="conv", "PD Hessian"="pdHess", "Runtime (min)"="runtime",
         "$\\rho_{R}$"="rho_R", "$\\rho_{SSB}$"="rho_SSB", "$\\rho_{\\overline{F}}$"="rho_Fbar") |>
  knitr::kable(escape = F) |>
  kableExtra::kable_styling(bootstrap_options = c("condensed","responsive")) |>
  kableExtra::row_spec(posdef, background = gray.colors(10,end=0.95)[10]) |>
  kableExtra::row_spec(thebest, bold=TRUE)
saveRDS(res_table, file = file.path(main.dir, "res_table.rds"))
saveRDS(res_table, file = file.path(to.dir, "res_table.rds"))
#res_table <- readRDS(file.path(main.dir, "res_table.rds"))
saveRDS(mods, file = file.path(main.dir, "mods.rds"))
saveRDS(mods_proj, file = file.path(main.dir, "mods_proj.rds"))
for(m in which(is_conv)) file.copy(from = file.path(write.dir, paste0("m",m)), to = main.dir, recursive = TRUE, overwrite = TRUE, copy.date = TRUE)
file.copy(from = file.path(write.dir, "MAA.png"), to = main.dir, overwrite = T)
file.copy(from = file.path(write.dir, "MAA.png"), to = to.dir, overwrite = T)

to.dir.m <- file.path(to.dir, paste0("m", c(1,2,12)))
from.dir <- file.path(main.dir, paste0("m", c(1,2,12)),"plots_png")
# plot output for models 1, 2, 8
# copy to vignette plots folder

file.copy(from=file.path(from.dir,"ref_points","FSPR_absolute.png"), 
          to=file.path(paste0(to.dir.m,"_FSPR_absolute.png")), overwrite = TRUE, copy.date = TRUE)
file.copy(from=file.path(from.dir,"ref_points","FSPR_relative.png"), 
          to=file.path(paste0(to.dir.m,"_FSPR_relative.png")), overwrite = TRUE)
file.copy(from=file.path(from.dir,"ref_points","Kobe_status.png"), 
          to=file.path(paste0(to.dir.m,"_Kobe_status.png")), overwrite = TRUE)
file.copy(from=file.path(from.dir,"results","SSB_F_trend.png"), 
          to=paste0(to.dir.m,"_SSB_F_trend.png"), overwrite = TRUE)
file.copy(from=file.path(from.dir,"retro","region_1_Fbar_retro_relative.png"), 
          to=paste0(to.dir.m,"_Fbar_retro_relative.png"), overwrite = TRUE)
file.copy(from=file.path(from.dir,"retro","stock_1_region_1_NAA_retro_relative.png"), 
          to=paste0(to.dir.m,"_NAA_retro_relative.png"), overwrite = TRUE)
file.copy(from=file.path(from.dir,"retro","stock_1_SSB_retro_relative.png"), 
          to=paste0(to.dir.m,"_SSB_retro_relative.png"), overwrite = TRUE)
file.copy(from=file.path(from.dir,"retro","stock_1_region_1_NAA_age_1_retro_relative.png"), 
      to=paste0(to.dir.m,"_NAA_age1_retro_relative.png"), overwrite = TRUE)

# vign5_res <- df.mods
# save(vign5_res, file=file.path(pkg.dir, "data","vign5_res.RData"))

# save results for testing
ex5_test_results <- list(pars=NULL, nll=NULL)
ex5_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
ex5_test_results$nll <- sapply(mods, function(x) x$opt$objective)
ex5_test_results$is_conv <- is_conv
dir.create(file.path(pkg.dir, "tests", "testthat", "fixtures"))
saveRDS(ex5_test_results, file=file.path(pkg.dir, "tests", "testthat", "fixtures","ex05_test_results.RDS"))
