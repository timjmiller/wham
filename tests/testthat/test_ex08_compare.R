# WHAM example 8: compare asap and wham models

# pkgbuild::compile_dll(debug = FALSE)
# pkgload::load_all()
# btime <- Sys.time(); devtools::test(filter = "ex08_compare"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~15 sec

context("Ex 8: Compare")

test_that("Ex 8 works",{

suppressWarnings(suppressMessages({

path_to_extdata <- system.file("extdata", package="wham")
write.dir <- tempdir(check=TRUE)
rds_files <- paste0(c("RDAT", "DAT", "STD", "COR", "PAR"),".RDS")
asap_files_rds <- lapply(file.path(path_to_extdata, "BASE_3", rds_files), readRDS)
names(asap_files_rds) <- rds_files
dir.create(file.path(write.dir, "BASE_3"))
dput(asap_files_rds[[1]], file.path(write.dir, "BASE_3", "BASE_3.RDAT"))
writeLines(asap_files_rds[[2]], file.path(write.dir,"BASE_3", "BASE_3.DAT"))
writeLines(asap_files_rds[[3]], file.path(write.dir,"BASE_3", "BASE_3.STD"))
writeLines(asap_files_rds[[4]], file.path(write.dir,"BASE_3", "BASE_3.COR"))
writeLines(asap_files_rds[[5]], file.path(write.dir, "BASE_3", "BASE_3.PAR"))

asap3 <- read_asap3_dat(file.path(write.dir, "BASE_3", "BASE_3.DAT"))

# just m1 = fixed effect recruitment
input <- prepare_wham_input(asap3, recruit_model = 2, # match asap model, which does not estimate stock-recruitment relationship
                            model_name = "WHAM-m1")
mod <- fit_wham(input, do.retro=F, do.osa=F, MakeADFun.silent=TRUE)

base <- read_asap3_fit(wd=file.path(write.dir,"BASE_3"), asap.name="BASE_3")
all.mods <- list(base, mod)
names(all.mods) <- c("ASAP","WHAM-m1")

tmp.dir <- tempdir(check=TRUE)

res <- compare_wham_models(all.mods, fdir=tmp.dir, do.table=T, table.opts = list(print = FALSE), plot.opts=list(kobe.prob=FALSE, browse = FALSE))

}))

})
