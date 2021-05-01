# WHAM example 8: compare asap and wham models

# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex8_compare.R"); etime <- Sys.time(); runtime = etime - btime;
# 0.25 min

context("Ex 8: Compare")

test_that("Ex 8 works",{
path_to_examples <- system.file("extdata", package="wham")
asap.dir <- file.path(path_to_examples,"BASE_3")

asap3 <- read_asap3_dat(file.path(asap.dir,"BASE_3.DAT"))

# just m1 = fixed effect recruitment
input <- suppressWarnings(prepare_wham_input(asap3, recruit_model = 2, # match asap model, which does not estimate stock-recruitment relationship
                            model_name = "WHAM-m1"))   
mod <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F, MakeADFun.silent=TRUE))

base <- read_asap3_fit(wd=asap.dir, asap.name="BASE_3")
mods <- list(base, mod)
names(mods) <- c("ASAP","WHAM-m1")

tmp.dir <- tempdir(check=TRUE)
res <- suppressWarnings(compare_wham_models(mods, fdir=tmp.dir, do.table=F, plot.opts=list(kobe.prob=FALSE)))

})
