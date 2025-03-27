#  Ecov and age-year effects on natural mortality
# To create test results see file.path(system.file("contribute", package="wham"), "copy_ex5.R")
# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all()
# btime <- Sys.time(); devtools::test(filter = "alt_ar1_loglike"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~20 sec

context("Alt. AR1 likelihoods")

test_that("Alt. AR1 likelihoods work",{

suppressWarnings(suppressMessages({

path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

cor_vals <- -1 + 2/(1 + exp(-c(1,2,1.5)))

NAA_re <- list(sigma="rec+1", cor="2dar1", cor_vals = array(cor_vals, dim = c(1,1,3)), decouple_recruitment = FALSE)
input <- prepare_wham_input(asap3, recruit_model = 2,
                            selectivity=list(model=rep("age-specific",3), re=c("none","none","none"), 
                              initial_pars=list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,1,0.5),c(0.5,1,1,1,1,1)), 
                              fix_pars=list(4:6,4:5,2:6)),
                            NAA_re = NAA_re)
NAA_re <- list(sigma="rec+1", cor="2dar1", cor_vals = array(cor_vals, dim = c(1,1,3)), decouple_recruitment = TRUE)
input_dr <- prepare_wham_input(asap3, recruit_model = 2,
                            selectivity=list(model=rep("age-specific",3), re=c("none","none","none"), 
                              initial_pars=list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,1,0.5),c(0.5,1,1,1,1,1)), 
                              fix_pars=list(4:6,4:5,2:6)),
                            NAA_re = NAA_re)

input_dr$par$logit_selpars[2,5] <- input$par$logit_selpars[2,5] <- 0 # original code started selpars at 0 (last 2 rows are fixed)

input_alt <- input 
input_alt_dr <- input_alt_dr
input_alt$data$use_alt_AR1 <- input_alt_dr$data$use_alt_AR1 <- 1

mod <- fit_wham(input, do.fit = FALSE, MakeADFun.silent=TRUE)
mod_alt <- fit_wham(input_alt, do.fit = FALSE, MakeADFun.silent=TRUE)
mod_dr <- fit_wham(input_dr, do.fit = FALSE, MakeADFun.silent=TRUE)
mod_alt_dr <- fit_wham(input_alt_dr, do.fit = FALSE, MakeADFun.silent=TRUE)

input$data$do_simulate_period[1] <- 0
fit <- fit_wham(input, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent=TRUE)
proj <- project_wham(fit, do.sdrep = FALSE, MakeADFun.silent=TRUE)
}))

expect_equal(as.numeric(mod$fn()), as.numeric(mod_alt$fn()), tolerance=1e-6) # nll
expect_equal(as.numeric(mod$fn() - mod_dr$fn()), -2.1916, tolerance = 1e-3) # nll
expect_equal(as.numeric(mod_dr$fn()), as.numeric(mod_alt_dr$fn()), tolerance=1e-6) # nll

})
