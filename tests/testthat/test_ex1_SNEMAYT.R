# Test WHAM example 1: Southern New England Mid-Atlantic Yellowtail Flounder
#   - Prepare wham, data, and model settings
#   - Fit several (slightly different) WHAM models:
#       m1: statistical catch-at-age (SCAA) model, but with recruitment estimated as random effects; multinomial age-compositions
#       m2: as m1, but with logistic normal age-compositions
#       m3: full state-space model (numbers at all ages are random effects), multinomial age-compositions
#       m4: full state-space model, logistic normal age-compositions
#   - Compare models by AIC and Mohn’s rho (retrospective analysis)
#   - Plots of input data, diagnostics, and results

# library(wham)
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex1_SNEMAYT.R"); etime <- Sys.time(); runtime = etime - btime;
# 5.5 min

context("Ex 1: SNEMA yellowtail")

test_that("Ex 1 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
ex1_test_results <- readRDS(file.path(path_to_examples,"ex1_test_results.rds"))

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

input1 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
	                            selectivity=list(model=rep("age-specific",3), 
                                	re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                	fix_pars=list(4:5,4,2:4)),
	                            NAA_re = list(sigma="rec", cor="iid"))
m1 <- suppressWarnings(fit_wham(input1, do.osa = F, MakeADFun.silent=TRUE)) # turn off OSA residuals to save time

#Like m1, but change age comp likelihoods to logistic normal
input2 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                                    selectivity=list(model=rep("age-specific",3), 
                                        re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                        fix_pars=list(4:5,4,2:4)),
                                    NAA_re = list(sigma="rec", cor="iid"),
                                    age_comp = "logistic-normal-miss0")
m2 <- suppressWarnings(fit_wham(input2, do.osa = F, MakeADFun.silent=TRUE)) # turn off OSA residuals to save time

#full state-space model, abundance is the state vector
input3 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
	                            selectivity=list(model=rep("age-specific",3), 
                                	re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                	fix_pars=list(4:5,4,2:4)),
	                            NAA_re = list(sigma="rec+1", cor="iid"))
m3 <- suppressWarnings(fit_wham(input3, do.osa = F, MakeADFun.silent=TRUE)) # turn off OSA residuals to save time

#Like m3, but change age comp likelihoods to logistic normal
input4 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                                    selectivity=list(model=rep("age-specific",3), 
                                        re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                        fix_pars=list(4:5,4,2:4)),
                                    NAA_re = list(sigma="rec+1", cor="iid"),
                                    age_comp = "logistic-normal-miss0")
m4 <- suppressWarnings(fit_wham(input4, do.osa = F, MakeADFun.silent=TRUE)) # turn off OSA residuals to save time

# Save list of all fit models
mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)

# Check neg-log-likelihoods are within 1e-6
nll <- sapply(mods, function(x) x$opt$obj)
expect_equal(nll, ex1_test_results$nll, tolerance=1e-6, scale=1)

# Compare models by AIC and Mohn's rho
tmp.dir <- tempdir(check=TRUE)
res <- compare_wham_models(mods, fdir=tmp.dir, table.opts=list(fname="ex1_table", sort=TRUE, print=FALSE))

# WHAM output plots for best model with projections
m4_proj <- project_wham(model=mods$m4, MakeADFun.silent=TRUE)
plot_wham_output(mod=m4_proj, out.type='html', dir.main=tmp.dir)

# hard to see which model fails bc they're indexed by m
# print out each one by one

# Check that m1 converged
m1_check <- check_convergence(m1, ret=TRUE)
expect_equal(m1_check$convergence, 0) # opt$convergence should be 0
expect_false(m1_check$na_sdrep) # sdrep should succeed
expect_lt(m1_check$maxgr, 1e-5) # maximum gradient should be < 1e-06
# Check m1 parameter values
# order of logit_selpars changed when modifying prepare_wham_input for time-varying selectivity
expect_equal(as.numeric(m1$opt$par), ex1_test_results$par[[1]], tolerance=1e-1)

# Check that m2 converged
m2_check <- check_convergence(m2, ret=TRUE)
expect_equal(m2_check$convergence, 0) # opt$convergence should be 0
expect_false(m2_check$na_sdrep) # sdrep should succeed
expect_lt(m2_check$maxgr, 1e-5) # maximum gradient should be < 1e-06
# Check m2 parameter values
expect_equal(as.numeric(m2$opt$par), ex1_test_results$par[[2]], tolerance=1e-1)

# Check that m3 converged
m3_check <- check_convergence(m3, ret=TRUE)
expect_equal(m3_check$convergence, 0) # opt$convergence should be 0
expect_false(m3_check$na_sdrep) # sdrep should succeed
expect_lt(m3_check$maxgr, 1e-4) # maximum gradient should be < 1e-06
# Check m3 parameter values
expect_equal(as.numeric(m3$opt$par), ex1_test_results$par[[3]], tolerance=1e-1)

# Check that m4 converged
m4_check <- check_convergence(m4, ret=TRUE)
expect_equal(m4_check$convergence, 0) # opt$convergence should be 0
expect_false(m4_check$na_sdrep) # sdrep should succeed
expect_lt(m4_check$maxgr, 1e-5) # maximum gradient should be < 1e-06
# Check m4 parameter values
expect_equal(as.numeric(m4$opt$par), ex1_test_results$par[[4]], tolerance=1e-1)

})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

