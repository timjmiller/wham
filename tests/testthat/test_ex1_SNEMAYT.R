# Test WHAM example 1: Southern New England Mid-Atlantic Yellowtail Flounder
#   - Prepare wham, data, and model settings
#   - Fit several (slightly different) WHAM models:
#       m1: statistical catch-at-age (SCAA) model, but with recruitment estimated as random effects; multinomial age-compositions
#       m2: as m1, but with logistic normal age-compositions
#       m3: full state-space model (numbers at all ages are random effects), multinomial age-compositions
#       m4: full state-space model, logistic normal age-compositions
#   - Compare models by AIC and Mohn’s rho (retrospective analysis)
#   - Plots of input data, diagnostics, and results
context("Ex 1: SNEMA yellowtail")

test_that("Ex 1 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
ex1_test_results <- readRDS(file.path(path_to_examples,"ex1_test_results.rds"))

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder")

# Make one or more selectivity blocks with age-specific parameters
age.specific = 1:3 # 3 age-specific blocks
not.age.specific = (1:input$data$n_selblocks)[-age.specific]
input = set_age_sel0(input, age.specific)
input$par$logit_selpars[not.age.specific,c(1:input$data$n_ages,input$data$n_ages + 3:6)] = Inf
input$par$logit_selpars[1,5] = Inf
input$par$logit_selpars[2,4] = Inf
input$par$logit_selpars[3,2] = Inf
# Now redefine the map argument for the selectivity parameters to estimate only selectivity parameters without initial values at lower and upper bounds.
input$map$logit_selpars = matrix(input$map$logit_selpars, input$data$n_selblocks, input$data$n_ages + 6)
input$map$logit_selpars[is.infinite(input$par$logit_selpars)] = NA
input$map$logit_selpars[!is.infinite(input$par$logit_selpars)] = 1:sum(!is.infinite(input$par$logit_selpars))
input$map$logit_selpars = factor(input$map$logit_selpars)
base = input

#SCAA, but with random effects for recruitment
temp = base
temp$random = c(temp$random, "log_R")
temp$map = temp$map[!(names(temp$map) %in% c("log_R_sigma", "mean_rec_pars"))]
temp$data$random_recruitment = 1
# m1 <- fit_wham(temp, do.retro=F, do.osa=F)
m1 <- suppressWarnings(fit_wham(temp))

# Check that m1 converged
m1_check <- check_convergence(m1, ret=TRUE)
expect_equal(m1_check$convergence, 0) # opt$convergence should be 0
expect_false(m1_check$na_sdrep) # sdrep should succeed
expect_lt(m1_check$maxgr, 1e-6) # maximum gradient should be < 1e-06

# Check m1 parameter values
expect_equal(as.numeric(m1$opt$par), ex1_test_results$m1par, tolerance=1e-3)

#Like m1, but change age comp likelihoods to logistic normal
temp = base
temp$data$age_comp_model_indices = rep(7, temp$data$n_indices)
temp$data$age_comp_model_fleets = rep(7, temp$data$n_fleets)
temp$data$n_age_comp_pars_indices = rep(1, temp$data$n_indices)
temp$data$n_age_comp_pars_fleets = rep(1, temp$data$n_fleets)
temp$par$index_paa_pars = rep(0, temp$data$n_indices)
temp$par$catch_paa_pars = rep(0, temp$data$n_fleets)
temp$map = temp$map[!(names(temp$map) %in% c("index_paa_pars", "catch_paa_pars"))]
temp$random = c(temp$random, "log_R")
temp$map = temp$map[!(names(temp$map) %in% c("log_R_sigma", "mean_rec_pars"))]
temp$data$random_recruitment = 1
m2 <- suppressWarnings(fit_wham(temp))

# Check that m2 converged
m2_check <- check_convergence(m2, ret=TRUE)
expect_equal(m2_check$convergence, 0) # opt$convergence should be 0
expect_false(m2_check$na_sdrep) # sdrep should succeed
expect_lt(m2_check$maxgr, 1e-6) # maximum gradient should be < 1e-06

# Check m2 parameter values
expect_equal(as.numeric(m2$opt$par), ex1_test_results$m2par, tolerance=1e-3)

#full state-space model, abundance is the state vector
temp = base
temp$data$use_NAA_re = 1
temp$data$random_recruitment = 0
temp$map = temp$map[!(names(temp$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
temp$map$log_R = factor(rep(NA, length(temp$par$log_R)))
temp$random = c(temp$random, "log_NAA")
m3 <- suppressWarnings(fit_wham(temp))

# Check that m3 converged
m3_check <- check_convergence(m3, ret=TRUE)
expect_equal(m3_check$convergence, 0) # opt$convergence should be 0
expect_false(m3_check$na_sdrep) # sdrep should succeed
expect_lt(m3_check$maxgr, 1e-6) # maximum gradient should be < 1e-06

# Check m3 parameter values
expect_equal(as.numeric(m3$opt$par), ex1_test_results$m3par, tolerance=1e-3)

#Like m3, but change age comp likelihoods to logistic normal
temp = base
temp$data$age_comp_model_indices = rep(7, temp$data$n_indices)
temp$data$age_comp_model_fleets = rep(7, temp$data$n_fleets)
temp$data$n_age_comp_pars_indices = rep(1, temp$data$n_indices)
temp$data$n_age_comp_pars_fleets = rep(1, temp$data$n_fleets)
temp$par$index_paa_pars = rep(0, temp$data$n_indices)
temp$par$catch_paa_pars = rep(0, temp$data$n_fleets)
temp$map = temp$map[!(names(temp$map) %in% c("index_paa_pars", "catch_paa_pars"))]
temp$data$use_NAA_re = 1
temp$data$random_recruitment = 0
temp$map = temp$map[!(names(temp$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
temp$map$log_R = factor(rep(NA, length(temp$par$log_R)))
temp$random = c(temp$random, "log_NAA")
m4 <- suppressWarnings(fit_wham(temp))

# Check that m4 converged
m4_check <- check_convergence(m4, ret=TRUE)
expect_equal(m4_check$convergence, 0) # opt$convergence should be 0
expect_false(m4_check$na_sdrep) # sdrep should succeed
expect_lt(m4_check$maxgr, 1e-6) # maximum gradient should be < 1e-06

# Check m4 parameter values
expect_equal(as.numeric(m4$opt$par), ex1_test_results$m4par, tolerance=1e-3)

# Save list of all fit models
mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)

# Check neg-log-likelihoods are within 1e-6
nll <- sapply(mods, function(x) x$opt$obj)
expect_equal(nll, ex1_test_results$nll, tolerance=1e-6, scale=1)

# Compare models by AIC and Mohn's rho
tmp.dir <- tempdir(check=TRUE)
res <- compare_wham_models(mods, fname="model_comparison", sort=TRUE, fdir=tmp.dir)

# WHAM output plots for best model with projections
plot_wham_output(mod=m4, out.type='html', dir.main=tmp.dir)

# remove files created during testing
unlink(tmp.dir, recursive=TRUE)

# # save objects to test against in future
# ex1_test_results <- list(nll=nll,
#                          m1par=as.numeric(m1$opt$par),
#                          m2par=as.numeric(m2$opt$par),
#                          m3par=as.numeric(m3$opt$par),
#                          m4par=as.numeric(m4$opt$par))
# saveRDS(ex1_test_results, file="/home/bstock/Documents/wham/inst/extdata/ex1_test_results.rds")
# save("mods", file="/home/bstock/Documents/wham/sandbox/ex1/ex1_models.RData")
})
