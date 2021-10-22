# Test WHAM example 10: Operating models and Management Strategy Evaluation
#   - Make a default wham input file without an ASAP3 dat file.
#   - Make an operating model, simulate data, and fit models
#   - Perform a simple management strategy evaluation using with an operating model with Beverton-Holt stock recruitment and a SCAA estimating model
#   - Make some plots

# library(wham)
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex10_om_mse.R"); etime <- Sys.time(); runtime = etime - btime;
# btime <- Sys.time(); testthat::test_file("~/work/wham/wham/tests/testthat/test_ex10_om_mse.R"); etime <- Sys.time(); runtime = etime - btime;
# ~4.2 min
# BCS I get 6.2 min

context("Ex 10: Operating models and Management Strategy Evaluation")

test_that("Ex 10 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
ex10_tests <- readRDS(file.path(path_to_examples,"ex10_tests.rds"))
tmp.dir <- tempdir(check=TRUE)

library(ggplot2)
library(tidyr)
library(dplyr)

#make a list of input components that prepare_wham_input can use to generate an input for fit_wham
make_digifish <- function(years = 1975:2014) {
    digifish = list()
    digifish$ages = 1:10
    digifish$years = years
    na = length(digifish$ages)
    ny = length(digifish$years)

    digifish$n_fleets = 1
    digifish$catch_cv = matrix(0.1, ny, digifish$n_fleets)
    digifish$catch_Neff = matrix(200, ny, digifish$n_fleets)
    digifish$n_indices = 1
    digifish$index_cv = matrix(0.3, ny, digifish$n_indices)
    digifish$index_Neff = matrix(100, ny, digifish$n_indices)
    digifish$fracyr_indices = matrix(0.5, ny, digifish$n_indices)
    digifish$index_units = rep(1, length(digifish$n_indices)) #biomass
    digifish$index_paa_units = rep(2, length(digifish$n_indices)) #abundance
    digifish$maturity = t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny))

    L = 100*(1-exp(-0.3*(1:na - 0)))
    W = exp(-11)*L^3
    nwaa = digifish$n_indices + digifish$n_fleets + 2
    digifish$waa = array(NA, dim = c(nwaa, ny, na))
    for(i in 1:nwaa) digifish$waa[i,,] = t(matrix(W, na, ny))

    digifish$fracyr_SSB = rep(0.25,ny)
    digifish$q = rep(0.3, digifish$n_indices)
    digifish$F = matrix(0.2,ny, digifish$n_fleets)

    digifish$selblock_pointer_fleets = t(matrix(1:digifish$n_fleets, digifish$n_fleets, ny))
    digifish$selblock_pointer_indices = t(matrix(digifish$n_fleets + 1:digifish$n_indices, digifish$n_indices, ny))
    return(digifish)
}
digifish = make_digifish()

selectivity = list(model = c(rep("logistic", digifish$n_fleets),rep("logistic", digifish$n_indices)),
    initial_pars = rep(list(c(5,1)), digifish$n_fleets + digifish$n_indices)) #fleet, index

M = list(initial_means = rep(0.2, length(digifish$ages)))

NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]))
NAA_re$sigma = "rec" #random about mean
NAA_re$use_steepness = 0
NAA_re$recruit_model = 2 #random effects with a constant mean
NAA_re$recruit_pars = exp(10)


#make input and operating model
input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M))
om = suppressWarnings(fit_wham(input, do.fit = FALSE))

#simulate data from operating model
set.seed(8675309)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#ex10_tests = list()
#fit estimating model that is the same as the operating model
fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
fit$mohns_rho = mohns_rho(fit) 
plot_wham_output(fit, dir.main=tmp.dir)
expect_equal(fit$opt$obj, ex10_tests$fit1$nll, tolerance=1e-6, scale=1)
expect_equal(fit$opt$par, ex10_tests$fit1$par, tolerance=1e-6, scale=1)
expect_equal(fit$mohns_rho, ex10_tests$fit1$mohns_rho, tolerance=1e-4, scale=1)
#ex10_tests$fit1 = list(nll = fit$opt$obj, par = fit$opt$par, mohns_rho = fit$mohns_rho)
#saveRDS(ex10_tests, "ex10_tests.rds")

#make input for a scaa estimation model (recruitment as fixed effects)
scaa_info = digifish
data_names = c("agg_catch", "catch_paa", "agg_indices","index_paa")
scaa_info[data_names] = newdata[data_names]

#the only thing that is different is how the parameterization of numbers at age is configured
scaa_NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]))
scaa_NAA_re$use_steepness = 0
scaa_NAA_re$recruit_model = 1 #recruitments as fixed effects

scaa_input = suppressWarnings(prepare_wham_input(basic_info = scaa_info, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M, recruit_model = 1))

#fit the scaa estimation model
scaa_fit = suppressWarnings(fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
scaa_fit$mohns_rho = mohns_rho(scaa_fit)
expect_equal(scaa_fit$opt$obj, ex10_tests$fit2$nll, tolerance=1e-6, scale=1)
expect_equal(scaa_fit$opt$par, ex10_tests$fit2$par, tolerance=1e-6, scale=1)
expect_equal(scaa_fit$mohns_rho, ex10_tests$fit2$mohns_rho, tolerance=1e-4, scale=1)
#ex10_tests$fit2 = list(nll = scaa_fit$opt$obj, par = scaa_fit$opt$par, mohns_rho = scaa_fit$mohns_rho)
#saveRDS(ex10_tests, "ex10_tests.rds")


##########################################################################################
# Make the operating model assume a Beverton-Holt stock recruit relationship
NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]))
NAA_re$sigma = "rec" #random about mean
NAA_re$use_steepness = 1 #ok because M, WAA, etc are constant
NAA_re$recruit_model = 3 #Beverton-Holt
NAA_re$recruit_pars = c(0.5, exp(10))

#make input object for operating model
bh_input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M))

#make the operating model
bh_om = suppressWarnings(fit_wham(bh_input, do.fit = FALSE))

set.seed(8675309)
sim_pop = bh_om$simulate(complete=TRUE)
temp = bh_input
temp$data = sim_pop

#fit estimating model that is the same as the operating model
bh_fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
bh_fit$mohns_rho = mohns_rho(bh_fit)
expect_equal(bh_fit$opt$obj, ex10_tests$fit3$nll, tolerance=1e-6, scale=1)
expect_equal(bh_fit$opt$par, ex10_tests$fit3$par, tolerance=1e-6, scale=1)
expect_equal(bh_fit$mohns_rho, ex10_tests$fit3$mohns_rho, tolerance=1e-4, scale=1)
#ex10_tests$fit3 = list(nll = bh_fit$opt$obj, par = bh_fit$opt$par, mohns_rho = bh_fit$mohns_rho)
#saveRDS(ex10_tests, "ex10_tests.rds")

#make input for a scaa estimation model (recruitment as fixed effects)
scaa_info = digifish
data_names = c("agg_catch", "catch_paa", "agg_indices","index_paa")
scaa_info[data_names] = sim_pop[data_names]

scaa_NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]))
scaa_NAA_re$use_steepness = 0
scaa_NAA_re$recruit_model = 1 #recruitments as fixed effects

scaa_input = suppressWarnings(prepare_wham_input(basic_info = scaa_info, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M, recruit_model = 1))

#fit the scaa estimation model
scaa_fit = suppressWarnings(fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
scaa_fit$mohns_rho = mohns_rho(scaa_fit)
expect_equal(scaa_fit$opt$obj, ex10_tests$fit4$nll, tolerance=1e-6, scale=1)
expect_equal(scaa_fit$opt$par, ex10_tests$fit4$par, tolerance=1e-6, scale=1)
expect_equal(scaa_fit$mohns_rho, ex10_tests$fit4$mohns_rho, tolerance=1e-4, scale=1)
#ex10_tests$fit4 = list(nll = scaa_fit$opt$obj, par = scaa_fit$opt$par, mohns_rho = scaa_fit$mohns_rho)
#saveRDS(ex10_tests, "ex10_tests.rds")


######################################
#now do a closed-loop simulation with the beverton-holt operating model where catch advice is set during the projection years.
#catch advice is set as the average catch over 5 years projected using F40

#first we use the prepare_projection function to set up an input that includes projection years.
#setting F in the projection years at the last value of the base period (0.2)
om_input <- suppressWarnings(prepare_projection(bh_om, proj.opts=list(n.yrs=39, use.last.F=FALSE, use.avg.F=FALSE,
                use.FXSPR=FALSE, proj.F=NULL, avg.yrs=NULL, proj.catch = rep(1,39),
                cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)))
#then we use fit_wham without fitting to set up a TMB model that can be used to simulate a population time series
# and associated catch and index observations with the stochastic assumptions made in the NAA_re options to generate bh_om
temp <- suppressWarnings(fit_wham(om_input, n.newton=n.newton, do.sdrep=F, do.retro=F, do.osa=F, do.check=F, do.proj=F, 
  MakeADFun.silent = TRUE, save.sdrep=FALSE, do.fit = F))

#First we generate a time series of the population with the stochastic assumptions made in the NAA_re options to generate bh_om
#a key component is to keep the same seed for each simulation of the MSE
set.seed(8675309)
pop_om_base_period = temp$simulate(complete=TRUE)
#now replace the realized abundance parameters in the operating model input
om_input$par$log_NAA = pop_om_base_period$log_NAA
# and reset the operating model
om <- suppressWarnings(fit_wham(om_input, n.newton=n.newton, do.sdrep=F, do.retro=F, do.osa=F, do.check=F, do.proj=F, 
  MakeADFun.silent = TRUE, save.sdrep=FALSE, do.fit = F))
plot(om$years_full, log(om$rep$NAA[,1]), type = 'l', xlab = "Year", ylab = "log(Recruits (1000s))")
lines(om$years_full, log(temp$rep$NAA[,1]), col = "red")
#the "temp" version of the operating model (red) just has the initial values of numbers at age (exp(10))

#Define how often an assessment (SCAA) and catch advice will be made
assess.interval = 3 #years.step = integer()
#catch advice
advice = rep(1,39)

#the first advice model will be completed at the end of the base period
#the SCAA model defined as above will be updated every 3 years (assess.interval) of the projection/evaluation period
scaa_step = suppressWarnings(fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
#the SCAA model with actual projections to make catch advice using F40
scaa_step_proj <- suppressWarnings(project_wham(scaa_step, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
            use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL, avg.rec.yrs = tail(scaa_step$years,30),
            cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE, MakeADFun.silent=TRUE))
NAAtrack = matrix(NA, 79, 13)
#making advice for the first 3 years of the evaluation period
advice[1:assess.interval] = rep(mean(scaa_step_proj$rep$pred_catch[scaa_step_proj$input$data$n_years_model + 1:5]), assess.interval)
for(y in seq(3,39,3)) {
    #set up projection of operating model
    om_input <- suppressWarnings(prepare_projection(bh_om, proj.opts=list(n.yrs=39, use.last.F=FALSE, use.avg.F=FALSE,
                use.FXSPR=FALSE, proj.F=NULL, proj.catch=advice, avg.yrs=NULL,
                cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)))
  temp <- suppressWarnings(fit_wham(om_input, n.newton=n.newton, do.sdrep=F, do.retro=F, do.osa=F, do.check=F, do.proj=F, 
    MakeADFun.silent = TRUE, save.sdrep=FALSE, do.fit = F))
    set.seed(8675309)
    updated_sim = temp$simulate(complete=TRUE)
    om_input$par$log_NAA = updated_sim$log_NAA
    NAAtrack[,y/3] = updated_sim$NAA[,1]
    lines(om$years_full, log(updated_sim$NAA[,1]), col = viridisLite::viridis(n= 39)[y], lty = 2)

    #increase terminal year
    scaa_info_step = make_digifish(min(digifish$years):(max(digifish$years)+y))
    scaa_info_step$agg_catch = rbind(updated_sim$agg_catch, updated_sim$agg_catch_proj[1:y,,drop=F])
    scaa_info_step$agg_indices = rbind(updated_sim$agg_indices, updated_sim$agg_indices_proj[1:y,,drop=F])
    scaa_info_step$catch_paa = abind::abind(updated_sim$catch_paa, updated_sim$catch_paa_proj[,1:y,,drop=F], along = 2)
    scaa_info_step$index_paa = abind::abind(updated_sim$index_paa, updated_sim$index_paa_proj[,1:y,,drop=F], along = 2)

    step_input = suppressWarnings(prepare_wham_input(basic_info = scaa_info_step, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M))
    scaa_step = suppressWarnings(fit_wham(step_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
    scaa_step_proj <- suppressWarnings(project_wham(scaa_step, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
                use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL, avg.rec.yrs = tail(scaa_step$years,30),
                cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE, MakeADFun.silent=TRUE))

    advice[y + 1:assess.interval] = rep(mean(scaa_step_proj$rep$pred_catch[scaa_step_proj$input$data$n_years_model + 1:5]), assess.interval)

}

expect_equal(advice, ex10_tests$advice, tolerance=1e-6, scale=1)
expect_equal(NAAtrack, ex10_tests$NAAtrack, tolerance=1e-6, scale=1)

#ex10_tests$advice = advice
#ex10_tests$NAAtrack = NAAtrack
#saveRDS(ex10_tests, "ex10_tests.rds")

#png()
par(mfrow = c(2,2))
pal = viridisLite::viridis(n=3)
plot(temp$years_full, log(pop_om_base_period$NAA[,1]), type = 'l', xlab = "Year", ylab = "log(recruits)", col = pal[1]) 
lines(temp$years_full, log(NAAtrack[,13]), col = pal[2]) 
abline(v = max(temp$years), col = pal[3])

plot(temp$years_full, pop_om_base_period$SSB, type = 'l', xlab = "Year", ylab = "SSB", ylim = c(0,450000), col  = pal[1]) 
lines(temp$years_full, updated_sim$SSB, col = pal[2]) 
abline(v = max(temp$years), col = pal[3])

plot(temp$years_full, pop_om_base_period$pred_catch, type = 'l', xlab = "Year", ylab = "Catch (mt)", ylim = c(0,80000), col = pal[1]) 
lines(temp$years_full, updated_sim$pred_catch, col = pal[2]) 
points(temp$years_full, scaa_info_step$agg_catch, col = pal[2])
abline(v = max(temp$years), col = pal[3])

pal = viridisLite::viridis(n=2)
plot(temp$years_full, updated_sim$FAA_tot[,10]/exp(updated_sim$log_FMSY), type = 'l', xlab = "Year", ylab = "F/FMSY", col = pal[1])
abline(v = max(temp$years), col = pal[2])
abline(h = 1, col = pal[2])
dev.off()

})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

