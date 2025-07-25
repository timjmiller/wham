# WHAM example 3: projections
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# To create test results see file.path(system.file("contribute", package="wham"), "copy_ex3.R"), note whether bias-correction used
# pkgbuild::compile_dll(debug = FALSE)
# pkgload::load_all()
# btime <- Sys.time(); devtools::test(filter = "ex03_projections"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~4 min

context("Ex 3: Projections")

test_that("Ex 3 works",{

suppressWarnings(suppressMessages({

# get results to check NLL and par estimates
tmp.dir <- tempdir(check=TRUE)
path_to_examples <- system.file("extdata", package="wham")
ex3_test_results <- readRDS(file.path(path_to_examples,"ex3_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"CPI.csv"), header=T)

# specify model: AR1 CPI, limiting effect on Bev-Holt
env <- list(
  label = "CPI",
  mean = as.matrix(env.dat$CPI), # CPI observations
  logsigma = as.matrix(log(env.dat$CPI_sigma)), # CPI standard error is given/fixed as data
  year = env.dat$Year,
  use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
  process_model = "ar1", # fit CPI as AR1 process
  recruitment_how = matrix("limiting-lag-1-linear",1,1)) # limiting (carrying capacity)

basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions
#basic_info <- list() #compare to previous versions

input <- prepare_wham_input(asap3, recruit_model = 3,
                            model_name = "Ex 3: Projections",
                            ecov = env,
                            NAA_re = list(sigma="rec+1", cor="iid"),
                            age_comp = "logistic-normal-pool0", # logistic normal pool 0 obs
                            basic_info = basic_info)

# selectivity = logistic, not age-specific
#   2 pars per block instead of n.ages
#   sel pars of indices 4/5 fixed at 1.5, 0.1 (neg phase in .dat file)
input$par$logit_selpars[1:4,7:8] <- 0 # original code started selpars at 0 (last 2 rows are fixed)

# ---------------------------------------------------------
## Fit model without projections
mod <- fit_wham(input, do.fit=F, MakeADFun.silent = TRUE)
expect_equal(length(mod$par), length(ex3_test_results$par), tolerance=1e-3) # parameter values
expect_equal(as.numeric(mod$fn(ex3_test_results$par)), ex3_test_results$nll, tolerance=1e-6) # nll
mod$par <- ex3_test_results$par
# mod$fn(mod$par)
input$par <- mod$env$parList(mod$par)

#mod <- fit_wham(input, do.proj=F, do.osa=F, do.retro=F)
mod <- fit_wham(input, do.proj=F, do.osa=F, do.retro=F, MakeADFun.silent = TRUE)

# Add projections to previously fit model
proj_opts <- list()
mod_proj <- list()
nll_proj <- numeric()

# default settings: 3 years, use last F, continue ecov
proj_opts[[1]] <- list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# 5 years, use last F, average ecov 1992-1996
proj_opts[[2]] <- list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=1992:1996, proj.ecov=NULL)

# 5 years, use last F, use last ecov
proj_opts[[3]] <- list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=TRUE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# 5 years, use last F, specify high CPI ~ 0.5
proj_opts[[4]] <- list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=matrix(c(0.5,0.7,0.4,0.5,0.55),ncol=1))

# 5 years, use last F, specify low CPI ~ -1.5
proj_opts[[5]] <- list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=matrix(c(-1.6,-1.3,-1,-1.2,-1.25),ncol=1))

# specify catch, 5 years
proj_opts[[6]] <- list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=c(10, 2000, 1000, 3000, 20), avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# specify F, 5 years
proj_opts[[7]] <- list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=c(0.001, 1, 0.5, .1, .2), proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# use FXSPR (avg.yrs defaults to last 5 years, 2007-2011), 5 years
proj_opts[[8]] <- list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# use avg F (avg.yrs defaults to last 5 years, 2007-2011), 3 years
proj_opts[[9]] <- list(n.yrs=3, use.last.F=FALSE, use.avg.F=TRUE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# use avg F 1992-1996, 10 years
proj_opts[[10]] <- list(n.yrs=10, use.last.F=FALSE, use.avg.F=TRUE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=1992:1996,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# use FMSY (avg.yrs defaults to last 5 years, 2007-2011), 5 years
proj_opts[[11]] <- list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FMSY=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL)

# default settings: 3 years, use last F, continue ecov, but projected recruitment deviations=0 and average NAA deviations (1992-1996) projected
proj_opts[[12]] <- list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, proj_R_opt = 4, proj_NAA_opt = 2, avg.yrs.NAA = list(list(1992:1996)))

# default settings: 3 years, use last F, continue ecov, but projected recruitment deviations are averged over 1992-1996 and no NAA deviations projected
proj_opts[[13]] <- list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, proj_R_opt = 3, proj_NAA_opt = 3, avg.yrs.R = list(1992:1996))

# default settings: 3 years, use last F, continue ecov, but projected recruitment converges to that used for 40%SPR BRPS and no NAA deviations projected
proj_opts[[14]] <- list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, proj_R_opt = 2, proj_NAA_opt = 3, avg.yrs.R = list(1992:1996))

# saveRDS(proj_opts, "c:/work/wham/wham/sandbox/proj_opts.RDS")
# saveRDS(mod, "c:/work/wham/wham/sandbox/ex3_mod.RDS")
# proj_opts <- readRDS("c:/work/wham/wham/sandbox/proj_opts.RDS")
# mod <- readRDS("c:/work/wham/wham/sandbox/ex3_mod.RDS")

# proj_opts[[13]]$avg.yrs.R[[1]]
# temp <- prepare_projection(mod,proj.opts = proj_opts[[12]])
# temp$data$proj_NAA_opt
# temp$data$proj_R_opt
# temp$data$NAA_re_model
# x <- temp$par$log_NAA
# x[] <- as.integer(temp$map$log_NAA)
# x[1,1,,]
#   x <- project_wham(mod, proj.opts=proj_opts[[1]], do.sdrep= F, MakeADFun.silent = FALSE)
# x$rep$NAA_devs[1,1,,]
# x$rep$NAA_devs[1,1,,]

###########################
# pkgbuild::compile_dll(debug = FALSE)
# pkgload::load_all()
# proj_opts <- readRDS("c:/work/wham/wham/sandbox/proj_opts.RDS")
# mod <- readRDS("c:/work/wham/wham/sandbox/ex3_mod.RDS")
#   m = 7
# mod_proj <- list()
#   # temp <- prepare_projection(mod, proj.opts=proj_opts[[m]])
#   mod_proj[[m]] <- project_wham(mod, proj.opts=proj_opts[[m]], do.sdrep= F, MakeADFun.silent = FALSE)
#   mod_proj[[m]]$input$data$proj_F_opt
#   # temp <- mod$simulate(complete = TRUE)
#   #test simulation works with projections included.
#   temp <- mod_proj[[m]]$simulate(complete=TRUE)

#   mod_proj[[m]]$rep$log_F_tot
#####################################

for(m in 1:length(proj_opts)) {

  mod_proj[[m]] <- project_wham(mod, proj.opts=proj_opts[[m]], do.sdrep= F, MakeADFun.silent = TRUE)
  
  # The !! allows the individual elements in the report to be be seen if there is an error. See ?testthat::quasi_label and example
  #  check length of fixed effects
  if(m<12){
    expect_equal(length(mod_proj[[!!m]]$par), length(mod$opt$par), tolerance=1e-6)
    #  check fixed effects are the same
    for(p in 1:length(mod_proj[[!!m]]$par)) expect_equal(mod_proj[[!!m]]$par[!!p], mod$opt$par[!!p], tolerance=1e-3) # parameter values
    #  check marginal nll is the same
    nll_proj[m] <-  mod_proj[[m]]$fn()
    expect_equal(as.numeric(nll_proj[!!m]), as.numeric(mod$opt$obj), tolerance=1e-6)
  }
  proj_yrs <-mod_proj[[m]]$input$data$n_years_model + 1:mod_proj[[m]]$input$data$n_years_proj 
  if(m == 12){
    no_devs <- mod_proj[[m]]$rep$NAA_devs[1,1,proj_yrs,1]
    avg_devs <- mod_proj[[m]]$rep$NAA_devs[1,1,proj_yrs,-1]
    for(i in 1:length(no_devs)) expect_equal(no_devs[!!i], 0.0, tolerance=1e-8)
    for(a in 1:NCOL(avg_devs)) for(i in 2:NROW(avg_devs)) expect_equal(avg_devs[!!i,!!a], avg_devs[1,a], tolerance=1e-8)
  }
  if(m == 13){
    no_devs <- mod_proj[[m]]$rep$NAA_devs[1,1,proj_yrs,-1]
    avg_devs <- mod_proj[[m]]$rep$NAA_devs[1,1,proj_yrs,1]
    for(i in 2:length(avg_devs)) expect_equal(avg_devs[!!i], avg_devs[1], tolerance=1e-8)
    for(a in 1:NCOL(no_devs)) for(i in 1:NROW(no_devs)) expect_equal(no_devs[!!i,!!a], 0.0, tolerance=1e-8)
  }
  if(m == 14){
    no_devs <- mod_proj[[m]]$rep$NAA_devs[1,1,proj_yrs,-1]
    avg_devs <- mod_proj[[m]]$rep$NAA_devs[1,1,proj_yrs,1]
    for(i in 1:length(avg_devs)) expect_equal(avg_devs[!!i], -0.2300036, tolerance=1e-6)
    for(a in 1:NCOL(no_devs)) for(i in 1:NROW(no_devs)) expect_equal(no_devs[!!i,!!a], 0.0, tolerance=1e-8)
  }
  
  #test simulation works with projections included.
  temp <- mod_proj[[m]]$simulate(complete=TRUE)

  # plot results
  plot_wham_output(mod_proj[[m]], dir.main=file.path(tmp.dir,paste0("proj_",m)), plot.opts = list(browse=FALSE))
}

}))
})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))
