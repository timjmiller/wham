# WHAM example 3: projections
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# btime <- Sys.time(); devtools::test("/home/bstock/Documents/wham"); etime <- Sys.time(); runtime = etime - btime;
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex3_projections.R"); etime <- Sys.time(); runtime = etime - btime;
# 5.7 min

context("Ex 3: Projections")

test_that("Ex 3 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
# ex3_test_results <- readRDS(file.path(path_to_examples,"ex3_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"CPI.csv"), header=T)

# specify model: AR1 CPI, limiting effect on Bev-Holt
env <- list(
  label = "CPI",
  mean = as.matrix(env.dat$CPI), # CPI observations
  logsigma = as.matrix(log(env.dat$CPI_sigma)), # CPI standard error is given/fixed as data
  year = env.dat$Year,
  use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
  lag = 1, # CPI in year t affects recruitment in year t+1
  process_model = "ar1", # fit CPI as AR1 process
  where = "recruit", # CPI affects recruitment
  how = 2) # limiting (carrying capacity)

input <- prepare_wham_input(asap3, recruit_model = 3,
                            model_name = "Ex 3: Projections",
                            ecov = env,
                            NAA_re = list(sigma="rec+1", cor="iid"),
                            age_comp = "logistic-normal-pool0") # logistic normal pool 0 obs

# selectivity = logistic, not age-specific
#   2 pars per block instead of n.ages
#   sel pars of indices 4/5 fixed at 1.5, 0.1 (neg phase in .dat file)
input$par$logit_selpars[1:4,7:8] <- 0 # original code started selpars at 0 (last 2 rows are fixed)

# ---------------------------------------------------------
## Fit model without projections
mod <- suppressWarnings(fit_wham(input, do.proj=F, do.osa=F, do.retro=F, MakeADFun.silent = TRUE))
# saveRDS(mod, file="m6.rds")
# mod <- readRDS("/home/bstock/Documents/wham/sandbox/ex3_projections/m6.rds")

# Add projections to previously fit model
mod_proj <- list()

# default settings: 3 years, use last F, continue ecov
mod_proj[[1]] <- project_wham(mod, proj.opts=list(n.yrs=3, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# 5 years, use last F, average ecov 1992-1996
mod_proj[[2]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=1992:1996, proj.ecov=NULL), MakeADFun.silent = TRUE)

# 5 years, use last F, use last ecov
mod_proj[[3]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=TRUE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# 5 years, use last F, specify high CPI ~ 0.5
mod_proj[[4]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=matrix(c(0.5,0.7,0.4,0.5),ncol=1)), MakeADFun.silent = TRUE)

# 5 years, use last F, specify low CPI ~ -1.5
mod_proj[[5]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=TRUE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=FALSE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=matrix(c(-1.6,-1.3,-1,-1.2),ncol=1)), MakeADFun.silent = TRUE)

# specify catch, 5 years
mod_proj[[6]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=c(10, 2000, 1000, 3000, 20), avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# specify F, 5 years
mod_proj[[7]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=FALSE, proj.F=c(0.001, 1, 0.5, .1, .2), proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# use FXSPR (avg.yrs defaults to last 5 years, 2007-2011), 5 years
mod_proj[[8]] <- project_wham(mod, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
              use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# use avg F (avg.yrs defaults to last 5 years, 2007-2011), 3 years
mod_proj[[9]] <- project_wham(mod, proj.opts=list(n.yrs=3, use.last.F=FALSE, use.avg.F=TRUE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# use avg F 1992-1996, 10 years
mod_proj[[10]] <- project_wham(mod, proj.opts=list(n.yrs=10, use.last.F=FALSE, use.avg.F=TRUE,
              use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=1992:1996,
              cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), MakeADFun.silent = TRUE)

# saveRDS(mod_proj, file="m6_proj.rds")
# mod_proj <- readRDS("/home/bstock/Documents/wham/sandbox/ex3_projections/m6_proj.rds")

#  check marginal nll is the same
nll_proj <-  sapply(mod_proj, function(x) x$opt$obj)
# mod$opt$obj


# plot results
tmp.dir <- tempdir(check=TRUE)
for(m in 1:length(mod_proj)){
  suppressWarnings(plot_wham_output(mod_proj[[m]], dir.main=tmp.dir))
  expect_equal(nll_proj[m], as.numeric(mod$opt$obj), tolerance=1e-6)
  expect_equal(as.numeric(mod_proj[[m]]$opt$par), as.numeric(mod$opt$par), tolerance=1e-3) # parameter values
}

})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))
