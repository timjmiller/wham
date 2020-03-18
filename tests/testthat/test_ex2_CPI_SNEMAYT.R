# Test WHAM example 2: Cold Pool Index effect on SNEMA Yellowtail Flounder recruitment
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# btime <- Sys.time(); devtools::test("/home/bstock/Documents/wham"); etime <- Sys.time(); runtime = etime - btime;
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex2_CPI_SNEMAYT.R"); etime <- Sys.time(); runtime = etime - btime;
# 12 min

context("Ex 2: CPI on yellowtail recruitment")

test_that("Ex 2 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
ex2_test_results <- readRDS(file.path(path_to_examples,"ex2_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"CPI.csv"), header=T)

df.mods <- data.frame(Recruitment = c(2,2,3,3,3,3,4),
                      Ecov_process = c(rep("rw",4),rep("ar1",3)),
                      Ecov_how = c(0,1,0,2,2,1,1), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- dplyr::select(df.mods, Model, tidyselect::everything()) # moves Model to first col

tmp.dir <- tempdir(check=TRUE)
mods <- list()
for(m in 1:n.mods){
  env <- list(
    label = "CPI",
    mean = as.matrix(env.dat$CPI),
    logsigma = as.matrix(log(env.dat$CPI_sigma)),
    year = env.dat$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # CPI in year t affects recruitment in year t+1
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    where = "recruit", # CPI affects recruitment
    how = df.mods$Ecov_how[m]) # 0 = no effect (but still fit Ecov to compare AIC), 1 = controlling (dens-indep mortality), 2 = limiting (carrying capacity), 3 = lethal (threshold), 4 = masking (metabolism/growth), 5 = directive (behavior)
  if(is.na(df.mods$Ecov_process[m])) env = NULL # if no ecov data

  # source("/home/bstock/Documents/wham/R/prepare_wham_input.R")
  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                              model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                              ecov = env,
                              selectivity=list(model=rep("logistic",6),
                                               initial_pars=c(rep(list(c(3,3)),4), list(c(1.5,0.1), c(1.5,0.1))),
                                               fix_pars=c(rep(list(NULL),4), list(1:2, 1:2))))

  # age comp logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # full state-space model, abundance is the state vector
  input$data$use_NAA_re = 1
  input$data$random_recruitment = 0
  input$map = input$map[!(names(input$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
  input$map$log_R = factor(rep(NA, length(input$par$log_R)))
  input$random = c(input$random, "log_NAA")

  # ---------------------------------------------------------
  ## Fit model
  # mods[[1]] <- fit_wham(input)
  mods[[m]] <- fit_wham(input, do.osa=F, do.retro=F)
  plot_wham_output(mod=mods[[m]], out.type='html', dir.main=tmp.dir)
}
# mod.list <- paste0("/home/bstock/Documents/wham/sandbox/ex2/",grep(".rds",list.files("/home/bstock/Documents/wham/sandbox/ex2"),value=TRUE))
# mods <- lapply(mod.list, readRDS)
# ex2_test_results <- list()
# ex2_test_results$pars <- lapply(mods, function(x) as.numeric(x$opt$par))
# ex2_test_results$nll <- sapply(mods, function(x) x$opt$obj)
# saveRDS(ex2_test_results, file="/home/bstock/Documents/wham/inst/extdata/ex2_test_results.rds")

# hard to see which model fails bc they're indexed by m
# print out each one by one
mcheck <- check_convergence(mods[[1]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[1]]$opt$par), ex2_test_results$pars[[1]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[1]]$opt$obj), ex2_test_results$nll[1], tolerance=1e-6) # nll

mcheck <- check_convergence(mods[[2]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[2]]$opt$par), ex2_test_results$pars[[2]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[2]]$opt$obj), ex2_test_results$nll[2], tolerance=1e-6) # nll

mcheck <- check_convergence(mods[[3]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[3]]$opt$par), ex2_test_results$pars[[3]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[3]]$opt$obj), ex2_test_results$nll[3], tolerance=1e-6) # nll

mcheck <- check_convergence(mods[[4]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[4]]$opt$par), ex2_test_results$pars[[4]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[4]]$opt$obj), ex2_test_results$nll[4], tolerance=1e-6) # nll

mcheck <- check_convergence(mods[[5]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[5]]$opt$par), ex2_test_results$pars[[5]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[5]]$opt$obj), ex2_test_results$nll[5], tolerance=1e-6) # nll

mcheck <- check_convergence(mods[[6]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[6]]$opt$par), ex2_test_results$pars[[6]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[6]]$opt$obj), ex2_test_results$nll[6], tolerance=1e-6) # nll

mcheck <- check_convergence(mods[[7]], ret=TRUE)
expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
expect_false(mcheck$na_sdrep) # sdrep should succeed
expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
expect_equal(as.numeric(mods[[7]]$opt$par), ex2_test_results$pars[[7]], tolerance=1e-3) # parameter values
expect_equal(as.numeric(mods[[7]]$opt$obj), ex2_test_results$nll[7], tolerance=1e-6) # nll
})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

