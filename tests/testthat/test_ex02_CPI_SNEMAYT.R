# Test WHAM example 2: Cold Pool Index effect on SNEMA Yellowtail Flounder recruitment
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# To create test results see file.path(system.file("contribute", package="wham"), "copy_ex2.R"), note whether bias-correction used
# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all()
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "ex02_CPI_SNEMAYT"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~15 sec

context("Ex 2: CPI on yellowtail recruitment")

test_that("Ex 2 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
ex2_test_results <- readRDS(file.path(path_to_examples,"ex2_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"CPI.csv"), header=T)

Ecov_how <- paste0(
  c("none", "controlling-", "none", "limiting-", "limiting-", "controlling-", "controlling-"), 
  c("", "lag-1-", "", rep("lag-1-",4)),
  c("", "linear", "", rep("linear", 4)))

df.mods <- data.frame(Recruitment = c(2,2,3,3,3,3,4),
                      Ecov_process = c(rep("rw",4),rep("ar1",3)),
                      Ecov_how = Ecov_how, stringsAsFactors=FALSE)

n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- dplyr::select(df.mods, Model, tidyselect::everything()) # moves Model to first col
basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

tmp.dir <- tempdir(check=TRUE)
mods <- mcheck <- list()
for(m in 1:n.mods){
  ecov <- list(
    label = "CPI",
    mean = as.matrix(env.dat$CPI),
    logsigma = as.matrix(log(env.dat$CPI_sigma)),
    year = env.dat$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    recruitment_how = matrix(df.mods$Ecov_how[m],1,1)) #

  # (not used in this vignette) can set Ecov = NULL to fit model without Ecov data
  if(is.na(df.mods$Ecov_process[m])) Ecov = NULL

  # Generate wham input from ASAP3 and Ecov data
  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                              model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                              ecov = ecov,
                              NAA_re = list(sigma="rec+1", cor="iid"),
                              age_comp = "logistic-normal-pool0", # logistic normal pool 0 obs
                              basic_info = basic_info)

  # Selectivity = logistic, not age-specific as in ex1
  #   2 pars per block instead of n.ages
  #   sel pars of indices 4/5 fixed at 1.5, 0.1 (specified via neg phase in ex2_SNEMAYT.dat)
  input$par$logit_selpars[1:4,7:8] <- 0 # last 2 rows will not be estimated (mapped to NA)

  # Fit model
  #mods[[m]] <- suppressWarnings(fit_wham(input, do.sdrep = F, do.retro=F, do.osa=F, MakeADFun.silent = TRUE, retro.silent = TRUE))
  mods[[m]] <- suppressWarnings(fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
  
  # mcheck[[m]] <- check_convergence(mods[[m]], ret=TRUE)
  # expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
  # expect_false(mcheck$na_sdrep) # sdrep should succeed
  # expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
  # The !! allows the individual elements in the report to be be seen if there is an error. See ?testthat::quasi_label and example
  expect_equal(length(mods[[!!m]]$par), length(ex2_test_results$par[[!!m]]), tolerance=1e-6) # nll
  expect_equal(as.numeric(mods[[!!m]]$fn(ex2_test_results$par[[!!m]])), ex2_test_results$nll[!!m], tolerance=1e-6) # nll
  # expect_equal(mods[[!!m]]$opt$par, ex2_test_results$par[[!!m]], tolerance=1e-3) # parameter values
  # expect_equal(as.numeric(mods[[!!m]]$opt$obj), ex2_test_results$nll[!!m], tolerance=1e-6) # nll

  # suppressWarnings(plot_wham_output(mod=mods[[m]], dir.main=tmp.dir, plot.opts = list(browse=FALSE)))
}
# ex2_test_results <- list()
# ex2_test_results$nll <- sapply(mods, function(x) x$opt$obj)
# ex2_test_results$par <- lapply(mods, function(x) x$opt$par)
# saveRDS(ex2_test_results, file.path(path_to_examples,"ex2_test_results.RDS"))

# hard to see which model fails bc they're indexed by m
# print out each one by one
# mcheck <- check_convergence(mods[[1]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[1]]$opt$par, ex2_test_results$par[[1]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[1]]$opt$obj), ex2_test_results$nll[1], tolerance=1e-6) # nll

# mcheck <- check_convergence(mods[[2]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[2]]$opt$par, ex2_test_results$par[[2]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[2]]$opt$obj), ex2_test_results$nll[2], tolerance=1e-6) # nll

# mcheck <- check_convergence(mods[[3]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[3]]$opt$par, ex2_test_results$par[[3]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[3]]$opt$obj), ex2_test_results$nll[3], tolerance=1e-6) # nll

# mcheck <- check_convergence(mods[[4]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[4]]$opt$par, ex2_test_results$par[[4]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[4]]$opt$obj), ex2_test_results$nll[4], tolerance=1e-6) # nll

# mcheck <- check_convergence(mods[[5]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[5]]$opt$par, ex2_test_results$par[[5]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[5]]$opt$obj), ex2_test_results$nll[5], tolerance=1e-6) # nll

# mcheck <- check_convergence(mods[[6]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[6]]$opt$par, ex2_test_results$par[[6]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[6]]$opt$obj), ex2_test_results$nll[6], tolerance=1e-6) # nll

# mcheck <- check_convergence(mods[[7]], ret=TRUE)
# expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# expect_false(mcheck$na_sdrep) # sdrep should succeed
# expect_lt(mcheck$maxgr, 1e-5) # maximum gradient should be < 1e-06
# expect_equal(mods[[7]]$opt$par, ex2_test_results$par[[7]], tolerance=1e-3) # parameter values
# expect_equal(as.numeric(mods[[7]]$opt$obj), ex2_test_results$nll[7], tolerance=1e-6) # nll
})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

