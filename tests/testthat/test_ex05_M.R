# WHAM example 5: Ecov and age-year effects on natural mortality
# To create test results see file.path(system.file("contribute", package="wham"), "copy_ex5.R")
# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all()
# btime <- Sys.time(); devtools::test(filter = "ex05_M"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~20 sec

context("Ex 5: Natural mortality")

test_that("Ex 5 works",{
path_to_examples <- system.file("extdata", package="wham")
ex5_test_results <- readRDS(file.path(path_to_examples,"ex5_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

Ecov_how <- paste0(
  c("none",rep("",2), rep("none", 9), rep("", 2)),
  c("", rep("lag-0-",2), rep("",9), rep("lag-0-",2)),
  c("", "linear", "poly-2", rep("",9), "linear", "poly-2"))

mean_model <- c(rep("fixed-M",6), "estimate-M", "weight-at-age", rep("estimate-M",6))
age_specific <- c(rep(NA,6),TRUE, NA, rep(FALSE, 6))

df.mods <- data.frame(M_model = c(rep("---",6),"age-specific","weight-at-age",rep("constant",6)),
                      mean_model = mean_model,
                      age_specific = age_specific,
                      M_re = c(rep("none",3),"ar1_a","ar1_y","ar1_ay",rep("none",3),"ar1_a", "ar1_y",rep("ar1_ay",3)),
                      Ecov_process = rep("ar1",14),
                      Ecov_how = Ecov_how, stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)

message(df.mods)
#to configure for bias correction of observation and/or process errors as in older versions of wham
#basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions
basic_info <- NULL

mods <- vector("list",n.mods)
mods_proj <- vector("list",n.mods)
tofit <- which(ex5_test_results$is_conv)  # should be 1, 2, 5, 8, 9, 11, 12
#tofit <- integer()
if(length(tofit)) for(m in tofit){
  # set up environmental covariate data and model options
  # see ?prepare_wham_input
  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    M_how = array(df.mods$Ecov_how[m],c(1,1,asap3[[1]]$dat$n_ages,1))) # n_Ecov x n_stocks x n_ages x n_regions

  mean_map <- NULL
  if(df.mods$mean_model[m] == "estimate-M"){
    if(df.mods$age_specific[m]) mean_map <- array(1:asap3[[1]]$dat$n_ages, dim = c(1,1,asap3[[1]]$dat$n_ages))
    else mean_map <- array(1, dim = c(1,1,asap3[[1]]$dat$n_ages))
  }
  M <- list(
    mean_model = df.mods$mean_model[m],
    re_model = matrix(df.mods$M_re[m], 1,1),
    means_map = mean_map
  )
  if(df.mods$mean_model[m] == "estimate-M" & !df.mods$age_specific[m]) M$initial_means = array(0.28, c(1,1,asap3[[1]]$dat$n_ages)) #n_stocks x n_regions x n_ages

  # Generate wham input from ASAP3 and Ecov data
  input <- suppressWarnings(prepare_wham_input(asap3, recruit_model = 2,
    model_name = paste0("m",m,": ", df.mods$mean_model[m]," + GSI link: ",df.mods$Ecov_how[m]," + M RE: ", df.mods$M_re[m]),
    ecov = ecov,
    selectivity=list(model=rep("logistic",6),
      initial_pars=c(rep(list(c(3,3)),4), list(c(1.5,0.1), c(1.5,0.1))),
      fix_pars=c(rep(list(NULL),4), list(1:2, 1:2))),
    NAA_re = list(sigma='rec+1',cor='iid'),
    M=M,
    age_comp = "logistic-normal-pool0",
    basic_info = basic_info))
  mods[[m]] <- suppressWarnings(fit_wham(input, do.fit=F, MakeADFun.silent = TRUE))
  # The !! allows the individual elements in the report to be be seen if there is an error. See ?testthat::quasi_label and example
  expect_equal(length(mods[[!!m]]$par), length(ex5_test_results$pars[[!!m]]), tolerance=1e-6) # nll
  expect_equal(as.numeric(mods[[!!m]]$fn(ex5_test_results$pars[[m]])), ex5_test_results$nll[!!m], tolerance=1e-6) # nll
  #mods[[m]] <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F, MakeADFun.silent = TRUE))
  # for(p in 1:length(mods[[!!m]]$opt$par)) expect_equal(as.numeric(mods[[!!m]]$opt$par[!!p]), ex5_test_results$pars[[!!m]][!!p], tolerance=1e-1) # parameter values
  # expect_equal(as.numeric(mods[[!!m]]$opt$obj), !!ex5_test_results$nll[m], tolerance=1e-6) # nll

  mods_proj[[m]] <- suppressWarnings(project_wham(mods[[m]], do.sdrep = F, MakeADFun.silent = TRUE))
  expect_equal(length(mods_proj[[!!m]]$par), length(ex5_test_results$pars[[!!m]]), tolerance=1e-6) # nll
  expect_equal(as.numeric(mods_proj[[!!m]]$fn(ex5_test_results$pars[[m]])), ex5_test_results$nll[!!m], tolerance=1e-6) # nll
}

})
