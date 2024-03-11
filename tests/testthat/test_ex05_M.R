# WHAM example 5: Ecov and age-year effects on natural mortality

# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all()
# btime <- Sys.time(); devtools::test(filter = "ex05_M"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~4.8 min

context("Ex 5: Natural mortality")

test_that("Ex 5 works",{
path_to_examples <- system.file("extdata", package="wham")
ex5_test_results <- readRDS(file.path(path_to_examples,"ex5_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

Ecov_how <- paste0(
  c("none",rep("",2), rep("none", 5), rep("", 2), "none"),
  c("", rep("lag-0-",2), rep("",5), rep("lag-0-",2), ""),
  c("", "linear", "poly-2", rep("",5), "linear", "poly-2", ""))

mean_model <- c(rep("fixed-M",3), "estimate-M", "weight-at-age", rep("estimate-M",5), "fixed-M")
age_specific <- c(rep(NA,3),TRUE, NA, rep(FALSE, 5), NA)

df.mods <- data.frame(M_model = c(rep("---",3),"age-specific","weight-at-age",rep("constant",5),"---"),
                      mean_model = mean_model,
                      age_specific = age_specific,
                      M_re = c(rep("none",6),"ar1_y","ar1_ay","none","none","ar1_ay"),
                      Ecov_process = rep("ar1",11),
                      Ecov_how = Ecov_how,
                      Ecov_link = c(0,1,2,rep(0,5),1,2,0), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)

basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions


mods <- vector("list",n.mods)
mods_proj <- vector("list",n.mods)
# tofit <- c(1:10,14:17)
tofit <- c(1:3,5:11)
for(m in tofit){
  # set up environmental covariate data and model options
  # see ?prepare_wham_input
  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    #lag = 0, # GSI in year t affects M in same year
    process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
    #where = "M", # GSI affects natural mortality
    M_how = array(df.mods$Ecov_how[m],c(1,1,asap3[[1]]$dat$n_ages,1))) # n_Ecov x n_stocks x n_ages x n_regions
    #link_model = c(NA,"linear","poly-2")[df.mods$Ecov_link[m]+1])

  m_model <- df.mods$M_model[m]
  if(df.mods$M_model[m] == '---') m_model = "age-specific"
  if(df.mods$M_model[m] %in% c("constant","weight-at-age")) est_ages = 1
  if(df.mods$M_model[m] == "age-specific") est_ages = 1:asap3$dat$n_ages
  if(df.mods$M_model[m] == '---') est_ages = NULL
  mean_map <- array(NA, c(1,1,asap3[[1]]$dat$n_ages)) #n_stocks x n_regions x n_ages
  if(df.mods$mean_model[m] == "estimate-M"){
    if(df.mods$age_specific[m]) mean_map[1,1,] <- 1:asap3[[1]]$dat$n_ages
    else mean_map[1,1,] <- 1
  }
  M <- list(
    mean_model = df.mods$mean_model[m],
    re_model = matrix(df.mods$M_re[m], 1,1),
    mean_map = mean_map
  )
  if(df.mods$mean_model[m] == "estimate-M" & !df.mods$age_specific[m]) M$initial_means = array(0.28, c(1,1,asap3[[1]]$dat$n_ages)) #n_stocks x n_regions x n_ages

  # Generate wham input from ASAP3 and Ecov data
  # input <- prepare_wham_input(asap3, recruit_model = 2,
  #                             model_name = "Ex 5: Yellowtail Flounder with GSI effects on M",
  #                             ecov = ecov,
  #                             M = M)
  print(paste("m:",m))
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
  # Fit model
  if(m == 11) {#convergence sensitive to starting values
    temp <- fit_wham(input, do.fit=F, MakeADFun.silent = TRUE)
    x <- temp$fn(ex5_test_results$pars[[11]])
    input$par <- temp$env$parList()
  }
  #mods[[m]] <- fit_wham(input, do.retro=F, do.osa=F)
  mods[[m]] <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F, MakeADFun.silent = TRUE))
  mods_proj[[m]] <- suppressWarnings(project_wham(mods[[m]], MakeADFun.silent = TRUE))
}

# ex5_test_results <- list()
# ex5_test_results$nll <- sapply(mods, function(x) x$opt$obj) 
# ex5_test_results$par <- lapply(mods, function(x) x$opt$par) #have to save with different order and names
# saveRDS(ex5_test_results, file.path(path_to_examples,"ex5_test_results.RDS"))
sapply(tofit, function(x) range(as.numeric(mods[[x]]$opt$obj) - ex5_test_results$nll[x]))
sapply(tofit, function(x) range(as.numeric(mods[[x]]$opt$par) - ex5_test_results$pars[[x]]))

expect_equal(as.numeric(mods[[1]]$opt$par), ex5_test_results$pars[[1]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[1]]$opt$obj), ex5_test_results$nll[1], tolerance=1e-6) # nll
expect_equal(mods[[1]]$opt$objective, mods_proj[[1]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[2]]$opt$par), ex5_test_results$pars[[2]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[2]]$opt$obj), ex5_test_results$nll[2], tolerance=1e-6) # nll
expect_equal(mods[[2]]$opt$objective, mods_proj[[2]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[3]]$opt$par), ex5_test_results$pars[[3]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[3]]$opt$obj), ex5_test_results$nll[3], tolerance=1e-6) # nll
expect_equal(mods[[3]]$opt$objective, mods_proj[[3]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[5]]$opt$par), ex5_test_results$pars[[5]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[5]]$opt$obj), ex5_test_results$nll[5], tolerance=1e-6) # nll
expect_equal(mods[[5]]$opt$objective, mods_proj[[5]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[6]]$opt$par), ex5_test_results$pars[[6]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[6]]$opt$obj), ex5_test_results$nll[6], tolerance=1e-6) # nll
expect_equal(mods[[6]]$opt$objective, mods_proj[[6]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[7]]$opt$par), ex5_test_results$pars[[7]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[7]]$opt$obj), ex5_test_results$nll[7], tolerance=1e-6) # nll
expect_equal(mods[[7]]$opt$objective, mods_proj[[7]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[8]]$opt$par), ex5_test_results$pars[[8]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[8]]$opt$obj), ex5_test_results$nll[8], tolerance=1e-6) # nll
expect_equal(mods[[8]]$opt$objective, mods_proj[[8]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[9]]$opt$par), ex5_test_results$pars[[9]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[9]]$opt$obj), ex5_test_results$nll[9], tolerance=1e-6) # nll
expect_equal(mods[[9]]$opt$objective, mods_proj[[9]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[10]]$opt$par), ex5_test_results$pars[[10]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[10]]$opt$obj), ex5_test_results$nll[10], tolerance=1e-6) # nll
expect_equal(mods[[10]]$opt$objective, mods_proj[[10]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

expect_equal(as.numeric(mods[[11]]$opt$par), ex5_test_results$pars[[11]], tolerance=1e-1) # parameter values
expect_equal(as.numeric(mods[[11]]$opt$obj), ex5_test_results$nll[11], tolerance=1e-6) # nll
expect_equal(mods[[11]]$opt$objective, mods_proj[[11]]$opt$objective, tolerance=1e-6) # projection shouldn't change nll

})
