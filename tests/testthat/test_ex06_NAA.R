# WHAM example 6: Numbers-at-age options

# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all()
# btime <- Sys.time(); devtools::test(filter = "ex06_NAA"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~30 sec

context("Ex 6: Numbers-at-age")

test_that("Ex 6 works",{
path_to_examples <- system.file("extdata", package="wham")
ex6_test_results <- readRDS(file.path(path_to_examples,"ex6_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

df.mods <- data.frame(NAA_cor = c('---','iid','ar1_y','iid','ar1_a','ar1_y','2dar1','iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('---',rep("rec",2),rep("rec+1",4),rep("rec",2),rep("rec+1",4)),
                      R_how = paste0(c(rep("none",7),rep("limiting-lag-1-linear",6))), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

#no more bias correction
# basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions
basic_info <- NULL

mods <- vector("list",n.mods)
mods_proj <- vector("list",n.mods)
#fit.mods <- c(1:2,4:8,10:13) # m3 and m9 don't converge
fit.mods <- which(ex6_test_results$is_conv)
for(m in fit.mods){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"], decouple_recruitment = FALSE)
  if(NAA_list$sigma == '---') NAA_list = NULL

  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    process_model = 'ar1', # "rw" or "ar1"
    recruitment_how = matrix(df.mods$R_how[m])) # n_Ecov x n_stocks

  input <- suppressWarnings(prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = "Ex 6: Numbers-at-age",
                              selectivity=list(model=rep("age-specific",3), re=c("none","none","none"), 
                                initial_pars=list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,0.5,1,1,1,1)), 
                                fix_pars=list(4:6,4,3:6)),
                              NAA_re = NAA_list,
                              ecov=ecov,
                              basic_info = basic_info,
                              age_comp = "logistic-normal-miss0")) # logistic normal, treat 0 obs as missing

  # Fit model
  print(m)
  # temp <- suppressWarnings(fit_wham(input, do.fit=F))
  # mods[[m]] <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F, MakeADFun.silent = TRUE))
  mods[[m]] <- suppressWarnings(fit_wham(input, do.fit=F, MakeADFun.silent = TRUE))
  expect_equal(length(mods[[!!m]]$par), length(ex6_test_results$par[[!!m]]), tolerance=1e-3)
  expect_equal(as.numeric(mods[[!!m]]$fn(ex6_test_results$par[[!!m]])), as.numeric(ex6_test_results$nll[!!m]), tolerance=1e-3)
  # if fitting the models...
  # expect_equal(as.numeric(mods[[!!m]]$opt$obj), ex6_test_results$nll[!!m], tolerance=1e-3)
  # expect_equal(as.numeric(mod$opt$par), ex6_test_results$pars[[m]], tolerance=1e-3) # parameter values
  # print(c(mods[[m]]$opt$obj, ex6_test_results$nll[m]))
  # mods_proj[[m]] <- suppressWarnings(project_wham(mods[[m]], MakeADFun.silent = TRUE))
}

})
