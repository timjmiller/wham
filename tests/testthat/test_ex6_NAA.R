# WHAM example 6: Numbers-at-age options

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref='naa')
# library(wham)
# devtools::load_all()
# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex6_NAA.R"); etime <- Sys.time(); runtime = etime - btime;
# 11 min

context("Ex 6: Numbers-at-age")

test_that("Ex 6 works",{
path_to_examples <- system.file("extdata", package="wham")
ex6_test_results <- readRDS(file.path(path_to_examples,"ex6_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)

df.mods <- data.frame(NAA_cor = c('---','iid','ar1_y','iid','ar1_a','ar1_y','2dar1','iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('---',rep("rec",2),rep("rec+1",4),rep("rec",2),rep("rec+1",4)),
                      GSI_how = c(rep(0,7),rep(2,6)), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
# df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

mods <- vector("list",n.mods)
mods_proj <- vector("list",n.mods)
fit.mods <- c(1:2,4:8,10:13) # m3 and m9 don't converge
for(m in fit.mods){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"])
  if(NAA_list$sigma == '---') NAA_list = NULL

  ecov <- list(
    label = "GSI",
    mean = as.matrix(env.dat$GSI),
    logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
    year = env.dat$year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    lag = 1, # GSI in year t affects Rec in year t + 1
    process_model = 'ar1', # "rw" or "ar1"
    where = c("none","recruit")[as.logical(df.mods$GSI_how[m])+1],
    how = df.mods$GSI_how[m], # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
    link_model = "linear")

  input <- suppressWarnings(prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = "Ex 6: Numbers-at-age",
                              selectivity=list(model=rep("age-specific",3), re=c("none","none","none"), 
                                initial_pars=list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,0.5,1,1,1,1)), 
                                fix_pars=list(4:6,4,3:6)),
                              NAA_re = NAA_list,
                              ecov=ecov,
                              age_comp = "logistic-normal-miss0")) # logistic normal, treat 0 obs as missing

  # Fit model
  mods[[m]] <- suppressWarnings(fit_wham(input, do.retro=F, do.osa=F, MakeADFun.silent = TRUE))
  mods_proj[[m]] <- suppressWarnings(project_wham(mods[[m]], MakeADFun.silent = TRUE))

  # expect_equal(as.numeric(mod$opt$par), ex6_test_results$pars[[m]], tolerance=1e-3) # parameter values
}

expect_equal(as.numeric(mods[[1]]$opt$obj), ex6_test_results$nll[1], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[2]]$opt$obj), ex6_test_results$nll[2], tolerance=1e-3) # nll
# expect_equal(as.numeric(mods[[3]]$opt$obj), ex6_test_results$nll[3], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[4]]$opt$obj), ex6_test_results$nll[4], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[5]]$opt$obj), ex6_test_results$nll[5], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[6]]$opt$obj), ex6_test_results$nll[6], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[7]]$opt$obj), ex6_test_results$nll[7], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[8]]$opt$obj), ex6_test_results$nll[8], tolerance=1e-3) # nll
# expect_equal(as.numeric(mods[[9]]$opt$obj), ex6_test_results$nll[9], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[10]]$opt$obj), ex6_test_results$nll[10], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[11]]$opt$obj), ex6_test_results$nll[11], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[12]]$opt$obj), ex6_test_results$nll[12], tolerance=1e-3) # nll
expect_equal(as.numeric(mods[[13]]$opt$obj), ex6_test_results$nll[13], tolerance=1e-3) # nll

})
