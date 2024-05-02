# WHAM example 4: Time-varying selectivity
# as in example 1
#   no environmental covariate
#   2 indices
#   fit to 1973-2016 data
#   age compositions = 7, logistic normal don't pool zero obs (ex 2 used 5, logistic normal pool` zero obs)
#   selectivity = age-specific
# as in example 2
#   selectivity = logistic

# To create test results see file.path(system.file("contribute", package="wham"), "copy_ex4.R"), note whether bias-correction used

# CURRENTLY: selectivity$initial_pars is different between ex4_selectivity and test_ex04_selectivity

# pkgbuild::compile_dll(debug = FALSE)
# pkgload::load_all()
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "ex04_selectivity"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~17 sec

context("Ex 4: Selectivity")

test_that("Ex 4 works",{
path_to_examples <- system.file("extdata", package="wham")
ex4_test_results <- readRDS(file.path(path_to_examples,"ex4_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
inv.logit <- function(x) exp(x)/(1+exp(x))

sel_model <- c(rep("logistic",3), rep("age-specific",5))
sel_re <- list(c("none","none","none"), # m1-m4 logistic
				c("iid","none","none"),
#				c("ar1","none","none"), #can't get a converged model for logistic with two re and two fe for
				c("2dar1","none","none"),
				c("none","none","none"), # m4-m8 age-specific
				c("iid","none","none"),
				c("ar1","none","none"), #now will map all age-specific (mean) parameters to one estimated value
				c("ar1_y","none","none"),
				c("2dar1","none","none"))
initial_pars <- c(rep(list(list(c(2,0.3),c(2,0.3),c(2,0.3))),3), rep(list(list(c(0.1,0.5,0.5,1,1,1),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,0.5,1,1,1,1))),5))
fix_pars <- c(list(NULL,NULL,NULL), rep(list(list(4:6,4,3:6)),5))
basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions

tmp.dir <- tempdir(check=TRUE)
n.mods <- length(sel_re)
mods <- mcheck <- selAA <- vector("list",n.mods)
for(m in 1:n.mods){
	# overwrite initial parameter values in ASAP data file (ex1_SNEMAYT.dat)
	input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2,
				selectivity=list(model=rep(sel_model[m],3), re=sel_re[[m]], initial_pars=initial_pars[[m]], fix_pars = fix_pars[[m]]),
				NAA_re = list(sigma='rec+1',cor='iid'),
				age_comp = "logistic-normal-miss0", # logistic normal, treat 0 obs as missing
				basic_info = basic_info)
 # age-specific selectivity
	# 	# # you can try not fixing any ages first
	# 	# input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2, 
	# 	# 			selectivity=list(model=rep("age-specific",3), re=sel_re[[m]], 
	# 	# 				initial_pars=list(rep(0.5,6), rep(0.5,6), rep(0.5,6))),
	# 	# 			NAA_re = list(sigma='rec+1',cor='iid'),
	#     #           age_comp = "logistic-normal-miss0", # logistic normal, treat 0 obs as missing
	#		#				basic_info = basic_info)
	# fit model
	#mods[[m]] <- suppressWarnings(fit_wham(input, do.osa=F, do.proj=F, do.retro=F, do.sdrep=F, MakeADFun.silent = TRUE)) 
	mods[[m]] <- suppressWarnings(fit_wham(input, do.fit=FALSE, MakeADFun.silent = TRUE))
	if(exists("err")) rm("err") # need to clean this up
  # The !! allows the individual elements in the report to be be seen if there is an error. See ?testthat::quasi_label and example
  expect_equal(length(mods[[!!m]]$par), length(ex4_test_results$par[[!!m]]), tolerance=1e-6) # nll
	# mcheck[[m]] <- check_convergence(mods[[m]], ret=TRUE)
	# expect_equal(mcheck[[!!m]]$convergence, 0) # opt$convergence should be 0
	# expect_false(mcheck[[!!m]]$na_sdrep) # sdrep should succeed
  expect_equal(as.numeric(mods[[!!m]]$fn(ex4_test_results$par[[!!m]])), ex4_test_results$nll[!!m], tolerance=1e-6) # nll
	# expect_equal(mods[[!!m]]$opt$par, ex4_test_results$par[[!!m]], tolerance=1e-1) # parameter values
	# expect_equal(as.numeric(mods[[!!m]]$opt$obj), ex4_test_results$nll[!!m], tolerance=1e-6) # nll	
}
# nll_nofit <- sapply(1:length(mods), function(x) mods[[x]]$fn(ex4_test_results$par[[x]]))
# nofit <- lapply(mods, function(x) {
# 	fit_wham(x$input, do.fit = F)
# })
# mods_fit <- mods
#nll <- sapply(mods, function(x) x$opt$obj)
#nll_nofit <- sapply(1:length(nofit), function(x) mods[[x]]$fn(ex4_test_results$par[[x]]))
# for(m in 1:length(mods)){
#   expect_equal(nll_nofit[!!m], ex4_test_results$nll[!!m], tolerance=1e-6) # nll
# }
# ex4_test_results <- list()
# ex4_test_results$nll <- sapply(mods, function(x) x$opt$obj) #have to remove old model 3 and old model 7 is done differently now.
# ex4_test_results$par <- lapply(mods, function(x) x$opt$par) #have to save with different order and names
# saveRDS(ex4_test_results, file.path(path_to_examples,"ex4_test_results.RDS"))

# for(m in 1:n.mods){
# 	mcheck <- check_convergence(mods[[m]], ret=TRUE)
# 	expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
# 	expect_false(mcheck$na_sdrep) # sdrep should succeed
# 	expect_equal(mods[[!!m]]$opt$par, ex4_test_results$par[[!!m]], tolerance=1e-1) # parameter values
# 	expect_equal(as.numeric(mods[[!!m]]$opt$obj), ex4_test_results$nll[!!m], tolerance=1e-6) # nll	
# }

})
