# WHAM example 4: Time-varying selectivity
# as in example 1
#   no environmental covariate
#   2 indices
#   fit to 1973-2016 data
#   age compositions = 7, logistic normal don't pool zero obs (ex 2 used 5, logistic normal pool` zero obs)
#   selectivity = age-specific
# as in example 2
#   selectivity = logistic

# btime <- Sys.time(); testthat::test_file("/home/bstock/Documents/wham/tests/testthat/test_ex4_selectivity.R"); etime <- Sys.time(); runtime = etime - btime;
# 4.5 min

context("Ex 4: Selectivity")

test_that("Ex 4 works",{
path_to_examples <- system.file("extdata", package="wham")
ex4_test_results <- readRDS(file.path(path_to_examples,"ex4_test_results.rds"))

asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
inv.logit <- function(x) exp(x)/(1+exp(x))

sel_model <- c(rep("logistic",5), rep("age-specific",4))
sel_re <- list(c("none","none","none"), # all logistic
				c("iid","none","none"),
				c("ar1","none","none"),
				c("ar1_y","none","none"),
				c("2dar1","none","none"),
				c("none","none","none"), # all age-specific
				c("iid","none","none"),
				c("ar1_y","none","none"),
				c("2dar1","none","none"))

tmp.dir <- tempdir(check=TRUE)
n.mods <- length(sel_re)
mods <- vector("list",n.mods)
selAA <- vector("list",n.mods)
# for(m in 1:n.mods){
for(m in c(1:3,5:6,8)){ # only models that converge
	if(sel_model[m] == "logistic"){ 
		input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2, 
					selectivity=list(model=rep("logistic",3), re=sel_re[[m]], initial_pars=list(c(1.5,0.2),c(2,0.2),c(2,0.2))),
					NAA_re = list(sigma='rec+1',cor='iid'))
	} else {
		# fix3: 1,4,5 / 4 / 2
		input <- prepare_wham_input(asap3, model_name=paste(paste0("Model ",m), sel_model[m], paste(sel_re[[m]], collapse="-"), sep=": "), recruit_model=2, 
					selectivity=list(model=rep("age-specific",3), re=sel_re[[m]], initial_pars=list(c(inv.logit(-4),0.5,0.5,1,0.5,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,0.5,0.5,0.5,0.5)), fix_pars=list(4,4,2)),
					NAA_re = list(sigma='rec+1',cor='iid'))
	}

	# overwrite age comp model (all models use logistic normal)
	input$data$age_comp_model_indices = rep(7, input$data$n_indices)
	input$data$age_comp_model_fleets = rep(7, input$data$n_fleets)
	input$data$n_age_comp_pars_indices = rep(1, input$data$n_indices)
	input$data$n_age_comp_pars_fleets = rep(1, input$data$n_fleets)
	input$par$index_paa_pars = rep(0, input$data$n_indices)
	input$par$catch_paa_pars = rep(0, input$data$n_fleets)
	input$map = input$map[!(names(input$map) %in% c("index_paa_pars", "catch_paa_pars"))]

	# fit model
	mods[[m]] <- fit_wham(input, do.osa=F, do.proj=F, do.retro=F) 
	if(exists("err")) rm("err") # need to clean this up
}

for(m in c(1:3,5:6,8)){
	plot_wham_output(mod=mods[[m]], out.type='html', dir.main=tmp.dir)

	mcheck <- check_convergence(mods[[m]], ret=TRUE)
	expect_equal(mcheck$convergence, 0) # opt$convergence should be 0
	expect_false(mcheck$na_sdrep) # sdrep should succeed
	expect_equal(as.numeric(mods[[m]]$opt$par), ex4_test_results$pars[[m]], tolerance=1e-3) # parameter values
	expect_equal(as.numeric(mods[[m]]$opt$obj), ex4_test_results$nll[m], tolerance=1e-6) # nll	
}

})
