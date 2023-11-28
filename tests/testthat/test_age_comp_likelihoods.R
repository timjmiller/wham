# Test WHAM composition likelihoods: 
#   - ensure likelhoods are evaluating correctly
#   - ensure simulation of Dirichlet, logistic-normal with pooling/missing values and mvtweedie is performing correctly

# pkgbuild::compile_dll("c:/work/wham/wham", debug = FALSE)
# pkgload::load_all("c:/work/wham/wham")
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "age_comp"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~19.6 sec

context("Test evaluation of various age composition likelihoods")

test_that("Age comp likelihoods evaluate correctly",{
# get results to check NLL
  path_to_examples <- system.file("extdata", package="wham")
  acomp_tests <- readRDS(file.path(path_to_examples,"age_comp_likelihood_test_values.rds"))
  tmp.dir <- tempdir(check=TRUE)

  asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

  models <- c(
    "multinomial",
    "dir-mult",
    "logistic-normal-miss0",
    "logistic-normal-ar1-miss0",
    "logistic-normal-pool0",
    "dirichlet-miss0",
    "dirichlet-pool0",
    "mvtweedie",
    "dir-mult-linear")

  inputs <- unfit <- list()
  for(i in 1:length(models)){
    inputs[[i]] <- suppressWarnings(prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
    selectivity=list(model=rep("age-specific",3), re=rep("none",3), 
          initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
        fix_pars=list(4:5,4,2:4)),
        age_comp = models[i],
      basic_info = list(bias_correct_process = TRUE, bias_correct_observation = TRUE)))
    unfit[[i]] <- suppressWarnings(fit_wham(inputs[[i]], do.osa = F, do.retro=F, do.fit = F, MakeADFun.silent = TRUE))
  }

  for(i in 1:length(models)){
    expect_equal(unfit[[i]]$rep$nll, acomp_tests[[i]]$nll, tolerance=1e-6, scale=1)
  }

  #check that missing and pooling is working for simulations of proportions
  sims <- siminputs <- list()
  for(i in 1:length(models)){ #check simulations of dirichlet and logistic-normal models
    index_Neff <- inputs[[i]]$data$index_Neff
    index_Neff[] <- 10000L
    siminputs[[i]] <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
      selectivity=list(model=rep("age-specific",3), re=rep("none",3), 
          initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
        fix_pars=list(4:5,4,2:4)),
        age_comp = models[i],
      basic_info = list(bias_correct_process = TRUE, bias_correct_observation = TRUE, index_Neff = index_Neff))
    if(i %in% c(3:5,8)) {
      siminputs[[i]]$par$index_paa_pars[,1] = -10
    } else siminputs[[i]]$par$index_paa_pars[,1] = 10
    if(i %in% c(1:2,8:9,11)) siminputs[[i]]$data$index_Neff[] <- 10000L
    set.seed(1234)
    temp = fit_wham(siminputs[[i]], do.osa = F, do.retro=FALSE, do.fit = FALSE, MakeADFun.silent = TRUE)
    temp$fn()
    sims[[i]] <- temp$simulate(complete=TRUE)
  }


  pred_pos = function(pred,input,pool0=NULL){
    obs = input$data$index_paa
    pred_pos = array(0, dim = dim(obs))
    for(i in 1:dim(obs)[1]) for(y in 1:dim(obs)[2]){
      pool0. = pool0
      if(is.null(pool0.)){ #no pooling or missing
        ind = 1:dim(obs)[3]
        pool0. = FALSE 
      } else ind = which(obs[i,y,]>0)
      if(pool0.){
        posind = 1
        for(a in 1:dim(obs)[3]){
          pred_pos[i,y,ind[posind]] = pred_pos[i,y,ind[posind]] + pred[i,y,a]
          if(obs[i,y,a] > 0 & posind < length(ind)) posind = posind + 1
        }
      } else pred_pos[i,y,ind] = pred[i,y,ind]/sum(pred[i,y,ind])
    }
    return(pred_pos)
  }

  pred_compare = list()
  for(i in 1:length(models)){
    pool0 <- NULL
    if(length(grep("pool", models[i]))) pool0 <- TRUE
    if(length(grep("miss", models[i]))) pool0 <- FALSE
    pred_compare[[i]] <- pred_pos(sims[[i]]$pred_index_paa, siminputs[[i]], pool0=pool0)
  }

  for(i in 3:7){
    expect_equal(sims[[i]]$index_paa, pred_compare[[i]], tolerance = 1e-2, scale=1)
  }

})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

