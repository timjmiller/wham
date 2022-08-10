# Test WHAM composition likelihoods: 
#   - make sure likelhoods are evaluating correctly

# library(wham)
# btime <- Sys.time(); devtools::test(filter = "age_comp"); etime <- Sys.time(); runtime = etime - btime;
# ~4.8 min

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
    "dirichlet-pool0")

  inputs <- unfit <- list()
  for(i in 1:length(models)){
    inputs[[i]] <- suppressWarnings(prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
    selectivity=list(model=rep("age-specific",3), re=rep("none",3), 
          initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
        fix_pars=list(4:5,4,2:4)),
        age_comp = models[i]))
    unfit[[i]] <- suppressWarnings(fit_wham(inputs[[i]], do.osa = F, do.retro=F, do.fit = F))
  }

  for(i in 1:length(models)){
    expect_equal(unfit[[i]]$rep$nll, acomp_tests$nll[i], tolerance=1e-6, scale=1)
  }

  #check that missing and pooling is working for simulations of proportions
  sims <- list()
  for(i in 3:7){ #check simulations of dirichlet and logistic-normal models

    inputs[[i]] <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
      selectivity=list(model=rep("age-specific",3), re=rep("none",3), 
          initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
        fix_pars=list(4:5,4,2:4)),
        age_comp = models[i])
    if(i %in% 3:5) {
      inputs[[i]]$par$index_paa_pars[,1] = -20
    } else inputs[[i]]$par$index_paa_pars[,1] = 20
    set.seed(1234)
    temp = fit_wham(inputs[[i]], do.osa = F, do.retro=FALSE, do.fit = FALSE)
    temp$fn()
    sims[[i]] <- temp$simulate(complete=TRUE)
  }


  pred_pos = function(pred,input,pool0=FALSE){
    obs = input$data$index_paa
    pred_pos = array(0, dim = dim(obs))
    for(i in 1:dim(obs)[1]) for(y in 1:dim(obs)[2]){
      ind = which(obs[i,y,]>0)
      if(pool0){
        posind = 1
        for(a in 1:dim(obs)[3]){
          pred_pos[i,y,ind[posind]] = pred_pos[i,y,ind[posind]] + pred[y,i,a]
          if(obs[i,y,a] > 0 & posind < length(ind)) posind = posind + 1
        }
      } else pred_pos[i,y,ind] = pred[ind]/sum(pred[ind])
    }
    return(pred_pos)
  }

  pred_compare = lapply(3:7, function(x) {
    pred_pos(sims[[x]]$pred_index_paa, inputs[[x]])
  })

})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

