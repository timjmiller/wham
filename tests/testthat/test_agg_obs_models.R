# Test WHAM aggregate observation likelihoods: 
#   - ensure likelhoods are evaluating correctly
#   - ensure simulations for each probability model is performing correctly

#dyn.unload(getLoadedDLLs()[["wham"]][["path"]])
pkgbuild::compile_dll(debug = FALSE)
pkgload::load_all()
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "agg_obs_model"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~35 sec

context("Test evaluation of various aggregate observation likelihoods")

test_that("Aggregate observation likelihoods evaluate correctly",{
# get results to check NLL

suppressWarnings(suppressMessages({

  path_to_examples <- system.file("extdata", package="wham")
  agg_obs_tests <- readRDS(file.path(path_to_examples,"agg_obs_likelihood_test_values.rds"))
  tmp.dir <- tempdir(check=TRUE)

  index_3_fn <- function(i, fit=fit0, seed = 123){

    pred <- fit$rep$pred_indices[,3]
    input <- fit$input
    input$data$index_model[3] <- i
    input$par <- fit$parList
    CV_min <- function(lambda) ((1-exp(-lambda))/lambda - exp(-lambda))^0.5 #as beta -> 0
    CVplus <- function(lambda,alpha){
      sqrt((1-exp(-lambda)) * (alpha+1)/(alpha*lambda) - exp(-lambda))
    }  
    if(i %in% c(5,13)) input$par$index_pars[3,1] <- log(3)

    if(i == 6){
      lambda <- 2.2
      alpha <- 0.1
      beta <- (1 - exp(-lambda))*pred/(lambda*alpha)
      input$data$agg_index_cv[,3] <- CVplus(lambda,alpha)
      input$par$index_pars[3,1] <- log(alpha/(1-alpha))
    }
    if(i == 7){
      CVplus <- function(lambda,alpha){
        sqrt((1-exp(-lambda)) * (alpha+1)/(alpha*lambda) - exp(-lambda))
      }
      CVplus_b <- function(lambda,beta,pred){
        sqrt((1-exp(-lambda))/lambda + beta/pred - exp(-lambda))
      }  
      lambda <- 2.2
      #beta <- 1000
      #alpha <- pred * (1 - exp(-lambda))/(lambda*beta)
      #beta <- pred * (1 - exp(-lambda))/(alpha * lambda)
      #CV_in <- CVplus_b(lambda,beta,pred)
      #alpha <- (pred/beta) * (CV_in^2 + exp(-lambda)) - 1
 
      CV_in <- 0.1 + CV_min(lambda)
      beta_max <- pred * CV_in^2
      true_beta_max <- pred * (CV_in^2 + exp(-lambda))
      beta <- pred*(CV_in^2 + exp(-lambda) - (1-exp(-lambda))/lambda)
      alpha <- (1 - exp(-lambda)) * pred/(lambda * beta)
      beta_frac <- beta/beta_max #constant
      input$par$index_pars[3,1] <- log(beta_frac[1]/(1-beta_frac[1]))
      # input$par$index_pars[3,1] <- log(beta)
      input$data$agg_index_cv[,3] <- CV_in
    }
    if(i == 8){
      lambda <- 2.2
      alpha <- 1
      beta <- pred * (1 - exp(-lambda))/(alpha * lambda)
      CV_in <- CVplus(lambda,alpha)
      input$par$index_pars[3,1] <- log(lambda)
      input$par$index_pars[3,2:3] <- -1e10
      input$data$agg_index_cv[,3] <- CV_in
    }
    if(i %in% 9:13){ #delta models
      lambda <- 2.2
      input$par$index_pars[3,2] <- log(lambda)
      input$par$index_pars[3,3:4] <- -1e10
    }

    nofit <- fit_wham(input, do.fit = FALSE)
    set.seed(seed)
    sims <- nofit$simulate(complete = TRUE)
    input$data <- sims
    out <- list(input = input)
    
    if(i %in% 6:8){
      parhat <- nofit$rep$index_annual_pars[3,,]
      muhat <- parhat[,1]
      phihat <- parhat[,2]
      phat <- parhat[,3]
      lambdahat <- muhat^(2-phat)/(phihat * (2 - phat)) # = lambda
      alphahat <- (2-phat)/(phat - 1) # = alpha
      betahat <- phihat * (phat-1)*muhat^(phat -1) # = beta

      out$pars <- list(lambda = lambda, alpha = alpha, beta = beta)
      out$diffs <-  list(lambda = lambda-lambdahat, alpha = alpha-alphahat, beta = beta-betahat)
    }
    return(out)
  }
  
  models <- c(
    "lognormal-1",
    "lognormal-2",
    "gamma-1",
    "gamma-2",
    "gengamma",
    "CPG-alpha",
    "CPG-beta",
    "CPG-lambda",
    "delta-lognormal-1",
    "delta-lognormal-2",
    "delta-gamma-1",
    "delta-gamma-2",
    "delta-gengamma")


  # input <- prepare_wham_input(asap3, recruit_model=2, model_name="SNEMA Yellowtail Flounder",
  #   selectivity=list(model=rep("age-specific",3), re=rep("none",3), 
  #         initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
  #       fix_pars=list(4:5,4,2:4)))
  asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
  
  #simplify initial abundance at age model
  input0 <- prepare_wham_input(asap3)
  N1_pars <- array(exp(10), dim = c(1,1,input0$data$n_ages)) 
  N1_pars[1,1,2] <- exp(-10)
  N1_map <- array(NA, dim = c(1,1,input0$data$n_ages)) 
  N1_map[1,1,1] <- 1
  NAA_info <- list(N1_model = "equilibrium", N1_pars = N1_pars, N1_map = N1_map)

  #stabilize selectivity estimation
  selectivity_info <- list(selectivity=list(model=rep("age-specific",3), re=rep("none",3), 
          initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
        fix_pars=list(4:5,4,2:4)))
  input <- prepare_wham_input(asap3, model_name="SNEMA Yellowtail Flounder", selectivity = selectivity_info,
    NAA_re = NAA_info)

  
  #add a third index where different observation models are used
  index_info <- list(
      n_indices = NCOL(input$data$agg_indices)+1,
      agg_indices = cbind(input$data$agg_indices[,c(1,2,1)]), #(n_years x n_indices)
      units_indices = input$data$units_indices[c(1,2,1)], #(n_indices) 1: biomass 2: numbers
      units_index_paa = input$data$units_index_paa[c(1,2,1)],
      agg_index_cv = cbind(input$data$agg_index_cv[,c(1,2,1)]), #(n_years x n_indices)
      index_Neff = cbind(input$data$index_Neff[,c(1,2,1)]), #(n_years x n_indices)
      fracyr_indices = cbind(input$data$fracyr_indices[,c(1,2,1)]), #(n_years x n_indices)
      use_indices = cbind(input$data$use_indices[,c(1,2,1)]), #(n_years x n_indices)
      use_index_paa = cbind(input$data$use_index_paa[,c(1,2,1)]),  #(n_years x n_indices)
      index_paa = input$data$index_paa[c(1,2,1),,],  #(n_indices x n_years x n_ages)
      selblock_pointer_indices = cbind(input$data$selblock_pointer_indices[,c(1,2,1)]), #(n_years x n_indices)
      waa_pointer_indices = input$data$waa_pointer_indices[c(1,2,1)], #(n_indices)
      index_regions = rep(1,input$data$n_indices+1)) #(n_indices)


  #base input
  input0 <- set_indices(input, index_info)

  inputs <- fits <- vector("list",  length(models))
  names(inputs) <- names(fits) <-models

  fit0 <- fit_wham(input0, do.osa = FALSE, do.retro=FALSE, do.fit = TRUE, MakeADFun.silent = FALSE)
  # inputs[[1]] <- input0
  # fits[[1]] <- fit0
  for(i in 1:5){
    agg_obs_model_info <- list(obs_model = list(fleet = "lognormal-1", index = c("lognormal-1","lognormal-1", models[i])))
    inputs[[i]] <- set_agg_obs_model(input0, agg_obs_model_info)
    fits[[i]] <- fit_wham(inputs[[i]], do.osa = FALSE, do.retro=FALSE, do.fit = TRUE, MakeADFun.silent = FALSE)
  }

  i <- 2
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)  
  #confidence interval for sigma inflation factor
  1 + exp((vals[1] + c(-1,1)*qnorm(0.975)*vals[2]))

  i <- 4
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)  
  #confidence interval for CV inflation factor
  1 + exp((vals[1] + c(-1,1)*qnorm(0.975)*vals[2]))

  i <- 5
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)  
  #confidence interval for alpha
  exp((vals[1] + c(-1,1)*qnorm(0.975)*vals[2]))


  #Simulate data for CPG models and do fits
  results <- list()# <- list(pars = list(), estimates = list())  

  for(i in 1:5){
    agg_obs_model_info <- list(obs_model = list(fleet = "lognormal-1", index = c("lognormal-1","lognormal-1", models[i])))
    temp <- index_3_fn(i, seed = 456)
    inputs[[i]] <- set_agg_obs_model(temp$input, agg_obs_model_info)
    fits[[i]] <- fit_wham(inputs[[i]], do.osa = FALSE, do.retro=FALSE, do.fit = TRUE, MakeADFun.silent = FALSE)    
  }

  i = 5
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)
  #confidence interval for alpha
  exp((vals[1] + c(-1,1)*qnorm(0.975)*vals[2]))
  #estimated and true value
  c(exp(fits[[i]]$parList$index_pars[3,1]), fits[[i]]$input$data$index_annual_pars[3,1,3])
  pars <- fits[[i]]$input$data$index_annual_pars[3,1,2:3]
  sigma <- pars[1]
  nu <- pars[2]
  theta <- 1/((sigma*nu)^2)
  cv <- sqrt(exp(lgamma(theta) + lgamma(theta+2.0/nu) - 2.0*lgamma(theta + 1.0/nu)) - 1)
  cv # 0.2

  for(i in 6:8){
    agg_obs_model_info <- list(obs_model = list(fleet = "lognormal-1", index = c("lognormal-1","lognormal-1", models[i])))
    temp <- index_3_fn(i, seed = 123)
    #nofit <- fit_wham(temp$input, do.fit = FALSE)    
    results[[i]] <- list(pars = temp$pars)
    inputs[[i]] <- set_agg_obs_model(temp$input, agg_obs_model_info)
    fits[[i]] <- fit_wham(inputs[[i]], do.osa = FALSE, do.retro=FALSE, do.fit = TRUE, MakeADFun.silent = FALSE)
    nofit <- fit_wham(inputs[[i]], do.fit = FALSE)    

    parhat <- fits[[i]]$rep$index_annual_pars[3,,]
    muhat <- parhat[,1]
    phihat <- parhat[,2]
    phat <- parhat[,3]

    lambdahat <- muhat^(2-phat)/(phihat * (2 - phat)) # = lambda
    alphahat <- (2-phat)/(phat - 1) # = alpha
    betahat <- phihat * (phat-1)*muhat^(phat -1) # = beta
    estimates <- list(lambda = lambdahat, alpha = alphahat, beta = betahat)
    results[[i]]$estimates <- estimates
  }



  i <- 6 #CPG-alpha
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)  
  #confidence interval for alpha
  1/(1+exp(-(vals[1] + c(-1,1)*qnorm(0.975)*vals[2])))
  #estimated and true value
  c(fits[[i]]$rep$alpha[1,3], fits[[i]]$input$data$alpha[1,3])

  i <- 7 #CPG-beta
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)
  #confidence interval for fraction of beta
  1/(1+exp(-(vals[1] + c(-1,1)*qnorm(0.975)*vals[2])))
  #estimated and true value
  c(fits[[i]]$rep$beta[1,3]/fits[[i]]$rep$beta_max[1,3], fits[[i]]$input$data$beta[1,3]/fits[[i]]$input$data$beta_max[1,3])

  i <- 8 #CPG-lambda
  vals <- c(fits[[i]]$parList$index_pars[3,1], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)
  #confidence interval for alpha
  fits[[i]]$rep$index_lambda_min[1,3] + exp((vals[1] + c(-1,1)*qnorm(0.975)*vals[2]))
  #estimated and true value
  c(fits[[i]]$rep$lambda[1,3], fits[[i]]$input$data$lambda[1,3])


  for(i in 9:13){
    agg_obs_model_info <- list(obs_model = list(fleet = "lognormal-1", index = c("lognormal-1","lognormal-1", models[i])))
    temp <- index_3_fn(i, seed = 456)
    inputs[[i]] <- set_agg_obs_model(temp$input, agg_obs_model_info)
    fits[[i]] <- fit_wham(inputs[[i]], do.osa = FALSE, do.retro=FALSE, do.fit = TRUE, MakeADFun.silent = FALSE)
  }
  i = 13
  vals <- c(fits[[i]]$parList$index_pars[3,2], fits[[i]]$sdrep$cov.fixed[names(fits[[i]]$opt$par)=="index_pars",names(fits[[i]]$opt$par)=="index_pars"]^0.5)
  #confidence interval for lambda
  exp((vals[1] + c(-1,1)*qnorm(0.975)*vals[2]))
  #estimated and true value
  c(exp(fits[[i]]$parList$index_pars[3,2]), -log(1-fits[[i]]$input$data$index_annual_pars[3,1,4]))

    agg_obs_model_info <- list(obs_model = list(fleet = "lognormal-1", index = c("lognormal-1","lognormal-1", models[i])))
    pred <- fit0$rep$pred_indices[,3]
    input <- fit0$input
    input$data$index_model[3] <- i
    input$par <- fit0$parList
    lambda <- 2.2
    input$par$index_pars[3,2] <- log(lambda)
    input$par$index_pars[3,3:4] <- -1e10
#    input$par$index_pars[3,1] <- fits[[5]]$parList$index_pars[3,1]
    nofit <- fit_wham(input, do.fit = FALSE)
    set.seed(678)
    sims <- nofit$simulate(complete = TRUE)
    # input$data$index_pos <- sims$index_pos
    # input$data$agg_indices[which(input$data$index_pos==0)] <- 0 
    input$data <- sims
    input <- set_agg_obs_model(input, agg_obs_model_info)
    # temp <- matrix(input$map$index_pars, 3,4)
    # temp[3,1] <- NA
    # input$map$index_pars <- factor(temp)
    # nofit <- fit_wham(input, do.osa = FALSE, do.retro=FALSE, do.fit = FALSE, MakeADFun.silent = FALSE) 
    fit <- fit_wham(input, do.osa = FALSE, do.retro=FALSE, do.fit = TRUE, MakeADFun.silent = FALSE)


}))
})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

