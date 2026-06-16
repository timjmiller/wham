# Test spasmaddock: 
#   - 2 temporal selectivity blocks for fishery, latter has RE

# pkgbuild::compile_dll(debug = FALSE)
# pkgload::load_all()
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "spasmaddock"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~1.5 min

context("Test Spasmaddock")

test_that("Spasmaddock selectivity and retros",{
# get results to check NLL

suppressWarnings(suppressMessages({

path_to_examples <- system.file("extdata", package="wham")
tests <- readRDS(file.path(path_to_examples,"spasmaddock_test_results.rds"))
tmp.dir <- tempdir(check=TRUE)
asap3 <- readRDS(file.path(path_to_examples,"spasmaddock_asap3.RDS"))


fix_fleet_sel = lapply(1:2, function(x) NA)
init_fleet_sel = lapply(1:2, function(x) c(2, 0.3))

## index age-specific selpars. Start them all at 0.5. 
fix_index_sel = lapply(1:3, function(x) NA)
init_index_sel = lapply(1:3, function(x) rep(0.5,9))

## fix a single age for each index selectivity ====
fix_index_sel[[1]] <- c(8)
init_index_sel[[1]] <- c(rep(0.50, 7), 1.00, 0.5)

fix_index_sel[[2]] <- c(1)
init_index_sel[[2]] <- c(1.0, rep(0.5,8))

fix_index_sel[[3]] <- c(9)
init_index_sel[[3]] <- c(rep(0.5,8),1.0)


sel_list <- list( model=c(rep("logistic", 2), rep("age-specific", 3) ),
  re=c( rep("none", 1), "2dar1", rep("none", 3) ) ,
  initial_pars=c(init_fleet_sel, init_index_sel) ,
  fix_pars=c(fix_fleet_sel, fix_index_sel)  )  


input <- prepare_wham_input(asap3, recruit_model = 2,
  NAA_re = list(cor='2dar1', sigma='rec+1', decouple_recruitment=TRUE),
  age_comp = "logistic-normal-miss0",
  selectivity = sel_list )
input$par <- tests$parList
  fit <- fit_wham(input, do.fit=TRUE, do.sdrep=FALSE, do.retro=FALSE, do.osa=FALSE, save.sdrep=FALSE, MakeADFun.silent = TRUE)
  fit<-do_retro_peels(fit, n.peels = 2, MakeADFun.silent=TRUE, retro.silent = TRUE)

  # test_results <- list()
  # test_results$nll <- c(fit$opt$obj, sapply(fit$peels[1:2], function(x) x$opt$obj)) 
  # test_results$parList <- fit$parList
  # test_results$par <- c(list(fit$opt$par), lapply(fit$peels[1:2], function(x) x$opt$par))
  # saveRDS(test_results, file.path(path_to_examples,"spasmaddock_test_results.RDS"))

  # Check neg-log-likelihoods are within 1e-6
  nll <- c(fit$opt$obj, sapply(fit$peels, function(x) x$opt$obj))
  parest <- c(list(fit$opt$par), lapply(fit$peels, function(x) x$opt$par))

  #nofit0$fn() == nofit1$fn()
  map_cor <- matrix(NA, input$data$n_selblocks, 2)
  map_cor[2,2] <- 1
  
  sel_list$map_re <- list(NULL, c(1,1), NULL, NULL, NULL)
  sel_list$re[2] <- "iid"
  sel_list$map_cor <- map_cor
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit0 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
  
  sel_list$map_re <- NULL
  sel_list$map_cor <- NULL
  sel_list$re[2] <- "ar1_y"
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit1 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
  
  #nofit2$fn() == nofit3$fn()
  sel_list$map_re <- list(NULL, c(1,2), NULL, NULL, NULL)
  sel_list$re[2] <- "iid"
  sel_list$map_cor <- map_cor
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit2 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
  
  sel_list$map_cor <- NULL
  sel_list$re[2] <- "ar1_y"
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit3 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

  #nofit4$fn() == nofit5$fn()
  sel_list$model[2] <- "age-specific"
  sel_list$initial_pars[[2]] <- c(rep(0.5,4), 1, rep(0.5,4))
  sel_list$fix_pars[[2]] <- 5
  sel_list$map_re <- NULL
  sel_list$re[2] <- "ar1_y"
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit4 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

  sel_list$re[2] <- "iid"
  sel_list$map_re <- list(NULL, c(rep(1,4), NA, rep(1,4)), NULL, NULL, NULL)
  sel_list$map_cor <- map_cor
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit5 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

  #nofit6$fn() == nofit7$fn()
  sel_list$re[2] <- "ar1_y"
  sel_list$map_cor <- NULL
  sel_list$map_re <- list(NULL, c(1:4, NA, 5:8), NULL, NULL, NULL)
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit6 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

  sel_list$re[2] <- "iid"
  sel_list$map_re <- NULL
  sel_list$map_cor <- map_cor
  input_rev <- set_selectivity(input, selectivity = sel_list)
  nofit7 <- fit_wham(input_rev, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

}))

  for(m in 1:length(nll)) expect_equal(as.numeric(nll[!!m]), as.numeric(tests$nll[!!m]), tolerance=1e-6, scale=1)
  for(m in 1:length(parest)) expect_equal(length(parest[[!!m]]), length(tests$par[[!!m]]), tolerance=1e-6, scale=1)
  expect_equal(nofit0$fn(),nofit1$fn(), tolerance=1e-6, scale=1)
  expect_equal(nofit2$fn(),nofit3$fn(), tolerance=1e-6, scale=1)
  expect_equal(nofit4$fn(),nofit5$fn(), tolerance=1e-6, scale=1)
  expect_equal(nofit6$fn(),nofit7$fn(), tolerance=1e-6, scale=1)
  expect_equal(length(nofit0$env$par),length(nofit1$env$par), tolerance=1e-6, scale=1)
  expect_equal(length(nofit2$env$par),length(nofit3$env$par), tolerance=1e-6, scale=1)
  expect_equal(length(nofit4$env$par),length(nofit5$env$par), tolerance=1e-6, scale=1)
  expect_equal(length(nofit6$env$par),length(nofit7$env$par), tolerance=1e-6, scale=1)
 
})

