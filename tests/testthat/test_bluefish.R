# Test Bluefish: 
#   - 2 fleets, and lots of indices, some with limited age ranges in age comp

# pkgbuild::compile_dll("c:/work/wham/wham", debug = FALSE)
# pkgload::load_all("c:/work/wham/wham")
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "bluefish"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~3 min

context("Test Bluefish")

test_that("Bluefish fitting performs correctly",{
# get results to check NLL

suppressWarnings(suppressMessages({

  path_to_examples <- system.file("extdata", package="wham")
  bluefish_tests <- readRDS(file.path(path_to_examples,"bluefish_test_results.rds"))
  tmp.dir <- tempdir(check=TRUE)

  asap3 <- read_asap3_dat(file.path(path_to_examples,"bluefish_23.dat"))


  input<-prepare_wham_input(asap3, recruit_model = 2, model_name = "multiwham_olddat", 
          NAA_re = list(sigma = "rec+1", cor = "2dar1", decouple_recruitment = FALSE), basic_info = list(percentSPR = 35,
                 bias_correct_process = TRUE, bias_correct_observation=TRUE, bias_correct_BRPs=FALSE))

  fit<-fit_wham(input, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE, MakeADFun.silent=TRUE) # run multi_WHAM
  fit<-do_retro_peels(fit, MakeADFun.silent=TRUE, retro.silent = TRUE)

  # bluefish_test_results <- list()
  # bluefish_test_results$nll <- c(fit$opt$obj, sapply(fit$peels, function(x) x$opt$obj)) 
  # bluefish_test_results$par <- c(list(fit$opt$par), lapply(fit$peels, function(x) x$opt$par))
  # saveRDS(bluefish_test_results, file.path(path_to_examples,"bluefish_test_results.RDS"))

  # Check neg-log-likelihoods are within 1e-6
  nll <- c(fit$opt$obj, sapply(fit$peels, function(x) x$opt$obj))
  parest <- c(list(fit$opt$par), lapply(fit$peels, function(x) x$opt$par))
}))

  for(m in 1:length(nll)) expect_equal(as.numeric(nll[!!m]), as.numeric(bluefish_tests$nll[!!m]), tolerance=1e-6, scale=1)
  for(m in 1:length(parest)) expect_equal(length(parest[[!!m]]), length(bluefish_tests$par[[!!m]]), tolerance=1e-6, scale=1)

})

