# Test WHAM example 2: Cold Pool Index effect on SNEMA Yellowtail Flounder recruitment
# Replicate Miller et al 2016 results
#   adds environmental covariate (CPI, treated as rw)
#   uses 5 indices (ex 1: only 2 indices)
#   only fit to 1973-2011 data (ex 1: 1973-2016)
#   age compositions = 5, logistic normal pool zero obs (ex 1: 7, logistic normal missing zero obs)
#   selectivity = logistic (ex 1: age-specific)

# To create test results see file.path(system.file("contribute", package="wham"), "copy_ex2.R"), note whether bias-correction used
# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all()
# library(wham)
# btime <- Sys.time(); devtools::test(filter = "ex02_CPI_SNEMAYT"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~15 sec

context("Ex 2: CPI on yellowtail recruitment")

test_that("Ex 2 works",{

  suppressWarnings(suppressMessages({

    # get results to check NLL and par estimates
    ex2_test_results <- readRDS(testthat::test_path("fixtures","ex02_test_results.RDS"))
    path_to_examples <- system.file("extdata", package="wham")
    #ex2_test_results <- readRDS(file.path(path_to_examples,"ex2_test_results.RDS"))

    asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
    env.dat <- read.csv(file.path(path_to_examples,"CPI.csv"), header=T)

    Ecov_how <- paste0(
      c("none", "controlling-", "none", "limiting-", "limiting-", "controlling-", "controlling-"), 
      c("", "lag-1-", "", rep("lag-1-",4)),
      c("", "linear", "", rep("linear", 4)))

    df.mods <- data.frame(Recruitment = c(2,2,3,3,3,3,4),
                          Ecov_process = c(rep("rw",4),rep("ar1",3)),
                          Ecov_how = Ecov_how, stringsAsFactors=FALSE)

    n.mods <- dim(df.mods)[1]
    df.mods$Model <- paste0("m",1:n.mods)
    df.mods <- df.mods |> dplyr::relocate(Model) # moves Model to first col
    #basic_info <- list(bias_correct_process=TRUE, bias_correct_observation=TRUE) #compare to previous versions
    basic_info <- NULL

    tmp.dir <- tempdir(check=TRUE)
    mods <- mcheck <- inputs <- list()
    for(m in 1:n.mods){
      ecov <- list(
        label = "CPI",
        mean = as.matrix(env.dat$CPI),
        logsigma = as.matrix(log(env.dat$CPI_sigma)),
        year = env.dat$Year,
        use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
        process_model = df.mods$Ecov_process[m], # "rw" or "ar1"
        recruitment_how = matrix(df.mods$Ecov_how[m],1,1)) #

      # (not used in this vignette) can set Ecov = NULL to fit model without Ecov data
      if(is.na(df.mods$Ecov_process[m])) Ecov = NULL

      # Generate wham input from ASAP3 and Ecov data
      inputs[[m]] <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[m],
                                  model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                                  ecov = ecov,
                                  NAA_re = list(sigma="rec+1", cor="iid"),
                                  age_comp = "logistic-normal-pool0", # logistic normal pool 0 obs
                                  basic_info = basic_info)

      # Selectivity = logistic, not age-specific as in ex1
      #   2 pars per block instead of n.ages
      #   sel pars of indices 4/5 fixed at 1.5, 0.1 (specified via neg phase in ex2_SNEMAYT.dat)
      inputs[[m]]$par$logit_selpars[1:4,7:8] <- 0 # last 2 rows will not be estimated (mapped to NA)
    }
    make_res <- FALSE
    for(m in 1:n.mods){
      # Fit model
      #mods[[m]] <- fit_wham(input, do.sdrep = F, do.retro=F, do.osa=F, MakeADFun.silent = TRUE, retro.silent = TRUE)
      mods[[m]] <- fit_wham(inputs[[m]], do.fit = make_res, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
    }
    if(make_res){
      ex2_test_results <- list()
      ex2_test_results$nll <- sapply(mods, function(x) x$opt$obj)
      ex2_test_results$par <- lapply(mods, function(x) x$opt$par)
      saveRDS(ex2_test_results, testthat::test_path("fixtures","ex02_test_results.RDS"))
    }  
      
    # The !! allows the individual elements in the report to be be seen if there is an error. See ?testthat::quasi_label and example
    for(m in 1:n.mods){
      expect_equal(length(mods[[!!m]]$par), length(ex2_test_results$par[[!!m]]), tolerance=1e-6) # nll
      expect_equal(as.numeric(mods[[!!m]]$fn(ex2_test_results$par[[!!m]])), ex2_test_results$nll[!!m], tolerance=1e-6) # nll
    }
  }))
})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

