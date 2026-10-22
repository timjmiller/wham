context("AIC function")

test_that("aic() calculates AIC from a WHAM model fit", {
  path_to_examples <- system.file("extdata", package="wham")
  ex2_test_results <- readRDS(file.path(path_to_examples,"ex2_test_results.rds"))

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
  df.mods <- dplyr::select(df.mods, Model, tidyselect::everything()) # moves Model to first col
  ecov <- list(
    label = "CPI",
    mean = as.matrix(env.dat$CPI),
    logsigma = as.matrix(log(env.dat$CPI_sigma)),
    year = env.dat$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
    process_model = df.mods$Ecov_process[1], # "rw" or "ar1"
    recruitment_how = matrix(df.mods$Ecov_how[1],1,1)) #

  input <- prepare_wham_input(asap3, recruit_model = df.mods$Recruitment[1],
                              model_name = "Ex 2: SNEMA Yellowtail Flounder with CPI effects on R",
                              ecov = ecov,
                              NAA_re = list(sigma="rec+1", cor="iid"),
                              age_comp = "logistic-normal-pool0") # logistic normal pool 0 obs
                              

  input$par$logit_selpars[1:4,7:8] <- 0 # last 2 rows will not be estimated (mapped to NA)

  mod <- fit_wham(input, do.osa = FALSE, do.brps = FALSE, do.retro = FALSE, do.sdrep = FALSE, MakeADFun.silent = TRUE)
  
  expect_equal(round(as.numeric(aic(mod)),3), -1472.279)
  expect_equal(round(as.numeric(aic(mod, conditional = TRUE)),3), -1800.114)

})
