# Test mvtweedie distribution

context("mvtweedie with SNEMAYT")

test_that("mvtweedie works",{

  #library(wham)
  path_to_examples <- system.file("extdata", package="wham")
  #wham.dir <- find.package("wham")
  write.dir <- tempdir(check=TRUE)

  # read asap3 data file and convert to input list for wham
  asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

  # Fit model
  input2 <- suppressWarnings(prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
  	                            selectivity=list(model=rep("age-specific",3),
                                  	re=rep("none",3),
                                  	initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),
                                  	fix_pars=list(4:5,4,2:4)),
  	                            NAA_re = list(sigma="rec", cor="iid"),
                                age_comp = "mvtweedie"))
  m2 <- suppressWarnings(fit_wham(input2, do.osa = F, do.sdrep = F, do.retro=F, n.newton = 1, MakeADFun.silent = TRUE)) # turn off OSA residuals to save time

  # Test against prior objective function
  expect_equal( as.numeric(m2$opt$objective), 2418.76, tol = 1e-1 )

})



