# Test WHAM example 12: Illustrate mvtweedie distribution

context("Ex 12: Priors and time-varying and environmental effects on catchability")

test_that("Ex 12 works",{

  #library(wham)
  wham.dir <- find.package("wham")
  write.dir <- tempdir(check=TRUE)

  # copy asap3 data file to working directory
  file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=TRUE)

  # read asap3 data file and convert to input list for wham
  asap3 <- read_asap3_dat( file.path(write.dir,"ex1_SNEMAYT.dat") )

  # Fit model
  input2 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
  	                            selectivity=list(model=rep("age-specific",3),
                                  	re=rep("none",3),
                                  	initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),
                                  	fix_pars=list(4:5,4,2:4)),
  	                            NAA_re = list(sigma="rec", cor="iid"),
                                age_comp = "mvtweedie")
  m2 <- fit_wham(input2, do.osa = F, n.newton = 1) # turn off OSA residuals to save time

  # Test against prior objective function
  expect_equal( as.numeric(m2$opt$objective), 2418.76, tol = 1e-1 )

})
