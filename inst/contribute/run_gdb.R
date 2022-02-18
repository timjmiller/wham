pkg.dir <- here::here()

# load wham R functions *without recompiling* cpp
devtools::load_all(path=pkg.dir, compile=FALSE)

dyn.load(file.path(pkg.dir, "src","wham_v0.so"))

# set up the problematic wham model, this is taken from ex1
asap3 <- read_asap3_dat(file.path(pkg.dir, "inst","extdata","ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                              selectivity=list(model=rep("age-specific",3), 
                                  re=rep("none",3), 
                                  initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                  fix_pars=list(4:5,4,2:4)),
                              NAA_re = list(sigma="rec", cor="iid"))

# don't use fit_wham, use TMB::MakeADFun
mod <- TMB::MakeADFun(input$data, input$par, DLL = "wham_v0", random = input$random, map = input$map)

# optimize
mod$opt <- stats::nlminb(mod$par, mod$fn, mod$gr, control = list(iter.max = 1000, eval.max = 1000))
for(i in 1:3) { # newton steps
  g <- as.numeric(mod$gr(mod$opt$par))
  h <- stats::optimHess(mod$opt$par, mod$fn, mod$gr)
  mod$opt$par <- mod$opt$par - solve(h, g)
  mod$opt$objective <- mod$fn(mod$opt$par)
}
mod$sdrep <- TMB::sdreport(mod)

# originally provided by https://github.com/kaskr/TMB_contrib_R/blob/master/TMBhelper/R/check_estimability.R    
test <- wham:::check_estimability(mod)
test$WhichBad

