# load wham
is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)
# read asap3 data file and convert to input list for wham
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

# ---------------------------------------------------------------
# model 1
#   recruitment expectation (recruit_model): random about mean (no S-R function)
#   recruitment deviations (NAA_re): independent random effects
#   selectivity: age-specific (fix sel=1 for ages 4-5 in fishery, age 4 in index1, and ages 2-4 in index2)
input1 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                                selectivity=list(model=rep("age-specific",3), 
                                    re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                    fix_pars=list(4:5,4,2:4)),
                                NAA_re = list(sigma="rec", cor="iid"), basic_info = basic_info)
m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time

# Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradiet should be < 1e-06)
check_convergence(m1)

# ---------------------------------------------------------------
# model 2
#   as m1, but change age comp likelihoods to logistic normal (treat 0 observations as missing)
input2 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                                    selectivity=list(model=rep("age-specific",3), 
                                        re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                        fix_pars=list(4:5,4,2:4)),
                                    NAA_re = list(sigma="rec", cor="iid"),
                                    age_comp = "logistic-normal-miss0", basic_info = basic_info)
m2 <- fit_wham(input2, do.osa = F) # turn off OSA residuals to save time

# Check that m2 converged
check_convergence(m2)

# ---------------------------------------------------------------
# model 3
#   full state-space model, numbers at all ages are random effects (NAA_re$sigma = "rec+1")
input3 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                                selectivity=list(model=rep("age-specific",3), 
                                    re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                    fix_pars=list(4:5,4,2:4)),
                                NAA_re = list(sigma="rec+1", cor="iid"), basic_info = basic_info)
m3 <- fit_wham(input3, do.osa = F) # turn off OSA residuals to save time

# Check that m3 converged
check_convergence(m3)

# ---------------------------------------------------------------
# model 4
#   as m3, but change age comp likelihoods to logistic normal
input4 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                                    selectivity=list(model=rep("age-specific",3), 
                                        re=rep("none",3), 
                                    initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)), 
                                        fix_pars=list(4:5,4,2:4)),
                                    NAA_re = list(sigma="rec+1", cor="iid"),
                                    age_comp = "logistic-normal-miss0", basic_info = basic_info)
m4 <- fit_wham(input4, do.retro = F, do.sdrep = F, do.osa = F) # turn off extras because it turns out convergence is not good.

# Check that m4 converged. Bad max absolute gradient
check_convergence(m4)

#retry starting at best estimates so far
input4$par <- m4$parList #best estimates so far
m4_rev <- fit_wham(input4, do.osa = F) 

# Now convergence is good
check_convergence(m4_rev)

#add OSA residuals
m4_rev <- make_osa_residuals(m4_rev)

# ------------------------------------------------------------
# Save list of all fit models
mods <- list(m1=m1, m2=m2, m3=m3, m4=m4_rev)
save("mods", file="ex1_models.RData")

# Compare models by AIC and Mohn's rho
res <- compare_wham_models(mods, table.opts=list(fname="ex1_table", sort=TRUE))
res$best

# Project best model, m4,
# Use default values: 3-year projection, use average selectivity, M, etc. from last 5 years
m4_proj <- project_wham(model=mods[[4]])

# WHAM output plots for best model with projections
plot_wham_output(mod=m4_proj)
