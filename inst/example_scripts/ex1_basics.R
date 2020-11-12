# load wham
library(wham)

# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

# copy asap3 data file to working directory
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)

# confirm you are in the working directory and it has the ASAP_SNEMAYT.dat file
list.files()

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat("ex1_SNEMAYT.dat")

# ---------------------------------------------------------------
# model 1
#   recruitment expectation (recruit_model): random about mean (no S-R function)
#   recruitment deviations (NAA_re): independent random effects
#   selectivity: age-specific (fix sel=1 for age 5 in fishery, age 4 in index1, and age 2 in index2)
input1 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
	                            selectivity=list(model=rep("age-specific",3), 
                                	re=rep("none",3), 
                                	initial_pars=list(c(0.5,0.5,0.5,0.5,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,0.5,0.5,0.5,0.5)), 
                                	fix_pars=list(5,4,2)),
	                            NAA_re = list(sigma="rec", cor="iid"))
m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time

# Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradiet should be < 1e-06)
check_convergence(m1)

# ---------------------------------------------------------------
# model 2
#   as m1, but change age comp likelihoods to logistic normal
input2 = input1
input2$data$age_comp_model_indices = rep(7, input2$data$n_indices)
input2$data$age_comp_model_fleets = rep(7, input2$data$n_fleets)
input2$data$n_age_comp_pars_indices = rep(1, input2$data$n_indices)
input2$data$n_age_comp_pars_fleets = rep(1, input2$data$n_fleets)
input2$par$index_paa_pars = rep(0, input2$data$n_indices)
input2$par$catch_paa_pars = rep(0, input2$data$n_fleets)
input2$map = input2$map[!(names(input2$map) %in% c("index_paa_pars", "catch_paa_pars"))]
m2 <- fit_wham(input2, do.osa = F) # turn off OSA residuals to save time

# Check that m2 converged
check_convergence(m2)

# ---------------------------------------------------------------
# model 3
#   full state-space model, numbers at all ages are random effects (NAA_re$sigma = "rec+1")
input3 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
	                            selectivity=list(model=rep("age-specific",3), 
                                	re=rep("none",3), 
                                	initial_pars=list(c(0.5,0.5,0.5,0.5,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,0.5,0.5,0.5,0.5)), 
                                	fix_pars=list(5,4,2)),
	                            NAA_re = list(sigma="rec+1", cor="iid"))
m3 <- fit_wham(input3, do.osa = F) # turn off OSA residuals to save time

# Check that m3 converged
check_convergence(m3)

# ---------------------------------------------------------------
# model 4
#   as m3, but change age comp likelihoods to logistic normal
input4 = input3
input4$data$age_comp_model_indices = rep(7, input4$data$n_indices)
input4$data$age_comp_model_fleets = rep(7, input4$data$n_fleets)
input4$data$n_age_comp_pars_indices = rep(1, input4$data$n_indices)
input4$data$n_age_comp_pars_fleets = rep(1, input4$data$n_fleets)
input4$par$index_paa_pars = rep(0, input4$data$n_indices)
input4$par$catch_paa_pars = rep(0, input4$data$n_fleets)
input4$map = input4$map[!(names(input4$map) %in% c("index_paa_pars", "catch_paa_pars"))]
m4 <- fit_wham(input4, do.osa = T) # do OSA residuals for m4 bc we'll show that output

# Check that m4 converged
check_convergence(m4)

# ------------------------------------------------------------
# Save list of all fit models
mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)
save("mods", file="ex1_models.RData")

# Compare models by AIC and Mohn's rho
res <- compare_wham_models(mods, fname="ex1_table", sort=TRUE)
res$best

# Project best model, m4,
# Use default values: 3-year projection, use average selectivity, M, etc. from last 5 years
m4_proj <- project_wham(model=mods$m4)

# WHAM output plots for best model with projections
# plot_wham_output(mod=m4, out.type='html')
plot_wham_output(mod=m4_proj, out.type='html')