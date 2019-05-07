# load wham
library(wham)

# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir = ""
dir.create(write.dir)
setwd(write.dir)

# copy asap3 data file to working directory
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ASAP_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)

# confirm you are in the working directory and it has the ASAP_SNEMAYT.dat file
list.files()

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat("ASAP_SNEMAYT.dat")
# asap3 <- read.asap3.dat.fn("ASAP_SNEMAYT.dat")
input <- prepare_wham_input(asap3, recruit_model=2, model_name="SNEMA Yellowtail Flounder")

# Make one or more selectivity blocks with age-specific parameters
age.specific = 1:3 # 3 age-specific blocks
not.age.specific = (1:input$data$n_selblocks)[-age.specific]
input = set_age_sel0(input, age.specific)
input$par$logit_selpars[not.age.specific,c(1:input$data$n_ages,input$data$n_ages + 3:6)] = Inf
input$par$logit_selpars[1,5] = Inf
input$par$logit_selpars[2,4] = Inf
input$par$logit_selpars[3,2] = Inf
# Now redefine the map argument for the selectivity parameters to estimate only selectivity parameters without initial values at lower and upper bounds.
input$map$logit_selpars = matrix(input$map$logit_selpars, input$data$n_selblocks, input$data$n_ages + 6)
input$map$logit_selpars[is.infinite(input$par$logit_selpars)] = NA
input$map$logit_selpars[!is.infinite(input$par$logit_selpars)] = 1:sum(!is.infinite(input$par$logit_selpars))
input$map$logit_selpars = factor(input$map$logit_selpars)
base = input

#SCAA, but with random effects for recruitment
temp = base
temp$random = "log_R"
temp$map = temp$map[!(names(temp$map) %in% c("log_R_sigma", "mean_rec_pars"))]
temp$data$random_recruitment = 1
m1 <- fit_wham(temp)

# Check that m1 converged (m1$opt$convergence should be 0, and the maximum gradiet should be < 1e-06)
check_convergence(m1)

#Like m1, but change age comp likelihoods to logistic normal
temp = base
temp$data$age_comp_model_indices = rep(7, temp$data$n_indices)
temp$data$age_comp_model_fleets = rep(7, temp$data$n_fleets)
temp$data$n_age_comp_pars_indices = rep(1, temp$data$n_indices)
temp$data$n_age_comp_pars_fleets = rep(1, temp$data$n_fleets)
temp$par$index_paa_pars = rep(0, temp$data$n_indices)
temp$par$catch_paa_pars = rep(0, temp$data$n_fleets)
temp$map = temp$map[!(names(temp$map) %in% c("index_paa_pars", "catch_paa_pars"))]
temp$random = "log_R"
temp$map = temp$map[!(names(temp$map) %in% c("log_R_sigma", "mean_rec_pars"))]
temp$data$random_recruitment = 1
m2 <- fit_wham(temp)

# Check that m2 converged
check_convergence(m2)

#full state-space model, abundance is the state vector
temp = base
temp$data$use_NAA_re = 1
temp$data$random_recruitment = 0
temp$map = temp$map[!(names(temp$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
temp$map$log_R = factor(rep(NA, length(temp$par$log_R)))
temp$random = "log_NAA"
m3 <- fit_wham(temp)

# Check that m3 converged
check_convergence(m3)

#Like m3, but change age comp likelihoods to logistic normal
temp = base
temp$data$age_comp_model_indices = rep(7, temp$data$n_indices)
temp$data$age_comp_model_fleets = rep(7, temp$data$n_fleets)
temp$data$n_age_comp_pars_indices = rep(1, temp$data$n_indices)
temp$data$n_age_comp_pars_fleets = rep(1, temp$data$n_fleets)
temp$par$index_paa_pars = rep(0, temp$data$n_indices)
temp$par$catch_paa_pars = rep(0, temp$data$n_fleets)
temp$map = temp$map[!(names(temp$map) %in% c("index_paa_pars", "catch_paa_pars"))]
temp$data$use_NAA_re = 1
temp$data$random_recruitment = 0
temp$map = temp$map[!(names(temp$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
temp$map$log_R = factor(rep(NA, length(temp$par$log_R)))
temp$random = "log_NAA"
m4 <- fit_wham(temp)

# Check that m4 converged
check_convergence(m4)

# Save list of all fit models
mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)
save("mods", file="ex1_models.RData")

# Compare models by AIC and Mohn's rho
res <- compare_wham_models(mods, fname="model_comparison", sort=TRUE)
res$best

# 3-year projection for best model
m4 <- project_wham(m4, n.years = 3)

# WHAM output plots for best model
plot_wham_output(mod=m4, out.type='html')
