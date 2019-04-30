# load wham
library(wham)

# create directory for analysis
write.dir <- "/media/brian/ExtraDrive1/brian/Documents/NRC/code/vign1"
dir.create(write.dir)
setwd(write.dir)

# copy asap3 data file to directory
wham.dir <- find.package("wham")
# file.copy(from=file.path(wham.dir,"extdata","SNEMA_ytl.dat"), to=write.dir, overwrite=FALSE)
file.copy(from="/media/brian/ExtraDrive1/brian/Documents/wham/inst/extdata/ASAP_SNEMAYT.dat", to=write.dir, overwrite=FALSE)

# confirm you are in the working directory and it has the ASAP_SNEMAYT.dat file
list.files()

# read asap3 data file and convert to input list for wham
asap3 <- read_asap3_dat("ASAP_SNEMAYT.dat")
# asap3 <- read.asap3.dat.fn("ASAP_SNEMAYT.dat")
input <- prepare_wham_input(asap3, recruit_model=2, model_name="SNEMA Yellowtail Flounder")

# Modify Fbar_ages
input$data$Fbar_ages = 4:5

# Make one or more selectivity blocks with age-specific parameters
age.specific = 1:3 # 3 age-specific blocks
not.age.specific = (1:input$data$n_selblocks)[-age.specific]
input = set_age_sel0(input, age.specific)
input$par$logit_selpars[not.age.specific,c(1:input$data$n_ages,input$data$n_ages + 3:6)] = Inf
input$par$logit_selpars[1,5] = Inf
input$par$logit_selpars[2,4] = Inf
input$par$logit_selpars[3,2] = Inf
input$map$logit_selpars = matrix(input$map$logit_selpars, input$data$n_selblocks, input$data$n_ages + 6)
input$map$logit_selpars[is.infinite(input$par$logit_selpars)] = NA
input$map$logit_selpars[!is.infinite(input$par$logit_selpars)] = 1:sum(!is.infinite(input$par$logit_selpars))
input$map$logit_selpars = factor(input$map$logit_selpars)
base = input

#SCAA, but with random effects for recruitment and index observation error variances fixed
temp = base
temp$random = "log_R"
temp$map = temp$map[!(names(temp$map) %in% c("log_R_sigma", "mean_rec_pars"))]
temp$data$random_recruitment = 1
m1 <- fit_wham(temp)
# m1 <- MakeADFun(temp$data,temp$par,DLL="wham_v1", random = temp$random, map = temp$map)
# m1 = fit.tmb.fn(m1, n.newton = 3)
# m1$rep$selblocks
# m1$sdrep

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
# m2 <- MakeADFun(temp$data,temp$par,DLL="wham_v1", random = temp$random, map = temp$map)
# m2 = fit.tmb.fn(m2, n.newton = 3)
# m2$sdrep

#full state-space model, abundance is the state vector
temp = base
temp$data$use_NAA_re = 1
temp$data$random_recruitment = 0
temp$map = temp$map[!(names(temp$map) %in% c("log_NAA", "log_NAA_sigma", "mean_rec_pars"))]
temp$map$log_R = factor(rep(NA, length(temp$par$log_R)))
temp$random = "log_NAA"
m3 <- fit_wham(temp)
# m3 <- MakeADFun(temp$data,temp$par,DLL="wham_v1", random = temp$random, map = temp$map)
# m3 = fit.tmb.fn(m3, n.newton = 3)
# m3$sdrep

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
# m4 <- MakeADFun(temp$data,temp$par,DLL="wham_v1", random = temp$random, map = temp$map)
# m4 = fit.tmb.fn(m4, n.newton = 3)
# m4$sdrep

# load("/media/brian/ExtraDrive1/brian/Documents/NRC/code/ex1_ICES_mine/ex1_models.RData")

# Compare models by AIC and Mohn's rho
mods <- list(m1=m1, m2=m2, m3=m3, m4=m4)
res <- compare_wham_models(mods, fname="model_comparison", sort=TRUE)
res$tab

# 3-year projection for best model
res$best
m4 <- project_wham(m4, n.years = 3)

# Save list of (all) fit models
save("mods", file="ex1_models.RData")

# is this really projecting?
# temp$data$use_indices[(temp$data$n_years_model-2):temp$data$n_years_model,] = 0

# WHAM output plots for best model
make_wham_plots(mod=m4, out.type='pdf')
