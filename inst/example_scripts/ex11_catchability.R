# WHAM example 11: Illustrate alternative configurations for catchability.
#   - using priors and random effects for catchability.
#   - estimating environmental effects on catchability (and recuitment at the same time)
#   - estimating time-varying random effects on catchability.

#library(wham)
# btime <- Sys.time()
# wham.dir <- find.package("wham")
# source(file.path(wham.dir, "example_scripts", "ex11_catchability.R"))
# etime <- Sys.time()
# runtime = etime - btime

# ~X min

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
# devtools::load_all("~/work/wham/wham", compile= TRUE, recompile=TRUE, reset = FALSE)

library(ggplot2)
library(tidyr)
library(dplyr)

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

#make a list of input components that prepare_wham_input can use to generate an input for fit_wham
#NOTE: specifying 2 indices
make_digifish <- function(years = 1975:2014) {
	digifish = list()
	digifish$ages = 1:10
	digifish$years = years
	na = length(digifish$ages)
	ny = length(digifish$years)

	digifish$n_fleets = 1
	digifish$catch_cv = matrix(0.1, ny, digifish$n_fleets)
	digifish$catch_Neff = matrix(200, ny, digifish$n_fleets)
	digifish$n_indices = 2
	digifish$index_cv = matrix(0.3, ny, digifish$n_indices)
	digifish$index_Neff = matrix(100, ny, digifish$n_indices)
	digifish$fracyr_indices = matrix(0.5, ny, digifish$n_indices)
	digifish$index_units = rep(1, length(digifish$n_indices)) #biomass
	digifish$index_paa_units = rep(2, length(digifish$n_indices)) #abundance
	digifish$maturity = t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny))

	L = 100*(1-exp(-0.3*(1:na - 0)))
	W = exp(-11)*L^3
	nwaa = digifish$n_indices + digifish$n_fleets + 2
	digifish$waa = array(NA, dim = c(nwaa, ny, na))
	for(i in 1:nwaa) digifish$waa[i,,] = t(matrix(W, na, ny))

	digifish$fracyr_SSB = rep(0.25,ny)
	digifish$q = rep(0.3, digifish$n_indices)
	digifish$F = matrix(0.2,ny, digifish$n_fleets)

	digifish$selblock_pointer_fleets = t(matrix(1:digifish$n_fleets, digifish$n_fleets, ny))
	digifish$selblock_pointer_indices = t(matrix(digifish$n_fleets + 1:digifish$n_indices, digifish$n_indices, ny))
	return(digifish)
}
digifish = make_digifish()

selectivity = list(model = c(rep("logistic", digifish$n_fleets),rep("logistic", digifish$n_indices)),
    initial_pars = rep(list(c(5,1)), digifish$n_fleets + digifish$n_indices)) #fleet, index

M = list(initial_means = rep(0.2, length(digifish$ages)))

NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]))
NAA_re$sigma = "rec" #random about mean
NAA_re$use_steepness = 0
NAA_re$recruit_model = 2 #random effects with a constant mean
NAA_re$recruit_pars = exp(10)

##########################################################################################
#1. set up one q to have a prior distribution
#use prior on q for second survey
catchability = list(prior_sd = c(NA, 0.3))

#make input and operating model
input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability)
om = fit_wham(input, do.fit = FALSE)

#simulate data from operating model
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
fit$mohns_rho = mohns_rho(fit) 
plot_wham_output(fit)

#This is provided in plot_wham_output()
png(file.path(write.dir,"prior_posterior_q.png"), width=7, height=7, res=300, units='in')
wham:::plot_q_prior_post(fit)
dev.off()

##########################################################################################
#2. add random effects on q for first index
catchability = list(prior_sd = c(NA, 0.3), re = c("iid", "none")) #time varying catchability on the first index, prior on the second

input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability)

#set value to simulate variation in q
input$par$q_repars[1] = log(0.2)
om = fit_wham(input, do.fit = FALSE)

#simulate data from operating model
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE)#, retro.silent = TRUE)
fit$mohns_rho = mohns_rho(fit) 

plot_wham_output(fit)

#This is similar to what is provided in plot_wham_output(), but with the true values from the simulation also plotted
png(file.path(write.dir,"q_time_series_1.png"), width=7, height=7, res=300, units='in')
pal = viridisLite::viridis(n=2)
plot(fit$years, fit$rep$q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,1), ylab = "q", xlab = "Year")
se = summary(fit$sdrep)
se = matrix(se[rownames(se) == "logit_q_mat",2], length(fit$years))
for( i in 1:input$data$n_indices){
	lines(fit$years, fit$rep$q[,i], lwd = 2, col = pal[i])
  polyy = c(fit$rep$q[,i]*exp(-1.96*se[,i]),rev(fit$rep$q[,i]*exp(1.96*se[,i])))
	polygon(c(fit$years,rev(fit$years)), polyy, col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
	lines(fit$years, newdata$q[,i], lwd = 2, col = pal[i], lty = 2)
}
legend("topright", legend = paste0("Index ", rep(1:input$data$n_indices, each = 2), c(" Est.", " True")), lwd = 2, col = rep(pal, each = 2), lty = c(1,2))
dev.off()


##########################################################################################
#3. add some Ecovs to the model. First fit ecov but no effects of population
ecov = list(
		label = c("Climate variable 1", "Climate variable 2"),
		process_model = c("ar1","ar1"),
		mean = cbind(rnorm(length(digifish$years)),rnorm(length(digifish$years))), 
		logsigma = log(c(0.01, 0.2)), 
		lag = c(0,0),
		years = digifish$years, 
		use_obs = matrix(1,length(digifish$years),2),  
		where = c("none","none"),
		indices = list(2, NULL),
		how = c(1,0))

#get rid of prior on second index. add AR1 random effects on q for first index
catchability = list(re = c("iid", "none"))

input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability, ecov = ecov)
#how[1] conflicts with where[1]. Use where[1].

#set sd and rho of ecov processes
input$par$Ecov_process_pars[2,] = log(c(0.1,0.2)) #sd
#cor pars are c(0.4,-0.3)
input$par$Ecov_process_pars[3,] = log((c(0.4,-0.3)-(-1))/(1-c(0.4,-0.3)))
#set value to simulate variation in q
input$par$q_repars[1] = log(0.2)
om = fit_wham(input, do.fit = FALSE)
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE)#, retro.silent = TRUE)
fit$mohns_rho = mohns_rho(fit) 

plot_wham_output(fit)

png(file.path(write.dir,"q_time_series_2.png"), width=7, height=7, res=300, units='in')
pal = viridisLite::viridis(n=2)
se = summary(fit$sdrep)
se = matrix(se[rownames(se) == "logit_q_mat",2], length(fit$years))
plot(fit$years, fit$rep$q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,1), ylab = "q", xlab = "Year")
for( i in 1:input$data$n_indices){
	lines(fit$years, fit$rep$q[,i], lwd = 2, col = pal[i])
  polyy = c(fit$rep$q[,i]*exp(-1.96*se[,i]),rev(fit$rep$q[,i]*exp(1.96*se[,i])))
	polygon(c(fit$years,rev(fit$years)), polyy, col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
	lines(fit$years, newdata$q[,i], lwd = 2, col = pal[i], lty = 2)
}
legend("topright", legend = paste0("Index ", rep(1:input$data$n_indices, each = 2), c(" Est.", " True")), lwd = 2, col = rep(pal, each = 2), lty = c(1,2))
dev.off()

#compare assumed and estimated ecov process pars
input$par$Ecov_process_pars
fit$parList$Ecov_process_pars


#compare true and estimated time-varying q
input$par$q_repars
fit$parList$q_repars
#estimated variability in q is lower than truth, but estimate has large SE
as.list(fit$sdrep, "Std")$q_repars

##########################################################################################
#4. add some Ecovs to the model and allow effects of first Ecov on q for second index
ecov = list(
		label = c("Climate variable 1", "Climate variable 2"),
		process_model = c("ar1","ar1"),
		mean = cbind(rnorm(length(digifish$years)),rnorm(length(digifish$years))), 
		logsigma = log(c(0.01, 0.2)), 
		lag = c(0,0),
		years = digifish$years, 
		use_obs = matrix(1,length(digifish$years),2),  
		where = c("q","none"),
		indices = list(2, NULL),
		how = c(1,0))

#get rid of prior on second index. add AR1 random effects on q for first index
catchability = list(re = c("iid", "none"))

input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability, ecov = ecov)
#how[1] no longer conflicts with where[1].

#set sd and rho of ecov processes
input$par$Ecov_process_pars[2,] = log(c(0.1,0.2)) #sd
#cor pars are c(0.4,-0.3)
input$par$Ecov_process_pars[3,] = log((c(0.4,-0.3)-(-1))/(1-c(0.4,-0.3)))
#set value to simulate variation in q
input$par$q_repars[1] = log(0.2)

#set value for Ecov_beta effect on q (dims are n_effects (2 + n_indices, max_n_poly, n_Ecov, n_ages)
input$par$Ecov_beta[4,1,1,] = 0.5

x = array(input$map$Ecov_beta, dim = dim(input$par$Ecov_beta))
x[4,1,1,] #all the ages mapped to use the same value. Only the first value is used for q or recruitment.

om = fit_wham(input, do.fit = FALSE)
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE)#, retro.silent = TRUE)
fit$mohns_rho = mohns_rho(fit) 

plot_wham_output(fit)

#compare assumed and estimated ecov process pars
input$par$Ecov_beta[4,1,1,1]
fit$parList$Ecov_beta[4,1,1,1]

#compare assumed and estimated ecov process pars
input$par$Ecov_process_pars
fit$parList$Ecov_process_pars


#compare true and estimated time-varying q
input$par$q_repars
fit$parList$q_repars
#estimated variability in q is lower than truth, but estimate has large SE
as.list(fit$sdrep, "Std")$q_repars

png(file.path(write.dir,"q_time_series_3.png"), width=7, height=7, res=300, units='in')
pal = viridisLite::viridis(n=2)
se = summary(fit$sdrep)
se = matrix(se[rownames(se) == "logit_q_mat",2], length(fit$years))
plot(fit$years, fit$rep$q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,1), ylab = "q", xlab = "Year")
for( i in 1:input$data$n_indices){
	lines(fit$years, fit$rep$q[,i], lwd = 2, col = pal[i])
  polyy = c(fit$rep$q[,i]*exp(-1.96*se[,i]),rev(fit$rep$q[,i]*exp(1.96*se[,i])))
	polygon(c(fit$years,rev(fit$years)), polyy, col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
	lines(fit$years, newdata$q[,i], lwd = 2, col = pal[i], lty = 2)
}
legend("topright", legend = paste0("Index ", rep(1:input$data$n_indices, each = 2), c(" Est.", " True")), lwd = 2, col = rep(pal, each = 2), lty = c(1,2))
dev.off()

##########################################################################################
#5. add some Ecovs to the model and allow effects of first Ecov on q for second index AND recruitment
ecov = list(
		label = c("Climate variable 1", "Climate variable 2"),
		process_model = c("ar1","ar1"),
		mean = cbind(rnorm(length(digifish$years)),rnorm(length(digifish$years))), 
		logsigma = log(c(0.01, 0.2)), 
		lag = c(0,0),
		years = digifish$years, 
		use_obs = matrix(1,length(digifish$years),2),  
		where = list(c("recruit","q"),"none"),
		indices = list(2, NULL),
		how = c(1,0))

#get rid of prior on second index. add AR1 random effects on q for first index
catchability = list(re = c("iid", "none"))

input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability, ecov = ecov)
#how[1] no longer conflicts with where[1].

#set sd and rho of ecov processes
input$par$Ecov_process_pars[2,] = log(c(0.1,0.2)) #sd
#cor pars are c(0.4,-0.3)
input$par$Ecov_process_pars[3,] = log((c(0.4,-0.3)-(-1))/(1-c(0.4,-0.3)))
#set value to simulate variation in q
input$par$q_repars[1] = log(0.2)

x = array(input$map$Ecov_beta, dim = dim(input$par$Ecov_beta))
x[1,1,1,] #all the ages mapped to use the same value. Only the first value is used for q or recruitment.
x[4,1,1,] #all the ages mapped to use the same value. Only the first value is used for q or recruitment.

#set value for Ecov_beta effect on q (dims are n_effects (2 + n_indices, max_n_poly, n_Ecov, n_ages)
input$par$Ecov_beta[4,1,1,] = 0.5
#set value for Ecov_beta effect on recruitment (dims are n_effects (2 + n_indices, max_n_poly, n_Ecov, n_ages)
input$par$Ecov_beta[1,1,1,] = -0.5


om = fit_wham(input, do.fit = FALSE)
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE)#, retro.silent = TRUE)
fit$mohns_rho = mohns_rho(fit) 

plot_wham_output(fit)

#compare assumed and estimated ecov effect on q for second index
input$par$Ecov_beta[4,1,1,1]
fit$parList$Ecov_beta[4,1,1,1]

#compare assumed and estimated ecov effect on recruitment
input$par$Ecov_beta[1,1,1,1]
fit$parList$Ecov_beta[1,1,1,1]

#SE for beta parameters is large, especially for recruitment effect
as.list(fit$sdrep, "Std")$Ecov_beta[c(4,1),1,1,1]

#compare assumed and estimated ecov process pars
input$par$Ecov_process_pars
fit$parList$Ecov_process_pars


#compare true and estimated time-varying q
input$par$q_repars
fit$parList$q_repars

#estimated variability in q is lower than truth, but estimate has large SE
as.list(fit$sdrep, "Std")$q_repars

png(file.path(write.dir,"q_time_series_4.png"), width=7, height=7, res=300, units='in')
pal = viridisLite::viridis(n=2)
se = summary(fit$sdrep)
se = matrix(se[rownames(se) == "logit_q_mat",2], length(fit$years))
plot(fit$years, fit$rep$q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,0.4), ylab = "q", xlab = "Year")
for( i in 1:input$data$n_indices){
	lines(fit$years, fit$rep$q[,i], lwd = 2, col = pal[i])
  polyy = c(fit$rep$q[,i]*exp(-1.96*se[,i]),rev(fit$rep$q[,i]*exp(1.96*se[,i])))
	polygon(c(fit$years,rev(fit$years)), polyy, col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
	lines(fit$years, newdata$q[,i], lwd = 2, col = pal[i], lty = 2)
}
legend("topright", legend = paste0("Index ", rep(1:input$data$n_indices, each = 2), c(" Est.", " True")), lwd = 2, col = rep(pal, each = 2), lty = c(1,2))
dev.off()

