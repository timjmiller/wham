# WHAM example 10: Operating models and Management Strategy Evaluation
#   - Make a default wham input file without an ASAP3 dat file.
#   - Make an operating model, simulate data, and fit models
#   - Perform a simple management strategy evaluation using with an operating model with Beverton-Holt stock recruitment and a SCAA estimating model
#   - Make some plots

# devtools::install_github("timjmiller/wham", dependencies=TRUE)
is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo
#by default do not perform bias-correction
if(!exists("basic_info")) basic_info <- NULL

library(ggplot2)
library(tidyr)
library(dplyr)

# create directory for analysis, e.g.
# write.dir <- "/path/to/save/ex2" on linux/mac
if(!exists("write.dir")) write.dir = tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

#2: self test

path_to_examples <- system.file("extdata", package="wham")
wham.dir <- find.package("wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))

#Make an input and fit a model in wham. This model was also fit in the first example vignette. 

#define selectivity model
selectivity=list(
  model=rep("age-specific",3), #define selectivity model
  re=rep("none",3), #define selectivity random effects model
  #define initial/fixed values for age-specific selectivity
  initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),
  fix_pars=list(4:5,4,2:4)) #which ages to not estimate age-specfic selectivity parameters

NAA_re = list(
  sigma="rec", #random effects for recruitment only
  cor="iid", #random effects are independent
  recruit_model = 2) #mean recruitment is the only fixed effect parameter

input <- prepare_wham_input(
  asap3, 
  selectivity = selectivity,
  NAA_re = NAA_re,
  model_name="Ex 1: SNEMA Yellowtail Flounder", basic_info = basic_info) 

mod <- fit_wham(input, do.osa = F, do.retro=F, MakeADFun.silent = T) # don't do retro peels and OSA residuals, don't show output during optimization

names(mod)

mod$sdrep #the standard errors of the fixed effects

#the simulate function
mod$simulate

#make a function to simulate and fit the matching estimating model
sim_fn <- function(om, self.fit = FALSE){
    input <- om$input
    input$data = om$simulate(complete=TRUE)
    if(self.fit) {
      fit <- fit_wham(input, do.osa = F, do.retro = F, MakeADFun.silent = T)
      return(fit)
    } else return(input) 
}
set.seed(123)
self_sim_fit <- sim_fn(mod, self.fit = TRUE)

png(file = "self_fit_1_ssb.png")
plot(mod$years, self_sim_fit$input$data$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod$years, self_sim_fit$rep$SSB, lty = 2, col = 'red')
dev.off()

# a way of saving results as you go to make sure you don't loose everything if there is a catastrophic error
set.seed(123)
n_sims <- 10 #more is better, but takes longer
sim_inputs <- replicate(n_sims, sim_fn(mod), simplify=F)
res = list(reps = list(), par.est = list(), par.se = list(), adrep.est = list(), adrep.se =list())
for(i in 1:length(sim_inputs)){
  cat(paste0("i: ",i, "\n"))
  tfit = fit_wham(sim_inputs[[i]], do.osa = F, do.retro = F, MakeADFun.silent = T)
  res$reps[[i]] = tfit$rep
  res$par.est[[i]] = as.list(tfit$sdrep, "Estimate")
  res$par.se[[i]] = as.list(tfit$sdrep, "Std. Error")
  res$adrep.est[[i]] = as.list(tfit$sdrep, "Estimate", report = TRUE)
  res$adrep.se[[i]] = as.list(tfit$sdrep, "Std. Error", report = TRUE)
  saveRDS(res, "self_test_res.RDS")
}

#plot results for bias of SSB
true = sapply(sim_inputs, function(x) x$data$SSB)
est = sapply(res$rep, function(x) return(x$SSB))
SSB_rel_resid = est/true - 1
resid_cis = apply(SSB_rel_resid,1,mean) + apply(SSB_rel_resid,1,sd)*qnorm(0.975)*t(matrix(c(-1,1),2,44))/sqrt(n_sims)

png(file = "ssb_rel_bias.png")
plot(mod$years, apply(SSB_rel_resid,1,mean), ylim = c(-0.1,0.05), ylab = "Estimated Relative Bias SSB", xlab = "Year")
lines(mod$years, resid_cis[,1])
lines(mod$years, resid_cis[,2])
abline(h = 0, lty = 2)
dev.off()


# 3. Fitting a model different from the operating model

input <- prepare_wham_input(
  asap3, 
  selectivity = selectivity,
  NAA_re = NAA_re,
  basic_info = basic_info,
  model_name="Ex 1: SNEMA Yellowtail Flounder") 
names(input)
input_3 <- input #for vignette


#3a: alternative age comp model assumption
input_em <- prepare_wham_input(
  asap3, 
  recruit_model=2, 
  model_name="Ex 1: SNEMA Yellowtail Flounder",
  selectivity=selectivity,
  NAA_re = NAA_re,
  basic_info = basic_info,
  age_comp = "logistic-normal-miss0") #the only difference from original input

#only replace the observations so that the other inputs for configuration are correct.
#IMPORTANT TO PROVIDE obsvec!
sim_fn = function(om, input, do.fit = FALSE){
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  om_data = om$simulate(complete = TRUE)
  input$data[obs_names] = om_data[obs_names]
  input <- set_osa_obs(input) #to transform catch and index paa appropriately
  if(do.fit) {
    fit = fit_wham(input, do.osa = F, do.retro = F, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(input,om_data))
}

#compare observations for multinomial and logistic normal
set.seed(123)
sim = sim_fn(mod, input_em, do.fit = FALSE)
sim[[1]]$data$obs[1:10,] #data.frame of observations and identifying information
sim[[2]]$obsvec[4:9] #showing obsvec that is simulated (multinomial frequencies)
sim[[1]]$data$obsvec[4:9] #showing obsvec that is used (MVN = logistic normal)
paa <- sim[[2]]$index_paa[1,1,]
paa #multinomial proportions : naa/sum(naa)
sim[[2]]$index_Neff[1] * paa # = simulated obsvec (multinomial frequencies)
log(paa[-6]) - log(paa[6]) # = inverse additive logistic transform for obsvec that is used
sim_3a <- sim #for vignette

set.seed(123)
simfit = sim_fn(mod, input_em, do.fit = TRUE)

png(file = "sim_fit_3a.png")
plot(mod$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
dev.off()


#3b: alternative effective sample size
sim_fn = function(om, Neff_om = 50, Neff_em = 200, do.fit = FALSE){
  om_input = om$input #the original input
  om_input$data$index_Neff[] = Neff_om #true effective sample size for the indices
  om_input$data$catch_Neff[] = Neff_om #true effective sample size for the indices
  om_input$par = as.list(om$sdrep, "Estimate") #assume parameter values estimated from original model
  om_alt = fit_wham(om_input, do.fit = F, MakeADFun.silent=T) #make unfitted model for simulating data
  om_data = om_alt$simulate(complete = TRUE)
  em_input = om_input #the same configuration other than Neff
  em_input$data$index_Neff[] = Neff_em #assumed effective sample size for the indices
  em_input$data$catch_Neff[] = Neff_em #assumed effective sample size for the indices
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  em_input$data[obs_names] = om_data[obs_names]
  em_input <- set_osa_obs(em_input) #to rescale the frequencies appropriately
  if(do.fit) {
    fit = fit_wham(input, do.osa = F, do.retro = F, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(em_input,om_data))
}
set.seed(123)
sim = sim_fn(mod, do.fit = FALSE)
sim[[1]]$data$obs[1:10,] #data.frame of observations and identifying information
sim[[2]]$obsvec[4:9] #showing obsvec that is simulated (multinomial frequencies)
sim[[1]]$data$obsvec[4:9] #showing obsvec that is used (MVN = logistic normal)
paa <- sim[[2]]$index_paa[1,1,] #simulated paa
paa #multinomial proportions : naa/sum(naa)
sim[[1]]$data$index_paa[1,1,] #paa in estimating model (the same)
sim[[2]]$index_Neff[1] * paa # = simulated obsvec (multinomial frequencies)
sim[[1]]$data$index_Neff[1] * paa # = em obsvec (multinomial frequencies with larger Neff)

sim_3b <- sim #for vignette

set.seed(123)
simfit = sim_fn(mod, do.fit = TRUE)

png(file = "sim_fit_3b.png")
plot(mod$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
dev.off()

#3c: Alternative population assumptions


NAA_re_om = list(
  sigma="rec+1", #random effects for recruitment and older age classes
  cor="iid", #random effects are independent
  recruit_model = 2) #mean recruitment is the only fixed effect parameter

NAA_re_em = list(
  sigma="rec", #random effects for recruitment only
  cor="iid", #random effects are independent
  recruit_model = 2) #mean recruitment is the only fixed effect parameter

input_om <- prepare_wham_input(
    asap3, 
    selectivity = selectivity,
    NAA_re = NAA_re_om,
    basic_info = basic_info,
    model_name="Ex 1: SNEMA Yellowtail Flounder") 

input_em <- prepare_wham_input(
    asap3, 
    selectivity = selectivity,
    NAA_re = NAA_re_em,
    basic_info = basic_info,
    model_name="Ex 1: SNEMA Yellowtail Flounder") 

om_ss <- fit_wham(input_om, do.osa = F, do.retro=F, MakeADFun.silent = T) # don't do retro peels and OSA residuals, don't show output during optimization

#no differences in treatment of observations
sim_fn = function(om, input, do.fit = FALSE){
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  om_data = om$simulate(complete = TRUE)
  input$data[obs_names] = om_data[obs_names]
  if(do.fit) {
    fit = fit_wham(input, do.osa = F, do.retro = F, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(input,om_data))
}
set.seed(123)
simfit = sim_fn(om_ss, input_em, do.fit = TRUE)

png(file = "sim_fit_3c.png")
plot(om_ss$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(om_ss$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
dev.off()

#3d: Alternative selectivity assumptions

selectivity_em = list(
  model = c(rep("logistic", input$data$n_fleets),rep("logistic", input$data$n_indices)),
  initial_pars = rep(list(c(5,1)), input$data$n_fleets + input$data$n_indices)) #fleet, index

input_em <- prepare_wham_input(
    asap3, 
    selectivity = selectivity_em,
    NAA_re = NAA_re,
    model_name="Ex 1: SNEMA Yellowtail Flounder") 
sim_fn = function(om, input, do.fit = FALSE){
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  om_data = om$simulate(complete = TRUE)
  input$data[obs_names] = om_data[obs_names]
  if(do.fit) {
    fit = fit_wham(input, do.osa = F, do.retro = F, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(input,om_data))
}
set.seed(123)
simfit = sim_fn(mod, input_em, do.fit = TRUE)


png(file = "sim_fit_3d_1.png")
plot(mod$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
dev.off()

png(file = "sim_fit_3d_2.png")
plot(mod$years, simfit[[2]]$SSB/exp(simfit[[2]]$log_SSB_FXSPR[,2]), type = 'l', ylab = "SSB/SSB(F40)", xlab = "Year")
lines(mod$years, simfit[[1]]$rep$SSB/exp(simfit[[1]]$rep$log_SSB_FXSPR[,2]), lty = 2, col = 'red')
dev.off()

# Part 4: Making an operating model without an ASAP file

make_info <- function(base_years = 1982:2021, ages = 1:10, Fhist = "updown", n_feedback_years = 0) { #changed years
    info <- list()
    info$ages <- ages
    info$years <- as.integer(base_years[1] - 1 + 1:(length(base_years) + n_feedback_years))
    na <- length(info$ages)
    ny <- length(info$years)
    info$n_regions <- 1L
    info$n_stocks <- 1L
    nby <- length(base_years)
    mid <- floor(nby/2)
    #up then down

    catch_info <- list()
    catch_info$n_fleets <- 1L
    catch_info$catch_cv <- matrix(0.1, ny, catch_info$n_fleets)
    catch_info$catch_Neff <- matrix(200, ny, catch_info$n_fleets)

    if(Fhist == "updown") catch_info$F <- matrix(0.2 + c(seq(0,0.4,length.out = mid),seq(0.4,0,length.out=nby-mid)),nby, catch_info$n_fleets)

    #down then up
    if(Fhist == "downup") catch_info$F <- matrix(0.2 + c(seq(0.4,0,length.out = mid),seq(0,0.4,length.out=nby-mid)),nby, catch_info$n_fleets)

    if(n_feedback_years>0) catch_info$F <- rbind(catch_info$F, catch_info$F[rep(nby, n_feedback_years),, drop = F]) #same F as terminal year for feedback period

    index_info <- list()
    index_info$n_indices <- 1L
    index_info$index_cv <- matrix(0.3, ny, index_info$n_indices)
    index_info$index_Neff <- matrix(100, ny, index_info$n_indices)
    index_info$fracyr_indices <- matrix(0.5, ny, index_info$n_indices)
    index_info$index_units <- rep(1, length(index_info$n_indices)) #biomass
    index_info$index_paa_units <- rep(2, length(index_info$n_indices)) #abundance
    index_info$q <- rep(0.3, index_info$n_indices)
    
    info$maturity <- array(t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny)), dim = c(1, ny, na))

    L <- 100*(1-exp(-0.3*(1:na - 0)))
    W <- exp(-11)*L^3
    nwaa <- index_info$n_indices + catch_info$n_fleets + 2
    info$waa <- array(NA, dim = c(nwaa, ny, na))
    for(i in 1:nwaa) info$waa[i,,] <- t(matrix(W, na, ny))

    info$fracyr_SSB <- cbind(rep(0.25,ny))

    catch_info$selblock_pointer_fleets <- t(matrix(1:catch_info$n_fleets, catch_info$n_fleets, ny))
    index_info$selblock_pointer_indices <- t(matrix(catch_info$n_fleets + 1:index_info$n_indices, index_info$n_indices, ny))
    return(list(basic_info = info, catch_info = catch_info, index_info = index_info))
}

stock_om_info <- make_info()
stock_basic_info <- c(basic_info, stock_om_info$basic_info)
catch_info <- stock_om_info$catch_info
index_info <- stock_om_info$index_info

selectivity_om <- list(
  model = c(rep("logistic", catch_info$n_fleets),rep("logistic", index_info$n_indices)),
  initial_pars = rep(list(c(5,1)), catch_info$n_fleets + index_info$n_indices)) #fleet, index

M_om <- list(initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))))

NAA_re <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = list(exp(10)),
  sigma_vals <- list(0.5) #sigma_R
)

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om)

stock_om <- fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)

stock_om_4 <- stock_om
save(stock_om_4, file = "stock_om_4.RData")


#4a. Setting NAA random effects variance parameters.


stock_om_info <- make_info()
stock_basic_info <- c(basic_info, stock_om_info$basic_info)
catch_info <- stock_om_info$catch_info
index_info <- stock_om_info$index_info

selectivity_om <- list(
  model = c(rep("logistic", catch_info$n_fleets),rep("logistic", index_info$n_indices)),
  initial_pars = rep(list(c(5,1)), catch_info$n_fleets + index_info$n_indices)) #fleet, index

M_om <- list(initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))))

NAA_re_om <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  cor="2dar1", #random effects correlated among age classes and across time (separably)
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = list(exp(10)),
  sigma_vals = list(c(0.6, rep(0.2,length(stock_basic_info$ages)-1))), #sigma_R
  cor_vals = list(c(0.7, 0.4)) #age, then year
)


stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om)

stock_om_input <- prepare_wham_input(basic_info = stock_om_info, selectivity = selectivity_om, NAA_re = NAA_re_om, M = M_om)

stock_om_input$par$trans_NAA_rho[1:2]
wham:::gen.logit(c(0.7,0.4), -1, 1,s=2) #age first, year second


# 4b. Setting M random effects variance parameters

stock_om_info <- make_info()
stock_basic_info <- c(basic_info, stock_om_info$basic_info)
catch_info <- stock_om_info$catch_info
index_info <- stock_om_info$index_info

selectivity_om <- list(
  model = c(rep("logistic", catch_info$n_fleets),rep("logistic", index_info$n_indices)),
  initial_pars = rep(list(c(5,1)), catch_info$n_fleets + index_info$n_indices)) #fleet, index

M_om <- list(
  mean_model = "estimate-M",
  re_model = matrix("ar1_ay", 1,1),
  initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))),
  sigma_vals = matrix(0.2, 1,1),
  cor_vals = array(c(0.7,0.4),dim = c(1,1,2)))
  
NAA_re_om <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  recruit_model = 2, #random effects with a constant mean
  sigma_vals = list(0.6), #sigma_R
  recruit_pars = list(exp(10))
)

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om)

#stock_om <- fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)


stock_om_input$par$log_NAA_sigma[1]
log(0.6)

stock_om_input$par$M_repars[1]
log(0.2)

stock_om_input$par$M_repars[2:3]
wham:::gen.logit(c(0.7,0.4), -1, 1) #age first, year second


# 4c. Setting selectivity random effects variance parameters

stock_om_info <- make_info()
stock_basic_info <- c(basic_info, stock_om_info$basic_info)
catch_info <- stock_om_info$catch_info
index_info <- stock_om_info$index_info

selectivity_om = list(
  model = c(rep("logistic", stock_om_info$n_fleets),rep("logistic", stock_om_info$n_indices)),
  initial_pars = rep(list(c(5,1)), stock_om_info$n_fleets + stock_om_info$n_indices),
  re = rep("2dar1", stock_om_info$n_fleets + stock_om_info$n_indices)) 

M_om <- list(
  mean_model = "estimate-M",
  initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))),
)
  
NAA_re_om <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  recruit_model = 2, #random effects with a constant mean
  sigma_vals = list(0.6), #sigma_R
  recruit_pars = list(exp(10))
)

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om)


















































#make input and operating model
input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M)
om = fit_wham(input, do.fit = FALSE)

#simulate data from operating model
set.seed(8675309)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
fit$mohns_rho = mohns_rho(fit) 
plot_wham_output(fit)

#make input for a scaa estimation model (recruitment as fixed effects)
scaa_info = digifish
data_names = c("agg_catch", "catch_paa", "agg_indices","index_paa")
scaa_info[data_names] = newdata[data_names]

#the only thing that is different is how the parameterization of numbers at age is configured
scaa_NAA_re = list(N1_pars = NAA_re$N1_pars)
#scaa_NAA_re$use_steepness = 0
scaa_NAA_re$recruit_model = 1 #recruitments as fixed effects

scaa_input = prepare_wham_input(basic_info = scaa_info, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M, recruit_model = 1)

#fit the scaa estimation model
scaa_fit = fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
scaa_fit$mohns_rho = mohns_rho(scaa_fit)

##########################################################################################
# Make the operating model assume a Beverton-Holt stock recruit relationship
NAA_re_bh <- NAA_re
NAA_re_bh$sigma <- "rec" #random about mean
#NAA_re$use_steepness = 1 #ok because M, WAA, etc are constant
NAA_re_bh$recruit_model <- 3 #Beverton-Holt
#NAA_re$recruit_pars = c(0.5, exp(10))
NAA_re_bh$recruit_pars <- list(c(0.5, 2e-4))
NAA_re_bh$sigma_vals <- list(0.5) #sigma_R

#make input object for operating model
bh_input = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re_bh, M = M)

#make the operating model
bh_om = fit_wham(bh_input, do.fit = FALSE)

set.seed(8675309)
sim_pop = bh_om$simulate(complete=TRUE)
temp = bh_input
temp$data = sim_pop

#fit estimating model that is the same as the operating model
bh_fit = fit_wham(temp, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
bh_fit$mohns_rho = mohns_rho(bh_fit)

#make input for a scaa estimation model (recruitment as fixed effects)
scaa_info = digifish
data_names = c("agg_catch", "catch_paa", "agg_indices","index_paa")
scaa_info[data_names] = sim_pop[data_names]

scaa_NAA_re = list(N1_pars = NAA_re$N1_pars)
#scaa_NAA_re$use_steepness = 0
scaa_NAA_re$recruit_model = 1 #recruitments as fixed effects

scaa_input = prepare_wham_input(basic_info = scaa_info, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M, recruit_model = 1)

#fit the scaa estimation model
scaa_fit = fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
scaa_fit$mohns_rho = mohns_rho(scaa_fit)

######################################
#now do a closed-loop simulation with the beverton-holt operating model where catch advice is set during the projection years.
#catch advice is set as the average catch over 5 years projected using F40

#first we use set up an input for the full length of base and management experiment years. Earlier we specified 40 base years in the model. 
# Now we will just create an input that extends this 39 more years for the management experiment.
digifish <- make_digifish(years = 1975:(2014+39))
bh_input_full = prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re_bh, M = M)
bh_input_full$data$n_years_model #79

# om_input <- prepare_projection(bh_om, proj.opts=list(n.yrs=39, use.last.F=FALSE, use.avg.F=FALSE,
#                 use.FXSPR=FALSE, proj.F=NULL, avg.yrs=NULL, proj.catch = rep(1,39),
#                 cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL))
#then we use fit_wham without fitting to set up a TMB model that can be used to simulate a population time series
# and associated catch and index observations with the stochastic assumptions made in the NAA_re options to generate bh_om
om <- fit_wham(bh_input_full, do.fit = F, MakeADFun.silent = TRUE)
# currently the om sets F in all years = 0.2
exp(om$rep$log_F_tot)

#First we generate a time series of the population with the stochastic assumptions made in the NAA_re options to generate bh_om
#a key component is to keep the same seed for to resimulate the population as F changes in the management experiment years to make sure the past population stays the same
set.seed(8675309)
pop_om_sim = list(y0 = om$simulate(complete=TRUE))


ab <- NAA_re_bh$recruit_pars[[1]]
ssb <- seq(1,400000, 100)
plot(pop_om_sim[[1]]$SSB[-bh_input_full$data$n_years_model], pop_om_sim[[1]]$NAA[1,1,-1,1], xlim = c(0,max(ssb)), ylab = "Recruits (1000s)", xlab = "SSB (mt)")
points(pop_om_sim[[1]]$SSB[-bh_input_full$data$n_years_model], pop_om_sim[[1]]$pred_NAA[1,1,-1,1], col = 'red')
lines(seq(1,400000, 100), ab[1]*ssb/(1 + ab[2]*ssb), col = 'red')


digifish_t <- make_digifish(years = 1975:(2014+0))
om_input_t = prepare_wham_input(basic_info = digifish_t, selectivity = selectivity, NAA_re = NAA_re_bh, M = M)
om_t <- fit_wham(om_input_t, do.fit = F)
om_t_proj <- project_wham(om_t)
dim(om_t_proj$env$parameters$log_NAA)

#create a scaa input with 2014 terminal year
scaa_info = digifish_t
n_y <- length(1975:(2014+0))

#data_names = c("agg_catch", "catch_paa", "agg_indices","index_paa")
catch_info <- list(agg_catch = pop_om_sim[[1]]$agg_catch[1:n_y,,drop = F], catch_paa = pop_om_sim[[1]]$catch_paa[,1:n_y,,drop = F])
index_info <- list(agg_indices = pop_om_sim[[1]]$agg_indices[1:n_y,,drop = F], index_paa = pop_om_sim[[1]]$index_paa[,1:n_y,,drop = F])
#scaa_info[data_names] = sim_pop[data_names]

scaa_NAA_re = list(N1_pars = NAA_re$N1_pars)
scaa_NAA_re$recruit_model = 1 #recruitments as fixed effects

scaa_input = prepare_wham_input(basic_info = scaa_info, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M, recruit_model = 1, catch_info = catch_info, index_info = index_info)

#fit the scaa estimation model
scaa_fit = fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
scaa_fit$mohns_rho <- mohns_rho(scaa_fit)

scaa_proj <- project_wham(scaa_fit)

t_input <- bh_input_full
t_input$par$log_NAA[] <- pop_om_sim[[1]]$log_NAA
t_om <- fit_wham(t_input, do.fit = F)
#now replace the realized abundance parameters in the operating model input
om_input$par$log_NAA[] = pop_om_base_period$log_NAA
# and reset the operating model
om <- fit_wham(om_input, do.sdrep=F, do.retro=F, do.osa=F, do.check=F, do.proj=F, do.fit = F,
  MakeADFun.silent = TRUE, save.sdrep=FALSE)
plot(om$years_full, log(om$rep$NAA[,1]), type = 'l', xlab = "Year", ylab = "log(Recruits (1000s))")
lines(om$years_full, log(temp$rep$NAA[,1]), col = "red")

#Define how often an assessment (SCAA) and catch advice will be made
assess.interval = 3 #years.step = integer()
#catch advice
advice = rep(1,39)

#the first advice model will be completed at the end of the base period
#the SCAA model defined as above will be updated every 3 years (assess.interval) of the projection/evaluation period
scaa_step = fit_wham(scaa_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
#the SCAA model with actual projections to make catch advice using F40
scaa_step_proj <- project_wham(scaa_step, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
            use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL, avg.rec.yrs = tail(scaa_step$years,30),
            cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE, MakeADFun.silent=TRUE)
NAAtrack = matrix(NA, 79, 13)
#making advice for the first 3 years of the evaluation period
advice[1:assess.interval] = rep(mean(scaa_step_proj$rep$pred_catch[scaa_step_proj$input$data$n_years_model + 1:5]), assess.interval)
for(y in seq(3,39,3)) {
    #set up projection of operating model
    om_input <- prepare_projection(bh_om, proj.opts=list(n.yrs=39, use.last.F=FALSE, use.avg.F=FALSE,
                use.FXSPR=FALSE, proj.F=NULL, proj.catch=advice, avg.yrs=NULL,
                cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL))
    temp <- fit_wham(om_input, n.newton=n.newton, do.sdrep=F, do.retro=F, do.osa=F, do.check=F, do.proj=F, 
    MakeADFun.silent = TRUE, save.sdrep=FALSE, do.fit = F)
    set.seed(8675309)
    updated_sim = temp$simulate(complete=TRUE)
    om_input$par$log_NAA = updated_sim$log_NAA
    NAAtrack[,y/3] = updated_sim$NAA[,1]
    lines(om$years_full, log(updated_sim$NAA[,1]), col = viridisLite::viridis(n= 39)[y], lty = 2)

    #increase terminal year
    scaa_info_step = make_digifish(min(digifish$years):(max(digifish$years)+y))
    scaa_info_step$agg_catch = rbind(updated_sim$agg_catch, updated_sim$agg_catch_proj[1:y,,drop=F])
    scaa_info_step$agg_indices = rbind(updated_sim$agg_indices, updated_sim$agg_indices_proj[1:y,,drop=F])
    scaa_info_step$catch_paa = abind::abind(updated_sim$catch_paa, updated_sim$catch_paa_proj[,1:y,,drop=F], along = 2)
    scaa_info_step$index_paa = abind::abind(updated_sim$index_paa, updated_sim$index_paa_proj[,1:y,,drop=F], along = 2)

    step_input = prepare_wham_input(basic_info = scaa_info_step, selectivity = selectivity, NAA_re = scaa_NAA_re, M = M)
    scaa_step = fit_wham(step_input, do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE)
    scaa_step_proj <- project_wham(scaa_step, proj.opts=list(n.yrs=5, use.last.F=FALSE, use.avg.F=FALSE,
                use.FXSPR=TRUE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL, avg.rec.yrs = tail(scaa_step$years,30),
                cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL), save.sdrep=FALSE, MakeADFun.silent=TRUE)

    advice[y + 1:assess.interval] = rep(mean(scaa_step_proj$rep$pred_catch[scaa_step_proj$input$data$n_years_model + 1:5]), assess.interval)

}

# purple: OM (ie simulated pop with F = 0)
# teal: EM with catch advice = F40
png(file.path(write.dir,"om_mse.png"), width=7, height=7, res=300, units='in')
par(mfrow = c(2,2))
pal = viridisLite::viridis(n=3)
plot(temp$years_full, log(pop_om_base_period$NAA[,1]), type = 'l', xlab = "Year", ylab = "log(recruits)", col = pal[1]) 
lines(temp$years_full, log(NAAtrack[,13]), col = pal[2]) 
abline(v = max(temp$years), lty=2)

plot(temp$years_full, pop_om_base_period$SSB, type = 'l', xlab = "Year", ylab = "SSB", ylim = c(0,450000), col  = pal[1]) 
lines(temp$years_full, updated_sim$SSB, col = pal[2]) 
abline(v = max(temp$years), lty=2)

plot(temp$years_full, pop_om_base_period$pred_catch, type = 'l', xlab = "Year", ylab = "Catch (mt)", ylim = c(0,80000), col = pal[1]) 
lines(temp$years_full, updated_sim$pred_catch, col = pal[2]) 
points(temp$years_full, scaa_info_step$agg_catch, col = pal[2])
abline(v = max(temp$years), lty=2)

y.max = max(1.1*updated_sim$FAA_tot[,10]/exp(updated_sim$log_FMSY), na.rm=T)
plot(temp$years_full, updated_sim$FAA_tot[,10]/exp(updated_sim$log_FMSY), type = 'l', xlab = "Year", ylim=c(0,y.max), ylab = "F/FMSY", col = pal[2])
abline(v = max(temp$years), lty=2)
abline(h = 1, col = 'red', lty=2)
dev.off()
