# WHAM example 10: Illustrate simulation capabilities.
# - Simulate data from a fitted model with input derived from an ASAP3 dat file.
# - Make an operating model, simulate data, and fit models.
# - Fit models to data simulated from a different model.
# - Make a wham input file **without an ASAP3 dat file**.
# - Configure different assumptions about the population for both the operating model and the estimating model.
# - Do a closed-loop simulation with operating model, fitting an estimating model, generating catch advice and incorporating it into the operating model.

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
  model_name="Ex 1: SNEMA Yellowtail Flounder") 

mod_1 <- fit_wham(input, do.osa = FALSE, do.retro=FALSE, MakeADFun.silent = TRUE) # don't do retro peels and OSA residuals, don't show output during optimization

saveRDS(mod_1, file = "mod_1.RDS")

names(mod_1)


mod_1$sdrep #the standard errors of the fixed effects

#the simulate function
mod_1$simulate

#make a function to simulate and fit the matching estimating model
sim_fn <- function(om, self.fit = FALSE, do.sdrep = FALSE){
  input <- om$input
  input$data = om$simulate(complete=TRUE)
  if(self.fit) {
    fit <- fit_wham(input, do.sdrep = do.sdrep, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
    return(fit)
  } else return(input) 
}
set.seed(123)
self_sim_fit <- sim_fn(mod_1, self.fit = TRUE)

png(file = "self_fit_1_ssb.png")
plot(mod_1$years, self_sim_fit$input$data$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod_1$years, self_sim_fit$rep$SSB, lty = 2, col = 'red')
dev.off()

# a way of saving results as you go to make sure you don't loose everything if there is a catastrophic error
set.seed(123)
n_sims <- 10 #more is better, but takes longer
sim_inputs <- replicate(n_sims, sim_fn(mod_1), simplify=FALSE)
res = list(reps = list(), par.est = list(), par.se = list(), adrep.est = list(), adrep.se =list())
for(i in 1:length(sim_inputs)){
  cat(paste0("i: ",i, "\n"))
  tfit = fit_wham(sim_inputs[[i]], do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
  res$reps[[i]] = tfit$rep
  # do.sdrep =TRUE to create these results
  # res$par.est[[i]] = as.list(tfit$sdrep, "Estimate")
  # res$par.se[[i]] = as.list(tfit$sdrep, "Std. Error")
  # res$adrep.est[[i]] = as.list(tfit$sdrep, "Estimate", report = TRUE)
  # res$adrep.se[[i]] = as.list(tfit$sdrep, "Std. Error", report = TRUE)
  saveRDS(res, "self_test_res.RDS")
}

#plot results for bias of SSB
true = sapply(sim_inputs, function(x) x$data$SSB)
est = sapply(res$reps, function(x) return(x$SSB))
SSB_rel_resid = est/true - 1
resid_cis = apply(SSB_rel_resid,1,mean) + apply(SSB_rel_resid,1,sd)*qnorm(0.975)*t(matrix(c(-1,1),2,44))/sqrt(n_sims)

png(file = "ssb_rel_bias.png")
plot(mod_1$years, apply(SSB_rel_resid,1,mean), ylim = c(-0.1,0.05), ylab = "Estimated Relative Bias SSB", xlab = "Year")
lines(mod_1$years, resid_cis[,1])
lines(mod_1$years, resid_cis[,2])
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
vign_10_3_input <- input #for vignette


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
    fit = fit_wham(input, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(input,om_data))
}

#compare observations for multinomial and logistic normal
set.seed(123)
sim = sim_fn(mod_1, input_em, do.fit = FALSE)
sim[[1]]$data$obs[1:10,] #data.frame of observations and identifying information
sim[[2]]$obsvec[4:9] #showing obsvec that is simulated (multinomial frequencies)
sim[[1]]$data$obsvec[4:9] #showing obsvec that is used (MVN = logistic normal)
paa <- sim[[2]]$index_paa[1,1,]
paa #multinomial proportions : naa/sum(naa)
sim[[2]]$index_Neff[1] * paa # = simulated obsvec (multinomial frequencies)
log(paa[-6]) - log(paa[6]) # = inverse additive logistic transform for obsvec that is used
vign_10_3a_sim <- sim #for vignette

set.seed(123)
simfit = sim_fn(mod_1, input_em, do.fit = TRUE)

png(file = "sim_fit_3a.png")
plot(mod_1$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod_1$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
dev.off()


#3b: alternative effective sample size
sim_fn = function(om, Neff_om = 50, Neff_em = 200, do.fit = FALSE){
  om_input = om$input #the original input
  om_input$data$index_Neff[] = Neff_om #true effective sample size for the indices
  om_input$data$catch_Neff[] = Neff_om #true effective sample size for the indices
  om_input$par = om$parList #as.list(om$sdrep, "Estimate") #assume parameter values estimated from original model
  om_alt = fit_wham(om_input, do.fit = FALSE, MakeADFun.silent=TRUE) #make unfitted model for simulating data
  om_data = om_alt$simulate(complete = TRUE)
  em_input = om_input #the same configuration other than Neff
  em_input$data$index_Neff[] = Neff_em #assumed effective sample size for the indices
  em_input$data$catch_Neff[] = Neff_em #assumed effective sample size for the indices
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  em_input$data[obs_names] = om_data[obs_names]
  em_input <- set_osa_obs(em_input) #to rescale the frequencies appropriately
  if(do.fit) {
    fit = fit_wham(em_input, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(em_input,om_data))
}
set.seed(123)
simfit = sim_fn(mod_1, do.fit = FALSE)
sim[[1]]$data$obs[1:10,] #data.frame of observations and identifying information
sim[[2]]$obsvec[4:9] #showing obsvec that is simulated (multinomial frequencies)
sim[[1]]$data$obsvec[4:9] #showing obsvec that is used (MVN = logistic normal)
paa <- sim[[2]]$index_paa[1,1,] #simulated paa
paa #multinomial proportions : naa/sum(naa)
sim[[1]]$data$index_paa[1,1,] #paa in estimating model (the same)
sim[[2]]$index_Neff[1] * paa # = simulated obsvec (multinomial frequencies)
sim[[1]]$data$index_Neff[1] * paa # = em obsvec (multinomial frequencies with larger Neff)

vign_10_3b_sim <- sim #for vignette

set.seed(123)
simfit = sim_fn(mod_1, do.fit = TRUE)

png(file = "sim_fit_3b.png")
plot(mod_1$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod_1$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
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

om_ss <- fit_wham(input_om, do.osa = FALSE, do.retro=FALSE, MakeADFun.silent = TRUE) # don't do retro peels and OSA residuals, don't show output during optimization

#no differences in treatment of observations
sim_fn = function(om, input, do.fit = FALSE){
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  om_data = om$simulate(complete = TRUE)
  input$data[obs_names] = om_data[obs_names]
  if(do.fit) {
    fit = fit_wham(input, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
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
    fit = fit_wham(input, do.sdrep = FALSE, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
    return(list(fit, om_data))
  } else return(list(input,om_data))
}
set.seed(123)
simfit = sim_fn(mod_1, input_em, do.fit = TRUE)


png(file = "sim_fit_3d_1.png")
plot(mod_1$years, simfit[[2]]$SSB, type = 'l', ylab = "SSB", xlab = "Year")
lines(mod_1$years, simfit[[1]]$rep$SSB, lty = 2, col = 'red')
dev.off()

png(file = "sim_fit_3d_2.png")
plot(mod_1$years, simfit[[2]]$SSB/exp(simfit[[2]]$log_SSB_FXSPR[,2]), type = 'l', ylab = "SSB/SSB(F40)", xlab = "Year")
lines(mod_1$years, simfit[[1]]$rep$SSB/exp(simfit[[1]]$rep$log_SSB_FXSPR[,2]), lty = 2, col = 'red')
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

    F_info <- list()
    if(Fhist == "updown") F_info$F <- matrix(0.2 + c(seq(0,0.4,length.out = mid),seq(0.4,0,length.out=nby-mid)),nby, catch_info$n_fleets)

    #down then up
    if(Fhist == "downup") F_info$F <- matrix(0.2 + c(seq(0.4,0,length.out = mid),seq(0,0.4,length.out=nby-mid)),nby, catch_info$n_fleets)

    if(n_feedback_years>0) F_info$F <- rbind(F_info$F, F_info$F[rep(nby, n_feedback_years),, drop = F]) #same F as terminal year for feedback period

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
    return(list(basic_info = info, catch_info = catch_info, index_info = index_info, F = F_info))
}

stock_om_info <- make_info()
stock_basic_info <- c(basic_info, stock_om_info$basic_info)
catch_info <- stock_om_info$catch_info
index_info <- stock_om_info$index_info
F_info <- stock_om_info$F

selectivity_om <- list(
  model = c(rep("logistic", catch_info$n_fleets),rep("logistic", index_info$n_indices)),
  initial_pars = rep(list(c(5,1)), catch_info$n_fleets + index_info$n_indices)) #fleet, index

M_om <- list(initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))))

NAA_re <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(stock_basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = list(exp(10)),
  sigma_vals = array(0.5, c(1,1,length(stock_basic_info$ages))) #sigma_R, only the value in the first age class is used.
)

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om,
  F = F_info
)

# stock_om <- fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)

# stock_om_4 <- stock_om
# save(stock_om_4, file = "stock_om_4.RData")


#4a. Setting NAA random effects variance parameters.

NAA_re_om <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(stock_basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  cor="2dar1", #random effects correlated among age classes and across time (separably)
  recruit_model = 2, #random effects with a constant mean
  recruit_pars = list(exp(10)),
  sigma_vals = array(c(0.6, rep(0.2,length(stock_basic_info$ages)-1)), c(1,1,length(stock_basic_info$ages))), #sigma_R and sigma_2+
  cor_vals = array(c(0.7, 0.4, 0.6), c(1,1,3)) #age, then year (2+, then recruitment)
)


stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om,
  F = F_info
)

stock_om_input$par$trans_NAA_rho[1:2]
wham:::gen.logit(c(0.7,0.4), -1, 1,s=2) #age first, year second


# 4b. Setting M random effects variance parameters

M_om <- list(
  #mean_model = "estimate-M", fixed, not estimated by default
  re_model = matrix("ar1_ay", 1,1),
  initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))),
  sigma_vals = matrix(0.2, 1,1),
  cor_vals = array(c(0.7,0.4),dim = c(1,1,2)))
  
NAA_re_om <- list(
  N1_pars = array(exp(10)*exp(-(0:(length(stock_basic_info$ages)-1))*M_om$initial_means[1]), dim = c(1,1,length(stock_basic_info$ages))),
  sigma = "rec", #random about mean
  recruit_model = 2, #random effects with a constant mean
  sigma_vals = array(0.6, c(1,1,length(stock_basic_info$ages))), #sigma_R
  recruit_pars = list(exp(10))
)

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om,
  F = F_info
)

#stock_om <- fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)

#all NA
stock_om_input$map$Mpars

exp(stock_om_input$par$log_NAA_sigma[1,1,1]) #Rsigma

exp(stock_om_input$par$M_repars[1,1,1]) #

stock_om_input$par$M_repars[1,1,2:3]
wham:::gen.logit(c(0.7,0.4), -1, 1) #age first, year second


M_om <- list(
  #mean_model = "estimate-M", fixed, not estimated by default
  re_model = matrix("ar1_ay", 1,1),
  initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages))),
  sigma_vals = matrix(0.2, 1,1),
  cor_vals = array(c(0.7,0.4),dim = c(1,1,2)))

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om,
  F = F_info
)

#all 1: one parameter estimated
stock_om_input$map$Mpars

# 4c. Setting selectivity random effects variance parameters

selectivity_om = list(
  model = c(rep("logistic", catch_info$n_fleets),rep("logistic", index_info$n_indices)),
  initial_pars = rep(list(c(5,1)), catch_info$n_fleets + index_info$n_indices),
  re = rep("2dar1", catch_info$n_fleets + index_info$n_indices),
  sigma_vals = c(0.4,0.1),
  cor_vals = cbind(rep(0.7,2),rep(0.4,2))
) 

#remove random effects on M
M_om <- list(
  initial_means = array(0.2, dim = c(1,1, length(stock_basic_info$ages)))
)
  
stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catch_info = catch_info,
  index_info = index_info,
  M = M_om,
  F = F_info
)

stock_om_input$par$sel_repars
c(log(0.4), wham:::gen.logit(c(0.7,0.4),-1,1))
c(log(0.1), wham:::gen.logit(c(0.7,0.4),-1,1))


# 4d. Setting catchability random effects variance parameters.

selectivity_om = list(
  model = c(rep("logistic", catch_info$n_fleets),rep("logistic", index_info$n_indices)),
  initial_pars = rep(list(c(5,1)), catch_info$n_fleets + index_info$n_indices)
) 

catchability_om = list(
  re = "ar1",
  sigma_vals = 0.1,
  cor_vals = 0.5) #could also define the mean q parameter here.

stock_om_input <- prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  catchability = catchability_om,
  catch_info = catch_info,
  index_info = index_info,
  M = M_om,
  F = F_info
)

stock_om_input$par$q_repars
c(log(0.1), wham:::gen.logit(0.5,-1,1))

stock_om = fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)

set.seed(123)
sim_q = stock_om$simulate(complete=TRUE)

png("sim_plot_4d.png")
plot(stock_om$years, sim_q$q[,1], type = 'l', ylab = "Catchability", xlab = "Year")
dev.off()


# 4e. Setting effects of Environmental covariates.

ecov_om = list(
    label = "Climate variable 1",
    process_model = "ar1",
    process_mean_vals = 0,
    process_sig_vals = 0.3,
    process_cor_vals = 0.5,
    mean = matrix(0, length(stock_basic_info$years), 1), #observations which will be simulated later
    logsigma = log(0.2), 
    year = stock_basic_info$years, 
    use_obs = matrix(1,length(stock_basic_info$years),1),  
    recruitment_how = matrix("controlling-lag-0-linear",1,1),
    beta_R_vals = array(0.3, dim = c(1,1,1))
 )

stock_om_input = prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  M = M_om, 
  catch_info = catch_info,
  index_info = index_info,
  ecov = ecov_om,
  F = F_info
)

#or equivalently
#stock_om_input <- set_ecov(stock_om_input, ecov_om)

png("sim_plot_4e_1.png")
plot(stock_om$years, sim_Ecov$Ecov_x, type = 'l', ylab = "Covariate", xlab = "Year")
dev.off()

png("sim_plot_4e_2.png")
plot(stock_om$years, sim_Ecov$NAA[1,1,,1], type = 'l', ylab = "Recruitment", xlab = "Year")
dev.off()


# 4f. Estimation model with incorrect M assumption

stock_om_input = prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  M = M_om, 
  catch_info = catch_info,
  index_info = index_info,
  F = F_info
)


stock_om = fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)

M_em = list(
  initial_means = array(0.3, dim = c(1,1, length(stock_basic_info$ages)))
)

em_input = prepare_wham_input(
  basic_info = stock_basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  M = M_em, 
  catch_info = catch_info,
  index_info = index_info,
  F = F_info
)

sim_fn = function(om, em_input, do.fit = FALSE){
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  om_sim = om$simulate(complete = TRUE)
  em_input$data[obs_names] = om_sim[obs_names]
  if(do.fit) {
    em_fit = fit_wham(em_input, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE)
    return(list(em_fit, om_sim))
  } else return(list(em_input,om_sim))
}
set.seed(123)
simfit = sim_fn(stock_om, em_input, do.fit = TRUE)

png("sim_fit_4f_1.png")
res <- cbind(simfit[[2]]$SSB,simfit[[1]]$rep$SSB)
plot(stock_om$years, res[,1], type = 'l', ylab = "SSB", xlab = "Year", ylim = c(0, max(res)))
lines(stock_om$years, res[,2], lty = 2, col = 'red')
dev.off()

png("sim_fit_4f_2.png")
res <- cbind(simfit[[2]]$SSB/exp(simfit[[2]]$log_SSB_FXSPR[,1]), simfit[[1]]$rep$SSB/exp(simfit[[1]]$rep$log_SSB_FXSPR[,1]))
plot(stock_om$years, res[,1], type = 'l', ylab = "SSB/SSB(F40)", xlab = "Year", ylim = c(0, max(res)))
lines(stock_om$years, res[,2], lty = 2, col = 'red')
dev.off()


# 5. Closed-loop simulation

#define useful functions first.

#First, a function to figure out F given catch and abundance
get_F_from_catch <- function(om, year, catch, Finit = 0.1, maxF = 10){
    rep = om$rep
    naa = rep$NAA[1,1,year,]
    Maa = rep$MAA[1,1,year,]
    sel_tot = rbind(rep$FAA[,year,]/max(exp(rep$log_FAA_tot[year,]))) #matrix nfleets x n_ages 
    waa = rbind(om$input$data$waa[om$input$data$waa_pointer_fleets, year,]) #matrix nfleets x n_ages 
    nfleets <- length(om$input$data$waa_pointer_fleets)
    get_catch = function(log_F, naa, sel, waa, Maa){
        Faa = exp(log_F) * sel_tot
        Zaa = Maa + apply(Faa,2,sum)
        Catch = 0
        for(a  in 1:length(naa)) for(f in 1:nfleets) Catch = Catch + waa[f,a] * naa[a] * Faa[f,a] *(1 - exp(-Zaa[a]))/Zaa[a];
        return(Catch)
    }
    obj = function(log_F) (catch - get_catch(log_F, naa, sel_tot, waa, Maa))^2
    opt = try(nlminb(log(Finit), obj))
    if(!is.character(opt)) Fsolve = exp(opt$par)[1] else Fsolve = maxF
    if(Fsolve>10) Fsolve = maxF
    print(paste0("Fsolve: ", Fsolve))
    return(Fsolve)
}

#A function to modify the population with the F given the catch
update_om_F = function(om, year, catch){
    rep = om$rep #generate the reported values given the parameters
    year_ind = which(om$years == year) #index corresponding to year
    Fsolve = get_F_from_catch(om, year_ind, catch) #find the F for the catch advice

    #have to be careful if more than one fleet
    FAA = rbind(rep$FAA[,year_ind,]) #n_fleets x n_ages
    FAA_tot <- apply(FAA,2,sum)
    age_ind <- which(FAA_tot == max(FAA_tot))[1]
    selAA = FAA/FAA_tot[age_ind] #sum(sel_all[i,]) = 1
    FAA_catch = Fsolve * selAA
    F_fleet = apply(FAA_catch, 1, max) #full F for each fleet
    if(om$input$data$F_config==1) {
      if(year_ind>1) om$input$par$F_pars[year_ind-1,] <- log(F_fleet) - log(apply(rbind(rep$FAA[,year_ind-1,]),1,max)) #change the F_dev to produce the right full F
      else om$input$par$log_F1[] <- log(F_fleet) #if year is the first year of the model, change F in year 1
    } else{ #alternative configuration of F_pars
      om$input$par$F_pars[year_ind,] <- log(F_fleet)
    }
    om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)
    return(om)
}


#A function for both creating a simulated population and updating the simulations during the operating model during the feedback period.
update_om_fn = function(om, seed = 123, interval.info = NULL, random = "log_NAA"){
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")

  if(!is.null(interval.info)){ #iteratively update F over assessment interval for the given catch advice
    for(y in interval.info$years){
      om = update_om_F(om, year = y, catch = interval.info$catch) #put in the right F values
      set.seed(seed)
      om_sim = om$simulate(complete=TRUE) #resimulate the population and observations
      om$input$data[obs_names] = om_sim[obs_names] #update any simulated data
      om$input$par[om$input$random] = om_sim[om$input$random] #update any simulated random effects
      # reset the om
      om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)
    }
  } else { #otherwise just (re)generate the population
    set.seed(seed)
    om_sim = om$simulate(complete=TRUE) #resimulate the population and observations
    om$input$data[obs_names] = om_sim[obs_names] #update any simulated data
    om$input$par[random] = om_sim[random] #update any simulated random effects
    # reset the om
    om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)
  }
  return(om)
}

# A function to make the estimation model input in a way such that the terminal year can be updated appropriately through the feedback period.
make_em_input = function(M_em, sel_em, NAA_em, om_data, em_years){
  info = make_info(base_years = em_years) #update the terminal year for the estimation model
  ind_em = 1:length(em_years) #year indices
  #fill in the data from the operating model simulation
  info$catch_info$agg_catch = om_data$agg_catch[ind_em,, drop = FALSE]
  info$index_info$agg_indices = om_data$agg_indices[ind_em,, drop = FALSE]
  info$catch_info$catch_paa = om_data$catch_paa[,ind_em,, drop = FALSE]
  info$index_info$index_paa = om_data$index_paa[,ind_em,, drop = FALSE]
  em_input <- prepare_wham_input(
    basic_info = info$basic_info, 
    selectivity = sel_em, 
    NAA_re = NAA_em, 
    M = M_em, 
    catch_info = info$catch_info,
    index_info = info$index_info)
  return(em_input)
}

#A function that will take the estimated model and generate catch advice from it.
advice_fn = function(em){
  #make 5 year projections using F40. Use average SSB/R and YPR inputs over most recent 5 years
  proj_opts = list(n.yrs=5, use.FXSPR=TRUE, avg.yrs=tail(em$years,5)) 
  em_proj = project_wham(em, do.sdrep = FALSE, proj.opts = proj_opts, MakeADFun.silent=TRUE) #projected version of the em
  advice = mean(apply(em_proj$rep$pred_catch[length(em_proj$years) + 1:5,,drop = FALSE],1,sum)) #mean of the projected catch over the next 5 years fishing at F40
  print(advice)
  return(advice)
}

#Finally, a function that uses all the other functions and steps through the feedback period updating the operating model, refitting the estimation model, 
# and generating the catch advice.
loop_through_fn = function(om, M_em, selectivity_em, NAA_re_em, assess_years, assess_interval = assess.interval, base_years){
  catches = numeric() #save the catch advice
  for(i in assess_years){
    print(i)
    #make the input for the estimation model

    em_input = make_em_input(M_em, selectivity_em, NAA_re_em, om_data = om$input$data, em_years = base_years[1]:i)
    #fit the estimation model
    em = fit_wham(em_input, do.sdrep = FALSE, do.retro = FALSE, do.osa=FALSE, MakeADFun.silent = TRUE) #no feedback period yet
    #make the catch advice
    advice = advice_fn(em)
    catches = c(catches, advice)
    #set the catch for the next assess_interval years
    interval.info = list(catch = advice, years = i + 1:assess_interval)
    #update the operating model with the right Fs and resimulate the data given those Fs
    om = update_om_fn(om, seed = 123, interval.info = interval.info)  
  }
  return(list(om = om, em  = em, catches = catches))
}

#Create the stock operating model. 
# We will specify a feedback period of 40 years beyond the 40 year base period. 
# The operating model is a generic stock with biological inputs defined in the `make_info` function. 
# We also remove the `random` element so that the inner optimization of random effects is not performed when calling `fit_wham` so that any random effects 
# are kept at simulated values.

stock_om_info <- make_info(base_years = 1982:2021, n_feedback_years = 40)
basic_info <- stock_om_info$basic_info
catch_info <- stock_om_info$catch_info
index_info <- stock_om_info$index_info
F_info <- stock_om_info$F

stock_om_input = prepare_wham_input(
  basic_info = basic_info, 
  selectivity = selectivity_om, 
  NAA_re = NAA_re_om, 
  M = M_om, 
  catch_info = catch_info,
  index_info = index_info,
  F = F_info
)
stock_om_input$random <- NULL #so inner optimization won't change simulated RE
stock_om <- fit_wham(stock_om_input, do.fit = FALSE, MakeADFun.silent = TRUE)

#Define how often to perform an assessment (`assess.interval`) and the years they occur (`assess.years`). 
# We also do a first call to `update_om_fn` to create a simulated population for the operating model using the default seed 123. 
# This same seed is used iteratively in the closed loop to keep the same simulated values throughout the base period and up to each time catch advice is defined. 

assess.interval = 4
base.years = make_info()$basic_info$years #no feedback period yet
first.year = head(base.years,1)
terminal.year = tail(base.years,1)
assess.years = seq(terminal.year, tail(stock_om$years,1)-assess.interval,by = assess.interval)

stock_om = update_om_fn(stock_om)

# Running the `loop_through_fn` function will take a little while.

looped_res = loop_through_fn(stock_om, M_em = M_em, selectivity_em = selectivity_om, NAA_re_em = NAA_re, assess_years = assess.years, base_years = base.years)
looped_rep <- looped_res$om$rep

png("looped_rep_5_1.png")
res <- cbind(stock_om$rep$SSB, looped_rep$SSB)
plot(stock_om$years, res[,1], type = "l", ylim = c(0,max(res)), ylab = "SSB", xlab = "Year")
lines(stock_om$years, res[,2], col = "red")
dev.off()

png("looped_rep_5_2.png")
res <- cbind(looped_rep$SSB[1:76]/exp(looped_rep$log_SSB_FXSPR[1:76,1]), looped_res$em$rep$SSB[1:76]/exp(looped_res$em$rep$log_SSB_FXSPR[1:76,1]))
plot(stock_om$years[1:76], res[,1], type = "l", ylim = c(0,max(res)), ylab = "SSB/SSB40", xlab = "Year")
lines(stock_om$years[1:76], res[,2], col = "red")
dev.off()

