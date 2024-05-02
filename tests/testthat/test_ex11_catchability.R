# Test WHAM example 11: Illustrate alternative configurations for catchability.
#   - using priors and random effects for catchability.
#   - estimating environmental effects on catchability (and recuitment at the same time)
#   - estimating time-varying random effects on catchability.
#   - tests of nll, pars, retro taken out because changes to simulation order will generally cause these tests to fail.

# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all(compile=FALSE)
# btime <- Sys.time(); devtools::test(filter = "ex11_catchability"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~11 min

context("Ex 11: Priors and time-varying and environmental effects on catchability")

test_that("Ex 11 works",{
# get results to check NLL and par estimates
path_to_examples <- system.file("extdata", package="wham")
tmp.dir <- tempdir(check=TRUE)

#make a list of input components that prepare_wham_input can use to generate an input for fit_wham
#NOTE: specifying 2 indices
# 3. A simple operating model

make_digifish <- function(years = 1975:2014) {
    digifish = list()
    digifish$ages <- 1:10
    digifish$years <- years
    digifish$n_fleets <- 1
    na = length(digifish$ages)
    ny = length(digifish$years)

    digifish$maturity = array(t(matrix(1/(1 + exp(-1*(1:na - na/2))), na, ny)), c(1,ny,na))

    L = 100*(1-exp(-0.3*(1:na - 0)))
    W = exp(-11)*L^3
    nwaa = 1
    digifish$waa = array(t(matrix(W, na, ny)), dim = c(1, ny, na))

    digifish$fracyr_SSB <- cbind(rep(0.25,ny))

    digifish$bias_correct_process <- TRUE
    digifish$bias_correct_observation <- TRUE
    return(digifish)
}
digifish = make_digifish()

catch_info <- list()
catch_info$n_fleets = 1
catch_info$catch_cv = matrix(0.1, length(digifish$years), digifish$n_fleets)
catch_info$catch_Neff = matrix(200, length(digifish$years), digifish$n_fleets)
catch_info$selblock_pointer_fleets = t(matrix(1:digifish$n_fleets, digifish$n_fleets, length(digifish$years)))
index_info <- list()
index_info$n_indices <- 2
index_info$index_cv = matrix(0.3, length(digifish$years), index_info$n_indices)
index_info$index_Neff = matrix(100, length(digifish$years), index_info$n_indices)
index_info$fracyr_indices = matrix(0.5, length(digifish$years), index_info$n_indices)
index_info$index_units = rep(1, length(index_info$n_indices)) #biomass
index_info$index_paa_units = rep(2, length(index_info$n_indices)) #abundance
index_info$selblock_pointer_indices = t(matrix(digifish$n_fleets + 1:index_info$n_indices, index_info$n_indices, length(digifish$years)))
F_info <- list(F = matrix(0.2,length(digifish$years), catch_info$n_fleets))

selectivity = list(model = c(rep("logistic", digifish$n_fleets),rep("logistic", index_info$n_indices)),
    initial_pars = rep(list(c(5,1)), digifish$n_fleets + index_info$n_indices)) #fleet, index

M = list(initial_means = array(0.2, c(1,1,length(digifish$ages))))

NAA_re = list(
  N1_model = "age-specific-fe", 
  N1_pars = array(exp(10)*exp(-(0:(length(digifish$ages)-1))*M$initial_means[1]), c(1,1,length(digifish$ages))))
NAA_re$sigma = "rec" #random about mean
#NAA_re$use_steepness = 0
NAA_re$recruit_model = 2 #random effects with a constant mean
NAA_re$recruit_pars = list(exp(10))

##########################################################################################
# 4. Setting up the q parameter for the second index to have a prior distribution.
catchability = list(initial_q = rep(0.3, index_info$n_indices), prior_sd = c(NA, 0.3))

#make input and operating model
input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability,
  index_info = index_info, catch_info = catch_info, F = F_info))
om_input <- input
om_input$random <- NULL
om = suppressWarnings(fit_wham(om_input, do.fit = FALSE, MakeADFun.silent = TRUE))

#simulate data from operating model
set.seed(0101010)
newdata = om$simulate(complete=TRUE)
#attr(newdata,"check.passed") <- NULL

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, do.retro=FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
suppressWarnings(plot_wham_output(fit, dir.main=tmp.dir, plot.opts = list(browse = FALSE)))

#This is provided in plot_wham_output(). The true simulated q is shown as the solid vertical line.
wham:::plot_q_prior_post(fit)
abline(v = newdata$q[1,2], lwd = 2)
dev.off()
#mods <- list(fit1 = fit)

##########################################################################################
# 5. Add random effects on q for first index
catchability = list(prior_sd = c(NA, 0.3), initial_q = rep(0.3, index_info$n_indices), re = c("iid", "none"), sigma_val = c(0.3,0.3))
#time varying catchability on the first index, prior on the second.
#set value to simulate variation in q to 0.2, only first sigma_val is used.

input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability,
  index_info = index_info, catch_info = catch_info, F = F_info))
om_input <- input
om_input$random <- NULL
om = suppressWarnings(fit_wham(om_input, do.fit = FALSE, MakeADFun.silent = TRUE))

#simulate data from operating model
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))
suppressWarnings(plot_wham_output(fit, dir.main=tmp.dir, plot.opts = list(browse = FALSE)))

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
#mods$fit2 <- fit

##########################################################################################
# 6. Add Environmental covariates to the model
# First fit ecov but no effects of population
ecov = list(
        label = c("Climate variable 1", "Climate variable 2"),
        process_model = c("ar1","ar1"),
        mean = cbind(rnorm(length(digifish$years)),rnorm(length(digifish$years))), 
        logsigma = log(c(0.01, 0.2)), 
        years = digifish$years, 
        use_obs = matrix(1,length(digifish$years),2),  
        q_how = matrix("none",2,2))
        #q_how = matrix(c("lag-0-linear",rep("none",3)),1,2))

#set mean, sd, and rho of ecov processes
ecov$process_mean_vals = c(0,0) #mean
ecov$process_sig_vals = c(0.1,0.2) #sd
ecov$process_cor_vals <- c(0.4,-0.3)  #cor

#get rid of prior on second index. add AR1 random effects on q for first index
catchability = list(initial_q = rep(0.3, index_info$n_indices), sigma_val = c(0.4,0.4), re = c("iid", "none"))
#time varying catchability on the first index, prior on the second.
#set value to simulate variation in q to 0.2, only first sigma_val is used.

input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability, ecov = ecov,
  index_info = index_info, catch_info = catch_info, F = F_info))
om_input <- input
om_input$random <- NULL
om = suppressWarnings(fit_wham(om_input, do.fit = FALSE, MakeADFun.silent = TRUE))
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE))

suppressWarnings(plot_wham_output(fit, dir.main=tmp.dir, plot.opts = list(browse = FALSE)))

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
#estimated variability in q is lower than truth, but the SE suggests the CI includes the true value
as.list(fit$sdrep, "Std")$q_repars
#mods$fit3 <- fit

##########################################################################################
# 7. Environmental effects on catchability
# add some Ecovs to the model and allow effects of first Ecov on q for second index
ecov$q_how <- matrix(c("none", "none","lag-0-linear","none"),2,2)
# ecov = list(
#         label = c("Climate variable 1", "Climate variable 2"),
#         process_model = c("ar1","ar1"),
#         mean = cbind(rnorm(length(digifish$years)),rnorm(length(digifish$years))), 
#         logsigma = log(c(0.01, 0.2)), 
#         years = digifish$years, 
#         use_obs = matrix(1,length(digifish$years),2),  
#         #q_how = matrix("none",2,2))
#         q_how = matrix(c("none", "none","lag-0-linear","none"),2,2)) #n_Ecov x n_indices

#set mean, sd, and rho of ecov processes
# ecov$process_mean_vals = c(0,0) #mean
# ecov$process_sig_vals = c(0.1,0.2) #sd
# ecov$process_cor_vals <- c(0.4,-0.3)  #cor

#set value for Ecov_beta effect on q (dims are n_effects (n_indices, n_Ecov, max_n_poly)
ecov$beta_q_vals <- array(0, dim = c(index_info$n_indices, length(ecov$label), 1))
ecov$beta_q_vals[2,1,1] <- 0.5

#add AR1 random effects on q for first index
#catchability = list(initial_q = rep(0.3, index_info$n_indices), sigma_val = c(0.4,0.4), re = c("iid", "none"))

input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability, ecov = ecov,
  index_info = index_info, catch_info = catch_info, F = F_info))


#check value for Ecov_beta effect on q (dims are n_effects (n_indices, n_Ecov, max_n_poly)
input$par$Ecov_beta_q

om = suppressWarnings(fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE))
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))

# expect_equal(fit$opt$obj, ex11_tests$fit4$nll, tolerance=1e-6, scale=1)
# expect_equal(fit$opt$par, ex11_tests$fit4$par, tolerance=1e-6, scale=1)
# expect_equal(fit$mohns_rho, ex11_tests$fit4$mohns_rho, tolerance=1e-4, scale=1)

suppressWarnings(plot_wham_output(fit, dir.main=tmp.dir, plot.opts = list(browse = FALSE)))

#compare assumed and estimated ecov effect on q
input$par$Ecov_beta_q[2,1,1]
fit$parList$Ecov_beta_q[2,1,1]

#compare assumed and estimated ecov process pars
input$par$Ecov_process_pars
fit$parList$Ecov_process_pars


#compare true and estimated AR1 (log) sigma for q
input$par$q_repars
fit$parList$q_repars
#estimated variability in q is lower than truth, but estimate has large SE
as.list(fit$sdrep, "Std")$q_repars

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
#mods$fit4 <- fit

##########################################################################################
#5. add some Ecovs to the model and allow effects of first Ecov on q for second index AND recruitment
ecov$recruitment_how <- matrix(c("controlling-lag-0-linear","none"),2,1)
# ecov = list(
#         label = c("Climate variable 1", "Climate variable 2"),
#         process_model = c("ar1","ar1"),
#         mean = cbind(rnorm(length(digifish$years)),rnorm(length(digifish$years))), 
#         logsigma = log(c(0.01, 0.2)), 
#         years = digifish$years, 
#         use_obs = matrix(1,length(digifish$years),2),  
#         #q_how = matrix("none",2,2))
#         q_how = matrix(c("none", "none","lag-0-linear","none"),2,2), #n_Ecov x n_indices
#         recruitment_how = matrix(c("controlling-lag-0-linear","none"),2,1)) #n_Ecov x n_stocks

#set mean, sd, and rho of ecov processes
# ecov$process_mean_vals = c(0,0) #mean
# ecov$process_sig_vals = c(0.1,0.2) #sd
# ecov$process_cor_vals <- c(0.4,-0.3)  #cor

#set value for Ecov_beta effect on q (dims are (n_indices, n_Ecov, max_n_poly)
# ecov$beta_q_vals <- array(0, dim = c(index_info$n_indices, length(ecov$label), 1))
# ecov$beta_q_vals[2,1,1] <- 0.5

#set value for Ecov_beta effect on q (dims are (n_stocks x n_Ecov x max_n_poly)
ecov$beta_R_vals <- array(0, dim = c(1, length(ecov$label), 1))
ecov$beta_R_vals[1,1,1] <- -0.5

#add AR1 random effects on q for first index
# catchability = list(initial_q = rep(0.3, index_info$n_indices), sigma_val = c(0.4,0.4), re = c("iid", "none"))

input = suppressWarnings(prepare_wham_input(basic_info = digifish, selectivity = selectivity, NAA_re = NAA_re, M = M, catchability = catchability, ecov = ecov,
  index_info = index_info, catch_info = catch_info, F = F_info))

om = suppressWarnings(fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE))
set.seed(0101010)
newdata = om$simulate(complete=TRUE)

#put the simulated data in an input file with all the same configuration as the operating model
temp = input
temp$data = newdata

#fit estimating model that is the same as the operating model
fit = suppressWarnings(fit_wham(temp, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE))

# expect_equal(fit$opt$obj, ex11_tests$fit5$nll, tolerance=1e-6, scale=1)
# expect_equal(fit$opt$par, ex11_tests$fit5$par, tolerance=1e-6, scale=1)
# expect_equal(fit$mohns_rho, ex11_tests$fit5$mohns_rho, tolerance=1e-4, scale=1)

suppressWarnings(plot_wham_output(fit, dir.main=tmp.dir, plot.opts = list(browse = FALSE)))

#compare assumed and estimated ecov effect on q for second index
input$par$Ecov_beta_q[2,1,1]
fit$parList$Ecov_beta_q[2,1,1]

#compare assumed and estimated ecov effect on recruitment
input$par$Ecov_beta_R[1,1,1]
fit$parList$Ecov_beta_R[1,1,1]

#SE for beta parameters is large, especially for recruitment effect
as.list(fit$sdrep, "Std")$Ecov_beta_q[2,1,1]
as.list(fit$sdrep, "Std")$Ecov_beta_R[1,1,1]

#compare assumed and estimated ecov process pars
input$par$Ecov_process_pars
fit$parList$Ecov_process_pars


#compare true and estimated time-varying q
input$par$q_repars
fit$parList$q_repars

#estimated variability in q is lower than truth, but estimate has large SE
as.list(fit$sdrep, "Std")$q_repars

pal = viridisLite::viridis(n=2)
se = summary(fit$sdrep)
se = matrix(se[rownames(se) == "logit_q_mat",2], length(fit$years))
plot(fit$years, fit$rep$q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,0.6), ylab = "q", xlab = "Year")
for( i in 1:input$data$n_indices){
    lines(fit$years, fit$rep$q[,i], lwd = 2, col = pal[i])
  polyy = c(fit$rep$q[,i]*exp(-1.96*se[,i]),rev(fit$rep$q[,i]*exp(1.96*se[,i])))
    polygon(c(fit$years,rev(fit$years)), polyy, col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
    lines(fit$years, newdata$q[,i], lwd = 2, col = pal[i], lty = 2)
}
legend("topright", legend = paste0("Index ", rep(1:input$data$n_indices, each = 2), c(" Est.", " True")), lwd = 2, col = rep(pal, each = 2), lty = c(1,2))
dev.off()

#mods$fit5 <- fit
# ex11_test_results <- list()
# ex11_test_results$nll <- sapply(mods, function(x) x$opt$obj) #have to remove old model 3 and old model 7 is done differently now.
# ex11_test_results$par <- lapply(mods, function(x) x$opt$par) #have to save with different order and names
# ex11_test_results$mohns_rho <- lapply(mods, function(x) {
#  mr <- mohns_rho(x)
#  c(mr$SSB, mr$Fbar, mr$naa)
# saveRDS(ex11_test_results, file.path(path_to_examples,"ex11_test_results.RDS"))

})

# # remove files created during testing
# teardown(unlink(tmp.dir, recursive=TRUE))

