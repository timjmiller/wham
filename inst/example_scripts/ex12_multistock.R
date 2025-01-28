# Vignette 12: Multiple regions, stocks, and movement

lib.loc <- NULL
#lib.loc <- "c:/work/wham/old_packages/lab"
library("wham", lib.loc = lib.loc)

path_to_examples <- system.file("extdata", package="wham")

#########################################################################################################
#fit the same model twice, simultaneously

two_stocks_asap <- read_asap3_dat(file.path(path_to_examples,c("ex1_SNEMAYT.dat","ex1_SNEMAYT.dat")))
ini_2_stocks <- prepare_wham_input(two_stocks_asap)
fit_2_stocks <- fit_wham(ini_2_stocks, do.osa = FALSE, do.retro= FALSE, do.sdrep = FALSE)
#saveRDS(fit_2_stocks, file.path(res_dir,"vign_5_fit_2_stocks.RDS"))

fit_2_stocks$rep$SSB


#########################################################
#different stocks, no movement
diff_stocks_asap <- read_asap3_dat(file.path(getwd(),"data",c("north.dat","south.dat")))
selectivity <- list(model = rep(c("logistic", "age-specific"),c(8,4)), n_selblocks = 12,
	fix_pars = c(rep(list(NULL),8), list(2:8,3:8,3:8,2:8)),
	initial_pars = c(rep(list(c(2,0.2)),8),list(rep(c(0.5,1),c(1,7)), rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(1,7)))))
diff_stocks_input <- prepare_wham_input(diff_stocks_asap, selectivity = selectivity)
fit_diff_stocks <- fit_wham(diff_stocks_input, do.osa = FALSE, do.retro= FALSE, do.sdrep = FALSE)

plot_wham_output(fit_diff_stocks)

###################
# change to 5 seasons, spawning is at 0.5 in middle of season 3
#movement for northern stock only, fish must move back prior to spawning

basic_info$n_seasons <- 5L
basic_info$fracyr_seasons <- rep(1/5,5)
basic_info$spawn_seasons <- c(3,3)
basic_info$fracyr_SSB <- fracyr_ssb-2/5 # this should be fixed in prepare_wham_input

n_ages <- length(basic_info$ages)
#each age other than 1 (recruitment) for north stock can be in either region on Jan 1 
basic_info$NAA_where <- array(1, dim = c(2,2,n_ages))
basic_info$NAA_where[1,2,1] <- 0 #stock 1, age 1 can't be in region 2 on Jan 1
basic_info$NAA_where[2,1,] <- 0 #stock 2, any age can't be in region 1 2 on Jan 1 (stock 2 doesn't move) 

n_seasons <- basic_info$n_seasons

move = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't
move$must_move <- array(0,dim = c(2,n_seasons,2))	
#if north stock in region 2 (south) must move back to region 1 (north) at the end of interval 2 before spawning
move$must_move[1,2,2] <- 1 #stock 1 must move at the end of season 2 from region 2
move$can_move <- array(0, dim = c(2,n_seasons,2,2))
move$can_move[1,c(1,4:5),1,2] <- 1 #only north stock can move in seasons after spawning
move$can_move[1,2,2,] <- 1 #north stock can (and must) move in last season prior to spawning back to north 
mus <- array(0, dim = c(2,n_seasons,2,1))
mus[1,1:n_seasons,1,1] <- 0.3 #initial value proportion moving to south = 0.3 (mean of prior)
mus[1,1:n_seasons,2,1] <- 0.3 #initial value proportion north to south = 0.3 (not intended to be used)
move$mean_vals <- mus 
move$mean_model <- matrix("stock_constant", 2,1)

input_move <- prepare_wham_input(
	basic_info = basic_info,
	NAA_re = NAA_list,
	selectivity = selectivity, 
	catch_info = catch_info, 
	index_info = index_info, 
	M = M_in, 
	F = F_in, 
	catchability = q_in,
	move = move)

nofit_move <- fit_wham(input_move, do.fit = FALSE, do.brps = FALSE)
saveRDS(nofit_move, file.path(res_dir,"vign_5_nofit_move.RDS"))

#first season
nofit_move$rep$seasonal_Ps_terminal[1,1,8,,]
#second season
nofit_move$rep$seasonal_Ps_terminal[1,2,8,,]
#3rd season
nofit_move$rep$seasonal_Ps_terminal[1,3,8,,]
#4th season
nofit_move$rep$seasonal_Ps_terminal[1,4,8,,]
#5th season
nofit_move$rep$seasonal_Ps_terminal[1,5,8,,]

x <- input_move$par$trans_mu
x[] <- as.integer(input_move$map$trans_mu)
x[1,,2,1] <- NA
input_move$map$trans_mu <- factor(x)

fit_move <- fit_wham(input_move, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
# fit_move <- do_reference_points(fit_move, do.sdrep = TRUE)
saveRDS(fit_move, file.path(res_dir,"vign_5_fit_2_stocks_move.RDS"))

#estimated movement
fit_move$rep$mu[1,8,1,1,1,2]

###################
#age-specific RE
move_age <- move
move_age$age_re <- matrix("none",2,1)
move_age$age_re[1,1] <- "iid"

input_move_age <- prepare_wham_input(
	basic_info = basic_info,
	NAA_re = NAA_list,
	selectivity = selectivity, 
	catch_info = catch_info, 
	index_info = index_info, 
	M = M_in, 
	F = F_in, 
	catchability = q_in,
	move = move_age)

length(unique(input_move_age$map$mu_re))

nofit_move_age <- fit_wham(input_move_age, do.fit = FALSE, do.brps = FALSE)
saveRDS(nofit_move_age, file.path(res_dir,"vign_5_nofit_move_age.RDS"))

#age-specific movement rates assuming initial values for fixed effects
nofit_move_age$rep$mu[1,1:8,1,1,1,2]

x <- input_move_age$par$trans_mu
x[] <- as.integer(input_move_age$map$trans_mu)
x[1,,2,1] <- NA
input_move_age$map$trans_mu <- factor(x)

input_move_age$par <- fit_move$parList
fit_move_age <- fit_wham(input_move_age, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
saveRDS(fit_move_age, file.path(res_dir,"vign_5_fit_2_stocks_move_age.RDS"))

fit_move$opt$obj
fit_move_age$opt$obj
#same

#the log(sd) of the RE
fit_move_age$parList$mu_repars[1,1,1,1,1]

###################
#year-specific RE
move_year <- move
move_year$year_re <- matrix("none",2,1)
move_year$year_re[1,1] <- "iid"

input_move_year <- prepare_wham_input(
	basic_info = basic_info,
	NAA_re = NAA_list,
	selectivity = selectivity, 
	catch_info = catch_info, 
	index_info = index_info, 
	M = M_in, 
	F = F_in, 
	catchability = q_in,
	move = move_year)

length(unique(input_move_year$map$mu_re))

nofit_move_year <- fit_wham(input_move_year, do.fit = FALSE, do.brps = FALSE)
saveRDS(nofit_move_year, file.path(res_dir,"vign_5_nofit_move_year.RDS"))
#age-specific movement rates assuming initial values for fixed effects
nofit_move_year$rep$mu[1,8,1,,1,2]

x <- input_move_year$par$trans_mu
x[] <- as.integer(input_move_year$map$trans_mu)
x[1,,2,1] <- NA
input_move_year$map$trans_mu <- factor(x)

input_move_year$par <- fit_move$parList
fit_move_year <- fit_wham(input_move_year, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
saveRDS(fit_move_year, file.path(res_dir,"vign_5_fit_2_stocks_move_year.RDS"))

fit_move$opt$obj
fit_move_year$opt$obj
#not the same

#the log(sd) of the RE
fit_move_year$parList$mu_repars[1,1,1,1,1]

plot(fit_move_year$years,fit_move_year$rep$mu[1,8,1,,1,2], ylab = "Move North to South", xlab = "Year")

###################
#prior distribution on movement parameters 


#just use prior once because it is constant over all seasons.
move_prior <- move
move_prior$use_prior <- array(0, dim = c(2,n_seasons,2,1))
#movement is not used in first season, but the parameter is constant across seasons. Could use it in any (single) season
move_prior$use_prior[1,1,1,1] <- 1
# sd on logit scale
move_prior$prior_sigma <- array(0, dim = c(2,n_seasons,2,1))
move_prior$prior_sigma[1,1,1,1] <- 0.2

input_move_prior <- prepare_wham_input(
	basic_info = basic_info,
	NAA_re = NAA_list,
	selectivity = selectivity, 
	catch_info = catch_info, 
	index_info = index_info, 
	M = M_in, 
	F = F_in, 
	catchability = q_in,
	move = move_prior)

nofit_move_prior <- fit_wham(input_move_prior, do.fit = FALSE, do.brps = FALSE)

x <- input_move_prior$par$trans_mu
x[] <- as.integer(input_move_prior$map$trans_mu)
x[1,,2,1] <- NA
input_move_prior$map$trans_mu <- factor(x)


#fixed effect estimated in fit_move
ind <- names(fit_move$parList)[!names(fit_move$parList) %in% "trans_mu"]
input_move_prior$par[ind] <- fit_move$parList[ind]


fit_move_prior <- fit_wham(input_move_prior, do.osa = FALSE, do.retro = FALSE, do.sdrep = FALSE, do.brps = FALSE)
saveRDS(fit_move_prior, file.path(res_dir,"vign_5_fit_2_stocks_move_prior.RDS"))

#estimated movement
fit_move$rep$mu[1,8,1,1,1,2]
fit_move_prior$rep$mu[1,8,1,1,1,2]


prop_AA <- t(sapply(1:33, function(y) sapply(2:8, function(x) fit_move$rep$NAA[1,2,y,x]/sum(fit_move$rep$NAA[1,,y,x]))))
matplot(fit_move$years, prop_AA, type = 'l', ylab = "Proportion of northern stock in the south", xlab = "Year", lty = 1, col = 1:7)
legend("topright", legend = paste0("Age ", 2:8), col = 1:7, lty  = 1)

prop_AA <- t(sapply(1:33, function(y) sapply(2:8, function(x) fit_move_year$rep$NAA[1,2,y,x]/sum(fit_move_year$rep$NAA[1,,y,x]))))
matplot(fit_move$years, prop_AA, type = 'l', ylab = "Proportion of northern stock in the south", xlab = "Year", lty = 1, col = 1:7)
legend("topright", legend = paste0("Age ", 2:8), col = 1:7, lty  = 1)
