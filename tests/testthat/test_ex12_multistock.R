# Test WHAM example 12: 
# 1. Load wham library and set up a directory for work

# 2. Create inputs and fitting model for a two stock model without movement
#     - Load asap files
#     - create arguments to `prepare_wham_input` and model component functions (e.g., `set_move`)
#     - examine components of unfitted and fitted model

# 3. Create inputs and fitting model with movement of 1 stock
#     - Load asap files
#     - create arguments to `prepare_wham_input` and model component functions
#     - examine components of unfitted and fitted model

# 4. Create inputs and fitting model with age-varying movement
#     - Load asap files
#     - create arguments to `prepare_wham_input` and model component functions
#     - examine components of unfitted and fitted model

# 5. Create inputs and fitting model with time-varying movement
#     - Load asap files
#     - create arguments to `prepare_wham_input` and model component functions
#     - examine components of unfitted and fitted model

# 6. Create inputs and fitting model with a prior distribution on movement 
#     - create arguments to `prepare_wham_input` and model component functions
#     - examine components of input and fitted model

# pkgbuild::compile_dll(debug = FALSE); pkgload::load_all(compile = FALSE)
# btime <- Sys.time(); devtools::test(filter = "ex12_multistock"); etime <- Sys.time(); runtime = etime - btime; runtime;
# ~1.5 min

context("Ex 12: Multiple stocks")

test_that("Ex 12 works",{

suppressWarnings(suppressMessages({

path_to_examples <- system.file("extdata", package="wham")
two_stocks_asap <- read_asap3_dat(file.path(path_to_examples,c("ex1_SNEMAYT.dat","ex1_SNEMAYT.dat")))
ini_2_stocks <- prepare_wham_input(two_stocks_asap)

nofit_2_stock <- fit_wham(ini_2_stocks, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

diff_stocks_asap <- read_asap3_dat(file.path(path_to_examples,c("north.dat","south.dat")))
selectivity <- list(model = rep(c("logistic", "age-specific"),c(8,4)), n_selblocks = 12,
    fix_pars = c(rep(list(NULL),8), list(2:8,3:8,3:8,2:8)),
    initial_pars = c(rep(list(c(2,0.2)),8),list(rep(c(0.5,1),c(1,7)), rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(1,7)))))
diff_stocks_input <- prepare_wham_input(diff_stocks_asap, selectivity = selectivity)
nofit_2_diff_stock <- fit_wham(diff_stocks_input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)


###################
# change to 5 seasons, spawning is at 0.5 in middle of season 3
#movement for northern stock only, fish must move back prior to spawning
input <- diff_stocks_input
basic_info <- list(
    n_stocks = input$data$n_stocks,
    n_regions = input$data$n_regions,
    region_names <- c("North_r", "South_r"),
    stock_names <- c("North_s", "South_s"),
    ages = 1:input$data$n_ages,
    n_fleets = input$data$n_fleets,
    maturity = input$data$mature,
    years = as.integer(input$years),
    waa = input$data$waa,
    waa_pointer_ssb = input$data$waa_pointer_ssb,
    spawn_regions = input$data$spawn_regions
)

basic_info$n_seasons <- 5L
basic_info$fracyr_seasons <- rep(1/5,5)
basic_info$spawn_seasons <- c(3,3)
basic_info$fracyr_SSB <- input$data$fracyr_SSB#-2/5 # this should be fixed in prepare_wham_input

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

MAA <- exp(input$par$M_re)
for(i in 1:2) for(j in 1:2) MAA[i,j,,] <- MAA[i,j,,]*exp(matrix(input$par$Mpars[i,j,], length(basic_info$years), length(basic_info$ages), byrow = TRUE))
M_in <- list(initial_MAA = MAA)

F_in <- list(F = matrix(0.3, length(basic_info$years), input$data$n_fleets))
q_in <- list(initial_q = rep(1e-6, input$data$n_indices))

NAA_list <- list(sigma = "rec+1", N1_model = rep("equilibrium",2))

input_move <- prepare_wham_input(diff_stocks_asap,
    basic_info = basic_info,
    NAA_re = NAA_list,
    selectivity = selectivity, 
    # catch_info = catch_info, 
    # index_info = index_info, 
    M = M_in, 
    F = F_in, 
    catchability = q_in,
    move = move)

nofit_move <- fit_wham(input_move, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)


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


###################
#age-specific RE
move_age <- move
move_age$age_re <- matrix("none",2,1)
move_age$age_re[1,1] <- "iid"

input_move_age <- prepare_wham_input(diff_stocks_asap,
    basic_info = basic_info,
    NAA_re = NAA_list,
    selectivity = selectivity, 
    # catch_info = catch_info, 
    # index_info = index_info, 
    M = M_in, 
    F = F_in, 
    catchability = q_in,
    move = move_age)

# length(unique(input_move_age$map$mu_re))
expect_equal(length(unique(input_move_age$map$mu_re)), input_move_age$data$n_ages+1)

nofit_move_age <- fit_wham(input_move_age, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

#age-specific movement rates assuming initial values for fixed effects
nofit_move_age$rep$mu[1,1:8,1,1,1,2]


###################
#year-specific RE
move_year <- move
move_year$year_re <- matrix("none",2,1)
move_year$year_re[1,1] <- "iid"

input_move_year <- prepare_wham_input(diff_stocks_asap,
    basic_info = basic_info,
    NAA_re = NAA_list,
    selectivity = selectivity, 
    # catch_info = catch_info, 
    # index_info = index_info, 
    M = M_in, 
    F = F_in, 
    catchability = q_in,
    move = move_year)

length(unique(input_move_year$map$mu_re)) #+1 for NAs

expect_equal(length(unique(input_move_year$map$mu_re)), input_move_year$data$n_years_model+1) 

nofit_move_year <- fit_wham(input_move_year, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

nofit_move_year$rep$mu[1,8,1,,1,2]

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
move_prior$mean_vals[1,1,1,1] #transform of mean of the prior distribution

input_move_prior <- prepare_wham_input(diff_stocks_asap,
    basic_info = basic_info,
    NAA_re = NAA_list,
    selectivity = selectivity, 
    # catch_info = catch_info, 
    # index_info = index_info, 
    M = M_in, 
    F = F_in, 
    catchability = q_in,
    move = move_prior)

expect_equal(length(unique(input_move_prior$map$mu_re)), 1) 
expect_equal(length(unique(input_move_prior$map$mu_prior_re)), 2) #+1 for NAs

nofit_move_prior <- fit_wham(input_move_prior, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

}))
})

