# load wham
is.repo <- try(pkgload::load_all(compile=FALSE)) #this is needed to run from repo without using installed version of wham
if(is.character(is.repo)) library(wham) #not using repo

# create directory for analysis, E.g.,
#write.dir <- "/path/to/save/output"
if(!exists("write.dir")) write.dir <- tempdir(check=TRUE)
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)
# read asap3 data file and convert to input list for wham
path_to_examples <- system.file("extdata", package="wham")
asap3 <- read_asap3_dat(file.path(path_to_examples,"ex2_SNEMAYT.dat"))
input <- prepare_wham_input(asap3) 

basic_info <- list(
    n_stocks = 1L,
    n_seasons = 1L,
    n_fleets = 1L,
    ages = 1:input$data$n_ages,
    n_fleets = input$data$n_fleets,
    fracyr_SSB = cbind(input$data$fracyr_SSB), #(n_years x n_stocks)
    maturity = input$data$mature, #(n_stocks x n_years x n_ages)
    years = as.integer(input$years), #(n_years)
    waa = input$data$waa, #(any no. x n_years x n_ages)
    waa_pointer_ssb = input$data$waa_pointer_ssb, #(n_stocks)
    spawn_regions = 1, #(n_stocks)
    spawn_seasons = 1 #(n_stocks)
)

catch_info <- list(
    n_fleets = NCOL(input$data$agg_catch), #(n_fleets)
    agg_catch = cbind(input$data$agg_catch), #(n_years x n_fleets)
    agg_catch_cv = cbind(sqrt(exp(input$data$agg_catch_sigma^2) - 1)), #(n_years x n_fleets)
    catch_paa = input$data$catch_paa, #(n_fleets x n_years x n_ages)
    use_catch_paa = cbind(input$data$use_catch_paa), #(n_years x n_fleets), 0: don't use, 1: use
    catch_Neff = cbind(input$data$catch_Neff), #(n_years x n_fleets)
    selblock_pointer_fleets = cbind(input$data$selblock_pointer_fleets), #(n_years x n_fleets)
    waa_pointer_fleets = input$data$waa_pointer_fleets, #(n_fleets)
    fleet_regions = rep(1,input$data$n_fleets) #(n_fleets)
)

index_info <- list(
    n_indices = NCOL(input$data$agg_indices),
    agg_indices = cbind(input$data$agg_indices), #(n_years x n_indices)
    units_indices = input$data$units_indices, #(n_indices) 1: biomass 2: numbers
    units_index_paa = input$data$units_index_paa, #(n_indices) 1: biomass 2: numbers
    agg_index_cv = cbind(sqrt(exp(input$data$agg_index_sigma^2) - 1)), #(n_years x n_indices)
    index_Neff = cbind(input$data$index_Neff), #(n_years x n_indices)
    fracyr_indices = cbind(input$data$fracyr_indices), #(n_years x n_indices)
    use_indices = cbind(input$data$use_indices), #(n_years x n_indices)
    use_index_paa = cbind(input$data$use_index_paa),  #(n_years x n_indices)
    index_paa = input$data$index_paa,  #(n_indices x n_years x n_ages)
    selblock_pointer_indices = cbind(input$data$selblock_pointer_indices), #(n_years x n_indices)
    waa_pointer_indices = input$data$waa_pointer_indices, #(n_indices)
    index_regions = rep(1,input$data$n_indices)) #(n_indices)

sel_info <- list(model = rep("logistic",6), n_selblocks = 6,
    fix_pars = list(NULL,NULL,NULL,NULL, 1:2, 1:2),
    initial_pars = list(c(2,0.2),c(2,0.2),c(2,0.2),c(2,0.2),c(1.5,0.1),c(1.5,0.1)))

#initial values for natural mortality, will not be estimated

MAA <- exp(matrix(c(input$par$Mpars), length(basic_info$years), length(basic_info$ages), byrow = TRUE) + 
    c(input$par$M_re[1,1,,])) #(nstocks x n_regions x n_years x n_ages)
M_info <- list(initial_MAA = array(MAA, dim = c(1,1,length(basic_info$years), length(basic_info$ages))))

#initial values for annual fully-selected fishing mortality
F_info <- list(F = matrix(0.3, length(basic_info$years), catch_info$n_fleets)) #(n_years x n_fleets)

#initial values for fully-selected catchability
q_info <- list(initial_q = rep(1e-6, index_info$n_indices))

#recruitment and "survival" random effects, initial N configuration is equilibrium (estimate initial recruitment and equilibrium fully-selected F)
NAA_info <- list(N1_pars = exp(input$par$log_N1))

#all at once 
input_all <- prepare_wham_input(
    basic_info = basic_info,
    selectivity = sel_info, 
    catch_info = catch_info, 
    index_info = index_info, 
    NAA_re = NAA_info,
    M = M_info, 
    F = F_info, 
    catchability = q_info)

#piece by piece
input_seq <- prepare_wham_input(basic_info = basic_info)
input_seq <- set_NAA(input_seq, NAA_re = NAA_info)
input_seq <- set_M(input_seq, M = M_info)
input_seq <- set_catch(input_seq, catch_info = catch_info)
input_seq <- set_F(input_seq, F = F_info)
input_seq <- set_indices(input_seq, index_info = index_info)
input_seq <- set_q(input_seq, catchability = q_info)
input_seq <- set_selectivity(input_seq, selectivity = sel_info)

input_asap <- prepare_wham_input(asap3, F = F_info, catchability = q_info, selectivity = sel_info) 

names(input_asap)

#compare 
nofit_all <- fit_wham(input_all, do.fit = FALSE)
nofit_seq <- fit_wham(input_seq, do.fit = FALSE)
nofit_asap <- fit_wham(input_asap, do.fit = FALSE)

c(nofit_all$fn(), nofit_seq$fn(),nofit_asap$fn()) #equal
c(sum(nofit_asap$par-nofit_all$par), sum(nofit_asap$par-nofit_seq$par)) #all equal

fit_seq <- fit_wham(input_seq, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
fit_all <- fit_wham(input_all, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
fit_asap <- fit_wham(input_asap, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE)
c(fit_seq$opt$obj, fit_all$opt$obj, fit_asap$opt$obj) # equal
c(sum(fit_asap$opt$par-fit_all$opt$par), sum(fit_asap$opt$par-fit_seq$opt$par)) #all equal

vign13_res <- list()
vign13_res$single_nofit = list(
        obj = c(nofit_all$fn(), nofit_seq$fn(),nofit_asap$fn()),
        par_diff = c(sum(nofit_asap$par-nofit_all$par), sum(nofit_asap$par-nofit_seq$par))
    )
vign13_res$single_fit = list(
        obj = c(fit_seq$opt$obj, fit_all$opt$obj, fit_asap$opt$obj),
        par_diff = c(sum(fit_asap$opt$par-fit_all$opt$par), sum(fit_asap$opt$par-fit_seq$opt$par))
    )

# ---------------------------------------------------------------
# two stocks, two regions, no movement
# create input using ASAP files, but then use individual elements to create an input by hand.
# without ASAP files, fill basic_info, catch_info, index_info, M_in, F_in, q_in, NAA_list components with appropriate values

diff_stocks_asap <- read_asap3_dat(file.path(path_to_examples,c("north.dat","south.dat")))
selectivity <- list(model = rep(c("logistic", "age-specific"),c(8,4)), n_selblocks = 12,
    fix_pars = c(rep(list(NULL),8), list(2:8,3:8,3:8,2:8)),
    initial_pars = c(rep(list(c(2,0.2)),8),list(rep(c(0.5,1),c(1,7)), rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(2,6)),rep(c(0.5,1),c(1,7)))))
input <- prepare_wham_input(diff_stocks_asap, selectivity = selectivity)


basic_info <- list(
    n_stocks = 2L,
    n_regions = 2L,
    region_names <- c("North_r", "South_r"),
    stock_names <- c("North_s", "South_s"),
    ages = 1:input$data$n_ages,
    n_seasons = 1L,
    n_fleets = input$data$n_fleets,
    fracyr_SSB = cbind(input$data$fracyr_SSB), #(n_years x n_stocks)
    maturity = input$data$mature, #(n_stocks x n_years x n_ages)
    years = as.integer(input$years),
    waa = input$data$waa, #(any no. x n_years x n_ages)
    waa_pointer_ssb = input$data$waa_pointer_ssb, #(n_stocks)
    spawn_regions = 1:2, #(n_stocks)
    spawn_seasons = c(1,1) #(n_stocks)
)

catch_info <- list(
    n_fleets = NCOL(input$data$agg_catch), #(n_fleets)
    agg_catch = cbind(input$data$agg_catch), #(n_years x n_fleets)
    agg_catch_cv = cbind(sqrt(exp(input$data$agg_catch_sigma^2) - 1)), #(n_years x n_fleets)
    catch_paa = input$data$catch_paa, #(n_fleets x n_years x n_ages)
    use_catch_paa = cbind(input$data$use_catch_paa), #(n_years x n_fleets), 0: don't use, 1: use
    catch_Neff = cbind(input$data$catch_Neff), #(n_years x n_fleets)
    selblock_pointer_fleets = cbind(input$data$selblock_pointer_fleets), #(n_years x n_fleets)
    waa_pointer_fleets = input$data$waa_pointer_fleets, #(n_fleets)
    fleet_regions = rep(c(1:2),c(2,2)) #(n_fleets)
)

index_info <- list(
    n_indices = NCOL(input$data$agg_indices),
    agg_indices = cbind(input$data$agg_indices), #(n_years x n_indices)
    units_indices = input$data$units_indices, #(n_indices) 1: biomass 2: numbers
    units_index_paa = input$data$units_index_paa, #(n_indices) 1: biomass 2: numbers
    agg_index_cv = cbind(sqrt(exp(input$data$agg_index_sigma^2) - 1)), #(n_years x n_indices)
    index_Neff = cbind(input$data$index_Neff), #(n_years x n_indices)
    fracyr_indices = cbind(input$data$fracyr_indices), #(n_years x n_indices)
    use_indices = cbind(input$data$use_indices), #(n_years x n_indices)
    use_index_paa = cbind(input$data$use_index_paa),  #(n_years x n_indices)
    index_paa = input$data$index_paa,  #(n_indices x n_years x n_ages)
    selblock_pointer_indices = cbind(input$data$selblock_pointer_indices), #(n_years x n_indices)
    waa_pointer_indices = input$data$waa_pointer_indices, #(n_indices)
    index_regions = rep(c(1:2),c(2,2))) #(n_indices)

#initial values for natural mortality, will not be estimated
MAA <- exp(input$par$M_re) #(nstocks x n_regions x n_years x n_ages)
for(i in 1:2) for(j in 1:2) MAA[i,j,,] <- MAA[i,j,,]*exp(matrix(input$par$Mpars[i,j,], length(basic_info$years), length(basic_info$ages), byrow = TRUE))
M_info <- list(initial_MAA = MAA)

#initial values for annual fully-selected fishing mortality
F_info <- list(F = matrix(0.3, length(basic_info$years), catch_info$n_fleets)) #(n_years x n_fleets)

#initial values for fully-selected catchability
q_info <- list(initial_q = rep(1e-6, index_info$n_indices))

#recruitment and "survival" random effects, initial N configuration is equilibrium (estimate initial recruitment and equilibrium fully-selected F), and inital median recruitment parameters
NAA_list <- list(sigma = "rec+1", N1_model = rep("equilibrium",2), recruit_pars = list(exp(10),exp(10)))

#all at once 
input_all <- prepare_wham_input(
    basic_info = basic_info,
    NAA_re = NAA_list,
    selectivity = selectivity, 
    catch_info = catch_info, 
    index_info = index_info, 
    M = M_info, 
    F = F_info, 
    catchability = q_info)

#piece by piece
input_seq <- prepare_wham_input(basic_info = basic_info)
input_seq <- set_NAA(input_seq, NAA_re = NAA_list)
input_seq <- set_M(input_seq, M = M_info)
input_seq <- set_catch(input_seq, catch_info = catch_info)
input_seq <- set_F(input_seq, F_info)
input_seq <- set_indices(input_seq, index_info = index_info)
input_seq <- set_q(input_seq, q_info)
input_seq <- set_selectivity(input_seq, selectivity = selectivity)

input_asap <- prepare_wham_input(diff_stocks_asap, selectivity = selectivity, NAA_re = NAA_list)

#compare 
nofit_all <- fit_wham(input_all, do.fit = FALSE, do.brps = FALSE)
nofit_seq <- fit_wham(input_seq, do.fit = FALSE, do.brps = FALSE)
nofit_asap <- fit_wham(input_asap, do.fit = FALSE, do.brps = FALSE)

c(nofit_all$fn(), nofit_seq$fn(),nofit_asap$fn()) #equal
c(sum(nofit_asap$par-nofit_all$par), sum(nofit_asap$par-nofit_seq$par)) #all equal


fit_seq <- fit_wham(input_seq, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE, do.brps = FALSE)
fit_all <- fit_wham(input_all, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE, do.brps = FALSE)
fit_asap <- fit_wham(input_asap, do.retro = FALSE, do.osa = FALSE, do.sdrep = FALSE, do.brps = FALSE)

c(fit_seq$opt$obj, fit_all$opt$obj, fit_asap$opt$obj) # equal
c(sum(fit_asap$opt$par-fit_all$opt$par), sum(fit_asap$opt$par-fit_seq$opt$par)) #all equal

res_dir <- file.path(getwd(),"temp")

vign13_res$multi_nofit = list(
        obj = c(nofit_all$fn(), nofit_seq$fn(),nofit_asap$fn()),
        par_diff = c(sum(nofit_asap$opt$par-nofit_all$opt$par), sum(nofit_asap$opt$par-nofit_seq$opt$par))
    )
vign13_res$multi_fit = list(
        obj = c(fit_seq$opt$obj, fit_all$opt$obj, fit_asap$opt$obj),
        par_diff = c(sum(fit_asap$opt$par-fit_all$opt$par), sum(fit_asap$opt$par-fit_seq$opt$par))
    )

saveRDS(vign13_res, file.path(write.dir,"vign13_res.RDS"))
