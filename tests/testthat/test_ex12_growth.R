# Test growth function

context("growth with SS simulated data")

test_that("growth works",{

  path_to_examples <- system.file("extdata", package="wham")

  data_file <- readRDS(file.path(path_to_examples,"ex12_data_file.rds"))
  SS_report <- readRDS(file.path(path_to_examples,"ex12_SS_report.rds"))
  envindex <- read.csv(file.path(path_to_examples,"PDO_annual_1970_2019.csv"))
  ex12_test_results <- readRDS(file.path(path_to_examples,"ex12_test_results.rds"))

  # Define some model info:
    min_year = SS_report$startyr # first year
    max_year = SS_report$endyr # last year
    n_ages = length(1:max(SS_report$agebins)) # number of ages
    fish_len = SS_report$lbins # length bins
    n_years = length(min_year:max_year) # number of years
    NAA_SS = SS_report$natage[SS_report$natage$`Beg/Mid` == 'B' & SS_report$natage$Yr >= min_year & SS_report$natage$Yr <= max_year, 14:ncol(SS_report$natage)] # abundance at age from SS
    LAA_SS = SS_report$growthseries[SS_report$growthseries$Yr >= min_year & SS_report$growthseries$Yr <= max_year, 6:ncol(SS_report$growthseries)] # mean length at age from SS
    LWpars = c(SS_report$Growth_Parameters$WtLen1, SS_report$Growth_Parameters$WtLen2) # LW parameters
    GWpars = c(SS_report$Growth_Parameters$K, SS_report$Growth_Parameters$Linf, 
               SS_report$endgrowth$Len_Beg[SS_report$endgrowth$Real_Age == 1], # This is L1
               SS_report$endgrowth$SD_Beg[2], SS_report$endgrowth$SD_Beg[11]) # growth parameters (vB function)
    Q_pars = c(exp(SS_report$parameters$Value[SS_report$parameters$Label == "LnQ_base_Fleet2(2)"])) # Q parameter
    # Selectivity parameters:
    int_selPos = grep(pattern =  "Size_inflection", x = SS_report$parameters$Label)[1]
    SelecParams = SS_report$parameters[int_selPos:nrow(SS_report$parameters), ]
    selpars1 = SelecParams$Value[1:2] # logistic
    selpars2 = SelecParams$Value[3:4] # logistic

  # Create input
    wham_data = list() # input data list
    # Basic information:
    wham_data$ages = 1:n_ages
    wham_data$lengths = fish_len
    wham_data$years = min_year:max_year
    wham_data$n_fleets = 1L
    wham_data$n_indices = 1L
    wham_data$Fbar_ages = 1L:10L
    wham_data$percentSPR = 40
    wham_data$percentFXSPR = 100
    wham_data$percentFMSY = 100
    wham_data$XSPR_R_avg_yrs = 1:n_years
    wham_data$XSPR_R_opt = 2
    wham_data$simulate_period = c(1,0)
    wham_data$bias_correct_process = 1
    wham_data$bias_correct_observation = 1
    wham_data$maturity = matrix(rep(SS_report$endgrowth[2:(n_ages+1),18], times = n_years),
                                ncol = n_ages, nrow = n_years, byrow = TRUE) 
    wham_data$fracyr_SSB = rep(0, times = n_years)
    # Agg catch information:
    wham_data$agg_catch = as.matrix(data_file$catch[2:(n_years+1), 4])
    wham_data$catch_cv = as.matrix(data_file$catch[2:(n_years+1), 5])
    # Length comps fishery:
    catch_pal_num = data_file$lencomp[data_file$lencomp$FltSvy == 1, 7:(7+length(wham_data$lengths)-1)]
    catch_pal_prop = t(apply(as.matrix(catch_pal_num),1, function(x) x/sum(x)))
    wham_data$catch_pal = catch_pal_prop
    wham_data$catch_NeffL = as.matrix(as.double(data_file$lencomp[data_file$lencomp$FltSvy == 1, 6]))
    wham_data$use_catch_pal = matrix(1, nrow = length(wham_data$years), ncol = wham_data$n_fleets)
    # Age comps fishery:
    wham_data$use_catch_paa = matrix(0, ncol = 1, nrow = n_years)
    # Agg indices information:
    wham_data$agg_indices = as.matrix(data_file$CPUE[,4])
    wham_data$index_cv = as.matrix(data_file$CPUE[,5])
    wham_data$units_indices = rep(0L, times = 1) # numbers
    wham_data$fracyr_indices = matrix(0, ncol = 1, nrow = n_years)
    # Length comps indices:
    index_pal_num = data_file$lencomp[data_file$lencomp$FltSvy == 2, 7:(7+length(wham_data$lengths)-1)]
    index_pal_prop = t(apply(as.matrix(index_pal_num),1, function(x) x/sum(x)))
    wham_data$index_pal = array(data = index_pal_prop, dim = c(1, dim(index_pal_prop)[1], dim(index_pal_prop)[2]))
    wham_data$index_NeffL =  as.matrix(as.double(data_file$lencomp[data_file$lencomp$FltSvy == 2, 6]))
    wham_data$use_index_pal = matrix(1, ncol = 1, nrow = n_years)
    # Age comps indices:
    wham_data$index_paa = array(NA, dim = c(wham_data$n_indices, n_years, n_ages))
    this_matrix = data_file$agecomp[,11:(10+n_ages)]
    wham_data$index_paa[1,,] = as.matrix(this_matrix/rowSums(this_matrix))
    wham_data$index_Neff = matrix(data_file$agecomp$Nsamp,
                                  nrow = n_years, ncol = 1)
    wham_data$use_index_paa = matrix(1L, nrow = n_years, ncol = 1)
    # Fishing mortality and selectivity:
    wham_data$selblock_pointer_fleets = matrix(1L, ncol = 1, nrow = n_years)
    wham_data$selblock_pointer_indices = matrix(2L, ncol = 1, nrow = n_years)
    wham_data$F = matrix(apply(X = SS_report$fatage[SS_report$fatage$Yr >= min_year &
                                                      SS_report$fatage$Yr <= max_year,
                                                    9:ncol(SS_report$fatage)], MARGIN = 1, FUN = max),
                         ncol = 1)
    # WAA pointers:
    wham_data$waa_pointer_indices = 1
    wham_data$waa_pointer_fleets = 2
    wham_data$waa_pointer_totcatch = 2
    wham_data$waa_pointer_ssb = 1
    wham_data$waa_pointer_jan1 = 1


  # Fit model 1
  my_input1 <- suppressWarnings(prepare_wham_input(model_name="whamGrowth1",
                               selectivity=list(model=rep("len-logistic",2), 
                                                re=rep("none",2), 
                                                initial_pars=list(selpars1,selpars2),
                                                fix_pars = list(2, 2),
                                                n_selblocks = 2),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                             recruit_model = 2,
                                             N1_pars = c(NAA_SS[1,1], 0),
                                             recruit_pars = mean(NAA_SS[,1])),
                               growth = list(model = 'vB_classic',
                                             re = c('ar1_y', 'none', 'none'),
                                             init_vals = GWpars[1:3],
                                             est_pars = 1:3,
                                             SD_vals = GWpars[4:5], SD_est = 1:2),
                               LW = list(init_vals = LWpars,
                                         re = c('none', 'none')),
                               catchability = list(re = 'none', initial_q = Q_pars, q_lower = 0,
                                                   q_upper = 10, prior_sd = NA),
                               basic_info = wham_data))

  # Extra change:
    my_input1$par$log_NAA = as.matrix(log(NAA_SS[-1,])) # set initial rec devs as OM
    my_input1$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
    my_input1$map$log_NAA_sigma = factor(NA) # fix sigma
    my_input1$map$log_N1_pars = factor(c(1,NA))
    my_input1$random = "growth_re"

  # Run model 1:  
  my_model_1 <- suppressWarnings(fit_wham(my_input1, do.osa = F, do.sdrep = F, do.retro=F, MakeADFun.silent = TRUE)) # turn off OSA residuals to save time


  # Fit model 2

  # ecov info:
    envindex = envindex[envindex$years >= min_year & envindex$years <= max_year, ]
    ecov <- list(
      label = "PDO",
      mean = as.matrix(envindex$index),
      logsigma = matrix(log(0.2), ncol = 1, nrow = n_years),
      year = envindex$years,
      use_obs = matrix(1, ncol=1, nrow=dim(envindex)[1]),
      lag = 0,
      ages = list(1:n_ages),
      process_model = 'ar1', # "rw" or "ar1"
      where = 'growth', 
      where_subindex = 1, # K
      how = 1)

  my_input3 <- suppressWarnings(prepare_wham_input(model_name="whamGrowth1",
                               selectivity=list(model=rep("len-logistic",2), 
                                                re=rep("none",2), 
                                                initial_pars=list(selpars1,selpars2),
                                                fix_pars = list(2, 2),
                                                n_selblocks = 2),
                               M = list(model = 'constant', re = 'none',
                                        initial_means = SS_report$Natural_Mortality[1,5],
                                        est_ages = 1),
                               NAA_re = list(sigma="rec", cor = 'iid', N1_model = 1,
                                             recruit_model = 2,
                                             N1_pars = c(NAA_SS[1,1], 0),
                                             recruit_pars = mean(NAA_SS[,1])),
                               growth = list(model = 'vB_classic',
                                             re = c('none', 'none', 'none'),
                                             init_vals = GWpars[1:3],
                                             est_pars = 1:3,
                                             SD_vals = GWpars[4:5], SD_est = 1:2),
                               ecov = ecov,
                               LW = list(init_vals = LWpars,
                                         re = c('none', 'none')),
                               catchability = list(re = 'none', initial_q = Q_pars, q_lower = 0,
                                                   q_upper = 10, prior_sd = NA),
                               basic_info = wham_data))

  # Extra change:
    my_input3$par$log_NAA = as.matrix(log(NAA_SS[-1,])) # set initial rec devs as OM
    my_input3$par$log_NAA_sigma = log(SS_report$sigma_R_in) # sigma as in SS
    my_input3$map$log_NAA_sigma = factor(NA) # fix sigma
    my_input3$map$log_N1_pars = factor(c(1,NA))
    my_input3$random = "Ecov_re"

  # Run model 1:  
  my_model_3 <- suppressWarnings(fit_wham(my_input3, do.osa = F, do.sdrep = F, do.retro=F, MakeADFun.silent = TRUE)) # turn off OSA residuals to save time

  # Test against prior objective function
  expect_equal( as.numeric(my_model_1$opt$objective), as.numeric(ex12_test_results$nll_mod1), tol = 1e-3 )
  expect_equal( as.numeric(my_model_3$opt$objective), as.numeric(ex12_test_results$nll_mod3), tol = 1e-3 )

})



