#' Prepare input data and parameters for WHAM model
#'
#' Prepares data and parameter settings for  \code{\link{fit_wham}}, optionally using an ASAP3 data file read into R by \code{\link{read_asap3_dat}}.
#' By default, this will set up a SCAA version like \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP} without any penalties or random effects.
#' If any asap3, basic_info, catch_info, and index_info arguments are not provided, some respective arbitrary
#' assumptions are made. The various options described below, such as \code{NAA_re} and \code{selectivity},
#' can still be used without the asap3 object.
#'
#' \code{recruit_model} specifies the stock-recruit model. See \code{wham.cpp} for implementation.
#'   \describe{
#'     \item{= 1}{SCAA (without NAA_re option specified) or Random walk (if NAA_re$sigma specified), i.e. predicted recruitment in year i = recruitment in year i-1}
#'     \item{= 2}{(default) Random about mean, i.e. steepness = 1}
#'     \item{= 3}{Beverton-Holt}
#'     \item{= 4}{Ricker}
#'   }
#' Note: \code{\link{fit_wham}} allows fitting a SCAA model (\code{NAA_re = NULL}), which estimates recruitment in every year as separate fixed effect parameters,
#' but in that case no stock-recruit function is estimated. A warning message is printed if \code{recruit_model > 2} and \code{NAA_re = NULL}.
#' If you wish to use a stock-recruit function for expected recruitment, choose recruitment deviations as random effects,
#' either only age-1 (\code{NAA_re = list(sigma="rec")}) or all ages (\code{NAA_re = list(sigma="rec+1")}, "full state-space" model).
#' See below for details on \code{NAA_re} specification.
#'
#' \code{ecov} specifies any environmental covariate data and models. See \code{\link{set_ecov}} for full details. 
#'
#' \code{selectivity} specifies options for selectivity, to overwrite existing options specified in the ASAP data file. See \code{\link{set_selectivity}} for full details. 
#'
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file. See \code{\link{set_M}} for full details. 
#'
#' \code{NAA_re} specifies options for random effects on numbers-at-age (NAA, i.e. state-space model or not). See \code{\link{set_NAA}} for full details. 
#' If \code{NULL}, a traditional statistical catch-at-age model is fit.
#'
#' \code{catchability} specifies options for catchability. See \code{\link{set_q}} for full details. If \code{NULL} and \code{asap3} is not NULL, a single catchability parameter for each index is used with initial values specified in ASAP file. If both are NULL, initial catchabilities for all indices = 0.3.
#'
#' \code{move} specifies options for movement if there are multiple regions. See \code{\link{set_move}} for full details. 
#'
#' \code{age_comp} specifies the age composition models for fleet(s) and indices. See \code{\link{set_age_comp}} for full details. 
#'
#' \code{catch_info} is an optional list of fishery catch information that can be used to set up these types of observations when there is no asap3 file given.  See \code{\link{set_catch}} for full details. Useful for setting
#' up an operating model to simulate population processes and observations. Also can be useful for setting up the structure of assessment model when asap3 has not been used.
#'
#' \code{index_info} is an optional list of survey/index information that can be used to set up these types of observations when there is no asap3 file given.  See \code{\link{set_indices}} for full details. Useful for setting
#' up an operating model to simulate population processes and observations. Also can be useful for setting up the structure of assessment model when asap3 has not been used.
#'
#' \code{basic_info} is an optional list of information that can be used to set up the population and types of observations when there is no asap3 file given. Particularly useful for setting
#' up an operating model to simulate population processes and observations. Also can be useful for setting up the structure of assessment model when asap3 has not been used.
#' The current options are:
#'   \describe{
#'     \item{$ages}{integer vector of ages (years) with the last being a plus group}
#'     \item{$years}{integer vector of years that the population model spans.}
#'     \item{$n_seasons}{number of seasons within year.}
#'     \item{$n_fleets}{number of fishing fleets.}
#'     \item{$fracyr_seasons}{proportions of year for each season within year (sums to 1).}
#'     \item{$F}{matrix (length(years) x n_fleets) of annual fishing mortality rates for each fleet to initialize the model.}
#'     \item{$waa}{array ((n_fleets + n_indices + n_stocks) x length(years) x length(ages)) of annual weight at at age for each fleet, each index, and spawning biomass for each stock.}
#'     \item{$maturity}{array (n_stocks x length(years) x length(ages)) of annual maturity at age for estimating spawning biomass for each stock.}
#'     \item{$fracyr_SSB}{matrix (n_years x n_stocks) (1 or length(years)) of yearly proportions (0-1) of the year at which to calculate spawning biomass.}
#'     \item{$spawn_seasons}{vector (n_stocks) of seasons in which each stock spawns.}
#'     \item{$spawn_regions}{vector (n_stocks) of regions in which each stock spawns.}
#'     \item{$NAA_where}{array (n_stocks x n_regions x n_ages) of 0/1 indicating where individuals of each stock may exist on January 1 of each year.}
#'     \item{$Fbar_ages}{integer vector of ages to use to average total F at age defining fully selected F for reference points. May not be clearly known until a model is fitted.}
#'     \item{$q}{vector (length(n_indices)) of catchabilities for each of the indices to initialize the model.}
#'     \item{$percentSPR}{(0-100) percentage of unfished spawning biomass per recruit for determining equilibrium fishing mortality reference point}
#'     \item{$percentFXSPR}{(0-100) percentage of SPR-based F to use in projections.}
#'     \item{$percentFMSY}{(0-100) percentage of Fmsy to use in projections.}
#'		 \item{$XSPR_input_average_years}{which years to average inputs to per recruit calculation (selectivity, M, WAA, maturity) for SPR-based reference points. Default is last 5 years (tail(1:length(years),5))}
#'     \item{$XSPR_R_avg_yrs}{which years to average recruitments for calculating SPR-based SSB reference points. Default is 1:length(years)}
#'     \item{$XSPR_R_opt}{1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). 5: use bias-corrected expected recruitment. For long-term projections, may be important to use certain years for XSPR_R_avg_yrs}
#'     \item{$simulate_process_error}{T/F vector (length = 9). When simulating from the model, whether to simulate any process errors for 
#'     (NAA, M, selectivity, q, movement, unidentified mortality, q priors, movement priors, Ecov). Only used for applicable random effects.}
#'     \item{$simulate_observation_error}{T/F vector (length = 3). When simulating from the model, whether to simulate  catch, index, and ecov observations.}
#'     \item{$simulate_period}{T/F vector (length = 2). When simulating from the model, whether to simulate base period (model years) and projection period.}
#'     \item{$bias_correct_process}{T/F. Perform bias correction of log-normal random effects for NAA.}
#'     \item{$bias_correct_BRPs}{T/F. Perform bias correction of analytic SSB/R and Y/R when there is bias correction of log-normal NAA. May want to use XSPR_R_opt = 5 for long-term projections.}
#'   }
#' If other arguments to \code{prepare_wham_input} are provided such as \code{selectivity}, \code{M}, and \code{age_comp}, the information provided there
#' must be consistent with \code{basic_info}. For example the dimensions for number of years, ages, fleets, and indices.
#'
#' @param asap3 (optional) list containing data and parameters (output from \code{\link{read_asap3_dat}})
#' @param recruit_model numeric, option to specify stock-recruit model (see details)
#' @param model_name character, name of stock/model
#' @param ecov (optional) named list of environmental covariate data and parameters (see \code{\link{set_ecov}} for full details)
#' @param selectivity (optional) list specifying selectivity options by block: models, initial values, parameters to fix, and random effects (see \code{\link{set_selectivity}} for full details)
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see \code{\link{set_M}} for full details)
#' @param NAA_re (optional) list specifying options for random effect on numbers-at-age, initial numbers at age, and recruitment model (see \code{\link{set_NAA}} for full details)
#' @param catchability (optional) list specifying options for priors and random effects on catchability (see \code{\link{set_q}} for full details)
#' @param age_comp (optional) character or named list, specifies age composition model for fleet(s) and indices (see \code{\link{set_age_comp}} for full details)
#' @param move (optional) list specifying movement/migration options for models with more than 1 region (see \code{\link{set_move}} for full details)
#' @param L (optional) list specifying "extra" mortality options (see \code{\link{set_L}} for full details)
#' @param F (optional) list specifying fishing mortality options (see \code{\link{set_F}} for full details)
#' @param catch_info (optional) list specifying catch information (see \code{\link{set_catch}} for full details)
#' @param index_info (optional) list specifying index informaiton (see \code{\link{set_indices}} for full details)
#' @param basic_info (optional) list of basic population information for use when asap3 is not provided (see details)
#'
#' @return a named list with the following components:
#'   \describe{
#'     \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{map}}{Named list defining how to optionally collect and fix parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{years}}{Numeric vector of years to fit WHAM model (specified in ASAP3 .dat file)}
#'     \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group (specified in ASAP3 .dat file)}
#'     \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
#'     \item{\code{log}}{list of character strings attempting to describe the input and what the model assumes.}
#'   }
#'
#' @seealso \code{\link{read_asap3_dat}}, \code{\link{fit_wham}}, \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}, \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ex1_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' mod = fit_wham(input)
#'
#' # no ASAP3 file, default parameter values and configuration
#' input = prepare_wham_input()
#' mod = fit_wham(input, fit = FALSE)
#' set.seed(8675309)
#' simdata = mod$simulate(complete=TRUE)
#' input$data = simdata
#' fit = fit_wham(input, do.osa = FALSE)
#' }
#'
#' @export
prepare_wham_input <- function(asap3 = NULL, model_name="WHAM for unnamed stock", recruit_model=2, ecov=NULL, selectivity=NULL, M=NULL, 
	NAA_re=NULL, catchability=NULL, age_comp=NULL, move=NULL, L=NULL, F=NULL, catch_info=NULL, index_info=NULL, basic_info = NULL){

	data = list()
	par = list()
	map = list()
	random = character()
	log <- list()
	input = list(
	  	data = data,
	  	par = par,
	  	map = map,
	  	random = random,
	  	years = NULL, years_full = NULL, ages.lab = NULL, model_name = model_name, asap3 = asap3, log = log)


	# input = list()
	# input$data = list()
	# input$par = list()
	# input$map = list()

	if(is.null(basic_info)) basic_info = list(recruit_model = recruit_model)
	else basic_info$recruit_model = recruit_model


	#things that cannot be known from asap files
	input$data$n_seasons = 1
	input$data$fracyr_seasons = 1
	if(!is.null(basic_info$fracyr_seasons)){
		input$data$n_seasons = length(basic_info$fracyr_seasons)
		input$data$fracyr_seasons = basic_info$fracyr_seasons
	}


	waa_opts <- NULL
	waa_names <- which(names(basic_info) %in% c("waa", "waa_pointer_M", "waa_pointer_ssb"))
	waa_opts <- basic_info[waa_names]
	
	#catch_info = catch
	#catch_names = c("n_fleets","agg_catch", "catch_paa", "catch_cv","catch_Neff", "use_catch_paa", "selblock_pointer_fleets")
	#if(any(names(basic_info) %in% catch_names)) catch_opts = basic_info[catch_names]

	#index_opts = indices
	#index_names = c("n_indices", "agg_indices", "index_paa", "fracyr_indices", "index_cv", "index_Neff", "units_indices",
	#	"units_index_paa", "use_indices", "use_index_paa", "selblock_pointer_indices")
	#if(any(names(basic_info) %in% index_names)) index_opts = basic_info[index_names]

	F_opts = F
	#F_names = c("F")
	#if(any(names(basic_info) %in% F_names)) F_opts = basic_info[F_names]
	#print("1")
	input$log$misc <- list("\n NOTE: WHAM version 1.2.0 makes major changes to the structure of some data, parameters, and reported objects. \n") 
	input$log$misc <- c(input$log$misc, 
	"NOTE: WHAM version 1.2.0 decouples random effects for recruitment and random effects for older ages by default.
	To obtain results from previous versions set NAA_re$decouple_recruitment = FALSE. \n") 

	if(!is.null(asap3)) {
		input$log$asap3 <- list()
		asap3_test = sapply(asap3, function(x) "dat" %in% names(x))
		if(!all(asap3_test)) stop("object passed to asap3 argument does not have the correct structure.\n")

		asap3 = lapply(asap3, function(x) return(x$dat))

		if(length(asap3)> 1) input$log$asap3 <- c(input$log$asap3, paste0(length(asap3), " asap3 dat files were processed. One stock per region without mixing will be assumed \n
				unless the movement argument is provided.\n"))

		if(length(asap3) == 1) input$log$asap3 <- c(input$log$asap3, paste0(length(asap3), " asap3 dat file was processed. A single stock and region will be assumed. \n"))
#print(2)
  	n_ages = sapply(asap3, function(x) x$n_ages)
  	if(length(unique(n_ages))!= 1) stop("differing numbers of age classes in the asap3 dat files. Make them equal before passing to wham.")

  	n_years = sapply(asap3, function(x) x$n_years)
  	if(length(unique(n_years))!= 1) stop("differing numbers of years in the asap3 dat files. Make them equal before passing to wham.")

  	input$asap3 = asap3
		input$data$n_stocks = length(asap3)
		input$data$n_regions = length(asap3)
		input$data$n_ages = n_ages[1]
  	input$data$n_years_model = n_years[1]
		input$data$years_use <- 1:input$data$n_years_model - 1
  	input$years <- asap3[[1]]$year1 + 1:asap3[[1]]$n_years - 1

		input$data$spawn_seasons <- rep(1, length(asap3))
		input$data$spawn_regions <- 1:length(asap3)
		input$data$NAA_where <- array(0, dim = c(input$data$n_stocks,input$data$n_regions,input$data$n_ages))
		for(s in 1:input$data$n_stocks){
			input$data$NAA_where[s,input$data$spawn_regions[s],] <- 1 #recruit only to spawn region on Jan 1
		}
		if(input$data$n_regions > 1 & is.null(basic_info$NAA_where)){
			input$log$misc <- c(input$log$misc, "\n There is more than 1 region and basic_info$NAA_where is not specified so assuming each stock only exists in spawning region in first year.\n")
			#input$data$NAA_where[s,-input$data$spawn_regions[s],] <- 1 #recruit only to spawn region on Jan 1
		}
		if(!is.null(basic_info$NAA_where)){
			input$data$NAA_where[] <- basic_info$NAA_where
		}
  	input$data$fracyr_SSB <- matrix(NA, input$data$n_years_model, input$data$n_stocks)
		for(i in 1:length(asap3)){
			if(input$data$n_seasons == 1){
				input$data$fracyr_SSB[,i] <- asap3[[i]]$fracyr_spawn
			} else {
				int_starts <- cumsum(c(0,input$data$fracyr_seasons))
				ind <- max(which(int_starts <= asap3[[i]]$fracyr_spawn))
				input$data$spawn_seasons[i] <- ind
				input$data$fracyr_SSB[,i] <- asap3[[i]]$fracyr_spawn - int_starts[ind] 
			}
		}
#print(4)
  
		input$data$mature = array(NA, dim = c(input$data$n_stocks, input$data$n_years_model, input$data$n_ages))
	  for(i in 1:length(asap3)) input$data$mature[i,,] = asap3[[i]]$maturity
	  input$data$Fbar_ages = seq(asap3[[1]]$Frep_ages[1], asap3[[1]]$Frep_ages[2])
	}
	else {
		#if no asap3 is provided, make some default values to
		input = initial_input_no_asap_fn(input, basic_info)
	}

  input$years_full = input$years

  input$options <- list()
	 print("start")
	#some basic input elements see the function code below
	input = set_basic_info(input, basic_info)
	print("basic_info")

	# Catch
	#input$data$n_seasons first defined here
	input = set_catch(input, catch_info)
	print("catch")

	# Indices/surveys
	input = set_indices(input, index_info)
	print("indices")
	#print(input$data$selblock_pointer_indices)

	# WAA in case we want to modify how weight-at age is handled
	input = set_WAA(input, waa_opts)
	print("WAA")

	# NAA and recruitment options
	input = set_NAA(input, NAA_re)
	print("NAA")

	q_opts = catchability
	if(any(names(basic_info) == "q") & !any(names(q_opts) == "initial_q")) q_opts$initial_q = basic_info$q

	input = set_q(input, q_opts)
	print("q")

	# Selectivity
	input = set_selectivity(input, selectivity)
	print("selectivity")
	#print(input$data$selblock_pointer_indices)

	# Age composition model
	input = set_age_comp(input, age_comp)
	print("age_comp")

	#in case we want to add alternative F options
	input = set_F(input, F_opts)
	print("F")

	#set up natural mortality
	input = set_M(input, M)
	print("M")

	#set up movement
	input = set_move(input, move)
	print("move")

	#set up "extra" mortality
	input = set_L(input, L)
	print("L")


	#set up ecov data and parameters. Probably want to make sure to do this after set_NAA.
	input = set_ecov(input, ecov)
	print("ecov")

	# add vector of all observations for one step ahead residuals ==========================
	input = set_osa_obs(input)
	print("osa_obs")

	# projection data will always be modified by 'prepare_projection'
	input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
	#print("proj")

	#set any parameters as random effects
	input = set_random(input)
	#print("random")
	cat(unlist(input$log, recursive=T))

	input$call <- match.call()

	return(input)
}


#s may be 2 for most rhos on cpp side which is unusual.

gen.logit <- function(x, low, upp, s=1) (log((x-low)/(upp-x)))/s

is_internal_call <- function(n_gen = 2, NameSpace = "wham"){
  te <- topenv(parent.frame(n_gen)) #2 because the call will be made inside of a package function
	return(isNamespace(te) && getNamespaceName(te) == NameSpace)	
}


set_basic_info <- function(input, basic_info){
	#this function adds basic_info to input
	input$ages.lab = paste0(1:input$data$n_ages, c(rep("",input$data$n_ages-1),"+"))
	if(!is.null(basic_info$ages)) {
		if(!is.integer(basic_info$ages) | length(basic_info$ages) != input$data$n_ages) stop("basic_info$ages has been specified, but it is not an integer vector or it is not = n_ages")
		else {
  		input$ages.lab = paste0(basic_info$ages, c(rep("",input$data$n_ages-1),"+"))
		}
	}
	input$stock_names <- paste0("stock_", 1:input$data$n_stocks)
	if(!is.null(basic_info$stock_names)) {
		if(length(basic_info$stock_names) == input$data$n_stocks) input$stock_names <- as.character(basic_info$stock_names)
	}
	input$region_names <- paste0("region_", 1:input$data$n_regions)
	if(!is.null(basic_info$region_names)) {
		if(length(basic_info$region_names) == input$data$n_regions) input$region_names <- as.character(basic_info$region_names)
	}

  input$data$n_years_model = length(input$years)
	input$data$years_use <- 1:input$data$n_years_model - 1
  #input$data$n_years_catch = length(input$years)
  #input$data$n_years_indices = length(input$years)
	input$data$recruit_model = rep(2,input$data$n_stocks)
  input$data$recruit_model[] = basic_info$recruit_model #this is made from argument of the same name to prepare_wham_input
	if(is.null(basic_info$bias_correct_process) | is.null(basic_info$bias_correct_observation)){
		input$log$misc <- c(input$log$misc, 
	"NOTE: WHAM version 1.2.0 forward by default does not bias correct any log-normal process or observation errors. To 
	configure these, set basic_info$bias_correct_process = TRUE and/or basic_info$bias_correct_observation = TRUE. \n")
	}
  input$data$bias_correct_pe = 0 #bias correct log-normal process errors?
  input$data$bias_correct_oe = 0 #bias correct log-normal observation errors?
  input$data$bias_correct_brps = 0 #bias correct SSB/R and Y/R when NAA re are bias-corrected?
  if(!is.null(basic_info$bias_correct_process)) input$data$bias_correct_pe = as.integer(basic_info$bias_correct_process)
  if(!is.null(basic_info$bias_correct_observation)) input$data$bias_correct_oe = as.integer(basic_info$bias_correct_observation)
  if(!is.null(basic_info$bias_correct_BRPs)) input$data$bias_correct_brps = as.integer(basic_info$bias_correct_BRPs)
  
	#(NAA, M, selectivity, q, movement, unidentified mortality, q priors, movement priors, Ecov). Only used for applicable random effects.
  sim_pe = rep(1,9)
  if(!is.null(basic_info$simulate_process_error)) sim_pe[] = as.integer(basic_info$simulate_process_error)
  sim_oe = rep(1,3)
  if(!is.null(basic_info$simulate_observation_error)) sim_oe[] = as.integer(basic_info$simulate_observation_error)

  input$data$do_simulate_N_re = sim_pe[1] #simulate state variable
  input$data$do_simulate_M_re = sim_pe[2] #simulate state variable
  input$data$do_simulate_sel_re = sim_pe[3] #simulate state variable
  input$data$do_simulate_q_re = sim_pe[4] #simulate state variable
  input$data$do_simulate_mu_re = sim_pe[5] #simulate state variable
  input$data$do_simulate_L_re = sim_pe[6] #simulate state variable
  input$data$do_simulate_q_prior_re = sim_pe[7] #simulate state variable
  input$data$do_simulate_mu_prior_re = sim_pe[8] #simulate state variable
  input$data$do_simulate_Ecov_re = sim_pe[9] #simulate state variable
  input$data$do_simulate_data =  sim_oe #simulate data types (catch, indices, Ecov)
  input$data$do_simulate_period = c(1,1) #simulate processes and/or observations in model, projection periods
	input$data$do_post_samp_N = 0 #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
	input$data$do_post_samp_M = 0 #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
	input$data$do_post_samp_mu = 0 #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
	input$data$do_post_samp_q = 0 #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
	input$data$do_post_samp_sel = 0 #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
	input$data$do_post_samp_Ecov = 0 #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated

  #input$data$simulate_period = rep(1,2) #simulate above items for (model years, projection years)
	input$data$do_SPR_BRPs = 0 #this will be changed when after model fit
	input$data$do_MSY_BRPs = 0 #this will be changed when after model fit
	input$data$SPR_weight_type = 0
	input$data$SPR_weights = rep(1/input$data$n_stocks, input$data$n_stocks)
	input$data$n_regions_is_small = 1
	input$data$use_alt_AR1 = 0

  input$data$percentSPR = 40 #percentage of unfished SSB/R to use for SPR-based reference points
  input$data$percentFXSPR = 100 # percent of F_XSPR to use for calculating catch in projections
  input$data$percentFMSY = 100 # percent of F_XSPR to use for calculating catch in projections
  # data$XSPR_R_opt = 3 #1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). See next line for years to average over.
  input$data$XSPR_R_opt = 2 # default = use average R estimates
  input$data$XSPR_R_avg_yrs = 1:input$data$n_years_model-1 #model year indices to use for averaging recruitment when defining SSB_XSPR (if XSPR_R_opt = 2,4)
	input$data$FXSPR_static_init = 0.5 #initial value for Newton search of static F (spr-based) reference point (inputs to spr are averages of annual values using avg_years_ind)
	input$data$FMSY_static_init = 0.5 #initial value for Newton search of static F (spr-based) reference point (inputs to spr are averages of annual values using avg_years_ind)
	input$data$avg_years_ind <- tail(1:input$data$n_years_model,5) - 1 #default values to average for BRPs and projections
  input$data$which_F_age = rep(input$data$n_ages,input$data$n_years_model) #plus group by default used to define full F (total) IN annual reference points for projections, only. prepare_projection changes it to properly define selectivity for projections.
  	#rep(1,input$data$n_years_model))
  input$data$which_F_age_static = input$data$n_ages #plus group, fleet 1 by default used to define full F (total) for static SPR-based ref points.

  #if(!is.null(basic_info$simulate_period)) input$data$simulate_period = basic_info$simulate_period

  if(!is.null(basic_info$percentSPR)) input$data$percentSPR = basic_info$percentSPR
  if(!is.null(basic_info$percentFXSPR)) input$data$percentFXSPR = basic_info$percentFXSPR
  if(!is.null(basic_info$percentFMSY)) input$data$percentFMSY = basic_info$percentFMSY
  if(!is.null(basic_info$XSPR_R_opt)) input$data$XSPR_R_opt = basic_info$XSPR_R_opt
	if(!is.null(basic_info$XSPR_input_average_years)) input$data$avg_years_ind = basic_info$XSPR_input_average_years - 1 #user input shifted to start @ 0  
  if(!is.null(basic_info$XSPR_R_avg_yrs)) input$data$XSPR_R_avg_yrs = basic_info$XSPR_R_avg_yrs - 1 #user input shifted to start @ 0
  if(input$data$XSPR_R_opt<5) {
  	weighting <- "corresponding annual"
  	weighting <- ifelse(input$data$XSPR_R_opt %in% c(2,4), "corresponding annual", "average of annual")
		annual_R_type <- ifelse(input$data$XSPR_R_opt==1, "estimated", "conditionally expected")
		input$log$misc <- c(input$log$misc, paste0("For annual SPR-based reference points, ", weighting, " ", annual_R_type, " recruitments are used. \n"))
	} else{
		input$log$misc <- c(input$log$misc, paste0("For annual SPR-based reference points, bias-corrected marginally (over time) expected recruitment is used. \n"))
	}
	# if(input$data$XSPR_R_opt %in% c(1,3)){
	# 	annual_R_type <- ifelse(input$data$XSPR_R_opt==1, "estimated", "conditionally expected")
	# 	input$log$misc <- c(input$log$misc, paste0("For annual SPR-based reference points, coresponding annual ", annual_R_type, " recruitments are used. \n"))
	# }
	# if(input$data$XSPR_R_opt %in% c(2,4)){
	# 	annual_R_type <- ifelse(input$data$XSPR_R_opt==1, "estimated", "conditionally expected")
	# 	input$log$misc <- c(input$log$misc, paste0("For annual SPR-based reference points, average of annual ", annual_R_type, " recruitments over years, ", 
	# 	paste0(input$years[input$data$XSPR_R_avg_yrs+1], collapse = ","), " are used. \n"))
	# }

	input$options$basic_info <- basic_info
  return(input)

}
