#' Prepare input data and parameters for WHAM model
#'
#' After the data file is read into R by \code{\link{read_asap3_dat}}, this function
#' prepares the data and parameter settings for \code{\link{fit_wham}}.
#' By default, this will set up a SCAA version like \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP3}.
#'
#' \code{recruit_model} specifies the stock-recruit model. See lines 442-476
#' in \code{wham.cpp} to see implementation.
#'   \describe{
#'     \item{= 1}{Random walk, i.e. predicted recruitment in year i = recruitment in year i-1}
#'     \item{= 2}{(default) Random about mean, i.e. steepness = 1}
#'     \item{= 3}{Beverton-Holt}
#'     \item{= 4}{Ricker}
#'   }
#'
#' \code{Ecov$how} specifies HOW the environmental covariate affects the \code{Ecov$where} process.
#" Options for recruitment are described in \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}:
#'   \describe{
#'     \item{= 1}{"controlling" (dens-indep mortality)}
#'     \item{= 2}{"limiting" (carrying capacity, e.g. Ecov determines amount of suitable habitat)}
#'     \item{= 3}{"lethal" (threshold, i.e. R --> 0 at some Ecov value)}
#'     \item{= 4}{"masking" (metabolic/growth, decreases dR/dS)}
#'     \item{= 5}{"directive" (e.g. behavioral)}
#'   }
#'
#' \code{Ecov} specifies any environmental covariate data and model. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$label}{Name(s) of the environmental covariate(s). Used in printing.}
#'     \item{$mean}{Mean observations (vector if 1D, matrix if multivariate). Missing values = NA.}
#'     \item{$sigma}{Observation standard error (vector if 1D, covariance matrix if multivariate). Same length as \code{mean}.}
#'     \item{$year}{Years corresponding to observations (vector of same length as \code{mean} and \code{sigma})}
#'     \item{$use_obs}{T/F (or 0/1) vector/matrix of the same dimension as \code{mean} and \code{sigma}.
#'     Use the observation? Can be used to ignore subsets of the environmental covariate without changing data files.}
#'     \item{$lag}{Offset between the environmental covariate observations and their affect on the stock.
#'     I.e. if SST in year \emph{t} affects recruitment in year \emph{t + 1}, set \code{lag = 1}.}
#'     \item{$process_model}{Process model for the environmental covariate. "rw" = random walk, "ar1" = 1st order autoregressive, NA = do not fit}
#'     \item{$where}{Where does the environmental covariate affect the population? "recuit" = recruitment,
#'     "growth" = growth, "mortality" = natural mortality.}
#'     \item{$how}{How does the environmental covariate affect the \code{where} process? These options are
#'     specific to the \code{where} process.}
#'   }
#'
#' @param asap3 list containing data and parameters (output from \code{\link{read_asap3_dat}})
#' @param recruit_model numeric, option to specify stock-recruit model (see details)
#' @param model_name character, name of stock/model
#' @param Ecov (optional) named list of environmental covariate data and parameters (see details)
#'
#' @return a named list with the following components:
#'   \describe{
#'     \item{\code{data}}{Named list of data, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{par}}{Named list of parameters, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{map}}{not sure what this does, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{random}}{Character vector of parameters to treat as random effects, passed to \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}}
#'     \item{\code{years}}{Numeric vector of years to fit WHAM model (specified in ASAP3 .dat file)}
#'     \item{\code{ages.lab}}{Character vector of age labels, ending with plus-group (specified in ASAP3 .dat file)}
#'     \item{\code{model_name}}{Character, name of stock/model (specified in call to \code{prepare_wham_input})}
#'   }
#'
#' @seealso \code{\link{read_asap3_dat}}, \code{\link{fit_wham}}, \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP3}, \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ASAP_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' mod = fit_wham(input)
#' }
#'
#' @export
prepare_wham_input <- function(asap3, recruit_model=2, model_name="WHAM for unnamed stock", Ecov=NULL){
  asap3 = asap3$dat
  which_indices <- which(asap3$use_index ==1)
  asap3$n_indices = length(which_indices)
  asap3$survey_index_units <- asap3$index_units[which_indices]
  asap3$survey_acomp_units <- asap3$index_acomp_units[which_indices]
  asap3$survey_WAA_pointers <- asap3$index_WAA_pointers[which_indices]
  asap3$survey_month <- matrix(asap3$index_month[which_indices], asap3$n_years, asap3$n_indices, byrow = TRUE)
  asap3$use_survey_acomp <- asap3$use_index_acomp[which_indices]
  asap3$index_sel_option <- asap3$index_sel_option[which_indices]
  asap3$index_sel_ini = asap3$index_sel_ini[which_indices]
  asap3$index_WAA_pointers = asap3$index_WAA_pointers[which_indices]
  asap3$use_catch_acomp <- rep(1,asap3$n_fleets) #default is to use age comp for catch
  asap3$IAA_mats <- asap3$IAA_mats[which_indices]
  asap3$use_survey <- asap3$use_index[which_indices]
  data = list(n_years_model = asap3$n_years)
  data$n_years_catch = asap3$n_years
  data$n_years_indices = asap3$n_years
  data$n_ages = asap3$n_ages
  data$n_fleets = asap3$n_fleets
  data$n_indices <- asap3$n_indices
  data$n_selblocks <- asap3$n_fleet_sel_blocks + asap3$n_indices
  data$selblock_models <- c(asap3$sel_block_option, asap3$index_sel_option)
  #data$selblock_types <- c(asap3$sel_block_option, asap3$index_sel_option)
  n_selpars <- sum(sapply(data$selblock_models, function(x) if(x==1) data$n_ages else if(x==2) 2 else if(x==3) 4))
  data$selblock_pointer_fleets = cbind(sapply(asap3$sel_block_assign, function(x) return(x)))
  data$selblock_pointer_indices = matrix(rep(asap3$n_fleet_sel_blocks + 1:data$n_indices, each = data$n_years_model), data$n_years_model, data$n_indices)
  data$age_comp_model_fleets = rep(1, data$n_fleets) #multinomial by default
  data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[data$age_comp_model_fleets]
  data$age_comp_model_indices = rep(1, data$n_indices) #multinomial by default
  data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[data$age_comp_model_indices]

  data$fracyr_SSB = rep(asap3$fracyr_spawn, data$n_years_model)
  data$mature = asap3$maturity
  i <- c(seq(1,(asap3$n_fleets+1)*2-1,2),(asap3$n_fleets+1)*2 + 1:2)
  WAA_pointers <- asap3$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
  data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
  data$waa_pointer_totcatch = asap3$WAA_pointers[data$n_fleets + 1]
  data$waa_pointer_indices = asap3$index_WAA_pointers
  data$waa_pointer_ssb = asap3$WAA_pointers[data$n_fleets + 2]
  data$waa_pointer_jan1 = asap3$WAA_pointers[data$n_fleets + 3]

  data$waa = array(NA, dim = c(length(asap3$WAA_mats), data$n_years_model, data$n_ages))
  for(i in 1:length(asap3$WAA_mats)) data$waa[i,,] = asap3$WAA_mats[[i]]

  data$agg_catch = matrix(NA, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets) data$agg_catch[,i] = asap3$CAA_mats[[i]][,data$n_ages + 1]
  data$agg_catch_sigma = asap3$catch_cv
  data$agg_catch_sigma[which(data$agg_catch_sigma < 1e-15)] = 100
  data$agg_catch_sigma = sqrt(log(data$agg_catch_sigma^2 + 1))
  data$catch_paa = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  for(i in 1:data$n_fleets)
  {
    temp = asap3$CAA_mats[[i]][,1:data$n_ages]
    temp[which(is.na(temp))] = 0
    temp[which(temp<0)] = 0
    data$catch_paa[i,,] = temp/apply(temp,1,sum)
  }
  data$catch_paa[is.na(data$catch_paa)] = 0
  data$use_catch_paa = matrix(1, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets)
  {
    if(asap3$use_catch_acomp[i] != 1) data$use_catch_paa[,i] = 0
    else for(y in 1:data$n_years_model) if(asap3$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2) data$use_catch_paa[y,i] = 0
  }

  data$catch_Neff = asap3$catch_Neff
  data$catch_aref = matrix(NA, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets) data$catch_aref[,i] = get_aref_fn(data$catch_paa[i,,])

  data$units_indices <- asap3$survey_index_units
  data$fracyr_indices = (asap3$survey_month-1)/12 #make sure that this is right
  data$agg_indices = matrix(NA, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices) data$agg_indices[,i] = asap3$IAA_mats[[i]][,2]
  data$use_indices = matrix(1, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices)
  {
    for(y in 1:data$n_years_model) if(asap3$IAA_mats[[i]][y,2] < 1e-15) data$use_indices[y,i] = 0
  }
  data$agg_index_sigma = matrix(NA, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices) data$agg_index_sigma[,i] = asap3$IAA_mats[[i]][,3]
  data$agg_index_sigma[which(data$agg_index_sigma < 1e-15)] = 100
  data$agg_index_sigma = sqrt(log(data$agg_index_sigma^2 + 1))

  data$units_index_paa <- asap3$survey_acomp_units

  data$index_paa = array(NA, dim = c(data$n_indices, data$n_years_model, data$n_ages))
  for(i in 1:data$n_indices)
  {
    temp = asap3$IAA_mats[[i]][,3 + 1:data$n_ages]
    temp[which(is.na(temp))] = 0
    temp[which(temp<0)] = 0
    data$index_paa[i,,] = temp/apply(temp,1,sum)
    #data$index_paa[i,,][sign(asap3$IAA_mats[[i]][,3 + 1:data$n_ages]) == -1] = 0
  }
  data$index_paa[is.na(data$index_paa)] = 0
  data$use_index_paa = matrix(1, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices)
  {
    if(asap3$use_survey_acomp[i] != 1) data$use_index_paa[,i] = 0
    else for(y in 1:data$n_years_model) if(asap3$IAA_mats[[i]][y,4 + data$n_ages] < 1e-15 | sum(data$index_paa[i,y,] > 1e-15)<2) data$use_index_paa[y,i] = 0
  }
  # print(data$use_index_paa[14,1])
  data$index_Neff = matrix(NA, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices) data$index_Neff[,i] = asap3$IAA_mats[[i]][,4 + data$n_ages]
  data$index_aref = matrix(NA, data$n_years_model, data$n_indices)
  for(i in 1:data$n_indices) data$index_aref[,i] = get_aref_fn(data$index_paa[i,,])
  data$q_lower <- rep(0,data$n_indices)
  data$q_upper <- rep(1000,data$n_indices)

  selpars_ini = matrix(NA, data$n_selblocks, data$n_ages + 6)
  for(i in 1:asap3$n_fleet_sel_blocks) selpars_ini[i,] = asap3$sel_ini[[i]][,1]
  for(i in (1:asap3$n_indices)) selpars_ini[i+asap3$n_fleet_sel_blocks,] = asap3$index_sel_ini[[i]][,1]
  phase_selpars = matrix(NA, data$n_selblocks, data$n_ages + 6)
  for(i in 1:asap3$n_fleet_sel_blocks) phase_selpars[i,] = asap3$sel_ini[[i]][,2]
  for(i in (1:asap3$n_indices)) phase_selpars[i+asap3$n_fleet_sel_blocks,] = asap3$index_sel_ini[[i]][,2]
  for(i in 1:data$n_selblocks)
  {
    if(data$selblock_model[i] == 1) phase_selpars[i,data$n_ages + 1:6] = -1
    if(data$selblock_model[i] %in% c(2,3)) phase_selpars[i,c(1:data$n_ages, data$n_ages + 3:6)] = -1
    if(data$selblock_model[i] == 4) phase_selpars[i,data$n_ages + 1:2] = -1
  }
  # print(selpars_ini)
  selpars_lo = selpars_hi = matrix(0, data$n_selblocks, data$n_ages + 6)
  selpars_hi[,1:data$n_ages] = 1
  selpars_hi[,data$n_ages + 1:6] = data$n_ages
  temp = matrix(NA, data$n_selblocks, data$n_ages + 6)
  temp[which(phase_selpars>0)] = 1:sum(phase_selpars>0)
  map = list(logit_selpars = factor(temp))

  data$selpars_lower = selpars_lo #only need these for estimated parameters
  data$selpars_upper = selpars_hi

  data$n_NAA_sigma <- 2 #by default
  data$NAA_sigma_pointers <- c(1,rep(2,data$n_ages-1))
  data$recruit_model = recruit_model
  data$n_M_re = data$n_ages
  data$MAA_pointer = 1:data$n_ages #rep(1,data$n_ages)
  data$M_sigma_par_pointer = rep(1,data$n_M_re)
  data$M_model = 0
  # data$recruit_model = 2 #random about mean
  data$N1_model = 0 #0: just age-specific numbers at age
  data$use_M_re = 0
  data$use_NAA_re = 0
  data$use_b_prior = 0
  data$random_recruitment = 0 #1 #make sure use_NAA_re = 0, recruitment is still a random effect.
  data$which_F_age = data$n_ages #plus group by default used to define full F for BRPs
  data$use_steepness = 0 #use regular SR parameterization by default, steepness still can be estimated as derived par.
  data$bias_correct_pe = 0 #bias correct log-normal process errors?
  data$bias_correct_oe = 0 #bias correct log-normal observation errors?
  data$Fbar_ages = 1:data$n_ages
  data$simulate_state = 1 #simulate any state variables
  data$percentSPR = 40 #percentage of unfished SSB/R to use for SPR-based reference points

  model_years <- asap3$year1 + 1:asap3$n_years - 1
  # add in environmental covariate data
  if(is.null(Ecov)){
    data$Ecov_obs <- matrix(1, nrow=1, ncol=1)
    data$Ecov_obs_sigma <- matrix(0, nrow=1, ncol=1)
    data$n_Ecov <- 1
    data$Ecov_use_obs <- matrix(0, nrow=1, ncol=1)
    data$Ecov_year <- matrix(0, nrow=1, ncol=1)
    data$year1_Ecov <- 0
    data$year1_model <- asap3$year1
    data$Ecov_lag <- 0
    data$Ecov_model <- 0
    data$Ecov_where <- 1
    data$Ecov_how <- 0
    data$Ecov_recruit <- 1
    data$Ecov_growth <- 1
    data$Ecov_mortality <- 1
    data$n_years_Ecov <- 1
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- 0
    data$Ecov_label <- "none"
  } else {
    if(class(Ecov$mean) == "matrix") {data$Ecov_obs <- Ecov$mean}
    else{
      warning("Ecov mean is not a matrix. Coercing to a matrix...")
      data$Ecov_obs <- as.matrix(Ecov$mean)
    }
    if(class(Ecov$sigma) == "matrix") {data$Ecov_obs_sigma <- Ecov$sigma}
    else{
      warning("Ecov sigma is not a matrix. Coercing to a matrix...")
      data$Ecov_obs_sigma <- as.matrix(Ecov$sigma)
    }
    if(!identical(dim(data$Ecov_obs_sigma), dim(data$Ecov_obs))) stop("Dimensions of Ecov mean != dimensions of Ecov sigma")
    if(class(Ecov$use_obs) == "matrix") {data$Ecov_use_obs <- Ecov$use_obs}
    else{
      warning("Ecov_use_obs is not a matrix with same dimensions as Ecov mean. Coercing to a matrix...")
      data$Ecov_use_obs <- as.matrix(as.integer(Ecov$use_obs))
    }
    if(!identical(dim(data$Ecov_use_obs), dim(data$Ecov_obs))) stop("Dimensions of Ecov_use_obs != dimensions of Ecov mean")
    if(length(Ecov$year) != dim(data$Ecov_obs)[1]) stop("Ecov year is not the same length as # rows in Ecov mean")
    data$Ecov_year <- as.numeric(Ecov$year)
    data$n_Ecov <- dim(data$Ecov_obs)[2] # num of covariates
    data$year1_Ecov <- Ecov$year[1]
    data$year1_model <- asap3$year1
    end_model <- tail(model_years,1)
    end_Ecov <- tail(Ecov$year,1)
    if(length(Ecov$label) == data$n_Ecov) {data$Ecov_label <- Ecov$label}
    else{
      warning("Number of Ecov labels not equal to number of Ecovs
              Setting Ecov labels = 'Ecov 1', 'Ecov 2', ...")
      data$Ecov_label = paste0("Ecov ",1:data$n_Ecov)
    }

    # check that Ecov year vector doesn't have missing gaps
    if(all(diff(model_years)!=1)) stop("Ecov years not continuous")

    # # pad Ecov if it starts after model year2 - max(lag)
    # if(data$year1_Ecov > data$year1_model - max(Ecov$lag) + 1){
    #   warning("Ecov does not start by model year 2 - max(lag). Padding Ecov...")
    #   data$Ecov_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(Ecov$lag)+1), ncol = data$n_Ecov), data$Ecov_obs)
    #   data$Ecov_obs_sigma <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(Ecov$lag)+1), ncol = data$n_Ecov), data$Ecov_obs_sigma)
    #   data$Ecov_use_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(Ecov$lag)+1), ncol = data$n_Ecov), data$Ecov_use_obs)
    #   data$Ecov_year <- c(seq(data$year1_model - max(Ecov$lag)+1, data$year1_Ecov-1), data$Ecov_year)
    #   data$year1_Ecov <- data$year1_model - max(Ecov$lag)+1
    # }
    # pad Ecov if it starts after model year1 - max(lag)
    if(data$year1_Ecov > data$year1_model - max(Ecov$lag)){
      warning("Ecov does not start by model year 1 - max(lag). Padding Ecov...")
      data$Ecov_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(Ecov$lag)), ncol = data$n_Ecov), data$Ecov_obs)
      data$Ecov_obs_sigma <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(Ecov$lag)), ncol = data$n_Ecov), data$Ecov_obs_sigma)
      data$Ecov_use_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(Ecov$lag)), ncol = data$n_Ecov), data$Ecov_use_obs)
      data$Ecov_year <- c(seq(data$year1_model - max(Ecov$lag), data$year1_Ecov-1), data$Ecov_year)
      data$year1_Ecov <- data$year1_model - max(Ecov$lag)
    }

    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      warning("Ecov last year is before model last year. Padding Ecov...")
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_obs_sigma <- rbind(data$Ecov_obs_sigma, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_year <- c(data$Ecov_year, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- dim(data$Ecov_obs)[1] # num years Ecov to model (padded)

    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- rep(NA, data$n_Ecov)
    for(i in 1:data$n_Ecov){
      data$ind_Ecov_out_start[i] <- which(data$Ecov_year==data$year1_model)-Ecov$lag[i]-1 # -1 is for cpp indexing
      data$ind_Ecov_out_end[i] <- which(data$Ecov_year==end_model)-Ecov$lag[i]-1 # -1 is for cpp indexing
    }

    if(!identical(length(Ecov$lag), length(Ecov$label), data$n_Ecov)) stop("Length of Ecov_lag and Ecov_label vectors not equal to # Ecov")
    data$Ecov_lag <- Ecov$lag
    if(!all(Ecov$process_model %in% c(NA,"rw", "ar1"))){
      stop("Ecov$process_model must be 'rw' (random walk), 'ar1', or NA (do not fit)")}
    if(is.na(Ecov$process_model) && Ecov$how !=0){
      stop("Ecov$process_model not chosen (NA) but Ecov$how specified.
       Either 1) choose an Ecov process model ('rw' or 'ar1'),
              2) turn off Ecov (set Ecov$how = 0 and Ecov$process_model = NA),
           or 3) fit Ecov but with no effect on population (Ecov$how = 0, Ecov$process_model 'rw' or 'ar1').")}
    data$Ecov_model <- sapply(Ecov$process_model, match, c("rw", "ar1"))

  #  data$n_Ecov_pars <- c(1,3)[data$Ecov_model] # rw: 1 par (sig), ar1: 3 par (phi, sig)
    if(!Ecov$where %in% c('recruit')){
      stop("Sorry, only Ecov effects on recruitment currently implemented.
      Set Ecov$where = 'recruit'.")}
    data$Ecov_where <- sapply(Ecov$where, match, c('recruit','growth','mortality'))
    data$Ecov_recruit <- ifelse(any(data$Ecov_where == 1), which(data$Ecov_where == 1), 0) # Ecov index to use for recruitment?
    data$Ecov_growth <- ifelse(any(data$Ecov_where == 2), which(data$Ecov_where == 2), 0) # Ecov index to use for growth?
    data$Ecov_mortality <- ifelse(any(data$Ecov_where == 3), which(data$Ecov_where == 3), 0) # Ecov index to use for mortality?

    if(!Ecov$how %in% c(0,1,2,4)){
      stop("Sorry, only Ecov effects on recruitment currently implemented.
      Set Ecov$how = 0 (no effect), 1 (controlling), 2 (limiting, Bev-Holt only), or 4 (masking).")}
    if(recruit_model == 1 & data$Ecov_recruit != 0){
      stop("Random walk recruitment cannot have an Ecov effect on recruitment.
      Either choose a different recruit_model (2, 3, or 4), or remove the Ecov effect.")
    }
    if(recruit_model == 4 & Ecov$how[data$Ecov_recruit] == 2){
      stop("'Limiting' Ecov effect on Ricker recruitment not implemented.
      Either set Ecov$how = 0 (no effect), 1 (controlling), or 4 (masking)...
      Or set recruit_model = 3 (Bev-Holt).")
    }
    data$Ecov_how <- Ecov$how
    # if(Ecov$where=="recruit") data$Ecov_how <- match(Ecov$how, c('type1','type2','type3'))
    # if(Ecov$where=='growth') data$Ecov_how <- match(Ecov$how, c('type1','type2','type3'))
    # if(Ecov$where=='mortality') data$Ecov_how <- match(Ecov$how, c('type1','type2','type3'))

    cat(paste0("Please check that the environmental covariates have been loaded
and interpreted correctly.

Model years: ", data$year1_model, " to ", end_model,"
Ecov years: ", data$year1_Ecov, " to ", end_Ecov,"

"))
    for(i in 1:data$n_Ecov){
      years <- data$Ecov_year[as.logical(data$Ecov_use_obs[,i])]
      cat(paste0("Ecov ",i,": ",Ecov$label[i],"
",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how+1]," effect on: ", c('recruitment','growth','mortality')[data$Ecov_where[i]],"

In model years:
"))
cat(years, fill=TRUE)
lastyr <- tail(years,1)
cat(paste0("Lag: ",data$Ecov_lag[i],"
Ex: ",Ecov$label[i]," in ",years[1]," affects ", c('recruitment','growth','mortality')[data$Ecov_where[i]]," in ",years[1+data$Ecov_lag[i]],"
    ",Ecov$label[i]," in ",lastyr," affects ", c('recruitment','growth','mortality')[data$Ecov_where[i]]," in ",lastyr+data$Ecov_lag[i],"
"))
      }
    } # end load Ecov

  # add vector of all observations for one step ahead residuals ==========================
  # 4 components: fleet catch (log), index catch (log), paa catch, paa index
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 1. log fleet catch
  x <- as.data.frame(data$agg_catch)
  colnames(x) <- paste0("fleet_", 1:data$n_fleets)
  x$year <- 1:data$n_years_catch
  tmp <- tidyr::gather(x, fleet, val, -year)
  tmp$val <- log(tmp$val) # shouldn't be any years with 0 catch... could make this robust later
  tmp$age <- NA
  tmp$type <- "logcatch"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 2. log index catch
  x <- as.data.frame(data$agg_indices)
  x[data$use_indices==0] <- NA # only include index data to fit in obsvec
  colnames(x) <- paste0("index_", 1:data$n_indices)
  x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
  tmp <- tidyr::gather(x, fleet, val, -year)
  tmp <- tmp[complete.cases(tmp),]
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_indices==0, turned to NA, and already removed
  tmp$age <- NA
  tmp$type <- "logindex"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 3. Ecov
  if(!all(data$Ecov_use_obs==0)){
    x <- as.data.frame(data$Ecov_obs)
    x[data$Ecov_use_obs==0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("Ecov_", 1:data$n_Ecov)
    # x$year <- 1:data$n_years_Ecov # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
    x$year <- seq(from=data$year1_Ecov-data$year1_model+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::gather(x, fleet, val, -year)
    tmp <- tmp[complete.cases(tmp),]
    tmp$age <- NA
    tmp$type <- "Ecov"
    obs <- rbind(obs, tmp[, obs.colnames])
}

  # # 4. paa catch
  # dimnames(data$catch_paa) <- list(fleet=paste0("fleet_", 1:data$n_fleets),
  #                                  year=1:data$n_years_catch,
  #                                  age=1:data$n_ages)
  # x <- as.data.frame(dplyr::as.tbl_cube(data$catch_paa, met_name = "val"))
  # x$type <- "paacatch"
  # obs <- rbind(obs, x[, obs.colnames])

  # # 5. paa index
  # dimnames(data$index_paa) <- list(fleet=paste0("index_", 1:data$n_indices),
  #                                  year=1:data$n_years_indices,
  #                                  age=1:data$n_ages)
  # x <- as.data.frame(dplyr::as.tbl_cube(data$index_paa, met_name = "val"))
  # x$type <- "paaindex"
  # obs <- rbind(obs, x[, obs.colnames])

  # order by year, fleet, age, type
  o <- order(as.numeric(obs$year), obs$fleet, as.numeric(obs$age), obs$type)
  obs <- obs[o,]

  # calculate obsvec indices in keep arrays
  obs$ind <- 1:dim(obs)[1]
  data$keep_C <- matrix(subset(obs, type=='logcatch')$ind, nrow=data$n_years_catch, ncol=data$n_fleets, byrow=TRUE)
 
  data$keep_I <- matrix(NA, nrow=data$n_years_indices, ncol=data$n_indices)
  # data$keep_I[data$use_indices==1] <- subset(obs, type=='logindex')$ind
  # xl <- apply(data$use_indices,1,function(r) which(r==1))
  xl <- lapply(seq_len(nrow(data$use_indices)), function(r) which(data$use_indices[r,]==1))
  Col <- unlist(xl)
  Row <- rep(1:data$n_years_indices, times=sapply(xl, length))
  data$keep_I[cbind(Row,Col)] <- subset(obs, type=='logindex')$ind
 
  data$keep_E <- matrix(NA, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
  # data$keep_E[data$Ecov_use_obs==1] <- subset(obs, type=='Ecov')$ind
  xl <- lapply(seq_len(nrow(data$Ecov_use_obs)), function(r) which(data$Ecov_use_obs[r,]==1))
  # xl <- apply(data$Ecov_use_obs,1,function(r) which(r==1))
  Col <- unlist(xl)
  Row <- rep(1:data$n_years_Ecov, times=sapply(xl, length))
  data$keep_E[cbind(Row,Col)] <- subset(obs, type=='Ecov')$ind  

  data$keep_Cpaa <- array(NA, dim=c(data$n_fleets, data$n_years_catch, data$n_ages))
  for(i in 1:data$n_fleets) data$keep_Cpaa[i,,] <- matrix(subset(obs, type=='paacatch' & fleet==paste0("fleet_",i))$ind, nrow=data$n_years_catch, ncol=data$n_ages, byrow=TRUE)
  data$keep_Ipaa <- array(NA, dim=c(data$n_indices, data$n_years_indices, data$n_ages))
  for(i in 1:data$n_indices) data$keep_Ipaa[i,,] <- matrix(subset(obs, type=='paaindex' & fleet==paste0("index_",i))$ind, nrow=data$n_years_indices, ncol=data$n_ages, byrow=TRUE)
  # subtract 1 bc TMB indexes from 0
  data$keep_C <- data$keep_C - 1
  data$keep_I <- data$keep_I - 1
  data$keep_E <- data$keep_E - 1
  data$keep_Cpaa <- data$keep_Cpaa - 1
  data$keep_Ipaa <- data$keep_Ipaa - 1

  data$obs <- obs
  data$obsvec <- obs$val

  # data$obsvec[data$keep_I[data$use_indices==1]+1] - log(data$agg_indices[data$use_indices==1])
  # data$obsvec[data$keep_E[data$Ecov_use_obs==1]+1] - data$Ecov_obs[data$Ecov_use_obs==1]

  par = list(mean_rec_pars = numeric(c(0,1,2,2)[recruit_model]))
  if(recruit_model==2) par$mean_rec_pars = 10
  if(recruit_model==4) par$mean_rec_pars[2] = -10
  par$logit_q = rep(-8, data$n_indices)
  par$log_F1 = rep(-2, data$n_fleets)
  par$F_devs = matrix(0, data$n_years_model-1, data$n_fleets)
  if(data$N1_model == 1) par$log_N1_pars = c(10,log(0.1))
  if(data$N1_model == 0) par$log_N1_pars = rep(10,data$n_ages)
  par$log_NAA_sigma = rep(0, data$n_NAA_sigma)
  par$logit_selpars = log(selpars_ini-selpars_lo) - log(selpars_hi - selpars_ini)
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] = -10
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] = 10
  # print(par$logit_selpars)
  n_catch_acomp_pars = c(0,1,1,3,1,2)[data$age_comp_model_fleets[which(apply(data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[data$age_comp_model_indices[which(apply(data$use_index_paa,2,sum)>0)]]
  par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  par$index_paa_pars = rep(0, sum(n_index_acomp_pars))
  par$log_NAA = matrix(10, data$n_years_model-1, data$n_ages)
  par$M_pars1 = log(asap3$M[1,]) #log(mean(asap3$M))
  par$M_sigma_pars = 0
  par$M_re = matrix(log(asap3$M[-1,]), data$n_years_model-1, data$n_M_re)
  par$log_b = log(0.305)
  par$log_R = rep(10, data$n_years_model-1) #/n_years_model-1, if used.
  par$log_R_sigma = 0
  par$log_catch_sig_scale = rep(0, data$n_fleets)
  par$log_index_sig_scale = rep(0, data$n_indices)

  # add environmental covariate parameters
  par$Ecov_re = matrix(0, data$n_years_Ecov, data$n_Ecov)
  par$Ecov_beta = rep(0, data$n_Ecov) # one for each Ecov, beta_R in eqns 4-5, Miller et al. (2016)
  par$Ecov_process_pars = matrix(0, 3, data$n_Ecov) # nrows = RW: 2 par (log_sig, Ecov1), AR1: 3 par (mu, phi, log_sig); ncol = N_ecov

  # turn off 3rd Ecov par if it's a RW
  tmp.pars <- par$Ecov_process_pars
  for(i in 1:data$n_Ecov) tmp.pars[3,i] <- ifelse(data$Ecov_model[i]==1, NA, 0)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)

  # turn off only Ecov_beta to fit Ecov process model without effect on population
  tmp <- par$Ecov_beta
  for(i in 1:data$n_Ecov) tmp[i] <- ifelse(data$Ecov_how[i]==0, NA, 0)
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_beta = factor(tmp)

  # turn off Ecov pars if no Ecov (re, process)
  # for any Ecov_model = NA, Ecov$how must be 0 and beta is already turned off
  data$Ecov_model[is.na(data$Ecov_model)] = 0 # turn any NA into 0
  tmp.re <- matrix(1:length(par$Ecov_re), data$n_years_Ecov, data$n_Ecov, byrow=FALSE)
  for(i in 1:data$n_Ecov){
    tmp.pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else tmp.pars[,i]
    tmp.re[,i] <- if(data$Ecov_model[i]==0) rep(NA,data$n_years_Ecov) else tmp.re[,i]
    if(data$Ecov_model[i]==1) tmp.re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
  }
  ind.notNA <- which(!is.na(tmp.re))
  tmp.re[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_re = factor(tmp.re)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_process_pars = factor(tmp.pars)

  map$log_catch_sig_scale = factor(rep(NA, data$n_fleets))
  map$log_index_sig_scale = factor(rep(NA, data$n_indices))
  map$M_pars1 = factor(rep(NA, length(par$M_pars1)))
  map$M_re = factor(rep(NA, length(par$M_re)))
  map$M_sigma_pars = factor(rep(NA, length(par$M_sigma_pars)))
  map$log_NAA = factor(rep(NA, length(par$log_NAA)))
  map$log_NAA_sigma = factor(rep(NA, length(par$log_NAA_sigma)))
  map$mean_rec_pars = factor(rep(NA, length(par$mean_rec_pars)))
  map$log_R_sigma = factor(rep(NA, length(par$log_R_sigma)))
  map$log_b = factor(rep(NA,length(par$log_b)))
  random = character()
  if(missing(model_name)) model_name = "WHAM for unnamed stock"
  return(list(data=data, par = par, map = map, random = random, years = model_years,
    ages.lab = paste0(1:data$n_ages, c(rep("",data$n_ages-1),"+")), model_name = model_name))
}

get_aref_fn = function(paa){
  n_years = NROW(paa)
  n_ages = NCOL(paa)
  aref = rep(-1, n_years)
  for(y in 1:n_years)
  {
    temp = paa[y,]
    for(a in 1:n_ages) if(temp[a] < 1.0e-15) temp[a] = 0.0
    if (sum(temp > 1.0e-15)>1)
    { #both requirements as well as total catch > 0 to include age comp in objective function
      paa[y,]=temp/sum(temp)
      for(a in 1:n_ages) if(paa[y,a] > 1.0e-15) aref[y] = a
      for(a in n_ages:1)
      {
        if(paa[y,a] > 1.0e-15)
        {
          aref[y] = a #last positive
          break
        }
      }
      #this part is necessary for logistic-normal age comp (type = 5)
      #note the aref can be associated with an observed 0 if needed.
      this_aref = aref[y] - 1 #start one less than last positive
      for(a in this_aref:1)
      {
        if(length(paa[y,a])== 0)
        {
          print(y)
          print(a)
          print(paa[y,])
          print(this_aref)
        }
        if(paa[y,a] > 1.0e-15) break #next to last is positive, don't change aref
        else aref[y] = a #move aref down one.
      }
    }
  }
  return(aref)
}
