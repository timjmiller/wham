#' Prepare input data and parameters for WHAM model
#'
#' After the data file is read into R by \code{\link{read_asap3_dat}}, this function
#' prepares the data and parameter settings for \code{\link{fit_wham}}.
#' By default, this will set up a SCAA version like \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}.
#'
#' \code{recruit_model} specifies the stock-recruit model. See \code{wham.cpp} to see implementation.
#'   \describe{
#'     \item{= 1}{Random walk, i.e. predicted recruitment in year i = recruitment in year i-1}
#'     \item{= 2}{(default) Random about mean, i.e. steepness = 1}
#'     \item{= 3}{Beverton-Holt}
#'     \item{= 4}{Ricker}
#'   }
#'
#' \code{ecov} specifies any environmental covariate data and model. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$label}{Name(s) of the environmental covariate(s). Used in printing.}
#'     \item{$mean}{Mean observations (matrix). Missing values = NA.}
#'     \item{$logsigma}{Observation standard error (log). Options:
#'       \describe{
#'         \item{Matrix with same dimensions as \code{$mean}}{Specified values (not estimated) for each time step }
#'         \item{Single value per ecov, numeric vector or matrix w/ dim 1 x n.ecov}{Specified value (not estimated) shared among time steps}
#'         \item{Vector with estimate options for each ecov, length = n.ecov}{
#'           \code{'est_1'}: Estimated, one value shared among time steps.
#'           \code{'est_re'}: Estimated value for each time step as random effects with two parameters (mean, var)}
#'       }
#'     }
#'     \item{$year}{Years corresponding to observations (vector of same length as \code{$mean} and \code{$logsigma})}
#'     \item{$use_obs}{T/F (or 0/1) vector/matrix of the same dimension as \code{$mean} and \code{$logsigma}.
#'     Use the observation? Can be used to ignore subsets of the ecov without changing data files.}
#'     \item{$lag}{Offset between the ecov observations and their affect on the stock.
#'     I.e. if SST in year \emph{t} affects recruitment in year \emph{t + 1}, set \code{lag = 1}.}
#'     \item{$process_model}{Process model for the ecov time-series. \code{"rw"} = random walk, \code{"ar1"} = 1st order autoregressive, \code{NA} = do not fit}
#'     \item{$where}{Where does the ecov affect the population? \code{"recuit"} = recruitment,
#'     \code{"M"} = natural mortality, \code{"growth"} = growth.}
#'     \item{$how}{How does the ecov affect the \code{$where} process? These options are
#'     specific to the \code{$where} process.}
#'     \item{$link_model}{Model describing ecov effect on the \code{$where} process. Options: 'linear' (default) or 'poly-x'
#'     where x = 2, ... (e.g. 'poly-2' specifies a quadratic model, \eqn{b0 + b1*ecov + b2*ecov^2 + ...}).}
#'   }
#'
#' \code{ecov$how} specifies HOW the ecov affects the \code{ecov$where} process.
#" Options for recruitment are described in \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}:
#'   \describe{
#'     \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'     \item{= 1}{"controlling" (dens-indep mortality)}
#'     \item{= 2}{"limiting" (carrying capacity, e.g. ecov determines amount of suitable habitat)}
#'     \item{= 3}{"lethal" (threshold, i.e. R --> 0 at some ecov value)}
#'     \item{= 4}{"masking" (metabolic/growth, decreases dR/dS)}
#'     \item{= 5}{"directive" (e.g. behavioral)}
#'   }
#' Options for M:
#'   \describe{
#'     \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'     \item{= 1}{effect on mean M (shared across ages)}
#'   }
#'
#' \code{selectivity} specifies options for selectivity, to overwrite existing options specified in the ASAP data file.
#' If \code{NULL}, selectivity options from the ASAP data file are used. \code{selectivity} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Selectivity model for each block. Vector with length = number of selectivity blocks. Each entry must be one of: "age-specific", "logistic", "double-logistic", or "decreasing-logistic".}
#'     \item{$re}{Time-varying (random effects) for each block. Vector with length = number of selectivity blocks.
#'                  If \code{NULL}, selectivity parameters in all blocks are constant over time and uncorrelated.
#'                  Each entry of \code{selectivity$re} must be one of the following options, where the selectivity parameters are:
#'                  \describe{
#'                    \item{"none"}{(default) are constant and uncorrelated}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by age/par (AR1), but not year}
#'                    \item{"ar1_y"}{correlated by year (AR1), but not age/par}
#'                    \item{"2dar1"}{correlated by year and age/par (2D AR1)}
#'                  }
#'                 }
#'     \item{$initial_pars}{Initial parameter values for each block. List of length = number of selectivity blocks. Each entry must be a vector of length # parameters in the block, i.e. \code{c(2,0.2)} for logistic or \code{c(0.5,0.5,0.5,1,1,0.5)} for age-specific with 6 ages.}
#'     \item{$fix_pars}{Which parameters to fix at initial values. List of length = number of selectivity blocks. E.g. model with 3 age-specific blocks and 6 ages, \code{list(c(4,5),4,c(2,3,4))} will fix ages 4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in block 3.}
#'   }
#'
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
#' If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
#' shared/mirrored among some but not all ages, modify \code{input$map$M_a} after calling \code{prepare_wham_input}
#' (see vignette for more details). \code{M} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Natural mortality model options are:
#'                    \describe{
#'                      \item{"constant"}{(default) estimate a single M shared across all ages}
#'                      \item{"age-specific"}{estimate M_a independent for each age}
#'                      \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
#'                        \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
#'                        \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.}
#'                    }
#'                  }
#'     \item{$re}{Time- and age-varying (random effects) on M. Note that random effects can only be estimated if
#'                mean M-at-age parameters are (\code{$est_ages} is not \code{NULL}).
#'                  \describe{
#'                    \item{"none"}{(default) M constant in time and uncorrelated across ages.}
#'                    \item{"iid"}{M varies by year and age, but uncorrelated.}
#'                    \item{"ar1_a"}{M correlated by age (AR1), but not year.}
#'                    \item{"ar1_y"}{M correlated by year (AR1), but not age.}
#'                    \item{"rw_y"}{Random walk on M by year, ages share variance but uncorrelated.}
#'                    \item{"2dar1"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan (2016)}.}
#'                  }
#'                 }
#'     \item{$initial_means}{Initial/mean M-at-age vector, with length = number of ages (if \code{$model = "age-specific"})
#'                          or 1 (if \code{$model = "weight-at-age" or "constant"}). If \code{NULL}, initial mean M-at-age values are taken
#'                          from the first row of the MAA matrix in the ASAP data file.}
#'     \item{$est_ages}{Vector of ages to estimate M (others will be fixed at initial values). E.g. in a model with 6 ages,
#'                      \code{$est_ages = 5:6} will estimate M for the 5th and 6th ages, and fix M for ages 1-4. If \code{NULL},
#'                      M at all ages is fixed at \code{M$initial_means} (if not \code{NULL}) or row 1 of the MAA matrix from the ASAP file (if \code{M$initial_means = NULL}).}
#'     \item{$logb_prior}{(Only if \code{$model = "weight-at-age"}) TRUE or FALSE (default), should a N(0.305, 0.08) prior be
#'                        used on log_b? Based on Fig. 1 and Table 1 (marine fish) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}.}
#'   }
#'
#' @param asap3 list containing data and parameters (output from \code{\link{read_asap3_dat}})
#' @param recruit_model numeric, option to specify stock-recruit model (see details)
#' @param model_name character, name of stock/model
#' @param ecov (optional) named list of environmental covariate data and parameters (see details)
#' @param selectivity (optional) list specifying selectivity options by block: models, initial values, parameters to fix, and random effects (see details)
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)
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
#'   }
#'
#' @seealso \code{\link{read_asap3_dat}}, \code{\link{fit_wham}}, \href{https://www.nefsc.noaa.gov/nft/ASAP.html}{ASAP}, \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ASAP_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' mod = fit_wham(input)
#' }
#'
#' @export
prepare_wham_input <- function(asap3, model_name="WHAM for unnamed stock", recruit_model=2, ecov=NULL, selectivity=NULL, M=NULL){
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

  # Selectivity
  data$n_selblocks <- asap3$n_fleet_sel_blocks + asap3$n_indices
  if(is.null(selectivity$model)) data$selblock_models <- c(asap3$sel_block_option, asap3$index_sel_option)
  if(!is.null(selectivity$model)){
    if(length(selectivity$model) != data$n_selblocks) stop("Length of selectivity$model must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$model %in% c("age-specific","logistic","double-logistic","decreasing-logistic"))) stop("Each model entry must be one of the following: 'age-specific','logistic','double-logistic','decreasing-logistic'")
    data$selblock_models <- match(selectivity$model, c("age-specific","logistic","double-logistic","decreasing-logistic"))
  }
  if(is.null(selectivity$re)) data$selblock_models_re <- rep(1, data$n_selblocks) # default: no RE on selectivity parameters
  if(!is.null(selectivity$re)){
    if(length(selectivity$re) != data$n_selblocks) stop("Length of selectivity$re must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$re %in% c("none","iid","ar1","ar1_y","2dar1"))) stop("Each selectivity$re entry must be one of the following: 'none','iid','ar1','ar1_y','2dar1'")
    data$selblock_models_re <- match(selectivity$re, c("none","iid","ar1","ar1_y","2dar1"))
  }
  data$n_selpars <- c(data$n_ages,2,4,2)[data$selblock_models] # num selpars per block
  data$selblock_pointer_fleets = cbind(sapply(asap3$sel_block_assign, function(x) return(x)))
  data$selblock_pointer_indices = matrix(rep(asap3$n_fleet_sel_blocks + 1:data$n_indices, each = data$n_years_model), data$n_years_model, data$n_indices)
  selblock_pointers <- cbind(data$selblock_pointer_fleets, data$selblock_pointer_indices)
  data$selblock_years <- matrix(0, nrow=data$n_years_model, ncol=data$n_selblocks)
  for(b in 1:data$n_selblocks) data$selblock_years[,b] <- apply(selblock_pointers, 1, function(x) b %in% x)
  data$n_years_selblocks <- apply(data$selblock_years, 2, sum)

  # Age composition model
  data$age_comp_model_fleets = rep(1, data$n_fleets) #multinomial by default
  data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[data$age_comp_model_fleets]
  data$age_comp_model_indices = rep(1, data$n_indices) #multinomial by default
  data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[data$age_comp_model_indices]

  data$fracyr_SSB = rep(asap3$fracyr_spawn, data$n_years_model)
  data$mature = asap3$maturity

  # Weight-at-age
  i <- c(seq(1,(asap3$n_fleets+1)*2-1,2),(asap3$n_fleets+1)*2 + 1:2)
  WAA_pointers <- asap3$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
  data$waa_pointer_fleets = WAA_pointers[1:data$n_fleets]
  data$waa_pointer_totcatch = asap3$WAA_pointers[data$n_fleets + 1]
  data$waa_pointer_indices = asap3$index_WAA_pointers
  data$waa_pointer_ssb = asap3$WAA_pointers[data$n_fleets + 2]
  data$waa_pointer_jan1 = asap3$WAA_pointers[data$n_fleets + 3]
  data$waa = array(NA, dim = c(length(asap3$WAA_mats), data$n_years_model, data$n_ages))
  for(i in 1:length(asap3$WAA_mats)) data$waa[i,,] = asap3$WAA_mats[[i]]

  # Catch
  data$agg_catch = matrix(NA, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets) data$agg_catch[,i] = asap3$CAA_mats[[i]][,data$n_ages + 1]
  data$agg_catch_sigma = asap3$catch_cv
  data$agg_catch_sigma[which(data$agg_catch_sigma < 1e-15)] = 100
  data$agg_catch_sigma = sqrt(log(data$agg_catch_sigma^2 + 1))
  data$catch_paa = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  data$use_agg_catch = matrix(1, data$n_years_model, data$n_fleets)
  for(i in 1:data$n_fleets)
  {
    temp = asap3$CAA_mats[[i]][,1:data$n_ages]
    temp[which(is.na(temp))] = 0
    temp[which(temp<0)] = 0
    data$catch_paa[i,,] = temp/apply(temp,1,sum)
    for(y in 1:data$n_years_model) if(asap3$CAA_mats[[i]][y,data$n_ages+1] < 1e-15) data$use_agg_catch[y,i] = 0
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

  # Indices/surveys
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

  # Prep selectivity initial values
  selpars_ini = matrix(NA, data$n_selblocks, data$n_ages + 6)
  if(is.null(selectivity$initial_pars)) {
    for(i in 1:asap3$n_fleet_sel_blocks) selpars_ini[i,] = asap3$sel_ini[[i]][,1]
    for(i in (1:asap3$n_indices)) selpars_ini[i+asap3$n_fleet_sel_blocks,] = asap3$index_sel_ini[[i]][,1]
  } else {
    if(!is.list(selectivity$initial_pars)) stop("selectivity$initial_pars must be a list")
    if(length(selectivity$initial_pars) != data$n_selblocks) stop("Length of selectivity$initial_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    for(b in 1:data$n_selblocks){
      if(length(selectivity$initial_pars[[b]]) != data$n_selpars[b]) stop(paste0("Length of vector ",b," in the selectivity$initial_pars list is not equal to the number of selectivity parameters for block ",b,": ",data$n_selpars[b]))
      if(data$selblock_models[b] == 1) selpars_ini[b,1:data$n_ages] = selectivity$initial_pars[[b]]
      if(data$selblock_models[b] %in% c(2,4)) selpars_ini[b,data$n_ages+1:2] = selectivity$initial_pars[[b]]
      if(data$selblock_models[b] == 3) selpars_ini[b,data$n_ages+3:6] = selectivity$initial_pars[[b]]
    }
  }

  # Prep selectivity map
  phase_selpars = matrix(1, data$n_selblocks, data$n_ages + 6)
  if(is.null(selectivity$fix_pars)){
    for(i in 1:asap3$n_fleet_sel_blocks) phase_selpars[i,] = asap3$sel_ini[[i]][,2]
    for(i in (1:asap3$n_indices)) phase_selpars[i+asap3$n_fleet_sel_blocks,] = asap3$index_sel_ini[[i]][,2]
  } else {
    if(!is.list(selectivity$fix_pars)) stop("selectivity$fix_pars must be a list")
    if(length(selectivity$fix_pars) != data$n_selblocks) stop("Length of selectivity$fix_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices).
      Use 'NULL' to not fix any parameters for a block, e.g. list(NULL,4,2) does not fix any pars in block 1")
    for(b in 1:data$n_selblocks){
      if(data$selblock_models[b] == 1) phase_selpars[b,selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] %in% c(2,4)) phase_selpars[b,data$n_ages+selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] == 3) selpars_ini[b,data$n_ages+2+selectivity$fix_pars[[b]]] = selectivity$initial_pars[[b]]
    }
  }
  for(i in 1:data$n_selblocks){
    if(data$selblock_models[i] == 1) phase_selpars[i,data$n_ages + 1:6] = -1
    if(data$selblock_models[i] %in% c(2,4)) phase_selpars[i,c(1:data$n_ages, data$n_ages + 3:6)] = -1
    if(data$selblock_models[i] == 3) phase_selpars[i,data$n_ages + 1:2] = -1
  }

  # For age-specific selectivity blocks, check for ages with ~zero catch and fix these at 0
  age_specific <- which(data$selblock_models==1)
  for(b in age_specific){
    ind = list(fleets = which(apply(data$selblock_pointer_fleets == b,2,sum) > 0))
    ind$indices = which(apply(data$selblock_pointer_indices == b,2,sum) > 0)
    paa = matrix(nrow = 0, ncol = data$n_ages)
    if(length(ind$fleets)) for(f in ind$fleets)
    {
      y = data$catch_paa[f,which(data$selblock_pointer_fleets[,f] == b & data$use_catch_paa[,f] == 1),]
      paa = rbind(paa,y)
    }
    if(length(ind$indices)) for(i in ind$indices)
    {
      y = data$index_paa[i,which(data$selblock_pointer_indices[,i] == b & data$use_index_paa[,i] == 1),]
      paa = rbind(paa,y)
    }
    y = apply(paa,2,sum)
    selpars_ini[b, which(y < 1e-5)] = 0
    phase_selpars[b, which(y < 1e-5)] = -1
  }
  data$selpars_est <- phase_selpars
  data$selpars_est[data$selpars_est == -1] = 0
  data$n_selpars_est <- apply(data$selpars_est > 0, 1, sum)
  selpars_lo = selpars_hi = matrix(0, data$n_selblocks, data$n_ages + 6)
  selpars_hi[,1:data$n_ages] = 1
  selpars_hi[,data$n_ages + 1:6] = data$n_ages
  temp = matrix(NA, data$n_selblocks, data$n_ages + 6)
  temp[which(phase_selpars>0)] = 1:sum(phase_selpars>0)
  data$selpars_lower = selpars_lo #only need these for estimated parameters
  data$selpars_upper = selpars_hi
  map = list(logit_selpars = factor(temp))

  data$n_NAA_sigma <- 2 #by default
  data$NAA_sigma_pointers <- c(1,rep(2,data$n_ages-1))
  data$recruit_model = recruit_model

  # natural mortality options, default = use values from ASAP file, no estimation
  data$n_M_a = data$n_ages
  data$M_model = 2
  data$use_b_prior = 0
  data$M_re_model = 1 # default = no RE / 'none'
  data$M_est <- rep(0, data$n_M_a) # default = don't estimate M
  M_first_est = 1
  M0_ini <- log(asap3$M[1,1])
  M_a_ini <- log(asap3$M[1,]) - M0_ini
  M_re_ini <- matrix(log(asap3$M[-1,])-matrix(M_a_ini + M0_ini,data$n_years_model-1,data$n_M_a,byrow=T), data$n_years_model-1, data$n_M_a)
  if(!is.null(M)){
    if(!is.null(M$model)){ # M model options
      if(!(M$model %in% c("constant","age-specific","weight-at-age"))) stop("M$model must be either 'constant', 'age-specific', or 'weight-at-age'")
      data$M_model <- match(M$model, c("constant","age-specific","weight-at-age"))
      if(M$model %in% c("constant","weight-at-age")){
        data$n_M_a = 1
        data$M_est = 0
        M_a_ini = 0
        M_re_ini = matrix(log(asap3$M[-1,1])-matrix(M0_ini,data$n_years_model-1,1,byrow=T), data$n_years_model-1, data$n_M_a)
        if(is.null(M$initial_means) & length(unique(asap$M[1,])) > 1) warning("Constant or weight-at-age M specified (so only 1 mean M parameter),
                                                                             but first row of MAA matrix has > 1 unique value.
                                                                             Initializing M at age-1 MAA values. To avoid this warning
                                                                             without changing ASAP file, specify M$initial_means.")
        if(!is.null(M$logb_prior)){
          if(!is.logical(M$logb_prior)) stop("M$logb_prior must be either TRUE or FALSE")
          if(M$logb_prior) data$use_b_prior = 1
        }
      }
    }
    if(!is.null(M$re)){
      if(length(M$re) != 1) stop("Length(M$re) must be 1")
      if(!(M$re %in% c("none","iid","ar1_a","ar1_y","rw_y","2dar1"))) stop("M$re must be one of the following: 'none','iid','ar1_a','ar1_y','rw_y','2dar1'")
      data$M_re_model <- match(M$re, c("none","iid","ar1_a","ar1_y","rw_y","2dar1"))
      # data$use_M_re <- 1
    }
    if(!is.null(M$initial_means)){
      if(length(M$initial_means) != data$n_M_a) stop("Length(M$initial_means) must be # ages (if age-specific M) or 1 (if weight-at-age M)")
      M0_ini <- log(M$initial_means[1])
      M_a_ini <- log(M$initial_means) - M0_ini
      # Use ASAP file M values
      # M_re_ini <- matrix(log(asap3$M[-1,])-matrix(M_a_ini,data$n_years_model-1,data$n_M_a,byrow=T), data$n_years_model-1, data$n_M_a)
      # overwrite ASAP file M values
      M_re_ini[] = 0
    }
    if(!is.null(M$est_ages)){
      if(!all(M$est_ages %in% 1:data$n_M_a)) stop("All M$est_ages must be in 1:n.ages (if age-specific M) or 1 (if constant or weight-at-age M)")
      data$M_est[M$est_ages] = 1 # turn on estimation for specified M-at-age
      M_first_est <- M$est_ages[1]
      if(is.null(M$initial_means)){
        M0_ini <- log(asap3$M[1,M_first_est])
        M_a_ini <- log(asap3$M[1,]) - M0_ini
      } else {
        M0_ini <- log(M$initial_means[M_first_est])
        M_a_ini <- log(M$initial_means) - M0_ini
      }
      M_re_ini[] = 0 # if estimating mean M for any ages, initialize yearly deviations at 0
    }
  }
  data$n_M_est <- sum(data$M_est)
  # if(data$n_M_est == 0 & data$M_re_model > 1) stop("At least one mean M must be estimated to estimate random effects on M")

  # data$recruit_model = 2 #random about mean
  data$N1_model = 0 #0: just age-specific numbers at age
  data$use_NAA_re = 0
  data$random_recruitment = 0 #1 #make sure use_NAA_re = 0, recruitment is still a random effect.
  data$which_F_age = data$n_ages #plus group by default used to define full F and F RP IN projections, only. prepare_projection changes it to properly define selectivity for projections.
  data$use_steepness = 0 #use regular SR parameterization by default, steepness still can be estimated as derived par.
  data$bias_correct_pe = 0 #bias correct log-normal process errors?
  data$bias_correct_oe = 0 #bias correct log-normal observation errors?
  data$Fbar_ages = 1:data$n_ages
  data$simulate_state = 1 #simulate any state variables
  data$percentSPR = 40 #percentage of unfished SSB/R to use for SPR-based reference points
  data$XSPR_R_opt = 3 #1(3): use annual R estimates(predictions) for annual SSB_XSPR, 2(4): use average R estimates(predictions). See next line for years to average over.
  data$XSPR_R_avg_yrs = 1:data$n_years_model #model year indices (TMB, starts @ 0) to use for averaging recruitment when defining SSB_XSPR (if XSPR_R_opt = 2,4)
  model_years <- asap3$year1 + 1:asap3$n_years - 1

  # --------------------------------------------------------------------------------
  # Environmental covariate data
  if(is.null(ecov)){
    data$Ecov_obs <- matrix(1, nrow=1, ncol=1)
    par.Ecov.obs.logsigma <- matrix(-2.3, nrow=1, ncol=1)
    map.Ecov.obs.logsigma <- factor(NA)
    par.Ecov.obs.sigma.par <- matrix(0, nrow=1, ncol=1)
    map.Ecov.obs.sigma.par <- factor(NA)
    data$Ecov_obs_sigma_opt <- 1
    data$n_Ecov <- 1
    data$Ecov_use_obs <- matrix(0, nrow=1, ncol=1)
    data$Ecov_year <- matrix(0, nrow=1, ncol=1)
    data$year1_Ecov <- 0
    data$year1_model <- asap3$year1
    data$Ecov_lag <- 0
    data$Ecov_model <- 0
    data$Ecov_where <- 1
    data$Ecov_how <- 0
    data$Ecov_poly <- 1
    data$Ecov_recruit <- 1
    data$Ecov_growth <- 1
    data$Ecov_mortality <- 1
    data$n_years_Ecov <- 1
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- 0
    data$Ecov_label <- "none"
    data$Ecov_use_re <- matrix(0, nrow=1, ncol=1)
  } else {
    if(class(ecov$mean) == "matrix") {data$Ecov_obs <- ecov$mean} else{
      warning("ecov mean is not a matrix. Coercing to a matrix...")
      data$Ecov_obs <- as.matrix(ecov$mean)
    }
    if(class(ecov$use_obs) == "matrix"){
      data$Ecov_use_obs <- ecov$use_obs
    } else{
      warning("ecov$use_obs is not a matrix with same dimensions as ecov$mean. Coercing to a matrix...")
      data$Ecov_use_obs <- as.matrix(as.integer(ecov$use_obs))
    }
    if(!identical(dim(data$Ecov_use_obs), dim(data$Ecov_obs))) stop("Dimensions of ecov$use_obs != dimensions of ecov$mean")

    # Handle Ecov sigma options
    data$n_Ecov <- dim(data$Ecov_obs)[2] # num Ecovs
    n_Ecov_obs <- dim(data$Ecov_obs)[1] # num Ecov obs
    if(class(ecov$logsigma) == "matrix"){
      # if(any(ecov$logsigma <= 0)) stop("ecov$logsigma must be positive")
      data$Ecov_obs_sigma_opt = 1
      par.Ecov.obs.logsigma <- ecov$logsigma
      if(!identical(dim(par.Ecov.obs.logsigma), dim(data$Ecov_obs))) stop("Dimensions of ecov$mean != dimensions of ecov$logsigma")
    }
    if(class(ecov$logsigma) == 'numeric'){
      # if(any(ecov$logsigma <= 0)) stop("ecov$logsigma must be positive")
      data$Ecov_obs_sigma_opt = 1
      print("ecov$logsigma is numeric. Coercing to a matrix...")
      # warning("Ecov sigma is not a matrix. Coercing to a matrix...")
      if(length(ecov$logsigma) == data$n_Ecov) par.Ecov.obs.logsigma <- matrix(rep(ecov$logsigma, each=n_Ecov_obs), ncol=data$n_Ecov)
      if(length(ecov$logsigma) == n_Ecov_obs && data$n_Ecov == 1) par.Ecov.obs.logsigma <- matrix(ecov$logsigma, ncol=1)
      if(length(ecov$logsigma) != data$n_Ecov && length(ecov$logsigma) != n_Ecov_obs) stop("ecov$logsigma is numeric but length is not equal to # of ecovs or ecov observations")
    }
    if(class(ecov$logsigma) == 'character'){
      if(ecov$logsigma == 'est_1'){ # estimate 1 Ecov obs sigma for each Ecov
        data$Ecov_obs_sigma_opt = 2
        par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov)
        map.Ecov.obs.logsigma <- matrix(rep(1:data$n_Ecov, each=n_Ecov_obs), ncol=data$n_Ecov)
        par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
        map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
      }
      if(ecov$logsigma == 'est_fe'){ # estimate Ecov obs sigma for each Ecov obs
        data$Ecov_obs_sigma_opt = 3
        par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov) # fixed effect inits
        map.Ecov.obs.logsigma <- matrix(1:(n_Ecov_obs*data$n_Ecov), nrow=n_Ecov_obs, ncol=data$n_Ecov) # est all
        par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
        map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
      }
      if(ecov$logsigma == 'est_re'){
        data$Ecov_obs_sigma_opt = 4
        par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov) # random effect inits
        map.Ecov.obs.logsigma <- matrix(1:(n_Ecov_obs*data$n_Ecov), nrow=n_Ecov_obs, ncol=data$n_Ecov) # est all
        par.Ecov.obs.sigma.par <- matrix(c(rep(-1.3, data$n_Ecov), rep(-0.5, data$n_Ecov)), ncol=data$n_Ecov, byrow=TRUE) # random effect pars
        map.Ecov.obs.sigma.par <- matrix(1:(2*data$n_Ecov), nrow=2, ncol=data$n_Ecov)
      }
    }
    if(data$Ecov_obs_sigma_opt == 1){ # Ecov sigma given, initialized at given values
      map.Ecov.obs.logsigma <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
      par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
      map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
    }

    if(length(ecov$year) != n_Ecov_obs) stop("ecov$year is not the same length as # rows in ecov$mean")
    data$Ecov_year <- as.numeric(ecov$year)
    data$year1_Ecov <- ecov$year[1]
    data$year1_model <- asap3$year1
    end_model <- tail(model_years,1)
    end_Ecov <- tail(ecov$year,1)
    if(length(ecov$label) == data$n_Ecov){
      data$Ecov_label <- ecov$label
    } else {
      warning("Number of Ecov labels not equal to number of Ecovs
              Setting Ecov labels = 'Ecov 1', 'Ecov 2', ...")
      data$Ecov_label = paste0("Ecov ",1:data$n_Ecov)
    }

    # # check that Ecov year vector doesn't have missing gaps
    # if(all(diff(model_years)!=1)) stop("Ecov years not continuous")

    # # pad Ecov if it starts after model year2 - max(lag)
    # if(data$year1_Ecov > data$year1_model - max(ecov$lag) + 1){
    #   warning("Ecov does not start by model year 2 - max(lag). Padding Ecov...")
    #   data$Ecov_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)+1), ncol = data$n_Ecov), data$Ecov_obs)
    #   data$Ecov_obs_sigma <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)+1), ncol = data$n_Ecov), data$Ecov_obs_sigma)
    #   data$Ecov_use_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)+1), ncol = data$n_Ecov), data$Ecov_use_obs)
    #   data$Ecov_year <- c(seq(data$year1_model - max(ecov$lag)+1, data$year1_Ecov-1), data$Ecov_year)
    #   data$year1_Ecov <- data$year1_model - max(ecov$lag)+1
    # }
    # pad Ecov if it starts after model year1 - max(lag)
    if(data$year1_Ecov > data$year1_model - max(ecov$lag)){
      print("ecov does not start by model year 1 - max(lag). Padding ecov...")
      # warning("Ecov does not start by model year 1 - max(lag). Padding Ecov...")
      data$Ecov_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), data$Ecov_obs)
      par.Ecov.obs.logsigma <- rbind(matrix(-1.3, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), par.Ecov.obs.logsigma)
      map.Ecov.obs.logsigma <- rbind(matrix(NA, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), map.Ecov.obs.logsigma)
      data$Ecov_use_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max(ecov$lag)), ncol = data$n_Ecov), data$Ecov_use_obs)
      data$Ecov_year <- c(seq(data$year1_model - max(ecov$lag), data$year1_Ecov-1), data$Ecov_year)
      data$year1_Ecov <- data$year1_model - max(ecov$lag)
    }

    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      print("Ecov last year is before model last year. Padding Ecov...")
      # warning("Ecov last year is before model last year. Padding Ecov...")
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      par.Ecov.obs.logsigma <- rbind(par.Ecov.obs.logsigma, matrix(-1.3, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      map.Ecov.obs.logsigma <- rbind(map.Ecov.obs.logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_year <- c(data$Ecov_year, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- dim(data$Ecov_obs)[1] # num years Ecov to model (padded)
    data$Ecov_use_re <- matrix(1, nrow=data$n_years_Ecov, ncol=data$n_Ecov)

    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- rep(NA, data$n_Ecov)
    for(i in 1:data$n_Ecov){
      data$ind_Ecov_out_start[i] <- which(data$Ecov_year==data$year1_model)-ecov$lag[i]-1 # -1 is for cpp indexing
      data$ind_Ecov_out_end[i] <- which(data$Ecov_year==end_model)-ecov$lag[i]-1 # -1 is for cpp indexing
    }

    if(!identical(length(ecov$lag), length(ecov$label), data$n_Ecov)) stop("Length of Ecov_lag and Ecov_label vectors not equal to # Ecov")
    data$Ecov_lag <- ecov$lag
    if(!all(ecov$process_model %in% c(NA,"rw", "ar1"))){
      stop("ecov$process_model must be 'rw' (random walk), 'ar1', or NA (do not fit)")
    }
    if(is.na(ecov$process_model) && ecov$how !=0){
      stop("ecov$process_model not chosen (NA) but ecov$how specified.
       Either 1) choose an ecov process model ('rw' or 'ar1'),
              2) turn off ecov (set ecov$how = 0 and ecov$process_model = NA),
           or 3) fit ecov but with no effect on population (ecov$how = 0, ecov$process_model 'rw' or 'ar1').")
    }
    data$Ecov_model <- sapply(ecov$process_model, match, c("rw", "ar1"))

  #  data$n_Ecov_pars <- c(1,3)[data$Ecov_model] # rw: 1 par (sig), ar1: 3 par (phi, sig)
    if(!ecov$where %in% c('recruit','M')){
      stop("Sorry, only ecov effects on recruitment and M currently implemented.
      Set ecov$where = 'recruit' or 'M'.")
    }
    data$Ecov_where <- sapply(ecov$where, match, c('recruit','M'))
    data$Ecov_recruit <- ifelse(any(data$Ecov_where == 1), which(data$Ecov_where == 1), 0) # Ecov index to use for recruitment?
    data$Ecov_growth <- ifelse(any(data$Ecov_where == 3), which(data$Ecov_where == 3), 0) # Ecov index to use for growth?
    data$Ecov_mortality <- ifelse(any(data$Ecov_where == 2), which(data$Ecov_where == 2), 0) # Ecov index to use for mortality?

    if(is.null(ecov$how)) stop("ecov$how must be specified")
    if(length(ecov$how) != data$n_Ecov) stop("ecov$how must be a vector of length(n.ecov)")
    if(data$Ecov_mortality > 0) if(!ecov$how[data$Ecov_mortality] %in% c(0,1)){
      stop("Sorry, only the following ecov effects on M are currently implemented.
      Set ecov$how = 0 (no effect) or 1 (effect on mean M, shared across ages).")
    }
    if(data$Ecov_recruit > 0) if(!ecov$how[data$Ecov_recruit] %in% c(0,1,2,4)){
      stop("Sorry, only the following ecov effects on recruitment are currently implemented.
      Set ecov$how = 0 (no effect), 1 (controlling), 2 (limiting, Bev-Holt only), or 4 (masking).")
    }
    if(recruit_model == 1 & data$Ecov_recruit != 0){
      stop("Random walk recruitment cannot have an ecov effect on recruitment.
      Either choose a different recruit_model (2, 3, or 4), or remove the Ecov effect.")
    }
    if(data$Ecov_recruit > 0) if(recruit_model == 4 & ecov$how[data$Ecov_recruit] == 2){
      stop("'Limiting' ecov effect on Ricker recruitment not implemented.
      Either set ecov$how = 0 (no effect), 1 (controlling), or 4 (masking)...
      Or set recruit_model = 3 (Bev-Holt).")
    }
    data$Ecov_how <- ecov$how
    data$Ecov_poly <- rep(1,data$n_Ecov)
    ecov_str <- as.list(rep('linear',data$n_Ecov))
    if(!is.null(ecov$link_model)){
      if(!is.na(ecov$link_model)){
        for(i in 1:data$n_Ecov){
          ecov_str[[i]] <- strsplit(ecov$link_model[i],"-")[[1]]
          if(!ecov_str[[i]][1] %in% c('linear','poly')) stop("Only 'linear' or 'poly-x' (x = 1, 2, ...) ecov link models implemented.")
          if(ecov_str[[i]]=='linear') data$Ecov_poly[i] <- 1
          if(ecov_str[[i]][1]=='poly') data$Ecov_poly[i] <- as.numeric(ecov_str[[i]][2])
          ecov_str[[i]] = ecov$link_model[i]
        }
      }
    }
    # if(ecov$where=="recruit") data$Ecov_how <- match(ecov$how, c('type1','type2','type3'))
    # if(ecov$where=='growth') data$Ecov_how <- match(ecov$how, c('type1','type2','type3'))
    # if(ecov$where=='mortality') data$Ecov_how <- match(ecov$how, c('type1','type2','type3'))

    cat(paste0("Please check that the environmental covariates have been loaded
and interpreted correctly.

Model years: ", data$year1_model, " to ", end_model,"
Ecov years: ", data$year1_Ecov, " to ", end_Ecov,"

"))
    for(i in 1:data$n_Ecov){
      years <- data$Ecov_year[as.logical(data$Ecov_use_obs[,i])]

      if(data$Ecov_where[i] == 1){ # recruitment
        cat(paste0("Ecov ",i,": ",ecov$label[i],"
",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how+1]," (",ecov_str[[i]],") effect on: ", c('recruitment','M')[data$Ecov_where[i]],"

In model years:
"))
      }
      if(data$Ecov_where[i] == 2){ # M
        cat(paste0("Ecov ",i,": ",ecov$label[i],"
",ecov_str[[i]]," effect on: ", c('recruitment','M')[data$Ecov_where[i]],"

In model years:
"))
      }

cat(years, fill=TRUE)
lastyr <- tail(years,1)
cat(paste0("Lag: ",data$Ecov_lag[i],"
Ex: ",ecov$label[i]," in ",years[1]," affects ", c('recruitment','M')[data$Ecov_where[i]]," in ",years[1+data$Ecov_lag[i]],"
    ",ecov$label[i]," in ",lastyr," affects ", c('recruitment','M')[data$Ecov_where[i]]," in ",lastyr+data$Ecov_lag[i],"
"))
    }
  } # end load Ecov

  # add vector of all observations for one step ahead residuals ==========================
  # 5 components: fleet catch (log), index catch (log), Ecov, paa catch, paa index
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 1. log fleet catch
  x <- as.data.frame(data$agg_catch)
  x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
  colnames(x) <- paste0("fleet_", 1:data$n_fleets)
  x$year <- 1:data$n_years_catch
  tmp <- tidyr::gather(x, fleet, val, -year)
  tmp <- tmp[complete.cases(tmp),]  
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_agg_catch==0, turned to NA, and removed
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

  # projection data will always be modified by 'prepare_projection'
  data$do_proj <- 0
  data$n_years_proj <- 0
  data$n_years_proj_Ecov <- 0
  data$avg_years_ind <- 0
  data$proj_F_opt <- 0
  data$proj_Fcatch <- 0

  # data$obsvec[data$keep_I[data$use_indices==1]+1] - log(data$agg_indices[data$use_indices==1])
  # data$obsvec[data$keep_E[data$Ecov_use_obs==1]+1] - data$Ecov_obs[data$Ecov_use_obs==1]

  # -------------------------------------------------------------------
  # Parameters
  par = list()
  par$mean_rec_pars = numeric(c(0,1,2,2)[recruit_model])
  if(recruit_model==2) par$mean_rec_pars = 10
  if(recruit_model==4) par$mean_rec_pars[2] = -10
  par$logit_q = rep(-8, data$n_indices)
  par$log_F1 = rep(-2, data$n_fleets)
  par$F_devs = matrix(0, data$n_years_model-1, data$n_fleets)
  if(data$N1_model == 1) par$log_N1_pars = c(10,log(0.1))
  if(data$N1_model == 0) par$log_N1_pars = rep(10,data$n_ages)
  par$log_NAA_sigma = rep(0, data$n_NAA_sigma)

  # selectivity pars
  par$logit_selpars = log(selpars_ini-selpars_lo) - log(selpars_hi - selpars_ini)
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] = -10
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] = 10
  # number of estimated selpars per block * number of years per block (only if that block has re)
  if(any(data$selblock_models_re > 1)){
    par$selpars_re <- rep(0, sum((data$selblock_models_re > 1)*data$n_selpars_est*data$n_years_selblocks))
  } else {
    par$selpars_re <- matrix(0)
    map$selpars_re <- factor(NA)
  }
  par$sel_repars <- matrix(0, nrow=data$n_selblocks, ncol=3)
  par$sel_repars[,1] <- log(0.1) # start sigma at 0.1, rho at 0
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 3) par$sel_repars[b,3] <- 0 # if ar1 over ages only, fix rho_y = 0
    if(data$selblock_models_re[b] == 4) par$sel_repars[b,2] <- 0 # if ar1 over years only, fix rho = 0
    # check if only 1 estimated sel par (e.g. because all but 1 age is fixed), can't estimate rho
    if(data$n_selpars_est[b] < 2) par$sel_repars[b,2] <- 0
  }

  # age comp pars
  n_catch_acomp_pars = c(0,1,1,3,1,2)[data$age_comp_model_fleets[which(apply(data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[data$age_comp_model_indices[which(apply(data$use_index_paa,2,sum)>0)]]
  par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  par$index_paa_pars = rep(0, sum(n_index_acomp_pars))
  par$log_NAA = matrix(10, data$n_years_model-1, data$n_ages)

  # natural mortality pars
  par$M0 <- M0_ini # mean M
  par$M_a <- M_a_ini # deviations by age
  par$M_re <- M_re_ini # deviations from mean M_a on log-scale, PARAMETER_ARRAY
  par$M_repars <- rep(0, 3)
  par$M_repars[1] <- log(0.1) # start sigma at 0.1, rho at 0
  if(data$M_re_model == 3) par$M_repars[3] <- 0 # if ar1 over ages only, fix rho_y = 0
  if(data$M_re_model == 4) par$M_repars[2] <- 0 # if ar1 over years only, fix rho_a = 0
  if(data$M_re_model == 5) {par$M_repars[2] <- 0; par$M_repars[3] <- Inf}  # if rw over years, fix rho_a = 0 and rho_y = 1
  # check if only 1 estimated mean M (e.g. because weight-at-age M or if all but 1 age is fixed), can't estimate rho_a
  if(data$n_M_est < 2) par$M_repars[2] <- 0
  par$log_b = log(0.305)

  par$log_R = rep(10, data$n_years_model-1) #/n_years_model-1, if used.
  par$log_R_sigma = 0
  par$log_catch_sig_scale = rep(0, data$n_fleets)
  par$log_index_sig_scale = rep(0, data$n_indices)

  # Ecov pars
  par$Ecov_re = matrix(0, data$n_years_Ecov, data$n_Ecov)
  max.poly <- max(data$Ecov_poly)
  par$Ecov_beta = matrix(0, nrow=max.poly, ncol=data$n_Ecov) # beta_R in eqns 4-5, Miller et al. (2016)
  par$Ecov_process_pars = matrix(0, 3, data$n_Ecov) # nrows = RW: 2 par (log_sig, Ecov1), AR1: 3 par (mu, phi, log_sig); ncol = N_ecov
  par$Ecov_obs_logsigma <- par.Ecov.obs.logsigma
  par$Ecov_obs_sigma_par <- par.Ecov.obs.sigma.par

  # turn off 3rd Ecov par if it's a RW
  tmp.pars <- par$Ecov_process_pars
  for(i in 1:data$n_Ecov) tmp.pars[3,i] <- ifelse(data$Ecov_model[i]==1, NA, 0)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)

  # turn off only Ecov_beta to fit Ecov process model without effect on population
  tmp <- par$Ecov_beta
  for(j in 1:data$n_Ecov){
    tmp[,j] <- ifelse(data$Ecov_how[j]==0, NA, 0)
    for(i in 1:max.poly){
      if(i > data$Ecov_poly) tmp[i,j] = NA
    }
  }
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_beta = factor(tmp)

  # turn off Ecov pars if no Ecov (re, process)
  # for any Ecov_model = NA, ecov$how must be 0 and beta is already turned off
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

  # map other pars
  map$log_catch_sig_scale = factor(rep(NA, data$n_fleets))
  map$log_index_sig_scale = factor(rep(NA, data$n_indices))

  if(data$n_M_est == 0) map$M0 = factor(NA)
  # map$M_a = factor(rep(NA, length(par$M_a)))
  tmp <- par$M_a
  tmp[M_first_est] = NA # M0 represents M for first estimated age, fix that M_a deviation = 0
  tmp[data$M_est==0] = NA
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$M_a <- factor(tmp)

  # map$M_re = factor(rep(NA, length(par$M_re)))
  tmp <- par$M_re
  # if(data$M_re_model == 1){
  #   tmp[] = NA
  # } else {
  #   tmp[,data$M_est==0] = NA # turn off RE for ages that aren't estimated
  # }
  if(data$M_re_model == 1) tmp[] = NA # either estimate RE for all ages or none at all
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$M_re <- factor(tmp)

  # M_repars: sigma_M, rho_M_a, rho_M_y
  # M re models: "none","iid","ar1_a","ar1_y","rw_y","2dar1"
  if(data$M_re_model == 1) tmp <- rep(NA,3) # no RE pars to estimate
  if(data$M_re_model == 2) tmp <- c(1,NA,NA) # estimate sigma
  if(data$M_re_model == 3) tmp <- c(1,2,NA) # ar1_a: estimate sigma, rho_a
  if(data$M_re_model == 4) tmp <- c(1,NA,2) # estimate sigma, rho_y
  if(data$M_re_model == 5) tmp <- c(1,NA,NA) # rw_y: estimate sigma only (fix rho_a = 0, rho_y = 1)
  if(data$M_re_model == 6) tmp <- 1:3 # 2dar1: estimate all
  if(data$n_M_est < 2) tmp[2] <- NA # can't estimate rho_a if M estimated for < 2 ages
  map$M_repars = factor(tmp)

  # map$M_sigma_pars = factor(rep(NA, length(par$M_sigma_pars)))
  map$log_NAA = factor(rep(NA, length(par$log_NAA)))
  map$log_NAA_sigma = factor(rep(NA, length(par$log_NAA_sigma)))
  map$mean_rec_pars = factor(rep(NA, length(par$mean_rec_pars)))
  map$log_R_sigma = factor(rep(NA, length(par$log_R_sigma)))
  map$log_b = factor(rep(NA,length(par$log_b)))
  map$Ecov_obs_logsigma <- factor(map.Ecov.obs.logsigma)
  map$Ecov_obs_sigma_par <- factor(map.Ecov.obs.sigma.par)

  # map selectivity RE
  tmp.sel.repars <- par$sel_repars
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 1) tmp.sel.repars[b,] <- rep(NA,3) # no RE pars to estimate
    if(data$selblock_models_re[b] == 2) tmp.sel.repars[b,2:3] <- rep(NA,2) # estimate sigma
    if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,3] <- NA # estimate sigma, rho
    if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,2] <- NA # estimate sigma, rho_y
    if(data$n_selpars_est[b] < 2) tmp.sel.repars[b,2] <- NA # can't estimate rho if only 1 selpar estimated
  }
  ind.notNA <- which(!is.na(tmp.sel.repars))
  tmp.sel.repars[ind.notNA] <- 1:length(ind.notNA)
  map$sel_repars = factor(tmp.sel.repars)

  # tmp.sel.re <- par$selpars_re
  # istart <- 1
  # for(b in 1:data$n_selblocks){
  #   iend <- istart + data$n_years_model*data$n_selpars[b] - 1
  #   if(data$selblock_models_re[b] == 1) tmp.sel.re[istart:iend] <- NA # no RE on selectivity
  #   istart <- istart + data$n_years_model*data$n_selpars[b]
  # }
  # ind.notNA <- which(!is.na(tmp.sel.re))
  # tmp.sel.re[ind.notNA] <- 1:length(ind.notNA)
  # map$selpars_re = factor(tmp.sel.re)

  random = character()
  if(data$Ecov_obs_sigma_opt == 4) random = "Ecov_obs_logsigma"
  if(any(data$selblock_models_re > 1)) random = c(random, "selpars_re")
  if(data$M_re_model > 1) random = c(random, "M_re")
  if(sum(data$Ecov_model) > 0) random = c(random, "Ecov_re")
  if(missing(model_name)) model_name = "WHAM for unnamed stock"
  return(list(data=data, par = par, map = map, random = random, years = model_years, years_full = model_years,
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
