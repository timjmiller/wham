#' Specify model and parameter configuration for movement when input$data$n_regions > 1
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param move (optional) list specifying movement options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{move} specifies estimation options for movement.
#' If \code{NULL}, no movement will occur. If there are multiple regions, each stock will be modeled separately in different regions without movement. 
#' \code{move} is a list with the following entries:
#'   \describe{
#'     \item{$stock_move}{length = n_stocks, T/F whether each stock can move. If not provided then movement will be defined below for all stocks.}
#'     \item{$separable}{length = n_stocks, T/F whether movement should be modeled separably from mortality or both occuring simultaneously.}
#'     \item{$mean_model}{"(mean) movement model options are:
#'                    \describe{
#'                      \item{"none"}{(default) no movement between regions.}
#'                      \item{"constant"}{estimate a single movement rate to each region shared across all stocks, seasons, ages, years}
#'                      \item{"season"}{estimate movement rates  to each region for each season shared across all stocks, ages, years}
#'                      \item{"stock_constant"}{estimate a movement rate for each stock to each region shared across all seasons, ages, years}
#'                      \item{"stock_season"}{estimate a movement rate for each stock each season to each region shared across all ages, years}
#'                    }
#'                  }
#'     \item{$re_model}{"random effects correlation structure options are:
#'                    \describe{
#'                      \item{"none"}{(default) no movement rate random effects.}
#' "age_iid","age_ar1","year_iid","year_ar1", "age_iid_year_iid", "age_ar1_year_iid", "age_iid_year_ar1", "age_ar1_year_ar1"
#'                      \item{"age_iid"}{movement rates for each age are random effects with means defined by mean_model}
#'                      \item{"season"}{estimate mortality rates  to each region for each season shared across all stocks, ages, years}
#'                      \item{"stock_constant"}{estimate a mortality rate for each stock to each region shared across all seasons, ages, years}
#'                      \item{"stock_season"}{estimate a mortality rate for each stock each season to each region shared across all ages, years}
#'                      \item{"iid_re"}{estimate independent random effects over years, for the region}
#'                      \item{"ar1_re"}{estimate random effect correlated over years, for the region}
#'                    }
#'                  }
#'     \item{$prior_sigma}{array (n_stocks x n_seasons x n_regions x n_regions - 1) of sd parameters for normal priors on movement parameters on transformed scale (-Inf,Inf)}
#'     \item{$use_prior}{array (n_stocks x n_seasons x n_regions x n_regions - 1) 0/1 indicator whether to use prior for movement parameters.}
#'     \item{$can_move}{array (n_stocks x n_ages x n_seasons x n_regions x n_regions - 1) 0/1 indicator whether movement can occur from one region to another.}
#'     \item{$must_move}{array (n_stocks x n_ages x n_seasons x n_regions) 0/1 indicator whether movement from region must occur.}
#'     \item{$initial_means}{Initial/mean M-at-region}
#'     \item{$sigma_vals}{Initial standard deviation by region value to use for the L random effects. Values are not used if \code{L$model} = "none".
#'     \item{$cor_vals}{Initial correlation values to use for the L random effects. If unspecified all initial values are 0}
#'   }

set_move = function(input, move)
{
  data = input$data
  par = input$par
  map = input$map
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mu_repars", "mu_re"))]
  
  data$mig_type = rep(1,data$n_stocks)
  if(!is.null(move$separable)) data$mig_type[] = as.integer(move$separable)

  data$trans_mu_prior_sigma = array(0.05, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  if(!is.null(move$prior_sigma)) data$trans_mu_prior_sigma[] = move$prior_sigma
  data$use_mu_prior = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  if(!is.null(move$use_prior)) data$use_mu_prior[] = move$use_prior
  data$can_move = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  if(!is.null(move$can_move)) data$can_move[] = move$can_move
  data$must_move = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions))
  if(!is.null(move$must_move)) data$must_move[] = move$must_move


############right here
  data$mu_model <- 1
  move_mods <- c("none","constant", "season")
  move_re_mods <- c("age_iid","age_ar1","year_iid","year_ar1", "age_iid_year_iid", "age_ar1_year_iid", "age_iid_year_ar1", "age_ar1_year_ar1")
  move_mods <- c(move_mods, paste0(move_mods[2:3], "_stock"))
  
  #L$model length is n_regions
  data$L_model = rep(0, data$n_regions)
  par$L_re = matrix(0, data$n_years_model, data$n_regions)
  par$L_repars = matrix(0,data$n_regions, 3) #mean, sig, rho
  map$L_re = matrix(NA, data$n_years_model, data$n_regions)
  map$L_repars = matrix(NA,data$n_regions, 3) #mean, sig, rho

  # natural mortality options, default = use values from ASAP file, no estimation
  if(!is.null(L)){
    if(!is.null(L$model)){ # L model options
      L_mods = c("none","constant","iid_re","ar1_re")
      if(!(L$model %in% L_mods)) stop(paste0("L$model must be one of these: ", paste0(L_mods, collapse=","))
      data$L_model[] = match([L$model, L_mods]) - 1
    }
  }
  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s
  k = 1
  for(r in 1:data$n_regions) {
    if(data$L_model[r] >0) {
      map$L_repars[r,1] <- k
      k <- k + 1
      if(!is.null(L$initial_means)){
        par$L_repars[r,1] = log(L$initial_means[r])
      }
    }
    if(data$L_model[r] > 1){
      map$L_repars[r,2] <- k
      k <- k + 1
      map$L_re[,r] <- 1
      if(!is.null(L$sigma_vals)){
        par$L_repars[r,2] = log(L$sigma_vals[r])
      }
    }
    if(data$L_model[r] > 2){
      map$L_repars[r,3] <- k
      k <- k + 1
      if(!is.null(L$cor_vals)){
        par$L_repars[r,3] = inv_trans_rho(L$cor_vals[r])
      }
    }
  }
  map$L_re[which(map$L_re==1)] <- 1:sum(map$L_re==1)
  map$L_re = factor(map$L_re)
  map$L_repars = factor(map$L_repars)

  input$data = data
  input$par = par
  input$map = map
  
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = wham:::set_random(input)
  return(input)

}