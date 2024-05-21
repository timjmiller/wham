#' Specify model and parameter configuration for catchability
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param catchability (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' 
#' \code{catchability} specifies options for catchability. If \code{NULL} and \code{asap3} is not NULL, a single catchability parameter for each index is used with initial values specified in ASAP file. If both are NULL, initial catchabilities for all indices = 0.3.
#' Otherwise, it is a list with the following optional entries:
#'   \describe{
#'     \item{$re}{Time-varying (random effects) for each index. Vector with length = number of indices.
#'                  Each entry of \code{catchability$re} must be one of the following options:
#'                  \describe{
#'                    \item{"none"}{(default) are constant}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by year (AR1)}
#'                  }
#'                 }
#'     \item{$initial_q}{Initial catchabilities for each index. vector length = number of indices. Will override values provided in \code{basic_info$q}.
#'        If \code{basic_info$q} and \code{asap3} are not provided, default q values are 0.3.}
#'     \item{$q_lower}{Lower bound for catchabilities for each index. vector length = number of indices. For indices with NULL components default lower values are 0.}
#'     \item{$q_upper}{Upper bound for catchabilities for each index. vector length = number of indices. For indices with NULL components default lower values are 1000.}
#'     \item{$prior_sd}{vector of NA and standard deviations to use for gaussian prior on logit transform of catchability parameter. Length = number of indices.
#'       Indices with non-NA values will have mean logit q as a random effect with mean determined by logit transform of \code{catchability$initial_q} and
#'       sigma as standard error.}
#'     \item{$sigma_val}{Vector of initial standard deviation values to use for annual random effects for each index. Values are not used if \code{q$re} = "none". Otherwise, a single value for all indices.}
#'     \item{$sigma_map}{Specify which sigma parameters to fix for the random effect deviations. Must be a vector with length = number of indices. 
#'                Use \code{NA} to fix a parameter and integers to estimate. Use the same integer for multiple indices to share the same sigma parameter.
#'                Not used if \code{re = 'none'} for all indices.}
#'     \item{$cor_vals}{Vector of initial correlation values to use for annual random effects for each index. If unspecified all initial values are 0. Only used if \code{q$re} = "ar1"}
#'     \item{$cor_map}{Specify which ar1 correlation parameters to fix for the random effect deviations. Must be a vector with length = number of indices. 
#'                If \code{re = 'ar1'}, each element (index) must be a single value. Use \code{NA} to fix a parameter and integers to estimate. 
#'                Use the same integer for multiple indices to share the same correlation parameter. Not used if \code{re = 'none'} or \code{re = 'iid'} for all indices.}
#'   }
#'
#' @return a named list with same elements as the input provided with catchability options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' catchability = list(re = c("iid", "none"))
#' input <- set_q(input, catchability = catchability) #independent time-varying random effects on q for first survey.
#' }
#'
#' @export
set_q = function(input, catchability = NULL){
	
  q_opts = catchability
  data = input$data
  par = input$par
  map = input$map
  if(is.null(data$n_indices)) stop("No index information has been added yet. Run set_indices() first. See ?set_indices.")

  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("q_prior_re", "q_re", "q_repars"))]
	input$log$q <- list()
	asap3 = input$asap3
	#if(is.null(input$asap3)) asap3 = NULL
  #else asap3 = input$asap3
  data$q_lower <- rep(0,data$n_indices)
  data$q_upper <- rep(1000,data$n_indices)
  if(!is.null(q_opts$q_lower)) {
		if(length(q_opts$q_lower) != data$n_indices) stop("the length of catchability$q_lower is not equal to the number of indices")
  	data$q_lower = q_opts$q_lower
  }
  if(!is.null(q_opts$q_upper)) {
		if(length(q_opts$q_upper) != data$n_indices) stop("the length of catchability$q_upper is not equal to the number of indices")
  	data$q_upper = q_opts$q_upper
  }
	par$logit_q = gen.logit(0.3, data$q_lower, data$q_upper)
  if(!is.null(asap3)) {
		k = 0
		for(i in 1:length(asap3)){
			ind = which(asap3[[i]]$use_index ==1)
			k_ind = k + 1:length(ind)
			par$logit_q[k_ind] = gen.logit(asap3[[i]]$q_ini[ind], data$q_lower[k_ind], data$q_upper[k_ind]) # use q_ini values from asap3 file
			k = max(k_ind)
		}
	}
	if(!is.null(q_opts$initial_q)) {
		if(length(q_opts$initial_q) != data$n_indices) stop("the length of catchability$initial_q is not equal to the number of indices")
		par$logit_q = gen.logit(q_opts$initial_q, data$q_lower, data$q_upper)
	}

	data$use_q_re = rep(0, data$n_indices)
	par$q_re = matrix(0, data$n_years_model, data$n_indices)
	map$q_re = matrix(NA, data$n_years_model, data$n_indices)
	par$q_repars = matrix(0, data$n_indices, 2)
	map$q_repars = matrix(NA, data$n_indices, 2)

	if(!is.null(q_opts$re)){
  	ind = which(q_opts$re %in% c("iid", "ar1"))
		map$q_re[,ind] = 1:(length(ind) * (data$n_years_model)) #turn on appropriate columns of q_re
		if(!is.null(q_opts$sigma_val)){
			if(length(q_opts$sigma_val) != data$n_indices) stop("the length of catchability$sigma_val provided is not equal to the number of indices")
			if(q_opts$sigma_val[ind]<0) stop("catchability$sigma_val must be greater than 0")
			par$q_repars[ind,1] = log(q_opts$sigma_val[ind])
		}
	  iids = q_opts$re == "iid"
	  ar1s = q_opts$re == "ar1"
		if(!is.null(q_opts$cor_val)){
			if(length(q_opts$cor_val) != data$n_indices) stop("the length of catchability$cor_val provided is not equal to the number of indices")
			if(abs(q_opts$cor_val[ind])>1) stop("it must be that -1 < catchability$cor_val < 1 ")
			if(any(iids & abs(q_opts$cor_val)>1e-10)) input$log$q <- c(input$log$q, "certain indices have re='iid' and cor_val not = 0. Those will values will be ignored. \n")
			par$q_repars[ind,2] = gen.logit(q_opts$cor_val[ind], -1, 1, 1) 
			par$q_repars[which(iids),2] = 0 #iids must be set to 0.

		}
	  if(!is.null(q_opts$sigma_map)){
	    if(length(q_opts$sigma_map) != data$n_indices) stop("catchability$sigma_map must be a vector of length = number of indices")
	    map$q_repars[which(!is.na(q_opts$sigma_map)),1] = 1:sum(!is.na(q_opts$sigma_map))
	  } else map$q_repars[ind,1] = 1:length(ind)

	  if(!is.null(q_opts$cor_map)){
	    if(length(q_opts$cor_map) != data$n_indices) stop("catchability$cor_map must be a list of length = number of indices")
	    if(any(iids & !is.na(q_opts$cor_map))) stop("certain indices have re='iid' and cor_map is not NA.")
	    map$q_repars[which(!is.na(q_opts$cor_map)),2] = sum(!is.na(map$q_repars[,1])) + 1:sum(!is.na(q_opts$cor_map))
	    #"iid" not necessary because already set to NA
	  } else {
			if(sum(ar1s)) map$q_repars[ind,2] = sum(!is.na(map$q_repars[,1])) + 1:sum(ar1s)	  	
	  }
		data$use_q_re[ind] = 1
	}
	map$q_repars = factor(map$q_repars)
  map$q_re = factor(map$q_re)

	data$use_q_prior = rep(0, data$n_indices)
	par$q_prior_re = rep(0, data$n_indices)
	map$q_prior_re = rep(NA, data$n_indices)
	data$logit_q_prior_sigma = rep(1, data$n_indices)
	map$logit_q = 1:data$n_indices
	if(!is.null(q_opts$prior)){
  	ind = which(!is.na(q_opts$prior_sd))
		data$use_q_prior[ind] = 1
		data$logit_q_prior_sigma[ind] = q_opts$prior_sd[ind]
		par$q_prior_re[ind] = par$logit_q[ind]
		map$q_prior_re[ind] = 1:length(ind)
		map$logit_q[ind] = NA #turn off logit q (mean of the prior) when random effects are used
	}
	map$q_prior_re = factor(map$q_prior_re)
	map$logit_q = factor(map$logit_q)

  input$data = data
  input$par = par
  input$map = map
	if(length(input$log$q))	input$log$q <- c("Catchability: \n", input$log$q)

  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
	input$options$q <- q_opts
  if(!is_internal_call()) cat(unlist(input$log$q, recursive=T))
	return(input)
}