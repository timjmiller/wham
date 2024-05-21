#'  Specify the age composition models for fleet(s) and indices.
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param age_comp specifies the age composition models for fleet(s) and indices. If \code{NULL}, the multinomial is used because this was the only option in ASAP. 
#'
#' The age composition models available are:
#'   \describe{
#'     \item{\code{"multinomial"}}{Multinomial. This is the default because it was the only option in ASAP. 0 parameters.}
#'     \item{\code{"dir-mult"}}{Saturating Dirichlet-multinomial, parameterized such that effective-sample-size is a nonlinear and saturating function with respect to input-sample-size. 1 parameter. Effective sample size is estimated by the model (\href{https://www.ccamlr.org/es/system/files/science_journal_papers/07candy.pdf}{Candy 2008})}
#'     \item{\code{"dirichlet-pool0"}}{Dirichlet, pooling zero observations with adjacent age classes. 1. parameter. See \href{https://www.sciencedirect.com/science/article/abs/pii/S0165783613003093}{Francis 2014} and \href{https://cdnsciencepub.com/doi/abs/10.1139/cjfas-2015-0532}{Albertsen et al. 2016}}
#'     \item{\code{"dirichlet-miss0"}}{}{Dirichlet, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-miss0"}}{Logistic normal, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-ar1-miss0"}}{Logistic normal, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-pool0"}}{Logistic normal, pooling zero observations with adjacent age classes. 1 parameter. See \href{https://doi.org/10.1093/icesjms/fsl024}{Schnute and Haigh (2007)} and \href{https://doi.org/10.1016/j.fishres.2013.12.015}{Francis (2014)}}.
#'     \item{\code{"logistic-normal-01-infl"}}{Zero-or-one inflated logistic normal. Inspired by zero-one inflated beta in \href{https://www.sciencedirect.com/science/article/abs/pii/S0167947311003628}{Ospina and Ferrari (2012)}. 3 parameters. . No OSA residuals.}
#'     \item{\code{"logistic-normal-01-infl-2par"}}{Zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters. No OSA residuals.}
#'     \item{\code{"mvtweedie"}}{Multivariate-tweedie, where the product of composition proportions and input sample sizes follows a distribution with mean equal to the product of predicted proportions and input sample size, and other parameters define the ratio of effective to input sample size (with is bounded 0 to Inf) and the probability of zeros. 2 parameters. No OSA residuals.}
#'     \item{\code{"dir-mult-linear"}}{Linear Dirichlet-multinomial, parameterized such that effective-sample-size is a linear function with respect to input-sample-size, estimating 1 parameter, \eqn{log(\theta)}, where the ratio of effective and input sample size is approximately \eqn{\theta / (1+\theta)}, i.e., the logistic transformation of the estimated parameter \eqn{log(\theta)}.  (\href{https://doi.org/10.1016/j.fishres.2016.06.005}{Thorson et al. 2017}) }
#'   }
#' The two Dirichlet-multinomial options will only differ when input-sample-size differs among years.  In these cases, the linear-Dirichlet multinomial is designed to decrease the effective sample size in each year by approximately the same proportion, while the saturating-Dirichlet multinomial will decrease the years with highest input-sample-size much more than those with lower input-sample-size.
#' One-step-ahead residuals will be calculated for all but options 8-10 when \code{do.osa=TRUE} (Nielsen et al. in prep.). An age composition model needs
#' to be specified for each fleet and index. If you would like all fleets and indices to use the same age composition likelihood, you 
#' can simply specify one of the strings above, i.e. \code{age_comp = "logistic-normal-miss0"}. If you do not want the same
#' age composition model to be used for all fleets and indices, you must specify a named list with the following entries:
#'   \describe{
#'     \item{$fleets}{A vector of the above strings with length = the number of fleets.}
#'     \item{$indices}{A vector of the above strings with length = the number of indices.}
#'   }
#'
#' @return a named list with same elements as the input provided with age composition likelihood options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' input <- set_age_comp(input, age_comp = "logistic-normal-miss0") #no longer multinomial
#' }
#'
#' @export
#'
set_age_comp = function(input, age_comp)
{
	data = input$data
	par = input$par
  map = input$map
  all_models <- c( "multinomial",
                   "dir-mult",
                   "dirichlet-miss0",
                   "dirichlet-pool0",
                   "logistic-normal-miss0",
                   "logistic-normal-ar1-miss0",
                   "logistic-normal-pool0",
                   "logistic-normal-01-infl",
                   "logistic-normal-01-infl-2par",
                   "mvtweedie",
                   "dir-mult-linear" )
  input$log$age_comp <- list()
  n_pars <- c(0,1,1,1,1,2,1,3,2,2,1)
  if(is.null(age_comp)){
    if(!is.null(data$n_fleets)) data$age_comp_model_fleets = rep(1, data$n_fleets) # multinomial by default
    if(!is.null(data$n_indices)) data$age_comp_model_indices = rep(1, data$n_indices) # multinomial by default
  } else {
    if(is.character(age_comp) & length(age_comp)==1){ # all use the same model
      name_change = age_comp == "dirichlet"
      if(any(name_change)){
        input$log$age_comp <- c(input$log$age_comp, "'dirichlet' is no longer an option and the old option is equivalent to 'dirichlet-pool0' so using that.\n")
        age_comp[which(name_change)] = "dirichlet-pool0"
      }
      themod <- match(age_comp, all_models)
      if(is.na(themod)) stop("age_comp option not recognized. See ?prepare_wham_input.")
      if(!is.null(data$n_fleets)) data$age_comp_model_fleets = rep(themod, data$n_fleets)
      if(!is.null(data$n_indices)) data$age_comp_model_indices = rep(themod, data$n_indices)
    } else {
      if(all(names(age_comp) %in% c("fleets","indices"))){
        name_change = list(fleets = age_comp$fleets == "dirichlet", indices = age_comp$indices == "dirichlet")
        if(any(unlist(name_change))){
          input$log$age_comp <- c(input$log$age_comp, "'dirichlet' is no longer an option and the old option is equivalent to 'dirichlet-pool0' so using that.\n")
          age_comp$fleets[which(name_change$fleets)] = "dirichlet-pool0"
          age_comp$indices[which(name_change$indices)] = "dirichlet-pool0"
        }
        themods <- match(age_comp$fleets, all_models)
        if(any(is.na(themods))) stop("age_comp$fleets option not recognized. See ?prepare_wham_input for available options.")
        if(!is.null(data$n_fleets)) if(length(themods) != data$n_fleets) stop("age_comp$fleets must have length = the number of fleets")
        data$age_comp_model_fleets = themods

        themods <- match(age_comp$indices, all_models)
        if(any(is.na(themods))) stop("age_comp$indices option not recognized. See ?prepare_wham_input for available options.")
        if(!is.null(data$n_indices)) if(length(themods) != data$n_indices) stop("age_comp$indices must have length = the number of indices")
        data$age_comp_model_indices = themods
      } else {
        stop("age_comp must either be a character or a named ('fleets','indices') list. See ?prepare_wham_input.")
      }
    }
  }

  # age comp pars
  if(!is.null(data$n_fleets)) { #catch data have been added
    par$catch_paa_pars = matrix(0,data$n_fleets, 3) 
    neff <- data$catch_Neff
    neff[neff <= 0] <- 1
    catch_neff <- apply(neff,2,mean, na.rm=TRUE)
    ind = which(data$age_comp_model_fleets %in% 5:7)
    par$catch_paa_pars[ind,1] <- 0.5*log(catch_neff[ind])
    map$catch_paa_pars = matrix(NA,data$n_fleets, 3)
    for(i in 1:data$n_fleets) if(sum(data$use_catch_paa[,i])){
      if(data$age_comp_model_fleets[i] %in% c(2:5,7,11))  map$catch_paa_pars[i,1] = 1
      if(data$age_comp_model_fleets[i] %in% c(6,9,10))  map$catch_paa_pars[i,1:2] = 1
      if(data$age_comp_model_fleets[i] %in% 8)  map$catch_paa_pars[i,1:3] = 1
    }
    nest = sum(map$catch_paa_pars,na.rm=TRUE)
    if(nest) map$catch_paa_pars[which(!is.na(map$catch_paa_pars))] = 1:nest
    map$catch_paa_pars = factor(map$catch_paa_pars)
  }
  
  if(!is.null(data$n_indices)) { #index data have been added
    par$index_paa_pars = matrix(0,data$n_indices, 3) 
    neff <- data$index_Neff
    neff[neff <= 0] <- 1
    index_neff <- apply(neff,2,mean, na.rm=TRUE)#[which(apply(data$use_index_paa,2,sum)>0)]
    ind = which(data$age_comp_model_indices %in% 5:7)
    par$index_paa_pars[ind,1] <- 0.5*log(index_neff[ind])
    map$index_paa_pars = matrix(NA,data$n_indices, 3)
    for(i in 1:data$n_indices) if(sum(data$use_index_paa[,i])){
      if(data$age_comp_model_indices[i] %in% c(2:5,7,11))  map$index_paa_pars[i,1] = 1
      if(data$age_comp_model_indices[i] %in% c(6,9,10))  map$index_paa_pars[i,1:2] = 1
      if(data$age_comp_model_indices[i] %in% 8)  map$index_paa_pars[i,1:3] = 1
    }
    nest = sum(map$index_paa_pars,na.rm=TRUE)
    if(nest) map$index_paa_pars[which(!is.na(map$index_paa_pars))] = 1:nest
    map$index_paa_pars = factor(map$index_paa_pars)
  }

	input$data = data
	input$par = par
  input$map = map
	if(length(input$log$age_comp))	input$log$age_comp <- c("Age composition: \n", input$log$age_comp)
  input$options$age_comp <- age_comp
  
  if(!is_internal_call()) { #check whether called by prepare_wham_input
    input <- set_osa_obs(input)
    cat(unlist(input$log$age_comp, recursive=T))
  }
 
	return(input)
}
