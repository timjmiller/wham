#' Specify the probability models for aggregate observation for fleets and indices.
#'
#' WHEN SIMULATING OBSERVATIONS THAT CAN BE ZERO: The number and location of the zero observations is random and the CV will be considered missing and treated as a random effect
#'   for those years. The estimated CVs for years that were 0 in the original data are needed if there is now a positive observations in that location. The CVs provided in the input are
#'   are overwritten with the estimated values (corrected for any estimated scalars for lognormal and gamma options). The necessary number of CV random effect parameters must be redetermined
#'   for each simulated data set (say, sim) before refitting the data by calling \code{sim <- set_agg_obs_model(sim, agg_obs_model)}. The presence of age composition observations may also be modified
#'   in that if the original data had both a positive aggregate observation and age comp observations, but now is a 0, then it is not possible to have age composition observations and WHAM 
#'   changes \code{data$use_catch_paa} or \code{data$use_index_paa} as necessary. The required call to \code{set_agg_obs_model} also calls \code{set_osa_obs} so the observation vector used 
#'   to calculate likelihods and OSA residuals is updated.
#'
#' @param input list containing data, parameters, map, and random elements (e.g., output from \code{\link{wham::prepare_wham_input}})
#' @param agg_obs_model list specifying the aggregate observation models, and other configuration options for each fleet and index. If \code{NULL}, the "lognormal-1" is used as it was previously 
#'   the only option. 
#' 
#'   \item{\code{obs_model}}{a single character string (applied to all fleets and indices) or a named list of two vectors ("fleet" and "index") of character strings of probability models 
#'     for annual aggregate observations:
#'     \describe{Below \eqn{\widehat x} is the predicted aggregate *positive* observation (for fleets) or index whereas \eqn{x} is the observed random variable that can = 0 for some models.
#'       \item{\code{"lognormal-1"}}{Log-normal distribution for observation: \eqn{E(\log(x)) = \log(\widehat{x})}, \eqn{\text{SD}(\log(x)) = \sqrt(\log(\text{CV}_+(x)^2 + 1)) * \phi}. 
#'         This will be used by default with \eqn{\phi = 1}.
#'       }
#'       \item{\code{"lognormal-2"}}{Log-normal distribution for observation: \eqn{E(\log(x)) = \log(\widehat{x})}, \eqn{\text{SD}(\log(x)) = \sqrt(\log(\text{CV}_+(x)^2 + 1)) * (\phi + 1)}. 
#'         \eqn{\phi} estimated by default. SD must be greater than that implied by input CV.
#'       }
#'       \item{\code{"gamma-1"}}{Gamma distribution for observation: \eqn{E(x) = \widehat x}, \eqn{SD(x) = \phi * \text{CV}_+(x)) * x}. \eqn{\phi} is set = 1 by default.}
#'       \item{\code{"gamma-2"}}{Gamma distribution for observation: \eqn{E(x) = \widehat x}, \eqn{SD(x) = (\phi + 1) * \text{CV}_+(x)) * x}; \eqn{\phi} is estimated by default. SD must be 
#'         greater than that implied by input CV.
#'       }
#'       \item{\code{"gengamma"}}{Generalized gamma distribution (Rigby et al. 2019 parameterization) for observation: \eqn{E(x) = \widehat x}, with parameter \eqn{\nu>0} estimated by 
#'         default, and parameter \eqn{\sigma>0} defined such that 
#'         \deqn{CV(x)^2 = \Gamma(\theta)\Gamma(\theta + 2/\nu)/\Gamma(\theta + 1/\nu)^2 - 1} where \eqn{\theta = (\sigma\nu)^{-2}}.
#'       } 
#'       \item{\code{"CPG-alpha"}}{Compound Poisson Gamma (Tweedie with power paraemter \eqn{1<p<2}) distribution for observation. For all versions:
#'         \deqn{E(x) = \widehat x p_+ = \mu = \lambda\alpha \beta}
#'         \deqn{V(x) = \mu^p * \phi  = (\text{CV}_+^2 + 1- p_+) * \mu^2/p_+ = \lambda\alpha(\alpha+1)\beta^2} with \eqn{1<p<2} where
#'         \eqn{p_+ = 1-\exp(-\lambda)} is the probability of a positive value (Poisson RV neq 0). For this version ("CPG-alpha"), \eqn{\alpha} is estimated by default and other parameters 
#'         are defined by \eqn{\alpha}, \eqn{\widehat x} and \eqn{CV_+}.
#'       }
#'       \item{\code{"CPG-beta"}}{Compound Poisson Gamma distribution with fraction of maximum \eqn{\beta} estimated by default (\eqn{0 < \beta < \mu_+ CV_+^2}).}
#'       \item{\code{"CPG-lambda"}}{Compound Poisson Gamma distribution with \eqn{\lambda} estimated by default with a minimum determined by \eqn{CV_+}.}
#'       \item{\code{"delta-X"}{delta distributions where X is one of the first five strictly positive distributions.}
#'     }
#'   }
#'   \item{\code{lambda_model}}{a single character string (applied to all fleets and indices) or a named list of two vectors ("fleet" and "index") of character strings. 
#'     For models that allow zeros, \eqn{p_+ = 1 - \exp{-lambda}} and the general model for lambda for \code{"CPG-lambda"} and \code{"delta-X"} models is 
#'     \deqn{\log\lambda = \beta_0 + \beta_1 \log\widehat P + \beta_2 \log\text{effort}}
#'     where \eqn{\beta_1>0} and \eqn{\beta_2>0} when used so that the probability of 0 decreases with the covariate. This option is specified by a character string that 
#'     may contain nothing (""), "0" (intercept or constant parameter), "P" (for effects of \eqn{\widehat P}, predicted abundance or catch), or "E" (for effort effects). For example, "0PE" 
#'     would indicate to estimate all parameters. When "0", "P", or "E" are not included the corresponding parameter is fixed at 0. These values are ignored for other aggregate observation model options.
#'   }
#'   \item{\code{effort}}{optional named list ("fleet" and "index") of two of matrices (n_years_model x n_fleets and n_years_model x n_indices) defining relative changes in effort 
#'     over time for each fleet and index that may be used to model the probability of observed 0s described above for the \code{"CPG-lambda"} or \code{"delta-X"} models. 
#'     See \code{\lambda_model} below.
#'   }
#'   \item{\code{obs_par_val}}{optional named list of two lists ("fleet" and "index") of vectors of initial/fixed parameters that may be used to configure models defined in 
#'     \code{obs_model}. A single value \eqn{\phi} for \code{"normal-1"}, \code{"normal-2"}, \code{"gamma-1"}, and \code{"gamma-2"} (set to 1 by default which will inflate SD/CV by 
#'     2 for "-2" models). A single value \eqn{\nu} for \code{"gen-gamma"} (set to 1 by default = gamma pdf). For \code{"CPG-alpha"} a single value \eqn{\alpha}. Set = 1 by default. For 
#'     \code{"CPG-beta"} a single value equal to the fraction of maximum \eqn{\beta} (set to 0.5 by default). For \code{"CPG-lambda"} models three values for the \code{beta} parameters in the order 
#'     defined for the  \eqn{\log\lambda} model. For "delta-" models 4 values: first for the CV scalar as described above, and the last three as described for \eqn{\log\lambda}}.
#'   }
#'   \item{\code{obs_par_map}}{Same structure as \code{obs_par_val}. Specifying whether to fix (use NA) or estimate (anything else) aggregate catch model parameters. For estimated parameters,
#'     each will be estimated uniquely regardless of common non-NA values.
#'   }
#'   \item{\code{obs_zero}}{optional named list ("fleet" and "index") of two of matrices (n_years_model x n_fleets and n_years_model x n_indices) of 0/1 integers defining years where
#'     fishing or sampling effort occurred but there was zero catch. Only used when models allowing observed 0s described above are used. E.g., \code{"CPG-X"} or \code{"delta-X"} models.
#'     the CV for positive observations will be assumed unavailable for these observations and the random effect model is used to estimate a CV for those years.
#'   }
#'   \item{\code{cv_model}}{Optional single character value or named list ("fleet" and "index") of two vectors of character values ("fixed-mean", iid" or "ar1"). "fixed-mean" is the default. For models
#'     that allow observed 0s, there will be missing zeros for those observations and an input CV is not available. The "fixed-mean" option will just use the mean CV of the non-zero observations in series. 
#'     The "iid" or "ar1" option will estimate a correspond CV as a random effect. The general model is
#'     \deqn{\log CV_y = \mu(1-\rho) + \rho \log CV_{y-1} + \epsilon_y}
#'     where \eqn{\epsilon_y \sim N(0,\sigma^2)}. For "ar1" \eqn{\rho} is estimated and set \eqn{\rho = 0} for the "iid" model. Specified values only used for \code{"CPG-X"} and 
#'     \code{"delta-X"} models.
#'   }
#'
#' @return a named list with same elements as the input provided with aggregate index options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' input <- set_agg_obs_model(input, age_comp = "gamma-2") #no longer lognormal
#' }
#'
#' @export
#'
set_agg_obs_model <- function(input, agg_obs_model)
{
	data <- input$data
	par <- input$par
  map <- input$map
  all_models <- c(
    "lognormal-1", #1
    "lognormal-2", #2
    "gamma-1", #3
    "gamma-2", #4
    "gengamma", #5
    "CPG-alpha", #6
    "CPG-beta", #7
    "CPG-lambda", #8
    "delta-lognormal-1", #9
    "delta-lognormal-2", #10
    "delta-gamma-1", #11
    "delta-gamma-2", #12
    "delta-gengamma" #13
    )
  if(is.null(data$n_fleets)) stop("the number of fleets has not been specified yet.")
  if(is.null(data$n_indices)) stop("the number of indices has not been specified yet.")

  lnms <- c("fleet", "index")
  n <- c(fleet = data$n_fleets, index = data$n_indices)
  agg_obs <- list(fleet = data$agg_catch, index = data$agg_indices)

  agg_obs_model = list(fleet = rep(1,n[1]), index = rep(1,n[2]))
  input$log$agg_obs_model <- list()
  if(!is.null(agg_obs_model$obs_model)){
    if(is.character(agg_obs_model$obs_model) & length(agg_obs_model$obs_model)==1){ # all use the same model
      themod <- match(agg_obs_model$obs_model, all_models)
      if(is.na(themod)) stop(paste0("agg_obs_model$obs_model option: ", themod, " not recognized. See ?set_agg_obs_model."))
      for(i in lnms) agg_obs_model[[i]] <- rep(themod, n[i])
    } else {
      if(!(is.list(agg_obs_model[["obs_model"]]) & (length(agg_obs_model[["obs_model"]]) == 2) & all(names(agg_obs_model[["obs_model"]]) %in% lnms))) stop("agg_obs_model$obs_model must be either be a character or a named ('fleet','index') list of length 2.")
      for(i in lnms){
        themods <- match(agg_obs_model[["obs_model"]][[i]], all_models)
        if(length(themods) != n[i]) stop(paste0("agg_obs_model$obs_model$",i, " must have length = the number of ", i, "."))
        if(any(is.na(themods))) stop(paste0("agg_obs_model$obs_model$", i, " option not recognized. See ?set_agg_obs_model for available options."))
        agg_obs_model[[i]] <- themods
      }
    }
  }
  data[c("catch_model","index_model")] <- agg_obs_model[lnms]

  input$log$agg_obs_model <- c(input$log$agg_obs_model, paste0("Fleet aggregate observation models are:\n",
    paste0("Fleet ", 1:data$n_fleets, ": ", all_models[data$catch_model], collapse = "\n"), "\n\n")
  )
  input$log$agg_obs_model <- c(input$log$agg_obs_model, paste0("Index aggregate observation models are:\n",
    paste0("Index ", 1:data$n_indices, ": ", all_models[data$index_model], collapse = "\n"), "\n\n")
  )

  # what is lambda (if used) a function of?
  lambda_mods <- list(fleet = rep("0", n[1]), index = rep("0", n[2]))
  lambda_opts <- c("0","P", "E")
  if(!is.null(agg_obs_model[["lambda_model"]])){
    if(is.character(agg_obs_model[["lambda_model"]]) & length(agg_obs_model[["lambda_model"]])==1){ # all use the same model
      if(!any(sapply(lambda_opts, grepl, x=agg_obs_model[["lambda_model"]]))) stop("agg_obs_model$lambda_model must include characters '0', 'P', or 'E'")
      lambda_mods <- list(fleet = rep(agg_obs_model[["lambda_model"]], n[1]), index = rep(agg_obs_model[["lambda_model"]], n[2]))
    } else lambda_mods <- agg_obs_model[["lambda_model"]]
    if(!(is.list(lambda_mods) & (length(lambda_mods) == 2) & all(names(lambda_mods) %in% lnms))) stop("agg_obs_model$lambda_model must be a list of length 2 with names 'fleet' and 'index'.")
    for(i in lnms) {
      if(length(lambda_mods[[i]]) != n[i]) stop(paste0("agg_obs_model$lambda_model$", i, "must have length equal to the number of fleets."))
      for(j in 1:n[i]){
        if(agg_obs_model[[i]][j] %in% 8:13) { #if specified models doesn't use lambda, then don't worry about it.
          if(!any(sapply(lambda_opts, grepl, x=lambda_mods[[i]][j]))) stop("agg_obs_model$lambda_model$", i, "[", j, "] must include characters '0', 'P', or 'E'")
        }
      }
    }
  }

    #columns of catch_pars or index_pars
    #  1: lognormal-1:        log(phi),NA,NA,NA
    #  2: lognormal-2:        log(phi),NA,NA,NA
    #  3: gamma-1:            log(phi),NA,NA,NA
    #  4: gamma-2:            log(phi),NA,NA,NA
    #  5: gengamma:           log(nu),NA,NA,NA
    #  6: CPG-alpha           logit(alpha),NA,NA,NA
    #  7: CPG-beta            logit(frac_beta),NA,NA,NA
    #  8: CPG-lambda          lambda_0,log(lambda_1),log(lambda_2),NA
    #  9: delta-lognormal-1   log(phi),lambda_0,log(lambda_1),log(lambda_2)
    # 10: delta-lognormal-2   log(phi),lambda_0,log(lambda_1),log(lambda_2)
    # 11: delta-gamma-1       log(phi),lambda_0,log(lambda_1),log(lambda_2)
    # 12: delta-gamma-2       log(phi),lambda_0,log(lambda_1),log(lambda_2)
    # 13: delta-gengamma      log(nu),lambda_0,log(lambda_1),log(lambda_2)

    #annual nll pars for each fleet and index or index_pars
    #  1: lognormal-1:        mu, sigma, NA, NA
    #  2: lognormal-2:        mu, sigma, NA, NA
    #  3: gamma-1:            alpha, beta, NA, NA
    #  4: gamma-2:            alpha, beta, NA, NA
    #  5: gengamma:           mu, sigma, nu, NA,
    #  6: CPG-alpha           mu, phi, power, NA
    #  7: CPG-beta            mu, phi, power, NA
    #  8: CPG-lambda          mu, phi, power, NA
    #  9: delta-lognormal-1   mu, sigma, NA, p(+)
    # 10: delta-lognormal-2   mu, sigma, NA, p(+)
    # 11: delta-gamma-1       alpha, beta, NA, p(+)
    # 12: delta-gamma-2       alpha, beta, NA, p(+)
    # 13: delta-gengamma      mu, sigma, nu, p(+)

  agg_par_vals <- list(fleet = matrix(0,n[1], 4), index = matrix(0,n[2], 4))
  
  #first, fix any unused lambda parameters to 0 aand then only use obs_par_val if lambda_model says to use them.
  for(i in lnms) for(j in 1:n[i]) {
    if(agg_obs_model[[i]][j] == 6) agg_par_vals[[i]][j,1] <- log(0.01/(0.99)) #Good to start at a very low alpha to allow lambda to be determined in internal Newton methods
    if(agg_obs_model[[i]][j] == 7) agg_par_vals[[i]][j,1] <- log(0.99/(0.01)) #Good to start at a very large beta fraction to allow lambda to be determined in internal Newton methods
    if(agg_obs_model[[i]][j] %in% 8:13){
      #which not specified to be estimated? Only beta_0 by default
      ind <- 1+which(!sapply(lambda_opts, grepl, x=lambda_mods[[i]][j])[2:3]) #TRUE if not used. 1 + integer(0) = integer(0)      
      #2nd and 3rd pars for CPG-lambda: model 8
      #3rd and 4th pars for delta-x: models 9:13
      if(agg_obs_model[[i]][j] %in% 9:13) ind <- ind + 1
      agg_par_vals[[i]][j,ind] <- -1e10  #exp(val) ~ 0, essentially no effect when not estimated
    }
  }

  #specify initial/fixed values
  if(!is.null(agg_obs_model[["obs_par_val"]])) {
    if(!(is.list(agg_obs_model[["obs_par_val"]]) & (length(agg_obs_model[["obs_par_val"]]) == 2) & all(names(agg_obs_model[["obs_par_val"]]) %in% lnms))) stop("agg_obs_model$obs_par_val must be a list of length 2 with names 'fleet' and 'index'.")
    for(i in lnms){
      if(!(is.list(agg_obs_model[["obs_par_val"]][[i]]) & (length(agg_obs_model[["obs_par_val"]][[i]]) != n[i]))) stop(paste0("agg_obs_model$obs_par_val$",i, " must be a list of length = number of ", i, "."))
      for(j in 1:n[i]) {
        par_vals <- agg_obs_model[["obs_par_val"]][[i]][[j]]
        if(agg_obs_model[[i]][j] %in% 1:7){
          if(length(par_vals) != 1) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]] must be a single value for model ", all_mods[agg_obs_model[[i]][j]],"."))
          if(agg_obs_model[[i]][j] < 7) {
            if(0>par_vals[1]) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]] must be a single value > 0 for model ", all_mods[agg_obs_model[[i]][j]],"."))
            agg_par_vals[[i]][j,1] <- log(par_vals[1])
          }
          if(agg_obs_model[[i]][j] == 7) {
            if(0>par_vals[1] | par_vals[1]>1) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]] must be a single value between 0 and 1 for model ", all_mods[agg_obs_model[[i]][j]],"."))
            agg_par_vals[[i]][j,1] <- log(par_vals[1]) - log(1-par_vals[1])
          }
        }
        if(agg_obs_model[[i]][j] == 8){
          if(length(par_vals) != 3) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]] must be 3 values for model ", all_mods[agg_obs_model[[i]][j]],"."))
          ind <- sapply(lambda_opts, grepl, x=lambda_mods[[i]][j])
          if(ind[1]) agg_par_vals[[i]][j,1] <- par_vals[1]
          if(ind[2]) {
            if(par_vals[2]<0) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]][2] defining effect of predicted catch or index must be greater than 0."))
            agg_par_vals[[i]][j,2] <- log(par_vals[2])
          }
          if(ind[3]) {
            if(par_vals[3]<0) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]][3] defining effect of predicted catch or index must be greater than 0."))
            agg_par_vals[[i]][j,3] <- log(par_vals[3])
          }
        }
        if(agg_obs_model[[i]][j] %in% 9:13){
          if(length(par_vals) != 4) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]] must be 4 values for model ", all_mods[agg_obs_model[[i]][j]],"."))
          if(0>par_vals[1]) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]][1] must be a single value > 0 for model ", all_mods[agg_obs_model[[i]][j]],"."))
          agg_par_vals[[i]][j,1] <- log(par_vals[1])
          ind <- sapply(lambda_opts, grepl, x=lambda_mods[[i]][j])
          if(ind[1]) agg_par_vals[[i]][j,2] <- par_vals[2]
          if(ind[2]) {
            if(par_vals[3]<0) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]][3] defining effect of predicted catch or index must be greater than 0."))
            agg_par_vals[[i]][j,3] <- log(par_vals[3])
          } 
          if(ind[3]) {
            if(par_vals[4]<0) stop(paste0("agg_obs_model$obs_par_val$",i, "[[", j, "]][4] defining effect of predicted catch or index must be greater than 0."))
            agg_par_vals[[i]][j,4] <- log(par_vals[4])
          } 
        }
      }
    }
  }
  par[c("catch_pars","index_pars")] <- agg_par_vals[lnms]
  
  #fix or estimate parameters
  agg_par_map <- list(fleet = matrix(NA,n[1], 4), index = matrix(NA,n[2], 4))
  for(i in lnms) {
    agg_par_map[[i]][which(agg_obs_model[[i]] %in% c(2,4:8, 8 + c(2,4:8))),1] <- 1 #estimate the first parameter for these models by default
    agg_par_map[[i]][which(agg_obs_model[[i]] %in% 9:13),2] <- 1 #estimate the the first parameter for probability of 0 for these models by default
  }
  #first, fix any unused lambda parameters to 0 aand then only use obs_par_val if lambda_model says to use them.
  for(i in lnms) for(j in 1:n[i]) if(agg_obs_model[[i]][j] %in% 8:13){
    #which not specified to be estimated? Only beta_0 by default
    ind <- 1+which(!sapply(lambda_opts, grepl, x=lambda_mods[[i]][j])[2:3]) #TRUE if not used. 1 + integer(0) = integer(0)
    #2nd and 3rd pars for CPG-lambda: model 8
    #3rd and 4th pars for delta-x: models 9:13
    if(agg_obs_model[[i]][j] %in% 9:13) ind <- ind + 1
    agg_par_map[[i]][j,ind] <- NA #turn off even if agg_obs_model$obs_par_map incorrectly says estimate
  }

  if(!is.null(agg_obs_model[["obs_par_map"]])) {
    if(!(is.list(agg_obs_model[["obs_par_map"]]) & (length(agg_obs_model[["obs_par_map"]]) == 2) & all(names(agg_obs_model[["obs_par_map"]]) %in% lnms))) stop("agg_obs_model$obs_par_map must be a list of length 2 with names 'fleet' and 'index'.")
    for(i in lnms){
      if(!(is.list(agg_obs_model[["obs_par_map"]][[i]]) & (length(agg_obs_model[["obs_par_map"]][[i]]) != n[i]))) stop(paste0("agg_obs_model$obs_par_map$",i, " must be a list of length = number of ", i, "."))
      for(j in 1:n[i]) {
        par_map <- agg_obs_model[["obs_par_map"]][[i]][[j]]
        if(agg_obs_model[[i]][j] %in% 1:7){
          if(length(par_map) != 1) stop(paste0("agg_obs_model$obs_par_map$",i, "[[", j, "]] must be a single value for model ", all_mods[agg_obs_model[[i]][j]],"."))
          agg_par_map[[i]][j,1] <- 1
          if(is.na(par_map)) agg_par_map[[i]][j,1] <- NA
        }
        if(agg_obs_model[[i]][j] == 8){
          if(length(par_map) != 3) stop(paste0("agg_obs_model$obs_par_map$",i, "[[", j, "]] must be 3 values for model ", all_mods[agg_obs_model[[i]][j]],"."))
          agg_par_map[[i]][j,1:3] <- 1
          agg_par_map[[i]][j,which(is.na(par_map))] <- NA
        }
        if(agg_obs_model[[i]][j] %in% 9:13){
          if(length(par_map) != 4) stop(paste0("agg_obs_model$obs_par_map$",i, "[[", j, "]] must be 4 values for model ", all_mods[agg_obs_model[[i]][j]],"."))
          agg_par_map[[i]][j,1:4] <- 1
          agg_par_map[[i]][j,which(is.na(par_map))] <- NA
        }
      }
    }
  }

  #map can now safely be defined
  for(i in lnms) agg_par_map[[i]][which(!is.na(agg_par_map[[i]]))] <- 1:sum(!is.na(agg_par_map[[i]]))
  map[c("catch_pars","index_pars")] <- lapply(agg_par_map[lnms],factor)

     
  #specify effort variables that may be needed for P(0) when modeling 0 observations
  effort <- list(fleet = matrix(1,data[["n_years_model"]],n[1]), index = matrix(1,data[["n_years_model"]],n[2]))
  if(any(unlist(agg_obs_model) %in% 8:13)) if(!is.null(agg_obs_model[["effort"]])) { #if specified model doesn't use lambda, then don't worry about it.
    if(!(is.list(agg_obs_model[["effort"]]) & (length(agg_obs_model[["effort"]]) == 2) & all(names(agg_obs_model[["effort"]]) %in% lnms))) stop("agg_obs_model$effort must be a list of length 2 with names 'fleet' and 'index'.")
    for(i in lnms){
      if(!(is.matrix(agg_obs_model[["effort"]][[i]]) & is.numeric(agg_obs_model[["effort"]][[i]]) & all(dim(agg_obs_model[["effort"]][[i]]) == c(data[["n_years_model"]],n[i])))) stop(paste0("agg_obs_model$effort$",i, 
        " must be a numeric matrix with dimensions: data$n_years_model x number of ", i, "."))
      effort[[i]] <- agg_obs_model[["effort"]][[i]]
    }
  }
  data[c("catch_effort","index_effort")] <- effort[lnms]

  obs_pos <- list(fleet = matrix(1L,data[["n_years_model"]],n[1]), index = matrix(1L,data[["n_years_model"]],n[2]))
  if(any(unlist(agg_obs_model) %in% 6:13)) {
    if(!is.null(agg_obs_model[["obs_zero"]])) { #if specified model doesn't allow 0s, then don't worry about it.
      if(!(is.list(agg_obs_model[["obs_zero"]]) & (length(agg_obs_model[["obs_zero"]]) == 2) & all(names(agg_obs_model[["obs_zero"]]) %in% lnms))) stop("agg_obs_model$obs_zero must be a list of length 2 with names 'fleet' and 'index'.")
      for(i in lnms){
        if(!(is.matrix(agg_obs_model[["obs_zero"]][[i]]) & is.integer(agg_obs_model[["obs_zero"]][[i]]) & all(dim(agg_obs_model[["obs_zero"]][[i]]) == c(data[["n_years_model"]],n[i])))) stop(paste0("agg_obs_model$obs_zero$",i, 
        " must be an integer (0/1) matrix with dimensions: data$n_years_model x number of ", i, "."))
        obs_pos[[i]][which(agg_obs_model[["obs_zero"]]==1)] <- 0L
      }
    } else{
      obs_pos[["fleet"]][which(data[["agg_catch"]]< 1e-10 & data[["use_agg_catch"]] == 1)] <- 0L
      obs_pos[["index"]][which(data[["agg_indices"]]< 1e-10 & data[["use_indices"]] == 1)] <- 0L
    }
  }
  data[c("catch_pos","index_pos")] <- obs_pos[lnms]
  obs_cv <- list(fleet = data[["agg_catch_cv"]], index = data[["agg_index_cv"]])
  cv_repars_val <- list(fleet = matrix(0,n[1],3), index = matrix(0,n[2],3))
  cv_repars_map <- list(fleet = matrix(NA,n[1],3), index = matrix(NA,n[2],3))
  n_cv_re <- c(fleet = 0, index = 0)
  cv_re <- list(fleet = numeric(), index = numeric())
  map_cv_re <- list(fleet = integer(), index = integer())
  use_cv_re <- list(fleet = rep(0L,n[1]), index = rep(0L,n[2]))
  #if no missing CVs, then don't do anything.
  #default for missing CVs, is just you the mean CV+ for the 0 years
  for(i in lnms) for(j in 1:n[i]) {
    n_ij <- sum(obs_pos[[i]][,j] == 0L)
    if((n_ij>0) & (agg_obs_model[[i]][j] %in% 6:13)) {
      cv_re[[i]][n_cv_re[i] + 1:n_ij] <- log(mean(obs_cv[[i]][,j][which(obs_pos[[i]][,j]==1)]))
      map_cv_re[[i]] <- c(map_cv_re[[i]], rep(NA, n_ij))
      n_cv_re[i] <- n_cv_re[i] + n_ij
      if(!is.null(agg_obs_model[["cv_model"]])){
        if(is.character(agg_obs_model[["cv_model"]]) && (length(agg_obs_model[["cv_model"]])==1)){
          if(agg_obs_model[["cv_model"]] != "fixed-mean"){
            use_cv_re[[i]][j] <- 1L
            map_cv_re[[i]][(n_cv_re[i]-n_ij):n_cv_re[i]] <- 1L
            cv_repars_map[[i]][j,1:2] <- 1L #default is iid: estimate mean and variances when
          }
          #only have to turn on correlation parameter if necessary
          if(agg_obs_model[["cv_model"]] == "ar1") cv_repars_map[[i]][j,3] <- 1L 
        } else{
          if(!(is.list(agg_obs_model[["cv_model"]]) && (length(agg_obs_model[["cv_model"]]) == 2) && all(names(agg_obs_model[["cv_model"]]) %in% lnms))) stop("agg_obs_model$cv_model must be either be a character or a named ('fleet','index') list of length 2.")
          if(agg_obs_model[["cv_model"]][[i]][j] != "fixed-mean") {
            use_cv_re[[i]][j] <- 1L
            map_cv_re[[i]][(n_cv_re[i]-n_ij):n_cv_re[i]] <- 1L
            cv_repars_map[[i]][j,1:2] <- 1L #default is iid: estimate mean and variances when
          }
          #only have to turn on correlation parameter if necessary
          if(agg_obs_model[["cv_model"]][[i]][j] == "ar1") cv_repars_map[[i]][j,3] <- 1L 
          n_cv_re[i] <- n_cv_re[i] + sum(obs_pos[[i]][,j] == 0L)
        }
      }
    }
  }
  for(i in lnms) map_cv_re[[i]][which(!is.na(map_cv_re[[i]]))] <- 1:sum(!is.na(map_cv_re[[i]]))
  for(i in lnms) cv_repars_map[[i]][which(!is.na(cv_repars_map[[i]]))] <- 1:sum(!is.na(cv_repars_map[[i]]))
  data[["do_simulate_cv_re"]] <- 0 #do not simulate CV random effects
  data[c("catch_cv_use_re", "index_cv_use_re")] <- use_cv_re[lnms]
  par[c("catch_cv_re", "index_cv_re")] <- cv_re[lnms]
  par[c("catch_cv_repars","index_cv_repars")] <- cv_repars_val[lnms]
  map[c("catch_cv_re","index_cv_re")] <- lapply(map_cv_re[lnms], factor)
  map[c("catch_cv_repars","index_cv_repars")] <- lapply(cv_repars_map[lnms], factor)
  input[["random"]] <- c(input[["random"]], c("catch_cv_re", "index_cv_re")[which(sapply(use_cv_re[lnms], sum)>0)])


	input$data <- data
	input$par <- par
  input$map <- map
  if(length(input$log$agg_obs_model))  input$log$agg_obs_model <- c(
    "--Aggregate Fleet and Index Observation Models---------------------------------------------------------------------------------------",
    "\n", 
    input$log$agg_obs_model,
    "-------------------------------------------------------------------------------------------------------------------------------------",
  "\n\n")
  input$options$agg_obs_model <- agg_obs_model
  
  if(is.null(input$by_pwi)) { #check whether called by prepare_wham_input
    input <- set_osa_obs(input)
    message(unlist(input$log$agg_obs_model, recursive=T))
  }
 
	return(input)
}

