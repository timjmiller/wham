#' Calculate Akaike's Information Criterion for a WHAM model
#'
#' Calculates AIC from a fitted WHAM model's negative log-likelihood
#' and number of estimated parameters.
#'
#' @param mod A fitted WHAM model object returned by \code{fit_wham()}.
#' @param conditional (TRUE/FALSE) When the model includes random effects, the default (\code{conditional = FALSE})
#' calculation uses the marginal likelihood and number of fixed effects parameters. If \code{conditional = TRUE},
#' the joint likelihood of the data conditional on the estimated random effects is used with an estimated effective degress of freedom
#' that is calculated using the approach descsribed by \href{https://doi.org/10.48550/arXiv.2411.14185}{Zhang et al. 2024} and code provided by Noel Cadigan.
#'
#' @return A numeric AIC value with attributes denoting the degrees of freedom, number of observations, and the type (marginal or conditional).
#' @export
#'
#' @examples
#' \dontrun{
#' mod <- fit_wham(input)
#' aic(mod)
#' aic(mod, conditional = TRUE)
#' }
aic <- function(mod, conditional = FALSE) {
  if(is.null(mod$TMB_version)) stop("mod$TMB_version does not exist. This must be a TMB model.")
  if(conditional) return(cAIC(mod))
  else{
    n.obs <- sum(!is.na(mod$env$data$obs$val))
    n.fe <- length(mod$par)
    if(is.null(mod[["opt"]][["objective"]])) {
      warning("The model appears to be unoptimized. The calculated AIC uses the likelihood at the current parameter values.")
      if(conditional) return(cAIC(mod))
      else{
        fe <- mod[["env"]][["last.par.best"]][mod[["env"]][["lfixed"]]()]
        AIC <- 2 * (mod[["fn"]](fe) + n.parm)
      }
    } else AIC <- 2 * (mod[["opt"]][["objective"]] + n.fe)
    AIC <- structure(AIC, df = n.fe, n = n.obs, type = "marginal")
    return(AIC)
  }
}


cAIC <- function(mod){

  # From Noel Cadigan using results from
  # Zhang, N., Cadigan, N. G., Thorson, J. T. 2024. A note on numerical evaluation of conditional Akaike information for nonlinear mixed-effects models.  
  #   arXiv:2411.14185 [stat.ME]. doi: 10.48550/arXiv.2411.14185
  
  if(is.null(mod[["opt"]][["objective"]])) {
    warning("The model appears to be unoptimized. The calculated AIC uses the likelihood at the current parameter values.")
  }
  n.parm <- length(mod$par)
  n.obs <- sum(!is.na(mod$env$data$obs$val))

  parDataMode <- mod$env$last.par.best
  nodata <- mod$env$data
  nodata$use_agg_catch[] <- 0 
  nodata$use_catch_paa[] <- 0 
  nodata$use_indices[] <- 0   
  nodata$use_index_paa[] <- 0
  nodata$Ecov_use_obs[] <- 0
  
  modND <- TMB::MakeADFun(nodata, mod$env$parList(par = mod$env$last.par.best), DLL = "wham", 
            random = mod$env$.random, map = mod$env$map, silent = TRUE)

  ## Marginal Precision matrix of random effects

  indx <- modND$env$lrandom() 
  ## use - for Hess because model returns negative loglikelihood
  Hess <- -Matrix::Matrix(modND$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE)
  cov_Psi_inv <- -Hess[indx,indx]; ## this is the marginal prec mat of REs;
  q <- nrow(cov_Psi_inv)             
  ddlr.r <- Hess[indx,indx]  
  
  ## Joint hessian etc  
  ## use - for Hess because model returns negative loglikelihood
  hess.m <- -Matrix::Matrix(mod$env$f(parDataMode,order=1,type="ADGrad"),sparse = TRUE) 
  ddlj.r <- hess.m[indx,indx]
  ddlj.r_inv <- Matrix::solve(ddlj.r)
  
  jnll <- mod$env$f(parDataMode) 
  cnll <- jnll - modND$env$f(parDataMode)
  
  edf <- (n.parm+q) - sum(diag(as.matrix(ddlj.r_inv %*% ddlr.r)))
  cAIC <- 2*(cnll + edf)
  #cAICc <- cAIC + 2*edf*(edf+1)/(n.obs - edf-1)
  AIC <- structure(cAIC, df = edf, n = n.obs, type = "conditional")
  return(AIC)
}
