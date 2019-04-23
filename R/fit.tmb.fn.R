#' Fit TMB model using nlminb
#'
#' Runs optimization on the TMB model using \code{\link[stats:nlminb]{stats::nlminb}}.
#' If specified, takes additional Newton steps and calculates standard deviations.
#' Internal function called by \code{\link{fit.wham.fn}}.
#'
#' @param model Output from \code{\link[TMB:MakeADFun]{TMB::MakeADFun}}.
#' @param n.newton Integer, number of additional Newton steps after optimization. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#'
#' @return \code{model}, appends the following:
#'   \describe{
#'     \item{\code{model$opt}}{Output from \code{\link[stats:nlminb]{stats::nlminb}}}
#'     \item{\code{model$date}}{System date}
#'     \item{\code{model$dir}}{Current working directory}
#'     \item{\code{model$rep}}{model$report()}
#'     \item{\code{model$TMB_version}}{Version of TMB installed}
#'     \item{\code{model$parList}}{List of parameters, \code{model$env$parList()}}
#'     \item{\code{model$final_gradient}}{Final gradient, \code{model$gr()}}
#'     \item{\code{model$sdrep}}{Estimated standard deviations for model parameters, \code{\link[TMB:sdreport]{TMB::sdreport}}}
#'   }
#'
#' @seealso \code{\link{fit.wham.fn}}, \code{\link{retro.fn}}
#'
fit.tmb.fn = function(model, n.newton=3, do.sdrep = TRUE)
{
  model$opt <- stats::nlminb(model$par, model$fn, model$gr, control = list(iter.max = 1000, eval.max = 1000))
  if(n.newton) for(i in 1:n.newton) { # Take a few extra newton steps
    g <- as.numeric(model$gr(model$opt$par))
    h <- stats::optimHess(model$opt$par, model$fn, model$gr)
    model$opt$par <- model$opt$par - solve(h, g)
    model$opt$objective <- model$fn(model$opt$par)
  }
  model$date = Sys.time()
  model$dir = getwd()
  model$rep <- model$report()
  model$TMB_version = packageVersion("TMB")
  model$parList = model$env$parList()
  model$final_gradient = model$gr()
  if(do.sdrep)
  {
    model$sdrep <- try(TMB::sdreport(model))
    model$is_sdrep = !is.character(model$sdrep)
    if(model$is_sdrep) model$na_sdrep = any(is.na(summary(model$sdrep)[,2]))
    else model$na_sdrep = NA
  }
  return(model)
}
