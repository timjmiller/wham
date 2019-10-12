#' Prepare input data and parameters to project an already fit WHAM model
#'
#' \code{prepare_projection} is an internal function called by \code{\link{project_wham}},
#' which in turn is called by \code{\link{fit_wham}} if \code{do.proj = TRUE}.
#'
#' @param model a previously fit wham model
#' @param proj.opts a named list with the following components:
#'   \describe{
#'     \item{\code{$n.yrs}}{integer, number of years to project/forecast. Default = \code{3}.}
#'     \item{\code{$use.lastF}}{T/F, use terminal year F for projections. Default = \code{TRUE}.}
#'     \item{\code{$use.avgF}}{T/F, use average F for projections.}
#'     \item{\code{$use.FXSPR}}{T/F, calculate F at X% SPR for projections.}
#'     \item{\code{$proj.F}}{vector, user-specified fishing mortality for projections. Length must equal \code{n.yrs}.}
#'     \item{\code{$proj.catch}}{vector, user-specified aggregate catch for projections. Length must equal \code{n.yrs}.}
#'     \item{\code{$avg.yrs}}{vector, specify which years to average over for calculating reference points. Default = last 5 model years, \code{tail(model$years, 5)}.}
#'   }
#'
#' @return same as \code{prepare_wham_input}, a list ready for \code{fit_wham}:
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
#' @seealso \code{\link{prepare_wham_input}}, \code{\link{project_wham}}
#'
prepare_projection = function(model, proj.opts)
{
write.dir <- "/home/bstock/Documents/wham/sandbox/ex1"
setwd(write.dir)
load("ex1_models_project.RData")
library(wham)
model=mods$m4
proj.opts=list(n.yrs=3, use.lastF=TRUE, use.FXSPR=FALSE, proj.F=NULL, proj.catch=NULL, avg.yrs=NULL)

  input1 <- model$input
  avg.yrs.ind <- match(proj.opts$avg.yrs, input1$years)
  avg_cols = function(x) apply(x, 2, mean, na.rm=TRUE)

  # modify data
  data <- input1$data; par <- input1$par; map <- input1$map; random <- input1$random;
  data$do_proj = 1
  data$n_years_proj = proj.opts$n.yrs
  data$mature <- rbind(data$mature, avg_cols(data$mature[avg.yrs.ind,]))



  return(list(data=data, par = par, map = map, random = random, years = model_years,
    ages.lab = input1$ages.lab, model_name = input1$model_name))
}
