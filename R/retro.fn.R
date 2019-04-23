#' Run retrospective analysis
#'
#' Internal function called by \code{\link{fit.wham.fn}}. Calls \code{\link{peel.fit.fn}}
#' to fit the model peeling off \code{1, 2, ..., n.peels} years of data.
#'
#' @param model Optimized TMB model, output from \code{\link{fit.tmb.fn}}.
#' @param n.peels Integer, number of peels to use in retrospective analysis. Default = \code{7}.
#' @param ran Character, specifies which parameters to treat as random effects. Default = \code{"log_NAA"}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters for each peel? Default = \code{FALSE}.
#' @param n.newton integer, number of additional Newton steps after optimization for each peel. Default = \code{0}.
#'
#' @return \code{peels}, a list of length \code{n.peels}, where entry \emph{i} is a model
#' fit by peeling off \emph{i} years of data.
#'
#' @seealso \code{\link{fit.wham.fn}}, \code{\link{peel.fit.fn}}
#'
retro.fn = function(model, n.peels = 7, ran = "log_NAA", do.sdrep = FALSE, n.newton = 0)
{
  temp = list(dat = model$env$data, par = model$parList, map = model$env$map, random = ran)
  if(n.peels>0) peels = list(peel.fit.fn(1, model = temp, do.sdrep = do.sdrep, n.newton = n.newton))
  if(n.peels>1) for(i in 2:n.peels) peels[[i]] = peel.fit.fn(i, model = temp, do.sdrep = do.sdrep, n.newton = n.newton)
  return(peels)
}
