#' Fit WHAM model
#'
#' Fits the compiled WHAM model using \code{\link[TMB:MakeADFun]{TMB::MakeADFun}} and
#' \code{\link[stats:nlminb]{stats::nlminb}}. Runs retrospective analysis if specified.
#'
#' @param input Named list with several components:
#'   \describe{
#'     \item{\code{input$data}}{Data to fit the assessment model to.}
#'     \item{\code{input$par}}{Parameters, a list of all parameter objects required by the user template (both random and fixed effects). See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{input$map}}{Map, a mechanism for collecting and fixing parameters. See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{input$random}}{Character vector defining the parameters to treat as random effect. See \code{\link[TMB]{MakeADFun}}.}
#'     \item{\code{input$years}}{Numeric vector of the years which the model spans. Not important for model fitting, but useful for plotting.}
#'     \item{\code{input$map.n}}{I don't know what this is.}
#'     \item{\code{input$model_name}}{Character, name of the model, e.g. \code{"Yellowtail flounder"}}
#'     \item{\code{input$ages.lab}}{Character vector of the age labels, e.g. \code{c("1","2","3","4+").}}
#'   }
#' @param n.newton integer, number of additional Newton steps after optimization. Passed to \code{\link{fit_tmb}}. Default = \code{3}.
#' @param do.sdrep T/F, calculate standard deviations of model parameters? See \code{\link[TMB]{sdreport}}. Default = \code{TRUE}.
#' @param do.retro T/F, do retrospective analysis? Default = \code{TRUE}.
#' @param n.peels integer, number of peels to use in retrospective analysis. Default = \code{7}.

#' @return \code{mod}, a fit TMB model
#'
#' @useDynLib wham
#' @export
#'
#' @seealso \code{\link{fit_tmb}}, \code{\link{retro}}
#'
#' @examples
#' \dontrun{
#' data("SNEMA_ytl") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input) # using default values
#' }
fit_wham = function(input, n.newton = 3, do.sdrep = TRUE, do.retro = TRUE, n.peels = 7)
{
  # wham.dir <- find.package("wham")
  # dyn.load( paste0(wham.dir,"/libs/", TMB::dynlib(version)) )
  mod <- TMB::MakeADFun(input$data,input$par, DLL = "wham", random = input$random, map = input$map)
  mod <- fit_tmb(mod, n.newton = n.newton, do.sdrep = do.sdrep)
  if(do.retro) mod$peels = retro(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels)
  mod$years = input$years
  mod$ages.lab = input$ages.lab
  mod$model_name = input$model_name
  return(mod)
}
