#' Extract retrospective results for plotting
#'
#' @param model A fit WHAM model, output from \code{\link{fit_wham}} with \code{do.retro = TRUE}.
#'
#' @return a named list with the components:
#'   \describe{
#'     \item{\code{SSB}}{Spawning stock biomass}
#'     \item{\code{Fbar}}{Fishing mortality}
#'     \item{\code{NAA}}{Numbers-at-age}
#'  }
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}
#'
#' @examples
#' \dontrun{
#' data("SNEMA_ytl") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input) # using default values: do.retro = T, n.peels = 7
#' x = retro_res(mod) # get retrospective results
#' }
retro_res = function(model) #get time series for retro plots
{
  npeels = length(model$peels)
  ny = model$env$data$n_years_model
  if(npeels)
  {
    SSB = lapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[1:(ny-x)]/model$rep$SSB[1:(ny-x)] - 1)
    Fbar = lapply(1:npeels, function(x) model$peels[[x]]$rep$Fbar[1:(ny-x)]/model$rep$Fbar[1:(ny-x)] - 1)
    na = model$env$data$n_ages
    NAA = lapply(1:na, function(a)
    {
      lapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[1:(ny-x),a]/model$rep$NAA[1:(ny-x),a] - 1)
    })
    res = list(SSB = SSB, Fbar = Fbar, NAA = NAA)
    return(res)
  }
  else stop("There are no peels in this model")
}
