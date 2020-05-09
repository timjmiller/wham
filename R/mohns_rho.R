#' Calculate Mohn's rho for a WHAM model with peels
#'
#' @param model A fit WHAM model, output from \code{\link{fit_wham}} with \code{do.retro = TRUE}.
#'
#' @return \code{rho}, a vector of Mohn's rho
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{retro}}
#'
#' @examples
#' \dontrun{
#' data("input4_SNEMAYT") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input4_SNEMAYT) # using default values: do.retro = T, n.peels = 7
#' mohns_rho(mod) # calculate Mohn's rho
#' }
mohns_rho = function(model)
{
  npeels = length(model$peels)
  ny = model$env$data$n_years_model
  na = model$env$data$n_ages
  if(npeels)
  {
    rho = c(
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[ny-x]/model$rep$SSB[ny-x] - 1)),
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$Fbar[ny-x]/model$rep$Fbar[ny-x] - 1)))#,
    #mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[ny-x,1]/model$rep$NAA[ny-x,1] - 1)))
    names(rho) = c("SSB","Fbar")#,"R")
    rho.naa = sapply(1:na, function(y)
    {
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[ny-x,y]/model$rep$NAA[ny-x,y] - 1))
    })
    names(rho.naa) = c("R", paste0("N", model$ages.lab[2:na]))
    rho = c(rho, rho.naa)
    return(rho)
  }
  else stop("There are no peels in this model")
}
