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
  data <- model$env$data
  ny = data$n_years_model
  nr = data$n_regions
  ns = data$n_stocks
  nf = data$n_fleets
  na = data$n_ages
  if(npeels)
  {
    rho  = list()
    rho$SSB <- rho$Fbar <- numeric()
    for(i in 1:ns) rho$SSB[i] <- mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[ny-x, i]/model$rep$SSB[ny-x, i] - 1))
    for(i in 1:nr) rho$Fbar[i] <- mean(sapply(1:npeels, function(x) {
      mean(model$peels[[x]]$rep$Fbar[ny-x,nf + i])/mean(model$rep$Fbar[ny-x,nf + i]) - 1
    }))
    rho$naa <- array(NA, c(ns, nf, na))
    for(s in 1:ns) for(r in 1:nr) for(a in 1:na) if(data$NAA_where[s,r,a]) {
      rho$naa[s,r,a] = mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[s,r,ny-x,a]/model$rep$NAA[s,r,ny-x,a] - 1))
    }
    dimnames(rho$naa)[[3]] = c("R", paste0("N", model$ages.lab[2:na]))
    return(rho)
  }
  else stop("There are no peels in this model")
}
