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
  na = data$n_ages
  if(npeels)
  {
    rho  = list()
    rho$SSB <- rho$Fbar <- numeric()
    for(i in 1:data$n_stocks) rho$SSB[i] <- mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[ny-x, i]/model$rep$SSB[ny-x, i] - 1))
    for(i in 1:data$n_regions) rho$Fbar[i] <- mean(sapply(1:npeels, function(x) {
      mean(model$peels[[x]]$rep$Fbar[ny-x,i])/mean(model$rep$Fbar[ny-x,i]) - 1
    }))
    #for(i in 1:data$n_regions) rho$Fbar[i] <- mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$Fbar[ny-x,i]/model$rep$Fbar[ny-x,i] - 1))
    #names(rho) = c("SSB","Fbar")#,"R")
    rho$naa <- array(NA, c(data$n_stocks, data$n_regions, data$n_ages))
    for(s in 1:data$n_stocks) for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]) {
      rho$naa[s,r,a] = mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[s,r,ny-x,a]/model$rep$NAA[s,r,ny-x,a] - 1))
    }
    dimnames(rho$naa)[[3]] = c("R", paste0("N", model$ages.lab[2:na]))
    #rho = c(rho, rho.naa)
    return(rho)
  }
  else stop("There are no peels in this model")
}
