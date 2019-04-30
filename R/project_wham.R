#' Project WHAM model
#'
#' After fitting a WHAM model, \code{project_wham} projects the model forward.
#'
#' @param mod list, output from \code{\link{fit_wham}}
#' @param n.years integer, number of years to project (default = 3)
#'
#' @return the input list \code{mod}, with an additional element \code{proj}.
#' \code{mod$proj} is a fit TMB model, same as output by \code{\link{fit_wham}}, except
#' with additional elements \code{mod$proj$proj_IAA} (projected index of abundance at-age, mean)
#' and \code{mod$proj$proj_IAA$sd} (projected index of abundance at-age, sd).
#'
#' @seealso \code{\link{fit_wham}}
#'
#' @examples
#' \dontrun{
#' m1 <- fit_wham(input1)
#' m1 <- project_wham(m1, n.years = 3)
#' }
#'
#' @export
project_wham <- function(mod, n.years = 3){
  temp <- list(data = mod$env$data, par = mod$parList, map = mod$env$map,
              random = unique(names(mod$env$par[mod$env$random])))
  temp$data$use_indices[(temp$data$n_years_model-(n.years-1)):temp$data$n_years_model,] = 0
  mod$proj <- fit_wham(temp, n.newton = 3, do.sdrep = TRUE, do.retro = FALSE)
  input <- summary(mod$proj$sdrep)
  input <- input[which(rownames(input) == "pred_IAA"),]
  mod$proj$proj_IAA <- list(pred = array(input[,1], dim = dim(mod$proj$rep$pred_IAA)))
  mod$proj$proj_IAA$sd <- array(input[,2], dim = dim(mod$proj$rep$pred_IAA))
  return(mod)
}
