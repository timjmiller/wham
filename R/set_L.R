#' Specify model and parameter configuration for "extra" mortality not directly attributed to natural mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param L (optional) list specifying "extra" mortality options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{L} specifies estimation options for "extra" mortality.
#' If \code{NULL}, This mortality source is not used. 
#' \code{L} is a list with the following entries:
#'   \describe{
#'     \item{$model}{length = n_regions. "extra" mortality model options are:
#'                    \describe{
#'                      \item{"none"}{(default) no extra mortality for this region.}
#'                      \item{"constant"}{estimate a single mean mortality for the region shared across all ages}
#'                      \item{"iid_re"}{estimate independent random effects over years, for the region}
#'                      \item{"ar1_re"}{estimate random effect correlated over years, for the region}
#'                    }
#'                  }
#'     \item{$initial_means}{Initial/mean L-at-region}
#'     \item{$sigma_vals}{Initial standard deviation by region value to use for the L random effects. Values are not used if \code{L$model} = "none".}
#'     \item{$cor_vals}{Initial correlation values to use for the L random effects. If unspecified all initial values are 0}
#'   }
#'
#' @export
set_L = function(input, L)
{
  data = input$data
  par = input$par
  map = input$map
  #asap3 = input$asap3
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("L_repars", "L_re"))]
  
  #L$model length is n_regions
  data$L_model = rep(0, data$n_regions)
  par$L_re = matrix(0, data$n_years_model, data$n_regions)
  par$L_repars = matrix(0,data$n_regions, 3) #mean, sig, rho
  map$L_re = matrix(NA, data$n_years_model, data$n_regions)
  map$L_repars = matrix(NA,data$n_regions, 3) #mean, sig, rho

  # natural mortality options, default = use values from ASAP file, no estimation
  if(!is.null(L)){
    if(!is.null(L$model)){ # L model options
      L_mods = c("none","constant","iid_re","ar1_re")
      if(!(L$model %in% L_mods)) stop(paste0("L$model must be one of these: ", paste0(L_mods, collapse=",")))
      data$L_model[] = match(L$model, L_mods) - 1
    }
  }
  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s
  k = 1
  for(r in 1:data$n_regions) {
    if(data$L_model[r] >0) {
      map$L_repars[r,1] <- k
      k <- k + 1
      if(!is.null(L$initial_means)){
        par$L_repars[r,1] = log(L$initial_means[r])
      }
    }
    if(data$L_model[r] > 1){
      map$L_repars[r,2] <- k
      k <- k + 1
      map$L_re[,r] <- 1
      if(!is.null(L$sigma_vals)){
        par$L_repars[r,2] = log(L$sigma_vals[r])
      }
    }
    if(data$L_model[r] > 2){
      map$L_repars[r,3] <- k
      k <- k + 1
      if(!is.null(L$cor_vals)){
        par$L_repars[r,3] = inv_trans_rho(L$cor_vals[r])
      }
    }
  }
  map$L_re[which(map$L_re==1)] <- 1:sum(map$L_re==1, na.rm = TRUE)
  map$L_re = factor(map$L_re)
  map$L_repars = factor(map$L_repars)

  input$data = data
  input$par = par
  input$map = map
  
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = set_random(input)
  input$options$L <- L
  return(input)

}