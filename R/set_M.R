#' Specify model and parameter configuration for natural mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
#' If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
#' shared/mirrored among some but not all ages, modify \code{input$map$M_a} after calling \code{prepare_wham_input}
#' (see vignette for more details). \code{M} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Natural mortality model options are:
#'                    \describe{
#'                      \item{"constant"}{(default) estimate a single mean M shared across all ages}
#'                      \item{"age-specific"}{estimate M_a independent for each age}
#'                      \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
#'                        \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
#'                        \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.}
#'                    }
#'                  }
#'     \item{$re}{Time- and age-varying (random effects) on M. Note that random effects can only be estimated if
#'                mean M-at-age parameters are (\code{$est_ages} is not \code{NULL}).
#'                 \describe{
#'                   \item{"none"}{(default) M constant in time and across ages.}
#'                   \item{"iid"}{M varies by year and age, but uncorrelated.}
#'                   \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
#'                   \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
#'                   \item{"2dar1"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan (2016)}.}
#'                 }
#'               }
#'     \item{$initial_means}{Initial/mean M-at-age vector, with length = number of ages (if \code{$model = "age-specific"})
#'                          or 1 (if \code{$model = "weight-at-age" or "constant"}). If \code{NULL}, initial mean M-at-age values are taken
#'                          from the first row of the MAA matrix in the ASAP data file.}
#'     \item{$est_ages}{Vector of ages to estimate M (others will be fixed at initial values). E.g. in a model with 6 ages,
#'                      \code{$est_ages = 5:6} will estimate M for the 5th and 6th ages, and fix M for ages 1-4. If \code{NULL},
#'                      M at all ages is fixed at \code{M$initial_means} (if not \code{NULL}) or row 1 of the MAA matrix from the ASAP file (if \code{M$initial_means = NULL}).}
#'     \item{$logb_prior}{(Only if \code{$model = "weight-at-age"}) TRUE or FALSE (default), should a N(0.305, 0.08) prior be
#'                        used on log_b? Based on Fig. 1 and Table 1 (marine fish) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}.}
#'     \item{$b_size}{if \code{$model = "weight-at-age"}, initial value for random effect with prior defined by logb_prior (intended for simulating data).
#'     \item{$sigma_val}{Initial standard deviation value to use for the M deviations. Values are not used if \code{M$re} = "none". Otherwise, a single value.
#'     \item{$cor_vals}{Initial correlation values to use for the M deviations. If unspecified all initial values are 0. When \code{M$re} = 
#'                  \describe{
#'                    \item{"iid" or "none"}{values are not used.}
#'                    \item{"ar1_a" or "ar1_y"}{cor_vals must be a single value.}
#'                    \item{"2dar1"}{2 values must be specified. First is for "age", second is for "year".}
#'                  }
#'                }
#'   }

set_M = function(input, M)
{
  data = input$data
  par = input$par
  map = input$map
  asap3 = input$asap3
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("log_b", "M_repars", "M_a","M_re"))]
  
  #data$n_M_a = data$n_ages

  #M_model length is n_regions; 1: age-constant, 2: age-specific,3: f(WAA), 4-6 as 1-3, but by stock-specific
  data$M_model = rep(4, data$n_regions)
  data$use_b_prior = 0
  #M_re_model length = n_regions; 1 = none, 2 = IID, 3 = ar1_a, 4 = ar1_y, 5 = 2dar1
  data$M_re_model = rep(1,data$n_regions) # default = no RE / 'none'
  #data$M_est <- rep(0, data$n_M_a) # default = don't estimate M
  #M_first_est = NA  
  M_re_ini = array(NA, dim = c(data$n_stocks, data$n_regions, data$n_years_model, data$n_ages))
  map$M_re = M_re_ini
  M_a_ini = array(NA, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  map$M_a = M_a_ini
  if(!is.null(asap3)) {
    for(i in 1:length(asap3)) for(r in 1:data$n_regions) {
      M_a_ini[i,r,] <- log(asap3[[i]]$M[1,])
      M_re_ini[i,r,,] <- matrix(log(asap3[[i]]$M), data$n_years_model, data$n_M_a) - matrix(M_a_ini[i,r,], data$n_years_model, data$n_ages, byrow=T)
    }
  }

  # natural mortality options, default = use values from ASAP file, no estimation
  if(is.null(M)) {
    if(is.null(asap3)){
      M_a_ini[] <- log(0.2)
      M_re_ini[] <- 0
    }
  }
  if(!is.null(M)){
    if(!is.null(M$model)){ # M model options
      M_mods = c("constant","age-specific","weight-at-age")
      M_mods = c(M_mods, paste0(M_mods,"-stock"))
      if(!(M$model %in% M_mods)) stop(paste0("M$model must be one of these: ", paste0(M_mods, collapse=","))
      if(!is.null(M$re)) {
        if(M$model == "age-specific" & M$re == "ar1_a") stop("Cannot estimate age-specific mean M and AR1 deviations M_a.
          If you want an AR1 process on M-at-age, set M$model = 'constant' and M$re = 'ar1_a'.")
        data$M_model <- match(M$model, M_mods)
        if(M$model %in% M_mods[c(1,3,4,6)]){ #c("constant","weight-at-age") or by stock
          if(!is.null(asap3)) for(i in 1:length(asap3)) for(r in 1:data$n_regions) M_a_ini[i,r,] = log(asap3[[i]]$M[1,1])
          else M_a_ini[] = log(0.2)
          if(!is.null(asap3)) {
            if(is.null(M$initial_means) & any(sapply(asap3,function(x) length(unique(x$M[1,])) > 1))) warning("Constant or weight-at-age M specified (so only 1 mean M parameter),
              but first row of MAA matrix has > 1 unique value.
              Initializing M at age-1 MAA values. To avoid this warning
              without changing ASAP file, specify M$initial_means.")
          }
          if(!is.null(M$logb_prior)){
            if(!is.logical(M$logb_prior)) stop("M$logb_prior must be either TRUE or FALSE")
            if(M$logb_prior) data$use_b_prior = 1
          }
        }
      }
    }
    if(!is.null(M$re)){
      if(length(M$re) != data$n_regions) stop("length(M$re) must be input$data$n_regions")
      if(any(!(M$re %in% c("none","iid","ar1_a","ar1_y","2dar1")))) stop("M$re must be one of the following: 'none','iid','ar1_a','ar1_y','2dar1'")
      data$M_re_model[] <- match(M$re, c("none","iid","ar1_a","ar1_y","2dar1"))
    }
    if(!is.null(M$initial_means)){
      if(!is.array(M$initial_means)) stop("M$initial_means must now be an array with dimensions = c(n_stocks,n_regions,n_ages)") 
      dimsM = dim(M$initial_means)
      if(length(dimsM) != 3) stop("dimensions of M$initial_means must be c(n_stocks,n_regions,n_ages)")
      if(!all(dimsM == dim(M_a_ini)) stop("dimensions of M$initial_means must be c(n_stocks,n_regions,n_ages)")
      M_a_ini[] <- log(M$initial_means)
      # overwrite ASAP file M values
      M_re_ini[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
    }
    if(!is.null(M$fixed)){
      if(!is.array(length(M$fixed)) stop("M$fixed must be an array with dimensions = c(n_stocks,n_regions,n_ages)")
      dimsM = dim(M$fixed)
      if(length(dimsM) != 3) stop("dimensions of M$fixed must be c(n_stocks,n_regions,n_ages)")
      if(!all(dimsM == dim(M_a_ini)) stop("dimensions of M$fixed must be c(n_stocks,n_regions,n_ages)")
      #if(!all(M$est_ages %in% 1:data$n_M_a)) stop("All M$est_ages must be in 1:n.ages (if age-specific M) or 1 (if constant or weight-at-age M)")
      #data$M_est[M$est_ages] = 1 # turn on estimation for specified M-at-age
      #M_first_est <- M$est_ages[1]
      M_re_ini[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
    }
  }
  #data$n_M_est <- sum(data$M_est)

  # natural mortality pars
  par$M_a <- M_a_ini # deviations by age
  for(i in data$M_models){
    if(!is.null(M$fixed)){
      if(i %in% c(1,3)) map$M_a[which(M$fixed==0)] <- 1 #constant value or M = f(WAA) across stock, region, and age
      if(i %in% c(4,6)) { #constant value or M = f(WAA) by stock, across, region, (and age)
        map$M_a[which(M$fixed==0)] <- 1
        k = 1
        for(s in 1:data$n_stocks) { 
          ind = which(map$M_a[s,r,]== 1)
          map$M_a[s,r,ind] = k
        }
        k = k + 1
      }
      if(i == 2) { #age-specific value across stock, region
        map$M_a[which(M$fixed==0)] <- 1 
        for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
          ind = which(map$M_a[s,r,]== 1)
          map$M_a[s,r,ind] = (1:data$n_ages)[ind]
        }
      }
      if(i == 5) { #age-specific value by stock across region
        map$M_a[which(M$fixed==0)] <- 1 
        k = 0
        for(s in 1:data$n_stocks) {
          for(r in 1:data$n_regions) {
            ind = which(map$M_a[s,r,]== 1)
            map$M_a[s,r,ind] = k + (1:data$n_ages)[ind]
          }
          k = k + data$n_ages
        }
      }
    } else{
      #do nothing because default is to not estimate any of the mean parameters
  }
  map$M_a = factor(map$M_a)

  par$M_re <- M_re_ini # deviations from mean M_a on log-scale, PARAMETER_ARRAY
  par$M_repars <- array(0, dim = c(data$n_stocks,data$n_regions,3))
  if(is.null(M$sigma_val)){
    par$M_repars[1] <- log(0.1) # start sigma at 0.1, rho at 0
  } else {
    if(length(M$sigma_val)==data$n_stocks) {
      for(s in 1:data$n_stocks) par$M_repars[s,,1] = log(M$sigma_val[s])
    }
    else stop("length of M$sigma_val must equal input$data$n_stocks.")
  }
  if(is.null(M$cor_vals)) {
    #already set to 0 above
    #if(data$M_re_model == 3) par$M_repars[3] <- 0 # if ar1 over ages only, fix rho_y = 0
    #if(data$M_re_model == 4) par$M_repars[2] <- 0 # if ar1 over years only, fix rho_a = 0
    } else{
    if(!is.matrix(M$cor_vals)) stop("M$cor_vals must be a matrix of values (n_regions x 2")
    if(NROW(M$cor_vals) != data$n_regions | NCOL(M$cor_vals) != 2) stop("M$cor_vals must be a matrix of values  with dimensions n_regions x 2")
    inv_trans_rho <- function(rho) 0.5 * (log(rho+1) - log(1-rho)) # 0.5 because needed transformation on cpp side is unusual.
    for(r in 1:n_regions){
      if(data$M_re_model[r] == 3) par$M_repars[,r,2] <- inv_trans_rho(M$cor_vals[r,1]) # if ar1 over ages only, fix rho_y = 0
      if(data$M_re_model[r] == 4) par$M_repars[,r,3] <- inv_trans_rho(M$cor_vals[r,2]) # if ar1 over years only, fix rho_a = 0
      if(data$M_re_model[r] == 5) par$M_repars[,r,2:3] <- inv_trans_rho(M$cor_vals[r,1:2]) # if 2dar1 years and ages
    }
  }
  # M_repars: sigma_M, rho_M_a, rho_M_y
  k = 1
  map$M_repars = array(NA, dim = dim(par$M_repars))
  for(r in 1:n_regions){
    if(data$M_re_model[r] == 1) #do nothing map$M_repars[,r,] <- NA # no RE pars to estimate
    if(data$M_re_model[r] == 2) { # estimate sigma
      map$M_repars[,r,1] <- k 
      k <- k + 1
    }
    if(data$M_re_model[r] == 3) {# ar1_a: estimate sigma, rho_a
      map$M_repars[,r,1] <- k 
      map$M_repars[,r,2] <- k + 1 
      k <- k + 2
    }
    if(data$M_re_model[r] == 4) {# ar1_y: estimate sigma, rho_y
      map$M_repars[,r,1] <- k 
      map$M_repars[,r,3] <- k + 1 
      k <- k + 2
    }
    if(data$M_re_model[r] == 5) { # 2dar1: estimate all
      map$M_repars[,r,1] <- k 
      map$M_repars[,r,2] <- k + 1 
      map$M_repars[,r,3] <- k + 2 
      k <- k + 3
    }
  }
  map$M_repars = factor(map$M_repars)

  # check if only 1 estimated mean M (e.g. because weight-at-age M or if all but 1 age is fixed), can't estimate rho_a
  # if(data$n_M_est < 2) par$M_repars[2] <- 0
  if(is.null(M$b_size)) par$log_b = log(0.305)
  else par$log_b = log(M$b_size)

  #tmp <- par$M_a
  #tmp[data$M_est==0] = NA
  #ind.notNA <- which(!is.na(tmp))
  #tmp[ind.notNA] <- 1:length(ind.notNA)
  #map$M_a <- factor(tmp)
  #if(data$M_model != 3) map$log_b = factor(rep(NA,length(par$log_b)))

  # M_re: "none","iid","ar1_a","ar1_y","2dar1"
  k = 0
  for(r in 1:data$n_regions)
    if(data$M_re_model[r] == 1) {}#do nothing NAs already defined no RE (either estimate RE for all ages or none at all)
    
    if(data$M_re_model[r] %in% c(2,5)){ # 2d ar1
      for(s in 1:data$n_stocks){
        ind = k + 1:prod(dim(map$M_re)[3:4]) # all y,a estimated
        map$M_re[s,r,,] <- ind
      }
      k = k + prod(dim(map$M_re)[3:4]) 
    }
    if(data$M_re_model[r] == 3){ # ar1_a (devs by age, constant by year)
      for(s in 1:data$n_stocks){
        ind = k + 1:dim(map$M_re)[4]
        map$M_re[s,r,,] <- matrix(ind, dim(map$M_re[3]), dim(map$M_re[4]), byrow = TRUE)
      }
      k = k + dim(map$M_re)[4] 
    }
    if(data$M_re_model[r] == 4){ # ar1_y (devs by year, constant by age)
      for(s in 1:data$n_stocks){
        ind = k + 1:dim(map$M_re)[3]
        map$M_re[s,r,,] <- matrix(ind, dim(map$M_re[3]), dim(map$M_re[4]))
      }
      k = k + dim(map$M_re)[3] 
    }
  }
  map$M_re <- factor(map$M_re)

  input$data = data
  input$par = par
  input$map = map
  
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = wham:::set_random(input)
  return(input)

}