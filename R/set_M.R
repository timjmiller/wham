#' Specify model and parameter configuration for natural mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
#' If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
#' shared/mirrored among some but not all ages, modify \code{M$means_map} (see vignette for more details). \code{M} is a list 
#' with the following entries:
#'   \describe{
#'     \item{$mean_model}{Character describing the type of model for M stock and regional models for natural mortality. Options are:
#'       \describe{
#'         \item{"fixed-M"}{Use initial values from ASAP3 dat files or \code{$initial_means} for (mean) M as fixed values. If no ASAP3 files
#'           and \code{$initial_means} is not provided, default is M = 0.2 for all stocks, regions and ages}
#'         \item{"estimate-M"}{estimate one or more (mean) M parameters. Default is to estimate a single M shared across all stocks and ages, but
#'           use \code{$means_map} to fix or estimate parameters for specific stocks, regions, ages.}
#'         \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
#'           \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
#'           \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.
#'           Default is to estimate a single model shared across all stocks and regions, but
#'           use \code{$means_map[s,r,1]} to fix or estimate the intercept for specific stocks, regions. See also \code{$logb_prior}
#'           and \code{$initial_b} configuring the slope on log scale.}
#'       }
#'     }
#'     \item{$initial_means}{array (n_stocks x n_regions x n_ages) of initial/mean M by stock, region and age. If \code{NULL}, initial 
#'       mean M-at-age values for a given stock and region are taken from the first row of the MAA matrix in the ASAP data file. If no
#'        ASAP data file, M = 0.2 is the default. If \code{$mean_model} is "weight-at-age" only 
#'       elements for the first age (\code{$initial_means[,,1]}) are used (for the intercept of log(M)).}
#'     \item{$means_map}{array (n_stocks x n_regions x n_ages) of NA or integers ( 0 <= max <= n_stocks * n_regions * n_ages) indicating 
#'       which ages to estimate (mean) M and whether to set any ages to be identical. E.g. in a model with 2 stock, 2 regions 
#'       and 6 ages with constant M estimated for each stock across regions and ages  \code{$M_ages_map[1,,] = 1} 
#'       and \code{$M_ages_map[2,,] = 2}. \code{$M_ages_map[1,1,] = c(NA,1,1,2,2,3)} will fix M for age 1 at the initial value, 
#'       and estimates for ages 2 and 3 are identical as are those for ages 4 and 5 and different from age 6+ for stock 1 and 
#'       region 1. If \code{NULL}, specifies all ages fixed at \code{M$initial_means}.  If \code{$mean_model} is "weight-at-age"
#'       these are used for all stocks and regions and only the elements for the first age (\code{$M_ages_map[,,1]}) 
#'       are used (for the intercept of log(M)).}
#'     \item{$intial_MAA}{array (n_stocks x n_regions x n_years x n_ages) of initial values for M at age. Intended to be uses when nothing pertaining to M estimated.}
#'     \item{$b_model}{"constant","stock","region", "stock-region" defining whether parameter is constant, stock-specific, region-specific, 
#'       stock- and region-specific. Only used if \code{M$mean_model} = "weight-at-age".}
#'     \item{$b_prior}{T/F, should a N(mu, 0.08) prior (where mu = log(0.305) by default) be used on log_b? Based on Fig. 1 and Table 1 
#'       (marine fish) in  \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}. (Only used if 
#'       \code{$mean_model} is "weight-at-age").}
#'     \item{$intial_b}{if any elements of \code{$mean_model} is "weight-at-age", initial value for mean b for weight-at-age
#'       model.}
#'     \item{$re_model}{Character matrix (n_stocks x n_regions) of options for time- and age-varying (random effects) on M by stock and region.
#'       Possible values are:
#'       \describe{
#'         \item{"none"}{(default) No random effects by age or year.}
#'         \item{"iid_a"}{uncorrelated M by age, constant in time.}
#'         \item{"iid_y"}{uncorrelated M by year, constant all ages.}
#'         \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
#'         \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
#'         \item{"iid_ay"}{M uncorrelated by year and age (2D).}
#'         \item{"ar1_ay"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan
#'           (2016)}.}
#'       }
#'     }
#'     \item{$re_map}{array (n_stocks x n_regions x n_ages) of NA and integers (1 <= max <= n_ages) indicating which ages, for a given 
#'       stock and region, have random effects (not NA) and whether to set RE for any ages to be identical. E.g. in a model with 2 stock, 
#'       2 regions and 6 ages, \code{$re_map[2,1,] = c(NA,1,1,2,2,3)} will not estimate RE for age 1, and those for ages 
#'       2 and 3 are identical as are those for ages 4 and 5 and different from age 6+ for stock 2 and region 1. If \code{NULL}, 
#'       and \code{$re_model} specifies M random effects at age, at least two ages must be 
#'       specified for correlation among ages to be estimated.}
#'     \item{$sigma_vals}{n_stocks x n_regions matrix Initial standard deviation value to use for the M random effects. Values are not used 
#'       if \code{M$re_model} = "none". Otherwise, a single value. If unspecified all values are 0.1.}
#'     \item{$cor_vals}{n_stocks x n_regions x 2 array of initial correlation values to use for the M deviations. If unspecified all initial 
#'       values are 0. When \code{M$re_model} = 
#'       \describe{
#'         \item{"iid_a", "iid_y", "iid_ay" or "none"}{values are not used.}
#'         \item{"ar1_a" }{first value cor_vals[s,r,1] is used.}
#'         \item{"ar1_y" }{second value cor_vals[s,r,2] is used.}
#'         \item{"ar1_ay"}{First is for "age", second is for "year".}
#'       }
#'     }
#'     \item{$sigma_map}{n_stocks x n_region matrix of NA or integers indicating which random effects sd is estimated and whether to 
#'       set any to be identical. If not supplied a single sd will be estimated for any stock and region where $re_model is other than "none".}
#'     \item{$cor_map}{n_stocks x n_region matrix x 2 array of NA or integers indicating which random effects correlation parameters are estimated
#'       and whether to set any to be identical. If not supplied a single value for age and/or year will be estimated for any stock and region where 
#'       $re_model is other than "none", "iid_a", "iid_y".}
#'   }
#'
#' @return a named list with same elements as the input provided with natural mortality options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' M = list(mean_model = "estimate-M")
#' input <- set_q(input, M = M) #estimate a constant M parameters
#' }
#'
#' @export
set_M = function(input, M)
{
  data = input$data
  par = input$par
  map = input$map
  asap3 = input$asap3
  input$log$M <- list()

  # elements of M: model, initial_means, means_map, logb_prior, intial_b, re_model, re_map, sigma_vals, cor_vals
  # data elements are: n_M_re, M_re_index, M_model, M_re_model, use_b_prior, log_b_model
  # data$waa_pointer_M is supplied in waa_info argument to set_WAA.
  # par elements are: Mpars, M_re, M_repars, log_b
  re_mods <- c("none", "iid_a", "iid_y", "ar1_a", "ar1_y", "iid_ay", "ar1_ay")
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("log_b", "M_repars", "Mpars","M_re"))]

  #M_model length is n_regions; 1: Estimate-M, 2: f(WAA), 3: f(WAA) by stock
  data$M_model = 1 #Fixed-M
  data$log_b_model = 1 #constant if estimated
  data$use_b_prior = 0
  # data$n_M_re <- matrix(1, data$n_stocks, data$n_regions)
  data$M_re_index <- array(1, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  data$n_M_re <- matrix(data$n_ages, data$n_stocks, data$n_regions)
  for(s in 1:data$n_stocks) for(r in 1:data$n_regions) data$M_re_index[s,r,] <- 1:data$n_ages
  #M_re_model length = n_regions; 1 = none, 2 = IID, 3 = ar1_a, 4 = ar1_y, 5 = 2dar1
  data$M_re_model = matrix(1,data$n_stocks,data$n_regions) # default = no RE / 'none'
  
  par$Mpars <- array(NA, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  map$Mpars <- par$Mpars
  par$Mpars[] <- log(0.2)
  
  par$M_re <- array(NA, dim = c(data$n_stocks, data$n_regions, data$n_years_model, data$n_ages))
  map$M_re <- par$M_re
  par$M_re[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
  
  par$M_repars <- array(NA, dim = c(data$n_stocks, data$n_regions, 3))
  map$M_repars <- par$M_repars
  par$M_repars[] <- 0
  par$M_repars[,,1] <- log(0.1)
 
  if(is.null(M$mean_model)) M$mean_model <- "fixed-M"
  if(is.null(M$re_model)) M$re_model <- matrix("none", data$n_stocks, data$n_regions)
  if(is.null(M$b_prior)) M$b_prior = FALSE

  if(is.null(M$intial_b)) M$intial_b = 0.305
  if(is.null(M$initial_mean)) {
    if(!is.null(asap3)){
      M$initial_means <- array(NA, dim = c(data$n_stocks, data$n_regions, data$n_ages))
      for(i in 1:length(asap3)) for(r in 1:data$n_regions) {
        M$initial_means[i,r,] <- exp(apply(log(asap3[[i]]$M),2,mean))
        # print(M$initial_means[i,r,])
        # print(log(asap3[[i]]$M))
        for(y in 1:data$n_years_model) {
          par$M_re[i,r,y,] <- log(asap3[[i]]$M[y,]) - log(M$initial_means[i,r,])
        }
      }
    }
    else M$intitial_means <- array(0.2, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  }

  par$log_b <- matrix(log(M$intial_b), data$n_stocks, data$n_regions)
  map$log_b <- matrix(NA, data$n_stocks, data$n_regions)
 
  if(length(M$mean_model) != 1 & !is.character(M$mean_model)) stop("M$mean_model must be 'fixed-M', 'estimate-M', or 'weight-at-age'")
  M_mods = c("fixed-M", "estimate-M", "weight-at-age")
  if(!(M$mean_model %in% M_mods)) stop(paste0("M$mean_model must be one of these: ", paste0(M_mods, collapse=",")))
  data$M_model[] <- ifelse(M$mean_model %in% M_mods[1:2], 1, 2) #only needs to know f(age) or f(waa)

  if(!is.null(M$initial_means)){
    if(!is.array(M$initial_means)) stop("M$initial_means must now be an array with dimensions = c(n_stocks,n_regions,n_ages)") 
    dimsM = dim(M$initial_means)
    if(length(dimsM) != 3) stop("dimensions of M$initial_means must be c(n_stocks,n_regions,n_ages)")
    if(!all(dimsM == dim(par$Mpars))) stop("dimensions of M$initial_means must be c(n_stocks,n_regions,n_ages)")  
    par$Mpars[] <- log(M$initial_means)
  }
  if(!is.null(M$initial_MAA)){
    if(!is.array(M$initial_MAA)) stop("M$initial_MAA must now be an array with dimensions = c(n_stocks,n_regions,n_years,n_ages)") 
    dimsM = dim(M$initial_MAA)
    if(length(dimsM) != 4) stop("dimensions of M$initial_MAA must be c(n_stocks,n_regions,n_years,n_ages)")
    if(!all(dimsM == dim(par$M_re))) stop("dimensions of M$initial_MAA must be c(n_stocks,n_regions,n_years,n_ages)")
    for(i in 1:data$n_stocks) for(r in 1:data$n_regions) {
      par$Mpars[i,r,] <- apply(log(M$initial_MAA[i,r,,]),2,mean)
      for(y in 1:data$n_years_model) par$M_re[i,r,y,] <- log(M$initial_MAA[i,r,y,]) - par$Mpars[i,r,]
    }
  }

  if(!is.null(M$re_model)){
    if(!is.matrix(M$re_model)) stop("M$re_model must be a character n_stocks x n_regions matrix.")
    dimsM = dim(M$re_model)
    if(length(dimsM) != 2) stop("dimensions of M$re_model must be n_stocks x n_regions.")
    if(!all(dimsM == dim(par$Mpars)[1:2])) stop("dimensions of M$re_model must be n_stocks x n_regions.")  
    if(any(!(M$re_model %in% re_mods))) {
      stop(paste0("Each M$re_model must be one of the following: ", paste(re_mods, collapse = ",")))
    }
  }

  if(!is.null(M$means_map)){
    if(!is.array(M$means_map)) stop("M$means_map must be an array with dimensions = c(n_stocks,n_regions,n_ages)")
    dimsM = dim(M$means_map)
    if(length(dimsM) != 3) stop("dimensions of M$means_map must be c(n_stocks,n_regions,n_ages)")
    if(!all(dimsM == dim(par$Mpars))) stop("dimensions of M$means_map must be c(n_stocks,n_regions,n_ages)")
  }

  if(!is.null(M$re_map)){
    if(!is.array(M$re_map)) stop("M$re_map must be an array with dimensions = c(n_stocks,n_regions,n_ages)")
    dimsM = dim(M$re_map)
    if(length(dimsM) != 3) stop("dimensions of M$re_map must be c(n_stocks,n_regions,n_ages)")
    if(!all(dimsM == dim(par$Mpars))) stop("dimensions of M$re_map must be c(n_stocks,n_regions,n_ages)")
    if(any(!(M$re_map %in% c(NA,1:data$n_ages)))) stop("Entries in M$re_map must be NA or in 1:n_ages.")
    k <- 0
    if(M$mean_model != "weight-at-age" ) for(s in 1:data$n_stocks) for(r in 1:data$n_regions) if(M$re_model[s,r] != "none"){
      temp <- M$re_map[s,r,]
      if(all(is.na(temp))) stop(paste0("Random effects have been specified for stock ", s, " in region ", r, ", but M$re_map[s,r,] are all NA."))
      ind <- unclass(factor(temp))
      if(M$re_model[s,r] %in% c("ar1_a","iid_a", "iid_ay","ar1_ay")){
        data$n_M_re[s,r] <- length(unique(temp[!is.na(temp)]))
        temp2 <- temp
        temp2[which(is.na(temp2))] <-0
        data$M_re_index[s,r,] <- unclass(factor(temp2))
        #temp[which(temp==0)] <- NA
        if(M$re_model[s,r] %in% c("ar1_a","iid_a")) {
          for(y in 1:dim(map$M_re)[3]) map$M_re[s,r,y,] <- ind + k 
          k <- k + data$n_M_re[s,r]
        }
        if(M$re_model[s,r] %in% c("ar1_ay","iid_ay")) {
          for(y in 1:dim(map$M_re)[3]) {
            map$M_re[s,r,y,] <- ind + k 
            k <- k + data$n_M_re[s,r]
          }
        }
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      if(M$re_model[s,r] %in% c("ar1_y","iid_y")){
        data$n_M_re[s,r] <- 1
        temp2 <- temp
        temp2[which(is.na(temp))] <- 0
        temp2[which(!is.na(temp))] <- 1
        data$M_re_index[s,r,] <- temp2
        for(a in 1:dim(map$M_re)[4]) map$M_re[s,r,,a] <- 1:dim(map$M_re)[3] + k 
        k <- k + dim(map$M_re)[3]
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      k <- k + 1
    }
  } else{
    k <- 0
    if(M$mean_model != "weight-at-age") for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("ar1_a","iid_a")){
        data$n_M_re[s,r] <- data$n_ages
        data$M_re_index[s,r,] <- 1:data$n_ages
        for(y in 1:dim(map$M_re)[3]) map$M_re[s,r,y,] <- 1:data$n_ages + k
        k <- k +  data$n_ages
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      if(M$re_model[s,r] %in% c("iid_ay","ar1_ay")){
        data$n_M_re[s,r] <- data$n_ages
        data$M_re_index[s,r,] <- 1:data$n_ages
        map$M_re[s,r,,] <- 1:(dim(map$M_re)[3]*data$n_ages) + k
        k <- k + dim(map$M_re)[3]*data$n_ages
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      if(M$re_model[s,r] %in% c("ar1_y","iid_y")){
        data$n_M_re[s,r] <- 1
        data$M_re_index[s,r,] <- 1
        for(a in 1:dim(map$M_re)[4]) map$M_re[s,r,,a] <- 1:dim(map$M_re)[3] + k
        k <- k + dim(map$M_re)[3]
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
    }

  }

  if(is.null(M$means_map)){ #constant M or weight at age being used.
    if(M$mean_model != "fixed-M") { #Estimating constant M or WAA intercept
      map$Mpars[] = 1
      #only use first value of M from asap files if estimating M for constant or f(WAA)
      for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
        if(length(unique(par$Mpars[s,r,])) > 1) {
          input$log$M <- c(input$log$M, 
          "M is estimated and no M$means_map is specified so only 1 mean M parameter), but for some stock/region, MAA has > 1 unique value.
          Initializing M at M$initial_means[1,1,1]. To avoid this warning without changing ASAP file, specify M$initial_means appropriately.\n")
        }
      }
      par$Mpars[] = par$Mpars[1,1,1]
    }
  }

  if((M$mean_model == "weight-at-age") & any(M$re_model == "ar1_a")) stop("Cannot estimate M as a function of weight at age and \n 
    random effects just on age. If you want random effects M parameters for on weight-at-age, set M$re_model = 'iid_y' or 'ar1_y'.")    
  if(any(M$re_model %in% c("ar1_a")) & !is.null(M$means_map)) { #this flag may need to be stock and region specific
    for(i in 1:data$n_stocks) for(j in 1:data$n_regions){
      if(length(unique(M$means_map[i,j,])) != 1) stop("If M$re_model is 'ar1_a', all M$means_map must be the same across ages \n
        (NA = fixed mean or an integer = estimated mean).")
    }
  }

  temp <- c("none", "ar1_a","ar1_y", "ar1_a","ar1_y", "ar1_ay", "ar1_ay")
  temp <- temp[match(M$re_model, re_mods)]
  data$M_re_model[] <- match(temp, c("none","ar1_a","ar1_y","ar1_ay"))

  if(M$mean_model == "weight-at-age"){
    if(is.null(M$b_model)){
      M$b_model <- "constant"
      input$log$M <- c(input$log$M, "M$b_model was not specified, so M as a function of weight at age will be used for all stocks and regions.\n")
    }
    if(length(M$b_model) != 1 | (!M$b_model %in% c("constant", "stock", "regions", "stock_region"))){
      stop("M$b_model must be a single value: 'constant', 'stock', 'region', or 'stock_region'.")
    }
    map$log_b[] <- 1 #log_b_model = 1: constant across stock, region
    map$Mpars[] <- 1 #the same "intercept": constant across stock, region
    if(M$b_model == "constant") par$Mpars[] <- par$Mpars[1,1,1]
    if(M$b_model == "stock") {
      data$log_b_model <- 2
      for(s in 1:data$n_stocks) {
        map$log_b[s,] <- s
        map$Mpars[s,,] <- s
        par$Mpars[s,,] <- par$Mpars[s,1,1]
      }
    }
    if(M$b_model == "region") {
      data$log_b_model <- 3
      for(r in 1:data$n_regions) {
        map$log_b[,r] <- r
        map$Mpars[,r,] <- r
        par$Mpars[,r,] <- par$Mpars[1,r,1]
      }
    }
    if(M$b_model == "stock_region") {
      data$log_b_model <- 4
      for(s in 1:data$n_stocks) for(r in 1:data$n_regions){
        map$log_b[s,r] <- r + (s-1) * data$n_regions
        map$Mpars[s,r,] <- r + (s-1) * data$n_regions
        par$Mpars[s,r,] <- par$Mpars[s,r,1]
      }
    }
    if(is.null(M$b_prior)){
      M$b_prior = FALSE
      input$log$M <- c(input$log$M, "M$b_prior was not specified, so prior for the b parameter of M_a = aW_a^b will not be used.\n")
    }
    if(length(M$b_prior) != 1 | !is.logical(M$b_prior)) stop("M$b_prior must be single value: TRUE or FALSE")
    if(M$b_prior) data$use_b_prior = 1
  }

  if(!is.null(M$sigma_vals)){
    if(!is.matrix(M$sigma_vals)) stop("M$sigma_vals must be a n_stocks x n_regions matrix.")
    dimsM = dim(M$sigma_vals)
    if(length(dimsM) != 2) stop("dimensions of M$sigma_vals must be n_stocks x n_regions.")
    if(!all(dimsM == dim(par$Mpars)[1:2])) stop("dimensions of M$sigma_vals must be n_stocks x n_regions.")
    if(any(M$sigma_vals<0)) stop("M$sigma_vals must be > 0.")
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) if(M$re_model[s,r] != "none") par$M_repars[s,r,1] <- log(M$sigma_vals[s,r])
  }
  if(!is.null(M$sigma_map)){
    if(!is.matrix(M$sigma_map)) stop("M$sigma_map must be a n_stocks x n_regions matrix.")
    dimsM = dim(M$sigma_map)
    if(length(dimsM) != 2) stop("dimensions of M$sigma_map must be n_stocks x n_regions.")
    if(!all(dimsM == dim(par$Mpars)[1:2])) stop("dimensions of M$sigma_map must be n_stocks x n_regions.")
  }

  #inv_trans_rho <- function(rho) 0.5 * (log(rho+1) - log(1-rho)) # 0.5 because needed transformation on cpp side is unusual.
  inv_trans_rho <- function(rho) (log(rho+1) - log(1-rho)) # 0.5 because needed transformation on cpp side is unusual.
  if(!is.null(M$cor_vals)){
    if(!is.array(M$cor_vals)) stop("M$cor_vals must be an array with dimensions = c(n_stocks,n_regions,2)")
    dimsM = dim(M$cor_vals)
    if(length(dimsM) != 3) stop("dimensions of M$cor_vals must be c(n_stocks,n_regions,2)")
    if(!all(dimsM == c(data$n_stocks,data$n_regions,2))) stop("dimensions of M$cor_vals must be c(n_stocks,n_regions,2)")
    if(any(abs(M$cor_vals) > 1)) stop("M$cor_vals must be > -1 and < 1.")
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("ar1_ay","ar1_y")) par$M_repars[s,r,3] <- inv_trans_rho(M$cor_vals[s,r,2])
      if(M$re_model[s,r] %in% c("ar1_ay","ar1_a")) par$M_repars[s,r,2] <- inv_trans_rho(M$cor_vals[s,r,1])
    }
  }
  if(!is.null(M$cor_map)){
    if(!is.array(M$cor_map)) stop("M$cor_map must be an array with dimensions = c(n_stocks,n_regions,2)")
    dimsM = dim(M$cor_map)
    if(length(dimsM) != 3) stop("dimensions of M$cor_map must be c(n_stocks,n_regions,2)")
    if(!all(dimsM == c(data$n_stocks,data$n_regions,2))) stop("dimensions of M$cor_map must be c(n_stocks,n_regions,2)")
  }

  # map mean M 
  if(is.null(M$means_map)){
    #estimate a single M for all ages or a single intercept for M = f(WAA) 
    if(M$mean_model != "fixed-M") map$Mpars[] = 1
  } else{ #use M$means_map
    if(M$mean_model == "weight-at-age"){ #only use the first age for the intercept
      map$Mpars[,,1] <- unclass(factor(M$means_map[,,1]))
    } else{
      #print(M$means_map)
      map$Mpars[] <- unclass(factor(M$means_map)) 
    }
  }

  #############
  #map M_repars
  if(!is.null(M$sigma_map)){
    M$sigma_map[which(M$re_model=="none")] <- NA
    M$sigma_map[] <- unclass(factor(M$sigma_map))
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) if(M$re_model[s,r] != "none") map$M_repars[s,r,1] <- M$sigma_map[s,r]
  } else{ #all the same value
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) if(M$re_model[s,r] != "none") map$M_repars[s,r,1] <- 1
  }
  k <- 0
  #print(map$M_repars[1,1,])
  if(any(!is.na(map$M_repars))) k <- max(map$M_repars, na.rm= T)
  #print(k)

  if(!is.null(M$cor_map)){
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("none","iid_y","ar1_y", "iid_ay")) M$cor_map[s,r,1] <- NA
      if(M$re_model[s,r] %in% c("none","iid_a","ar1_a", "iid_ay")) M$cor_map[s,r,2] <- NA
    }
    M$cor_map[] <- k + unclass(factor(M$cor_map))
  } else { #same values for age and same values for year
    M$cor_map <- array(NA, dim = c(data$n_stocks, data$n_regions, 2))
    if(any(M$re_model %in% c("ar1_a", "ar1_ay"))) {
      M$cor_map[,,1] <- k + 1
      k <- k + 1
    }
    if(any(M$re_model %in% c("ar1_y", "ar1_ay"))) {
      M$cor_map[,,2] <- k + 1
      k <- k + 1
    }
  }
  #print(M$cor_map[1,1,])
  #print(map$M_repars[1,1,])
  for(s in 1:data$n_stocks)for(r in 1:data$n_regions) map$M_repars[s,r,2:3] <- M$cor_map[s,r,]
  #print(map$M_repars[1,1,])
  #############

  map$M_repars = factor(map$M_repars)
  map$M_re <- factor(map$M_re)
  map$Mpars = factor(map$Mpars)
  map$log_b = factor(map$log_b)

  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$M)) input$log$M <- c("Natural Mortality: \n", input$log$M)
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = set_random(input)
  input$options$M <- M
  if(!is_internal_call()) cat(unlist(input$log$M, recursive=T))
  return(input)

}

