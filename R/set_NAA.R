#' Specify model and parameter configuration for numbers at age
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param NAA_re (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' 
#' If \code{NAA_re = NULL}, a traditional statistical catch-at-age model is fit (NAA = pred_NAA for all ages, deterministic). Otherwise,
#' \code{NAA_re} specifies numbers-at-age configuration. It is a list with the following possible entries:
#'   \describe{
#'     \item{$decouple_recruitment}{T/F determining whether correlation structure of recruitment is independent of RE deviations for older ages 
#'        (default = TRUE). Only applicable for \code{NAA_re$sigma = "rec+1"} and correlation across ages is specified. If TRUE and \code{NAA_re$cor = "ar1_a"}, only deviations for ages>1 
#'        have the correlation structure. If TRUE and NAA_re$cor is "ar1_y" or "2dar1" separate year correlation parameters are estimated for recruitment and older
#'        ages.
#'     }
#'     \item{$sigma}{Which ages allow deviations from the predicted NAA given NAA from previous time step? Must be a single character string described below or a vector
#'                    of length n_stocks. If length = 1, assumptions will be applied to all stocks. Options are specified with the strings:
#'                    \describe{
#'                      \item{"rec"}{Random effects on recruitment (deviations), all other ages deterministic}
#'                      \item{"rec+1"}{"Full state space" model with 2 estimated \code{sigma_a}, one for recruitment and one shared among other ages}
#'                    }
#'                  }
#'     \item{$sigma_vals}{Array (n_stocks x n_regions x n_ages) of initial standard deviation values to use for the NAA deviations. Values are not used if recruit_model = 1 and \code{NAA_re$sigma} is
#'                  not specified. Only those for age 1 in the spawning region are used if \code{NAA_re$sigma} = "rec".
#'                  If \code{NAA_re$sigma_map} is defined, the user must ensure that the configuration is compatible with \code{NAA_re$sigma_vals}
#'                }
#'     \item{$sigma_map}{You can specify a more complex parameter configuration by entering an integer array (nstocks x n_regions x n_ages), where each entry points to the
#'                   NAA_sigma to use for that stock and age. E.g., for 2 stocks, 2 regions, and 6 ages, array(rep(c(1,2,2,3,3,3), each = 4),c(2,2,6)) will estimate 3 
#'                   sigmas, with recruitment (age-1) deviations having their own \code{sigma_R}, ages 2-3 sharing \code{sigma_2}, and ages 4-6 sharing \code{sigma_3}. 
#'                   All parameters being the same for both stocks and across regions. The user must be sure that a compatible \code{NAA_re$sigma} configuration is defined.
#'                   Values are not used if recruit_model = 1 and \code{NAA_re$sigma} is not specified.
#'                  }
#'     \item{$cor}{Correlation structure for the NAA deviations. Must be a single character string described below or a vector
#'                  of length n_stocks. If length = 1, assumptions will be applied to all stocks. Options are:
#'                  \describe{
#'                    \item{"iid"}{NAA deviations vary by year and age, but uncorrelated.}
#'                    \item{"ar1_a"}{NAA deviations correlated by age (AR1).}
#'                    \item{"ar1_y"}{NAA deviations correlated by year (AR1).}
#'                    \item{"2dar1"}{NAA deviations correlated by year and age (2D AR1).}
#'                  }
#'                }
#'     \item{$cor_vals}{Array (n_stocks x n_regions x 3) of initial correlation values to use for the NAA deviations. The first value correspond to correlation across age.
#'                  The second is for yearly correlation for NAA deviations for all ages if recruitment is not decoupled, or otherwise just ages after recruitemnt. The third correlation is for annual 
#'                  recruitment deviations and used only when recruitment is decoupled. If unspecified all initial values are 0. If \code{NAA_re$cor} = "iid", values are not used. 
#'                  If \code{NAA_re$cor} = "ar1_a", those for yearly correlation are not used, and vice versa for "ar1_y". 
#'                }
#'     \item{$cor_map}{n_stocks x n_region matrix x 3 array of NA or integers indicating which random effects correlation parameters are estimated
#'       and whether to set any to be estimated identically. The 3 values for a given stock and region are for age, year ("survival"), and year ("recruitment"). 
#'       Both age and years values are used for \code{NAA_re$cor = "2dar1"}, but only the appropriate values are used for \code{NAA_re$cor = "ar1_a"}, or \code{"ar1_y"}.
#'       When \code{NAA_re$decouple_recruitment = TRUE}, the third value is used for both \code{NAA_re$sigma = "rec"} or \code{"rec+1"}. For 
#'       \code{NAA_re$decouple_recruitment = FALSE}, only the second value is used and applies to both recruitment and "survival" when \code{NAA_re$sigma = "rec+1"}, 
#'       and just recruitment when \code{NAA_re$sigma = "rec"}. If not supplied stock-specific values for age and/or year will be estimated for all regions where 
#'       \code{NAA_re$cor} is other than "none", "iid".
#'     }
#'     \item{$N1_model}{Character vector (n_stocks) determining which way to model the initial numbers at age:
#'       \describe{
#'          \item{"age-specific-fe"}{(default) age- and region-specific fixed effects parameters}
#'          \item{"equilibrium"}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'          \item{"iid-re"}{(default) age- and region-specific iid random effects parameters. 2 parameters: mean and sd for log NAA}
#'          \item{"ar1-re"}{(default) age- and region-specific random effects parameters. 3 parameters: mean and sd, and cor for log NAA}
#'       }
#'     }
#'     \item{$N1_pars}{An (n_stocks x n_regions x n_ages) array. If N1_model = 0, then this should be filled with the initial values to use for abundance at age by stock and region in the first year. 
#'        If N1_model = 1 (equilibrium assumption), only the first two values in the ages dimension are used: the (s,r,1) value is recruitment for stock (and region) and (s,r,2) is the fully-selected 
#'        equilibrium fishing mortality rate generating the rest of the numbers at age in the first year.
#'     }
#'     \item{$recruit_model}{Integer vector (n_stocks) determining how to model recruitment. Overrides \code{recruit_model} argument to \code{prepare_wham_input}. Must make sure \code{NAA_re$sigma}, \code{NAA_re$cor}
#'        and \code{ecov} are properly specified.
#'       \describe{
#'           \item{1}{SCAA, estimating annual recruitments as fixed effects.}
#'           \item{2}{estimating a mean recruitment with annual recruitments as random effects}
#'           \item{3}{Beverton-Holt stock-recruitment with annual recruitments as random effects}
#'           \item{4}{Ricker stock-recruitment with annual recruitments as random effects}
#'       }
#'     }
#'     \item{$recruit_pars}{list (length = n_stocks) of vectors of initial parameters for recruitment model. If $recruit_model is 3 or 4, parameters are "alpha" and "beta".
#'     }
#'   }
#'
#' @return a named list with same elements as the input provided with abundance modeling options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' NAA = list(sigma = "rec")
#' input <- set_q(input, NAA_re = NAA) #estimate recruitment as random effects
#' }
#'
#' @export
set_NAA = function(input, NAA_re=NULL)
{

  data = input$data
  par = input$par
  map = input$map
  asap3 = input$asap3
  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s 
  input$log$NAA <- list()
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mean_rec_pars", "log_N1", "log_NAA_sigma", "trans_NAA_rho", "log_NAA"))]
  
  #if(is.null(input$asap3)) asap3 = NULL
  #else asap3 = input$asap3


  #set up initial NAA
  #0: just age-specific numbers at age
  data$N1_model = rep(0, data$n_stocks)
  
  # data$decouple_recruitment <- 0 #until all examples, tests, vignettes are changed
  data$decouple_recruitment <- 1 #decouple is default now!
  if(!is.null(NAA_re$decouple_recruitment)) data$decouple_recruitment <- as.integer(NAA_re$decouple_recruitment)

  par$log_N1 = array(0,dim = c(data$n_stocks,data$n_regions,data$n_ages))
  map$log_N1 = array(NA,dim = c(data$n_stocks,data$n_regions,data$n_ages))
  par$N1_repars = array(0,dim = c(data$n_stocks,data$n_regions,3))
  map$N1_repars = array(NA,dim = c(data$n_stocks,data$n_regions,3))
  init_NAA = log(exp(10)*exp(-(0:(data$n_ages-1))*0.2))
  if(!is.null(NAA_re$N1_model)) {
    options = c("age-specific-fe", "equilibrium","iid-re", "ar1-re")
    k = 1
    for(s in 1:data$n_stocks) {
      if(!(NAA_re$N1_model[s] %in% options)) stop("NAA_re$N1_model must all be 'age-specific-fe', 'equilibrium', 'iid-re' or 'ar1-re'.")
      if(NAA_re$N1_model[s] == options[1]) data$N1_model[s] = 0
      if(NAA_re$N1_model[s] == options[2]) data$N1_model[s] = 1
      if(NAA_re$N1_model[s] %in% options[3:4]) {
        data$N1_model[s] = 2
        map$N1_repars[s,,1] = k
        map$N1_repars[s,,2] = k + 1
        k = k + 2
        if(NAA_re$N1_model[s] == options[4]){ #ar1 with age
          map$N1_repars[s,,3] = k
          k = k + 1
        }
      }
    }
  }
  map$N1_repars = factor(map$N1_repars)
  k = 1
  if(any(data$N1_model == 2) & any(data$N1_model != 2)) stop("If any initial numbers at age are treated as RE, then all must.")
  for(s in 1:data$n_stocks) {
    if(data$N1_model[s] == 0){
      for(r in 1:data$n_regions) for(a in 1:data$n_ages) {
        if(data$NAA_where[s,r,a] == 1) {
          if(!is.null(asap3)) {
            par$log_N1[s,r,a] = log(asap3[[s]]$N1_ini[a]) # use N1_ini values from asap3 file
          } else{
            par$log_N1[s,r,a] = init_NAA[a]
          }
          map$log_N1[s,r,a] = k
          k <- k + 1
        }
      }
    }
    if(data$N1_model[s] == 1) { #equilibrium assumption, 2 pars per stock
      par$log_N1[s,data$spawn_regions[s],1:2] = c(10,log(0.1)) # allowed in wham.cpp but no option to set here (must be hard-coded after calling prepare_wham_input)
      map$log_N1[s,data$spawn_regions[s],1:2] = k + 0:1
      k = k + 2
    }
    if(data$N1_model[s] == 2) { #RE
      for(r in 1:data$n_regions) for(a in 1:data$n_ages) {
        if(data$NAA_where[s,r,a] ==1) {
          par$log_N1[s,r,a] = init_NAA[a]
          map$log_N1[s,r,a] = k
          k <- k + 1
        }
      }
    }
  }
  if(!is.null(NAA_re[["N1_pars"]])){
    par$log_N1[] = log(NAA_re$N1_pars)
  }
  map$log_N1 = factor(map$log_N1)

  # NAA_re options for beyond year 1
  # default = SCAA for each stock
  data$NAA_re_model = rep(0, data$n_stocks)
  par$log_NAA_sigma = array(0, c(data$n_stocks, data$n_regions, data$n_ages))
  map$log_NAA_sigma = array(NA, c(data$n_stocks, data$n_regions, data$n_ages))
  par$trans_NAA_rho = array(0,c(data$n_stocks, data$n_regions, 3))
  map$trans_NAA_rho = array(NA,c(data$n_stocks, data$n_regions, 3))
  par$log_NAA = array(10,dim = c(data$n_stocks, data$n_regions, data$n_years_model-1, data$n_ages))
  map$log_NAA = array(NA,dim = c(data$n_stocks, data$n_regions, data$n_years_model-1, data$n_ages))
  for(s in 1:data$n_stocks){
    map$log_NAA[s,data$spawn_regions[s],,1] <- 1 #change to unique values later
  }

  #if(is.null(NAA_re$sigma)){ 
  #  data$n_NAA_sigma <- 0
  #  data$NAA_sigma_pointers <- rep(1,data$n_ages)
  #}
  if(!is.null(NAA_re$sigma)){
    k <- 1
    if(!length(NAA_re$sigma) %in% c(1,data$n_stocks)) stop("NAA_re$sigma length must be 1 or equal to the number of stocks.")
    if(length(NAA_re$sigma) == 1) {
      input$log$NAA <- c(input$log$NAA, paste0("\n Same NAA_re$sigma being used for all stocks (", NAA_re$sigma[[1]][1], ").\n"))
      #NAA_re$sigma = rep(list(NAA_re$sigma), data$n_stocks)
      map$log_NAA_sigma[,,1] <- k
      data$NAA_re_model[] <- 1
      #for(s in 1:data$n_stocks) for(r in 1:data$n_regions) if(data$NAA_where[s,r,1]==1) map$log_NAA[s,r,,1] <- 1 #change to unique values later
      if(NAA_re$sigma[[1]][1] == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
        map$log_NAA_sigma[,,-1] <- k+1
        data$NAA_re_model[] <- 2
        #also RE for ages>1
        for(s in 1:data$n_stocks) for(r in 1:data$n_regions) for(a in 2:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] <- 1 #change to unique values later
      }
    } else {
      for(s in 1:data$n_stocks) {
        if(NAA_re$sigma[[s]][1] == "rec"){
          data$NAA_re_model[s] <- 1
          map$log_NAA_sigma[s,,1] <- k
          #below is already done above for SCAA
          #map$log_NAA[s,data$spawn_regions[s],,1] = 1 #change to unique values later
          k <- k + 1
          #data$n_NAA_sigma <- 1
          #data$NAA_sigma_pointers <- rep(1,data$n_ages)
        } else {
          if(NAA_re$sigma[[s]][1] == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
            data$NAA_re_model[s] <- 2
            for(r in 1:data$n_regions) map$log_NAA_sigma[s,r,] <- c(k, rep(k+1, data$n_ages-1))
            for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] <- 1 #change to unique values later
            k <- k + 2
            #data$n_NAA_sigma <- 2
            #data$NAA_sigma_pointers <- c(1,rep(2,data$n_ages-1))
          } else {
            if(length(NAA_re$sigma[[s]]) != data$n_ages) stop("each element of NAA_re$sigma must either be 'rec' (random effects on recruitment only), 
    'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
    or a vector with length == n.ages specifying which sigma_a to use for each age.")
            #if(length(NAA_re$sigma[[s]]) == data$n_ages){
              if(any(is.na(unique(NAA_re$sigma[[s]])))){
                #use SCAA because of na values
                #data$n_NAA_sigma <- 0
                #data$NAA_sigma_pointers <- rep(1,data$n_ages)            
              } else {
                tmp <- unique(NAA_re_sigma[[s]])
                ind <- 1:length(tmp)
                ind <- ind[match(NAA_re_sigma[[s]],tmp)] - 1
                for(r in 1:data$n_regions) map$log_NAA_sigma[s,r,] <- k + ind
                for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] <- 1 #change to unique values later
                k <- max(k + ind, na.rm=T)
                #data$n_NAA_sigma <- max(unique(NAA_re$sigma), na.rm=T)
                #data$NAA_sigma_pointers <- NAA_re$sigma
              }
            #}
          }
        }
      }
    }
    if(!is.null(NAA_re$sigma_vals)) {
      if(!is.array(NAA_re$sigma_vals)) stop("NAA_re$sigma_vals must be an array with dimensions: n_stocks x n_regions x n_ages.")
      else {
        if(any(dim(NAA_re$sigma_vals) != c(data$n_stocks, data$n_regions, data$n_ages))) stop("NAA_re$sigma_vals must be an array with dimensions: n_stocks x n_regions x n_ages.")
        par$log_NAA_sigma[] <- log(NAA_re$sigma_vals)
      }
    }

    if(is.null(NAA_re$cor)) NAA_re$cor <- "iid"
    if(!(length(NAA_re$cor) %in% c(1,data$n_stocks))) stop("NAA_r$cor must have length 1 or n_stocks")
    k <- 0
    constant <- length(NAA_re$cor)==1
    for(s in 1:data$n_stocks) {
      ind <- ifelse(constant,1,s)
      if(!constant) k <- max(c(0, map$trans_NAA_rho), na.rm = TRUE)
      if(!is.null(NAA_re$cor[ind])){
        if(!NAA_re$cor[ind] %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor must be one of 'iid','ar1_a','ar1_y','2dar1'")
        if(NAA_re$cor[ind] %in% c("ar1_a")) map$trans_NAA_rho[s,,1] <- k+1
        if(NAA_re$cor[ind] %in% c("ar1_y")) {
          if(data$NAA_re_model[s] == 2) map$trans_NAA_rho[s,,2] <- k+1
          if(data$NAA_re_model[s] == 1 & !data$decouple_recruitment) map$trans_NAA_rho[s,,2] <- k+1
          if(data$decouple_recruitment) map$trans_NAA_rho[s,,3] <- k+2
        }
        if(NAA_re$cor[ind] == "2dar1") for(r in 1:data$n_regions) {
          map$trans_NAA_rho[s,r,1:2] <- k + 1:2
          if(data$decouple_recruitment) map$trans_NAA_rho[s,r,3] <- k + 3
        } else {
          # NAA_re$cor[s] <- 'iid'
        }
      }
    }
    if(!is.null(NAA_re$cor_vals)) {
      if(!is.array(NAA_re$cor_vals)) stop("NAA_re$cor_vals must be an array with dimensions: n_stocks x n_regions x 3.")
      else {
        if(any(dim(NAA_re$cor_vals) != c(data$n_stocks, data$n_regions, 3))) stop("NAA_re$cor_vals must be an array with dimensions: n_stocks x n_regions x 3.")
        for(s in 1:data$n_stocks) {
          if(NAA_re$cor[s] %in% c("2dar1","ar1_a")) par$trans_NAA_rho[s,,1] <- inv_trans_rho(NAA_re$cor_vals[s,,1])
          if(NAA_re$cor[s] %in% c("2dar1","ar1_y")){
            par$trans_NAA_rho[s,,2] <- inv_trans_rho(NAA_re$cor_vals[s,,2])
            if(data$decouple_recruitment) par$trans_NAA_rho[s,,3] <- inv_trans_rho(NAA_re$cor_vals[s,,3])
          }
        }
      }
    }
    if(!is.null(NAA_re$sigma_map)) {
      if(!is.array(NAA_re$sigma_map)) stop("NAA_re$sigma_map must be an array with dimensions = nstocks x nregions x nages")
      if(any(dim(NAA_re$sigma_map) != c(data$n_stocks, data$n_regions, data$n_ages))) stop("dimensions of NAA_re$sigma_map array are not c(nstocks,nregions,nages)")
      map$log_NAA_sigma[] <- NAA_re$sigma_map
    }
    if(!is.null(NAA_re$cor_map)) {
      if(!is.array(NAA_re$cor_map)) stop("NAA_re$cor_map must be an array with dimensions = nstocks x nregions x 3")
      if(any(dim(NAA_re$cor_map) != c(data$n_stocks, data$n_regions, 3))) stop("dimensions of NAA_re$cor_map array are not c(nstocks,nregions,3)")
      map$trans_NAA_rho[] <- NAA_re$cor_map
    }
  }
  #map$trans_NAA_rho[which(!is.na(map$trans_NAA_rho))] <- 1:sum(!is.na(map$trans_NAA_rho))
  map$trans_NAA_rho <- factor(map$trans_NAA_rho)
  map$log_NAA[which(!is.na(map$log_NAA))] <- 1:sum(!is.na(map$log_NAA))
  map$log_NAA <- factor(map$log_NAA)
  map$log_NAA_sigma <- factor(map$log_NAA_sigma)

  if(any(data$recruit_model > 2 & data$NAA_re_model == 0)) input$log$NAA <- c(input$log$NAA, "NOTE: SCAA model specified, yearly recruitment deviations estimated as fixed effects. Stock-recruit function also specified. WHAM will fit the SCAA model but without estimating a stock-recruit function.
    This message will not appear if you set recruit_model = 2 (random about mean).")

  #set up recruitment
  if(!is.null(NAA_re$recruit_model)) {
    data$recruit_model[] = NAA_re$recruit_model #overrides recruit_model argument to wham::prepare_wham_input
    for(s in 1:data$n_stocks) if(data$recruit_model[s] > 1 & data$NAA_re_model[s] == 0) {
      stop("NAA_re$recruit_model[s] > 1 has been specified, but NAA_re$sigma[[s]] must either be 'rec' (random effects on recruitment only), 
      'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
      or a vector with length == n.ages specifying which sigma_a to use for each age.")
    }
  }
  par$mean_rec_pars = matrix(0, data$n_stocks, 2)
  map$mean_rec_pars = matrix(NA, data$n_stocks, 2)
  
  for(s in 1:data$n_stocks) {
    if(!is.null(NAA_re$recruit_pars[[s]])){
      if(data$recruit_model[s] == 2) par$mean_rec_pars[s,1] = log(NAA_re$recruit_pars[[s]][1])
      if(data$recruit_model[s] %in% 3:4) par$mean_rec_pars[s,1:2] = log(NAA_re$recruit_pars[[s]][1:2])
    } else{
      if(data$recruit_model[s]==2) {
        if(!is.null(asap3)) par$mean_rec_pars[s,1] = log(asap3[[s]]$N1_ini[1]) # initialize R0 at initial age-1
        else par$mean_rec_pars[s,1] = 10
      }  
      if(data$recruit_model[s]==4) par$mean_rec_pars[s,2] = -10
    }
    if(data$NAA_re_model[s] > 0){
      if(data$recruit_model[s]==2) map$mean_rec_pars[s,1] = s
      if(data$recruit_model[s]>2) map$mean_rec_pars[s,] = 1:2 + (s-1) * 2
    }
  }
  map$mean_rec_pars = factor(map$mean_rec_pars)

  input$data = data
  input$par = par
  input$map = map
	if(length(input$log$NAA))	input$log$NAA <- c("NAA: \n", input$log$NAA)
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	#input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

  if(!is_internal_call()) cat(unlist(input$log$NAA, recursive=T))
	#set any parameters as random effects
  #print(sort(names(input$data)))
	input$random = NULL
	input = set_random(input)
  input$options$NAA_re <- NAA_re
  return(input)
}