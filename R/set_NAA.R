#' Specify model and parameter configuration for numbers at age
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param NAA_re (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' 
#' If \code{NAA_re = NULL}, a traditional statistical catch-at-age model is fit (NAA = pred_NAA for all ages, deterministic). Otherwise,
#' \code{NAA_re} specifies numbers-at-age configuration. It is a list with the following possible entries:
#'   \describe{
#'     \item{$sigma}{Which ages allow deviations from pred_NAA? Common options are specified with the strings:
#'                    \describe{
#'                      \item{"rec"}{Random effects on recruitment (deviations), all other ages deterministic}
#'                      \item{"rec+1"}{"Full state space" model with 2 estimated \code{sigma_a}, one for recruitment and one shared among other ages}
#'                    }
#'                   Alternatively, you can specify a more complex structure by entering a vector with length = n.ages, where each entry points to the
#'                   NAA_sigma to use for that age. E.g. c(1,2,2,3,3,3) will estimate 3 \code{sigma_a}, with recruitment (age-1) deviations having their
#'                   own \code{sigma_R}, ages 2-3 sharing \code{sigma_2}, and ages 4-6 sharing \code{sigma_3}. If length = 1, assumptions will be applied to all stocks.
#'                  }
#'     \item{$sigma_vals}{Initial standard deviation values to use for the NAA deviations. Values are not used if recruit_model = 1 and NAA_re$sigma is
#'                  not specifed. Otherwise when \code{NAA_re$sigma} =
#'                  \describe{
#'                    \item{"rec"}{must be a single value.}
#'                    \item{"rec+1"}{2 values must be specified. First is for the first age class (recruits), second is for all other ages.}
#'                    \item{vector of values (length = number of age classes)}{either 1 value or the number of values is equal to the number of unique values provided to \code{NAA_re$sigma}.}
#'                  }
#'                }
#'     \item{$cor}{Correlation structure for the NAA deviations. Options are:
#'                  \describe{
#'                    \item{"iid"}{NAA deviations vary by year and age, but uncorrelated.}
#'                    \item{"ar1_a"}{NAA deviations correlated by age (AR1).}
#'                    \item{"ar1_y"}{NAA deviations correlated by year (AR1).}
#'                    \item{"2dar1"}{NAA deviations correlated by year and age (2D AR1).}
#'                  }
#'                }
#'     \item{$cor_vals}{Initial correlation values to use for the NAA deviations. If unspecified all initial values are 0. When \code{NAA_re$cor} = 
#'                  \describe{
#'                    \item{"iid"}{values are not used.}
#'                    \item{"ar1_a" or "ar1_y"}{cor_vals must be a single value.}
#'                    \item{"2dar1"}{2 values must be specified. First is for "age", second is for "year".}
#'                  }
#'                }
#'     \item{$N1_model}{Integer vector (n_stocks) determining which way to model the initial numbers at age:
#'       \describe{
#'          \item{"age-specific-fe"}{(default) age- and region-specific fixed effects parameters}
#'          \item{"equilibrium"}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'          \item{"iid-re"}{(default) age- and region-specific iid random effects parameters. 2 parameters: mean and sd for log NAA}
#'          \item{"ar1-re"}{(default) age- and region-specific random effects parameters. 3 parameters: mean and sd, and cor for log NAA}
#'       }
#'     }
#'     \item{$N1_pars}{if N1_model = 0, then these would be the initial values to use for abundance at age in the first year. If N1_model = 1, This would be the
#'        initial numbers in the first age class and the equilibrium fishing mortality rate generating the rest of the numbers at age in the first year.
#'     }
#'     \item{$recruit_model}{Integer determining how to model recruitment. Overrides \code{recruit_model} argument to \code{prepare_wham_input}. Must make sure \code{NAA_re$sigma}, \code{NAA_re$cor}
#'        and \code{ecov} are properly specified.
#'       \describe{
#'           \item{1}{SCAA, estimating all recruitements as fixed effects or a random walk if NAA_re$sigma specified}
#'           \item{2}{estimating a mean recruitment with yearly recruitements as random effects}
#'           \item{3}{Beverton-Holt stock-recruitment with yearly recruitements as random effects}
#'           \item{4}{Ricker stock-recruitment with yearly recruitements as random effects}
#'       }
#'     }
#'     \item{$use_steepness}{T/F determining whether to use a steepness parameterization for a stock-recruit relationship. Only used if recruit_model>2}.
#'     \item{$recruit_pars}{vector of initial parameters for recruitment model. If use_steepness=F, parameters are "alpha" and "beta"
#'        otherwise they are steepness and R0.
#'     }
#'   }

set_NAA = function(input, NAA_re=NULL)
{

  data = input$data
  par = input$par
  map = input$map
  asap3 = input$asap3
  inv_trans_rho <- function(rho, s = 2) (log(rho+1) - log(1-rho))/s # 0.5 because needed transformation on cpp side is unusual.
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mean_rec_pars", "log_N1", "log_NAA_sigma", "trans_NAA_rho","logR_proj", "log_NAA"))]
  
  #if(is.null(input$asap3)) asap3 = NULL
  #else asap3 = input$asap3


  #set up initial NAA
  #0: just age-specific numbers at age
  data$N1_model = rep(0, data$n_stocks)
  par$log_N1 = array(0,dim = c(data$n_stocks,data$n_regions,data$n_ages))
  map$log_N1 = array(NA,dim = c(data$n_stocks,data$n_regions,data$n_ages))
  par$N1_repars = array(0,dim = c(data$n_stocks,data$n_regions,3))
  map$N1_repars = array(NA,dim = c(data$n_stocks,data$n_regions,3))
  if(!is.null(NAA_re$N1_model)) {
    options = c("age-specific-fe", "equilibrium","iid-re", "ar1-re")
    k = 1
    for(s in 1:data$n_stocks) {
      if(!(NAA_re$N1_model[s] %in% options)) stop("NAA_re$N1_model must all be 'age-specific-fe', 'equilibrium', 'iid-fe' or 'ar1-fe'.")
      if(NAA_re$N1_model[s] == options[1]) data$N1_model[s] = 0
      if(NAA_re$N1_model[s] == options[2]) data$N1_model[s] = 1
      if(NAA_re$N1_model[s] %in% options[3:4]) {
        data$N1_model[s] = 2
        map$N1_repars[s,,1] = k
        map$N1_repars[s,,2] = k + 1
        k = k + 2
        if(NAA_re$N1_model[s] == options[4]){ #ar1 with age
          map$N1_repars[s,,2] = k
          k = k + 1
        }
      }
    }
  }
  map$N1_repars = factor(map$N1_repars)
  k = 1
  init_NAA = log(exp(10)*exp(-(0:(data$n_ages-1))*0.2))
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
  if(!is.null(NAA_re[["N1"]])){
    par$log_N1[] = log(NAA_re$N1)
  }
  map$log_N1 = factor(map$log_N1)

  # NAA_re options for beyond year 1
  # default = SCAA for each stock
  data$NAA_re_model = rep(0, data$n_stocks)
  par$log_NAA_sigma = matrix(0,data$n_stocks, data$n_ages)
  map$log_NAA_sigma = matrix(NA, data$n_stocks, data$n_ages)
  par$trans_NAA_rho = matrix(0,data$n_stocks, 2)
  map$trans_NAA_rho = matrix(NA, data$n_stocks, 2)
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
    k = 0
    if(!length(NAA_re$sigma) %in% c(1,data$n_stocks)) stop("NAA_re$sigma length must be 1 or equal to the number of stocks.")
    if(length(NAA_re$sigma == 1) & data$n_stocks>1) {
      cat("\n Same NAA_re$sigma being used for all stocks.\n")
      NAA_re$sigma = rep(list(NAA_re$sigma), data$n_stocks)
    }
    for(s in 1:data$n_stocks) if(NAA_re$sigma[[s]][1] == "rec"){
      data$NAA_re_model[s] = 1
      map$log_NAA_sigma[s,1] = k
      #below is already done above for SCAA
      #map$log_NAA[s,data$spawn_regions[s],,1] = 1 #change to unique values later
      k = k + 1
      #data$n_NAA_sigma <- 1
      #data$NAA_sigma_pointers <- rep(1,data$n_ages)
    } else {
      if(NAA_re$sigma[[s]][1] == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
        data$NAA_re_model[s] = 2
        map$log_NAA_sigma[s,] = c(k, rep(k+1, data$n_ages-1))
        for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] = 1 #change to unique values later
        k = k + 2
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
            tmp = unique(NAA_re_sigma[[s]])
            ind = 1:length(tmp)
            ind = ind[match(NAA_re_sigma[[s]],tmp)] - 1
            map$log_NAA_sigma[s,] <- k + ind
            for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] = 1 #change to unique values later
            k <- max(k + ind, na.rm=T)
            #data$n_NAA_sigma <- max(unique(NAA_re$sigma), na.rm=T)
            #data$NAA_sigma_pointers <- NAA_re$sigma
          }
        #}
      }
    }
    if(!is.null(NAA_re$sigma_vals[[s]])) {
      if(!(length(NAA_re$sigma_vals[[s]]) %in% c(1,data$n_ages))) stop(paste0("length of NAA_re$sigma_vals[[s]] must be 1 or ", data$n_ages, "."))
      par$log_NAA_sigma[s,] <- log(NAA_re$sigma_vals[[s]])
    }
    
    if(length(NAA_re$cor == 1) & data$n_stocks>1) {
      cat("\n Same NAA_re$cor being used for all stocks.\n")
      NAA_re$cor = rep(list(NAA_re$cor), data$n_stocks)
    }
    if(!is.null(NAA_re$cor[[s]])){
      if(!NAA_re$cor[[s]] %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor[[s]] must be one of 'iid','ar1_a','ar1_y','2dar1'")
      if(NAA_re$cor[[s]] == "ar1_a") map$trans_NAA_rho[s,1] <- 1
      if(NAA_re$cor[[s]] == "ar1_y") map$trans_NAA_rho[s,2] <- 1
      if(NAA_re$cor[[s]] == "2dar1") map$trans_NAA_rho[s,] <- 1
    } else {
      NAA_re$cor[[s]] <- 'iid'
    }
    if(!is.null(NAA_re$cor_vals[[s]])) {
      if(!length(NAA_re$cor_vals[[s]]) %in% 1:2) stop(paste0("length of NAA_re$cor_vals[[s]] is not consistent with other elements of NAA_re$cor."))
      if(length(NAA_re$cor_vals[[s]]) == 2) par$trans_NAA_rho[s,] <- inv_trans_rho(NAA_re$cor_vals[[s]])
      if(length(NAA_re$cor_vals[[s]]) == 1) {
        if(NAA_re$cor[[s]] == "ar1_a") {
          par$trans_NAA_rho[s,1] <- inv_trans_rho(NAA_re$cor_vals[[s]])
        }
        if(NAA_re$cor[[s]] == "ar1_y") {
          par$trans_NAA_rho[s,2] <- inv_trans_rho(NAA_re$cor_vals[[s]])
        }
      }
    }
  }
  map$trans_NAA_rho[which(!is.na(map$trans_NAA_rho))] <- 1:sum(!is.na(map$trans_NAA_rho))
  map$trans_NAA_rho <- factor(map$trans_NAA_rho)
  map$log_NAA[which(!is.na(map$log_NAA))] <- 1:sum(!is.na(map$log_NAA))
  map$log_NAA <- factor(map$log_NAA)
  map$log_NAA_sigma <- factor(map$log_NAA_sigma)

  if(any(data$recruit_model > 2 && data$NAA_re_model == 0)) warning("SCAA model specified, yearly recruitment deviations estimated as fixed effects. Stock-recruit function also specified. WHAM will fit the SCAA model but without estimating a stock-recruit function.
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

  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	#input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
  #print(sort(names(input$data)))
	input$random = NULL
	input = set_random(input)

  return(input)
}