#' Specify model and parameter configuration for movement when input$data$n_regions > 1
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param move (optional) list specifying movement options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{move} specifies estimation options for movement.
#' If \code{NULL}, no movement will occur. If there are multiple regions, each stock will be modeled separately in different regions without movement. 
#' \code{move} is a list with the following entries:
#'   \describe{
#'     \item{$stock_move}{length = n_stocks, T/F whether each stock can move. If not provided then movement will be defined below for all stocks.}
#'     \item{$separable}{length = n_stocks, T/F whether movement should be modeled separably from mortality or both occuring simultaneously.}
#'     \item{$mean_model}{matrix (n_regions x (n_regions-1)): model options for fixed effects (mean and possibly variance) for each movement parameter are:
#'                    \describe{
#'                      \item{"none"}{(default) no movement between regions.}
#'                      \item{"constant"}{estimate a single movement rate to each region shared across all stocks, seasons, ages, years}
#'                      \item{"season"}{estimate movement rates  to each region for each season shared across all stocks, ages, years}
#'                      \item{"stock_constant"}{estimate a movement rate for each stock to each region shared across all seasons, ages, years}
#'                      \item{"stock_season"}{estimate a movement rate for each stock each season to each region shared across all ages, years}
#'                    }
#'     }
#'     \item{$age_re}{matrix (n_regions x (n_regions-1)): options for age random effects (for each mean parameter defined in \code{move$mean_model}):
#'                    \describe{
#'                      \item{"none"}{(default) no movement rate random effects by age.}
#'                      \item{"iid"}{independent movement rate random effects by age.}
#'                      \item{"ar1"}{allow first order autoregressive correlation of movement rate random effects by age.}
#'                    }
#'     }
#'     \item{$year_re}{matrix (n_regions x (n_regions-1)): options for yearly random effects (for each mean parameter defined in \code{move$mean_model}):
#'                    \describe{
#'                      \item{"none"}{(default) no movement rate random effects by year.}
#'                      \item{"iid"}{independent movement rate random effects by year.}
#'                      \item{"ar1"}{allow first order autoregressive correlation of movement rate random effects by year.}
#'                    }
#'     }
#'     \item{$prior_sigma}{array (n_stocks x n_seasons x n_regions x n_regions - 1) of sd parameters for normal priors on mean movement parameters on transformed scale (-Inf,Inf)}
#'     \item{$use_prior}{array (n_stocks x n_seasons x n_regions x n_regions - 1) 0/1 indicator whether to include prior for mean movement parameters in joint log-likelihood.}
#'     \item{$can_move}{array (n_stocks x n_seasons x n_regions x n_regions) 0/1 indicator whether movement can occur from one region to another.}
#'     \item{$must_move}{array (n_stocks x n_seasons x n_regions) 0/1 indicator whether movement from region must occur.}
#'     \item{$mean_vals}{array (n_stocks x n_seasons x n_regions x n_regions-1) of initial movement rate parameters *from* each region. Usage depends on \code{move$mean_model}.}
#'     \item{$sigma_vals}{array (n_stocks x n_seasons x n_regions x n_regions -1) of initial standard deviations to use for random effects. Usage depends on \code{move$age_re} and \code{move$year_re}.}
#'     \item{$cor_vals}{array (n_stocks x n_seasons x n_regions x n_regions - 1x 2) of initial correlation values to use for random effects. Usage depends on \code{move$age_re} and \code{move$year_re}.
#'        cor_vals[,,,,1] is for correlation with age, and cor_vals[,,,,2] is for correlation with year.}
#'   }
#'
#' @export
set_move = function(input, move)
{
  data = input$data
  par = input$par
  map = input$map
  input$log$move <- list()

  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mu_repars", "mu_re"))]

  data$can_move = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions))
  if(!is.null(move$can_move)) {
    data$can_move[] = move$can_move
    if(sum(data$can_move) == 0) input$log$move <- c(input$log$move, "All of move$can_move = 0, so model assumes no movement for any stocks.\n")
  }

  data$use_mu_prior = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  data$mig_type = rep(0,data$n_stocks)
  data$mu_model <- matrix(1, data$n_regions, data$n_regions-1)
  data$trans_mu_prior_sigma = array(0.1, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  data$must_move = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions))
  par$mu_prior_re = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  par$trans_mu = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  par$mu_re = array(0, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_years_model, data$n_regions, data$n_regions-1))
  par$mu_repars = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1, 3))
  map$mu_prior_re = array(NA, dim = dim(par$mu_prior_re))
  map$trans_mu = array(NA, dim = dim(par$trans_mu))
  map$mu_re = array(NA, dim = dim(par$mu_re))
  map$mu_repars = array(NA, dim = dim(par$mu_repars))

  if(!is.null(move$must_move)) data$must_move[] = move$must_move
  if(data$n_stocks>1) for(s in 1:data$n_stocks) if(sum(data$can_move[s,,,]) == 0) input$log$move <- c(input$log$move, paste0("Model assumes no movement for stock ", s, ".\n"))

  if(data$n_regions ==1 | sum(data$can_move)==0){
    map$trans_mu = factor(map$trans_mu)
    map$mu_repars = factor(map$mu_repars)
    map$mu_prior_re = factor(map$mu_prior_re)
    map$mu_re = factor(map$mu_re)
    input$data = data
    input$par = par
    input$map = map
    return(input)
  }

  if(!is.null(move$use_prior)) data$use_mu_prior[] = move$use_prior
  mean_mods <- c("none","constant", "season")
  mean_mods <- c(mean_mods, paste0("stock_", mean_mods[2:3]))
  if(is.null(move$mean_model)){
    if(is.null(move$can_move) | data$n_regions == 1) {
      move$mean_model <- matrix("none", data$n_regions, data$n_regions-1)
    } else{
      if(sum(move$can_move)>0) move$mean_model <- matrix("constant", data$n_regions, data$n_regions-1)
      else {
        if(data$n_regions==1) move$mean_model <- matrix("none", data$n_regions, data$n_regions-1)
        else move$mean_model <- matrix("constant", data$n_regions, data$n_regions-1)
      }
    }
    input$log$move <- c(input$log$move, paste0("\n move$mean_model was not specified and set to ", move$mean_model[1], " for all movement parameters based on data$n_regions and move$can_move if provided. \n"))
  }
  if(!is.matrix(move$mean_model)) stop(paste0("move$mean_model must be a n_regions x n_region-1 matrix filled with one of ", paste0(mean_mods, collapse = ", "), " when provided."))
  if(any(!move$mean_model %in% mean_mods)) stop(paste0("elements of move$mean_model must each be one of ", paste0(mean_mods, collapse = ", "), " when provided."))
  if(all(move$mean_model == "none")) {
    if(data$n_regions>1) input$log$move <- c(input$log$move, "all move$mean_model = 'none', so model assumes no movement for any stocks.\n")
    data$can_move[] <- 0
    map$trans_mu = factor(map$trans_mu)
    map$mu_repars = factor(map$mu_repars)
    map$mu_prior_re = factor(map$mu_prior_re)
    map$mu_re = factor(map$mu_re)
    input$data = data
    input$par = par
    input$map = map
    return(input)
  } else {
    if(sum(data$can_move)==0){ 
      input$log$move <- c(input$log$move, "move$model is not 'none', but all data$can_move = 0, so model assumes no movement for any stock.\n")
      map$trans_mu = factor(map$trans_mu)
      map$mu_repars = factor(map$mu_repars)
      map$mu_prior_re = factor(map$mu_prior_re)
      map$mu_re = factor(map$mu_re)
      input$data = data
      input$par = par
      input$map = map
      return(input)
    }
  }
  
  #if we've gotten this far then some parameters will be estimated
  #define data$mu_model and map for mu_re and mu_repars
  data$mu_model[] = 1 #this value needs to be between 1-8 after below
  data$mu_model[which(move$mean_model == "constant")] = 1 #already defined
  data$mu_model[which(move$mean_model == "stock_constant")] = 5 #already defined
  data$mu_model[which(move$mean_model == "season")] = 9 #already defined
  data$mu_model[which(move$mean_model == "stock_season")] = 13 #already defined
  re_mods <- c("none","iid", "ar1")
  if(is.null(move$age_re)) move$age_re <- matrix("none", data$n_regions, data$n_regions -1)
  if(!is.matrix(move$age_re)) stop(paste0("move$age_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  if(!all(dim(move$age_re) == data$n_regions - 0:1)) stop(paste0("move$age_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  if(is.null(move$year_re)) move$year_re = matrix("none", data$n_regions, data$n_regions -1)
  if(!is.matrix(move$year_re)) stop(paste0("move$year_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  if(!all(dim(move$year_re) == data$n_regions - 0:1)) stop(paste0("move$year_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  ind <- which(move$age_re != "none" & move$year_re == "none")
  data$mu_model[ind] <- data$mu_model[ind] + 1
  ind <- which(move$age_re == "none" & move$year_re != "none")
  data$mu_model[ind] <- data$mu_model[ind] + 2
  ind <- which(move$age_re != "none" & move$year_re != "none")
  data$mu_model[ind] <- data$mu_model[ind] + 3
  
  if(is.null(move$can_move)) for(r in 1:data$n_regions) {
    for(rr in 1:data$n_regions)
      rrless <- rr
      if(rr>r) rrless <- rr-1
      if(move$mean_model[r,rr] != "none") if(r!=rr){
        data$can_move[,,r,rr] <- 1
      }
  }

  use_mu_prior <- data$use_mu_prior
  map$trans_mu = array(NA, dim = dim(par$trans_mu)) 
  if(any(data$use_mu_prior>0)) {
    map$mu_prior_re <- array(NA, dim = dim(par$mu_prior_re))
    use_mu_prior[] = 0 #need to set use_mu_prior based on mean_model
  }
  
  #define map for mu_re and mu_repars
  i <- re_i <- 1
  for(r in 1:data$n_regions) {
    k <- 1
    for(rr in 1:data$n_regions) if(rr!= r) {
      if(data$mu_model[r,k] %in% 2:4) {
        if(sum(data$can_move[,,r,rr])>0) {
          map$mu_repars[,,r,k,1] <- i
          i <- i + 1
          if(data$mu_model[r,k] %in% c(2,4) & move$age_re[r,k] == "ar1") {
            map$mu_repars[,,r,k,2] <- i
            i <- i + 1
          }
          if(data$mu_model[r,k] %in% c(3,4) & move$year_re[r,k] == "ar1") {
            map$mu_repars[,,r,k,3] <- i
            i <- i + 1
          }
          for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) { #mu_re are constant across stock, season
            if(data$mu_model[r,k] == 2) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
            if(data$mu_model[r,k] == 3) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
            if(data$mu_model[r,k] == 4) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
          }
          re_i <- re_i + 1
        }
      }
      if(data$mu_model[r,k] %in% 6:8) {
        for(s in 1:data$n_stocks) {
          if(sum(data$can_move[s,,r,rr])>0) {
            map$mu_repars[s,,r,k,1] <- i
            i <- i + 1
            if(data$mu_model[r,k] %in% c(6,8) & move$age_re[r,k] == "ar1") {
              map$mu_repars[s,,r,k,2] <- i
              i <- i + 1
            }
            if(data$mu_model[r,k] %in% c(7,8) & move$year_re[r,k] == "ar1") {
              map$mu_repars[s,,r,k,3] <- i
              i <- i + 1
            }
            for(t in 1:data$n_seasons) { #mu_re are constant across season
              if(data$mu_model[r,k] == 6) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
              if(data$mu_model[r,k] == 7) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
              if(data$mu_model[r,k] == 8) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
            }
            re_i <- re_i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 10:12){ #season
        for(t in 1:data$n_seasons) {
          if(sum(data$can_move[,t,r,rr])>0) {
            map$mu_repars[,t,r,k,1] <- i
            i <- i + 1
            if(data$mu_model[r,k] %in% c(10,12) & move$age_re[r,k] == "ar1") {
              map$mu_repars[,t,r,k,2] <- i
              i <- i + 1
            }
            if(data$mu_model[r,k] %in% c(10,12) & move$year_re[r,k] == "ar1") {
              map$mu_repars[,t,r,k,3] <- i
              i <- i + 1
            }
            for(s in 1:data$n_stocks) { #mu_re are constant across stock
              if(data$mu_model[r,k] == 10) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
              if(data$mu_model[r,k] == 11) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
              if(data$mu_model[r,k] == 12) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
            }
            re_i <- re_i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 14:16){ #stock,season
        for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) {
          if(data$can_move[s,t,r,rr]) {
            map$mu_repars[s,t,r,k,1] <- i
            i <- i + 1
            if(data$mu_model[r,k] %in% c(14,16) & move$age_re[r,k] == "ar1") {
              map$mu_repars[s,t,r,k,2] <- i
              i <- i + 1
            }
            if(data$mu_model[r,k] %in% c(15,16) & move$year_re[r,k] == "ar1") {
              map$mu_repars[s,t,r,k,3] <- i
              i <- i + 1
            } #mu_re are different for each stock, season
            if(data$mu_model[r,k] == 14) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
            if(data$mu_model[r,k] == 15) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
            if(data$mu_model[r,k] == 16) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
            re_i <- re_i + 1
          }
        }
      }
      k <- k + 1
    }
  }
  # if(data$mu_model %in% 6:8){ #stock
  #   i <- re_i <- 1
  #   for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[s,,r,rr])>0) {
  #         map$mu_repars[s,,r,k,1] <- i
  #         i <- i + 1
  #         if(data$mu_model %in% c(6,8)) {
  #           map$mu_repars[s,,r,k,2] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model %in% c(7,8)) {
  #           map$mu_repars[s,,r,k,3] <- i
  #           i <- i + 1
  #         }
  #         for(t in 1:data$n_seasons) {
  #           if(data$mu_model == 6) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
  #           if(data$mu_model == 7) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
  #           if(data$mu_model == 8) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
  #         }
  #         re_i <- re_i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 10:12){ #season
  #   i <- re_i <- 1
  #   for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[,t,r,rr])>0) {
  #         map$mu_repars[,t,r,k,1] <- i
  #         i <- i + 1
  #         if(data$mu_model %in% c(10,12)) {
  #           map$mu_repars[,t,r,k,2] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model %in% c(10,12)) {
  #           map$mu_repars[,t,r,k,3] <- i
  #           i <- i + 1
  #         }
  #         for(s in 1:data$n_stocks) {
  #           if(data$mu_model == 10) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
  #           if(data$mu_model == 11) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
  #           if(data$mu_model == 12) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
  #         }
  #         re_i <- re_i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 14:16){ #stock,season
  #   i <- re_i <- 1
  #   for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(data$can_move[s,t,r,rr]) {
  #         map$mu_repars[s,t,r,k,1] <- i
  #         i <- i + 1
  #         if(data$mu_model %in% c(14,16)) {
  #           map$mu_repars[s,t,r,k,2] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model %in% c(15,16)) {
  #           map$mu_repars[s,t,r,k,3] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model == 14) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
  #         if(data$mu_model == 15) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
  #         if(data$mu_model == 16) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
  #         re_i <- re_i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }

  #define map and "use" elements for mean mu (trans_mu) parameters and prior-based RE (mu_prior_re) parameters 
  i <- 1
  for(r in 1:data$n_regions) {
    k <- 1
    for(rr in 1:data$n_regions) if(rr!= r) {
      if(data$mu_model[r,k] %in% 1:4) { #constant
        if(sum(data$can_move[,,r,rr])>0) {
          if(sum(data$use_mu_prior[,,r,k])>0) { #mean mu is a random effect with defined prior
            use_mu_prior[1,1,r,k] <- 1 #only evaluate the likelihood once
            map$mu_prior_re[,,r,k] <- i #all values the same
          } else {
            map$trans_mu[,,r,k] <- i #all values the same
          }
          i <- i + 1
        }
      }
      if(data$mu_model[r,k] %in% 5:8) { #stock
        for(s in 1:data$n_stocks) {
          if(sum(data$can_move[s,,r,rr])>0) {
            #mean mu is a random effect with defined prior
            if(sum(data$use_mu_prior[s,,r,k])>0) {
              use_mu_prior[s,1,r,k] <- 1 #only evaluate the likelihood once
              map$mu_prior_re[s,,r,k] <- i #all values the same
            } else {
              map$trans_mu[s,,r,k] <- i #all values the same
            }
            i <- i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 9:12) { #season
        for(t in 1:data$n_seasons) {
          if(sum(data$can_move[,t,r,rr])>0) {
            if(sum(data$use_mu_prior[,t,r,k])>0) { #mean mu is a random effect with defined prior
              use_mu_prior[1,t,r,k] <- 1 #only evaluate the likelihood once
              map$mu_prior_re[,t,r,k] <- i #all values the same
            } else {
              map$trans_mu[,t,r,k] <- i #all values the same
            }
            i <- i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 13:16) { #stock_season
        for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) {
          if(sum(data$can_move[s,t,r,rr])>0) {
            if(sum(data$use_mu_prior[s,t,r,k])>0) { #mean mu is a random effect with defined prior
              use_mu_prior[s,t,r,k] <- 1 #only evaluate the likelihood once
              map$mu_prior_re[s,t,r,k] <- i #all values the same
            } else {
              map$trans_mu[s,t,r,k] <- i #all values the same
            }
            i <- i + 1
          }
        }
      }
      k <- k + 1
    }
  }
  # if(data$mu_model %in% 5:8) { #stock
  #   i <- 1
  #   for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[s,,r,rr])>0) {
  #         #mean mu is a random effect with defined prior
  #         if(sum(data$use_mu_prior[s,,r,k])>0) {
  #           use_mu_prior[s,1,r,k] <- 1 #only evaluate the likelihood once
  #           map$mu_prior_re[s,,r,k] <- i #all values the same
  #         } else {
  #           map$trans_mu[s,,r,k] <- i #all values the same
  #         }
  #         i <- i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 9:12) { #season
  #   i <- 1
  #   for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[,t,r,rr])>0) {
  #         if(sum(data$use_mu_prior[,t,r,k])>0) { #mean mu is a random effect with defined prior
  #           use_mu_prior[1,t,r,k] <- 1 #only evaluate the likelihood once
  #           map$mu_prior_re[,t,r,k] <- i #all values the same
  #         } else {
  #           map$trans_mu[,t,r,k] <- i #all values the same
  #         }
  #         i <- i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 13:16) { #stock_season
  #   i <- 1
  #   for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[s,t,r,rr])>0) {
  #         if(sum(data$use_mu_prior[s,t,r,k])>0) { #mean mu is a random effect with defined prior
  #           use_mu_prior[s,t,r,k] <- 1 #only evaluate the likelihood once
  #           map$mu_prior_re[s,t,r,k] <- i #all values the same
  #         } else {
  #           map$trans_mu[s,t,r,k] <- i #all values the same
  #         }
  #         i <- i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  data$use_mu_prior <- use_mu_prior 

  if(!is.null(move$separable)) {
    if(any(!move$separable)) {
      input$log$move <- c(input$log$move, "NOTE: movement and mortality must be assumed to occur seequentially for now.")
      move$separable[] <- TRUE
    }
    data$mig_type[] = as.integer(!move$separable)
  }
  if(data$n_regions>1) {
    if(length(unique(move$separable))==1) {
      if(move$separable[1] == 0) input$log$move <- c(input$log$move, "movement and mortality will be assumed to occur simultaneously within each season.\n")
      if(move$separable[1] == 1) input$log$move <- c(input$log$move, "movement and mortality will be assumed to occur sequentially within each season.\n")
    }
    else{
      for(s in 1:data$n_stocks){
        if(data$mig_type[s] == 1) input$log$move <- c(input$log$move, paste0("for stock ", s, ", movement and mortality will be assumed to occur simultaneously within each season.\n"))
        if(data$mig_type[s] == 0) input$log$move <- c(input$log$move, paste0("for stock ", s, ", movement and mortality will be assumed to occur sequentially within each season.\n"))
      }
    }
  }
  if(!is.null(move$prior_sigma)) data$trans_mu_prior_sigma[] = move$prior_sigma

  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s
  k <- 1
  if(!is.null(move$mean_vals)){
    for(s in 1:data$n_stocks) {
      if(data$mig_type[s]) { #simultaneous
        for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
          #ind <- which(data$can_move[s,t,r,-r]==1) #not necessary because trans_mu only used if can_move == 1
          par$trans_mu[s,,t,r,]<- log(move$mean_vals[s,t,r,])
        }
      } else { #separable
        for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
          #print(data$can_move[s,t,r,-r])
          #ind <- which(data$can_move[s,t,r,-r]==1) #not necessary because trans_mu only used if can_move == 1
          p <- move$mean_vals[s,t,r,] #n_r - 1 long!
          par$trans_mu[s,t,r,] <- log(p) - log(1-sum(p)) #additive
        }
      }
    }
  }
  if(!is.null(move$sigma_vals)){
    par$mu_repars[,,,,1] = log(move$sigma_vals)
  }
  if(!is.null(move$cor_vals)){
    par$mu_repars[,,,,2] = inv_trans_rho(move$cor_vals[,,,,1])
    par$mu_repars[,,,,3] = inv_trans_rho(move$cor_vals[,,,,2])
  }
  map$trans_mu = factor(map$trans_mu)
  map$mu_repars = factor(map$mu_repars)
  map$mu_prior_re = factor(map$mu_prior_re)
  map$mu_re = factor(map$mu_re)

  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$move)) input$log$move <- c("Movement: \n", input$log$move)
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = set_random(input)
  input$options$move <- move
  if(!is_internal_call()) cat(unlist(input$log$move, recursive=T))
  return(input)

}