#' Specify configuration for environmental covariates, effects on the population, and parameter values
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param ecov (optional) named list of environmental covariate data and parameters (see details)
#' 
#' \code{ecov} specifies any environmental covariate data and model as well as the effect on the population. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$label}{Name(s) of the environmental covariate(s). Used in printing.}
#'     \item{$mean}{Mean observations (matrix). number of years x number of covariates. Missing values = NA.}
#'     \item{$logsigma}{Configure observation standard errors. Options:
#'       \describe{
#'         \item{Matrix of \eqn{log} standard errors with same dimensions as \code{$mean}}{Specified values for each time step }
#'         \item{log standard errors for each covariate, numeric vector or matrix w/ dim 1 x n.ecov}{Specified value the same for all time steps}
#'         \item{estimation option (for all covariates). character string:}{
#'           \code{"est_1"}: Estimated, one value shared among time steps.
#'           \code{"est_re"}: Estimated value for each time step as random effects with two parameters (mean, var)}
#'         \item{list of two elements.}{
#'           First is the matrix of log standard errors or the vector of single values for each covariate as above. 
#'           Second is a character vector of estimation options (\code{NA}, \code{"est_1"},\code{"est_re"}) for each covariate. 
#'           For covariates with non-NA values, values in the first element are ignored.}
#'       }
#'     }
#'     \item{$year}{Years corresponding to observations (vector of same length as \code{$mean} and \code{$logsigma})}
#'     \item{$use_obs}{T/F (or 0/1) vector/matrix of the same dimension as \code{$mean} and \code{$logsigma}.
#'     Use the observation? Can be used to ignore subsets of the ecov without changing data files.}
#'     \item{$lag}{integer vector of offsets between the ecov observation/process and their affect on the stock.
#'     I.e. if SST in year \emph{t} affects recruitment in year \emph{t + 1}, set \code{lag = 1}. May also be a list (length=n_Ecov) of vectors (length = 2+n_indices) if multiple effects of one or more Ecovs are modeled.}
#'     \item{$process_model}{Process model for the ecov time-series. \code{"rw"} = random walk, \code{"ar1"} = 1st order autoregressive, \code{NA} = do not fit}
#'     \item{$process_mean_vals}{vector of (initial) mean values for the ecov time-series.}
#'     \item{$process_sig_vals}{vector of (initial) standard deviation values for the ecov time-series.}
#'     \item{$process_cor_vals}{vector of (initial) correlation values for the ecov time-series.}
#'     \item{$where}{Character vector for where each ecov affects the population. \code{"recruit"} = recruitment,
#'     \code{"M"} = natural mortality, \code{"q"} = catchability for indices, \code{"none"} = fit ecov process model(s) but without an effect on the
#'     population. May also be a list (element for each ecov) of character vectors ("none", "recruit", "M", and/or "q") so each ecov can multiple effects.}
#'     \item{$indices}{indices that each ecov affects. Must be a list (length = n_Ecov), where each element is a vector of indices (1:n_indices). Must be provided when any of \code{where} = "q"}
#'     \item{$link_model}{vector of (orthogonal polynomial order) models for effects of each ecov on the \code{$where} process. Options: "none", "linear" (default) or "poly-x"
#'     where x = 2, ... (e.g. "poly-2" specifies a quadratic model, \eqn{b_0 + b_1*ecov + b_2*ecov^2 + ...}). Or a list (length = n_Ecov) of character vectors (same options) for modeling
#'      effects on (1) recruitment, (2) M, and catchabilities for (3) index 1,...,n_indices (length of the vector is 2 + n_indices).}
#'     \item{$beta_vals}{a list (length = n_Ecov) of lists (length is 2 + n_indices) of matrices of beta (effect size) values (rows determined by \code{link_model},cols is n_ages) for 
#'      effects on (1) recruitment, (2) M, and catchabilities for (3) index 1,...,n_indices. Age-specific values only used for M.}
#'     \item{$ages}{Ages that each ecov affects. Must be a list of length n.ecov, where each element is a vector of ages. Default = list(1:n.ages) to affect all ages.}
#'     \item{$how}{How does the ecov affect the \code{$where} process? An integer vector (length = n_Ecov). If corresponding \code{$where} = "none", then this is ignored.
#'     These options are specific to the \code{$where} process.
#'       \describe{
#'         \item{Recruitment options (see \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}):}{
#'           \describe{
#'             \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'             \item{= 1}{"controlling" (dens-indep mortality)}
#'             \item{= 2}{"limiting" (carrying capacity, e.g. ecov determines amount of suitable habitat)}
#'             \item{= 3}{"lethal" (threshold, i.e. R --> 0 at some ecov value)}
#'             \item{= 4}{"masking" (metabolic/growth, decreases dR/dS)}
#'             \item{= 5}{"directive" (e.g. behavioral)}
#'           }}
#'         \item{M options:}{
#'           \describe{
#'             \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'             \item{= 1}{effect on mean M (shared across ages)}
#'           }}
#'         \item{Catchability options:}{
#'           \describe{
#'             \item{= 0}{none (but fit ecov time-series model in order to compare AIC)}
#'             \item{= 1}{effect on one or more catchabilities (see \code{$indices)})}
#'           }}
#'       }
#'     }
#'   }

set_ecov = function(input, ecov)
{
  data = input$data
  par = input$par
  map = input$map

  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("Ecov_obs_logsigma_re", "Ecov_obs_sigma_par", "Ecov_obs_log_sigma","Ecov_process_pars", "Ecov_re", "Ecov_beta"))]
  
  #define new dimensions for all effects a given ecov can have on assessment model
  #currently the order is recruitment, M, index 1, ..., n_indices
  n_effects = 2 + data$n_indices
  index_effects = 2+1:data$n_indices

  # --------------------------------------------------------------------------------
  # Environmental covariate data
  if(is.null(ecov)){
    data$Ecov_obs <- matrix(1, nrow=1, ncol=1)
    par.Ecov.obs.logsigma <- matrix(-2.3, nrow=1, ncol=1)
    map.Ecov.obs.logsigma <- factor(NA)
    par.Ecov.obs.logsigma.re <- matrix(-2.3, nrow=1, ncol=1)
    map.Ecov.obs.logsigma.re <- factor(NA)
    par.Ecov.obs.sigma.par <- matrix(0, nrow=1, ncol=1)
    map.Ecov.obs.sigma.par <- factor(NA)
    data$Ecov_obs_sigma_opt <- 1
    data$n_Ecov <- 1
    data$Ecov_use_obs <- matrix(0, nrow=1, ncol=1)
    data$Ecov_year <- matrix(0, nrow=1, ncol=1)
    data$year1_Ecov <- 0
    data$year1_model <- input$years[1]
    #data$Ecov_lag <- 0 #This is not used anywhere
    data$Ecov_model <- 0
    data$Ecov_where = matrix(0, data$n_Ecov, n_effects)
    data$Ecov_how <- 0
    #not used
    #data$Ecov_poly <- 1
    data$n_years_Ecov <- 1
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- matrix(0, data$n_Ecov, n_effects)
    data$Ecov_label <- "none"
    data$Ecov_use_re <- matrix(0, nrow=1, ncol=1)
    Ecov_poly <- matrix(1,data$n_Ecov, n_effects)
  } else {
    if(class(ecov$mean)[1] == "matrix") {data$Ecov_obs <- ecov$mean} else{
      warning("ecov$mean is not a matrix. Coercing to a matrix...")
      data$Ecov_obs <- as.matrix(ecov$mean)
    }
    if(class(ecov$use_obs)[1] == "matrix"){
      data$Ecov_use_obs <- ecov$use_obs
    } else{
      warning("ecov$use_obs is not a matrix with same dimensions as ecov$mean. Coercing to a matrix...")
      data$Ecov_use_obs <- as.matrix(as.integer(ecov$use_obs))
    }
    if(!identical(dim(data$Ecov_use_obs), dim(data$Ecov_obs))) stop("Dimensions of ecov$use_obs != dimensions of ecov$mean")

    # Handle Ecov sigma options
    data$n_Ecov <- dim(data$Ecov_obs)[2] # num Ecovs
    n_Ecov_obs <- dim(data$Ecov_obs)[1] # num Ecov obs
    data$Ecov_obs_sigma_opt = rep(1, data$n_Ecov) # Ecov sigma given, initialized at given values, not estimated by default

    if(length(ecov$year) != n_Ecov_obs) stop("ecov$year is not the same length as # rows in ecov$mean")
    data$Ecov_year <- as.numeric(ecov$year)
    data$year1_Ecov <- ecov$year[1]
    data$year1_model <- input$years[1]
    end_model <- tail(input$years,1)
    end_Ecov <- tail(ecov$year,1)
    if(length(ecov$label) == data$n_Ecov){
      data$Ecov_label <- ecov$label
    } else {
      warning("Number of Ecov labels not equal to number of Ecovs
              Setting Ecov labels = 'Ecov 1', 'Ecov 2', ...")
      data$Ecov_label = paste0("Ecov ",1:data$n_Ecov)
    }
    
    #default NAs for parameter matrix of observation standard deviations
    par.Ecov.obs.logsigma = matrix(NA, n_Ecov_obs, data$n_Ecov)

    logsig_more = list()
    if(is.list(ecov$logsigma)){
      n = length(ecov$logsigma)
      if(n>1) logsig_more = ecov$logsigma[[2]] #further elements than the second are ignored.
      ecov$logsigma = ecov$logsigma[[1]]
    }
    if(class(ecov$logsigma)[1] == "matrix"){
      par.Ecov.obs.logsigma <- ecov$logsigma
      if(!identical(dim(par.Ecov.obs.logsigma), dim(data$Ecov_obs))) stop("Dimensions of ecov$mean != dimensions of ecov$logsigma")
    }
    if(class(ecov$logsigma)[1] == 'numeric'){
      #data$Ecov_obs_sigma_opt[] = 1 #defined above
      print("ecov$logsigma is numeric. Coercing to a matrix...")
      if(length(ecov$logsigma) == data$n_Ecov) par.Ecov.obs.logsigma <- matrix(rep(ecov$logsigma, each=n_Ecov_obs), ncol=data$n_Ecov)
      if(length(ecov$logsigma) == n_Ecov_obs && data$n_Ecov == 1) par.Ecov.obs.logsigma <- matrix(ecov$logsigma, ncol=1)
      if(length(ecov$logsigma) != data$n_Ecov && length(ecov$logsigma) != n_Ecov_obs) stop("ecov$logsigma is numeric but length is not equal to # of ecovs or ecov observations")
    }

    #set up and check length of ecov$process_model
    if(length(ecov$process_model) == 1) ecov$process_model = rep(ecov$process_model, data$n_Ecov) #use the single value for all Ecovs
    if(length(ecov$process_model) != data$n_Ecov) stop("length of ecov$process_model must be either 1 or the number of Ecovs")
    
    for(i in 1:data$n_Ecov) {
      if(is.na(ecov$process_model[i]) & any(ecov$where[[i]] %in% c("recruit", "M", "q"))){ #ecov$how !=0){
      stop(paste0("ecov$process_model ", i, " is turned off (NA) but ecov$where includes 'M', 'recruit', or 'q'.
       Either 1) choose an ecov process model ('rw' or 'ar1'),
              2) turn off ecov (set ecov$where[i] = 'none' and ecov$process_model = NA),
           or 3) fit ecov but with no effect on population (ecov$where[i] = 'none', ecov$process_model[i] = 'rw' or 'ar1')."))
      }
    }

    #now over write ecov$logsigma with logsig_more if available because the fixed obs var matrices have been defined
    if(length(logsig_more)) ecov$logsigma = logsig_more
    
    #default values for: 
    #map of random effects
    map.Ecov.obs.logsigma.re <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
    #map of observation standard deviation
    map.Ecov.obs.logsigma <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
    #initial values variance parameter (fixed effects) for random effects 
    par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
    #map of variance parameter (fixed effects) for random effects
    map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars

    if(class(ecov$logsigma)[1] == 'character'){
      #check that estimation options are right
      if(!all(ecov$logsigma %in% c("est_1", "est_re", NA))){
        stop("ecov$logsigma or ecov$logsigma[[2]] is character and must be NA (do not estimate), 'est_1' (single variance parameter), or 'est_re' (iid re annual variance parameters)")
      }
      if(length(ecov$logsigma) == 1) ecov$logsigma = rep(ecov$logsigma, data$n_Ecov) #use the single value for all Ecovs
      #check length of estimation options
      if(length(ecov$logsigma) != data$n_Ecov) stop("length of ecov$logsigma when character must be either 1 or the number of Ecovs")

      for(i in 1:data$n_Ecov) {
        if(is.na(ecov$process_model[[i]])){
          #data$Ecov_obs_sigma_opt[i] == 1 # already defined above        
        }
        if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_1'){ # estimate 1 Ecov obs sigma for each Ecov
          data$Ecov_obs_sigma_opt[i] = 2
          par.Ecov.obs.logsigma[,i] <- -1.3 # matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov)
          map.Ecov.obs.logsigma[,i] <- i #matrix(rep(1:data$n_Ecov, each=n_Ecov_obs), ncol=data$n_Ecov)
        }
        #This option is not discussed anywhere and probably not useful. map for "est_1" could be modified to have estimate blocks of fixed effects
        # if(ecov$logsigma == 'est_fe'){ # estimate Ecov obs sigma for each Ecov obs
        #   data$Ecov_obs_sigma_opt[i] = 3
        #   par.Ecov.obs.logsigma <- matrix(-1.3, nrow=n_Ecov_obs, ncol=data$n_Ecov) # fixed effect inits
        #   map.Ecov.obs.logsigma <- matrix(1:(n_Ecov_obs*data$n_Ecov), nrow=n_Ecov_obs, ncol=data$n_Ecov) # est all
        #   par.Ecov.obs.sigma.par <- matrix(-1.3, nrow=2, ncol=data$n_Ecov)
        #   map.Ecov.obs.sigma.par <- matrix(NA, nrow=2, ncol=data$n_Ecov) # turn off RE pars
        # }
        if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_re'){
          data$Ecov_obs_sigma_opt[i] = 4
          map.Ecov.obs.logsigma[,i] <- NA #matrix(1:(n_Ecov_obs*data$n_Ecov), nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation of fixed effects
          par.Ecov.obs.sigma.par[,i] <- c(-1.3, -2.3) #matrix(c(rep(-1.3, data$n_Ecov), rep(-2.3, data$n_Ecov)), ncol=data$n_Ecov, byrow=TRUE) # random effect pars
          map.Ecov.obs.sigma.par[,i] <- max(0, map.Ecov.obs.sigma.par, na.rm =T) + 1:2 #matrix(1:(2*data$n_Ecov), nrow=2, ncol=data$n_Ecov)
        }
      }
    }

    # # check that Ecov year vector doesn't have missing gaps
    # pad Ecov if it starts after model year1 - max(lag)
    #make ecov$lag compatible with multiple effects
    if(is.null(ecov$lag)) stop("ecov$lag needs to be provided for each ecov")
    if(!is.list(ecov$lag)) ecov$lag = lapply(ecov$lag, function(x) rep(x,n_effects))
    max.lag = max(unlist(ecov$lag))
    if(data$year1_Ecov > data$year1_model - max.lag){
      print("one or more ecov does not start by model year 1 - max(lag). Padding ecov...")
      data$Ecov_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max.lag), ncol = data$n_Ecov), data$Ecov_obs)
      par.Ecov.obs.logsigma <- rbind(matrix(par.Ecov.obs.logsigma[1,], nrow = data$year1_Ecov-(data$year1_model-max.lag), ncol = data$n_Ecov, byrow=T), par.Ecov.obs.logsigma)
      map.Ecov.obs.logsigma <- rbind(matrix(NA, nrow = data$year1_Ecov-(data$year1_model-max.lag), ncol = data$n_Ecov), map.Ecov.obs.logsigma)
      data$Ecov_use_obs <- rbind(matrix(0, nrow = data$year1_Ecov-(data$year1_model-max.lag), ncol = data$n_Ecov), data$Ecov_use_obs)
      data$Ecov_year <- c(seq(data$year1_model - max.lag, data$year1_Ecov-1), data$Ecov_year)
      data$year1_Ecov <- data$year1_model - max.lag
    }

    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      print("Ecov last year is before model last year. Padding Ecov...")
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      par.Ecov.obs.logsigma <- rbind(par.Ecov.obs.logsigma, matrix(par.Ecov.obs.logsigma[NROW(par.Ecov.obs.logsigma),], nrow = end_model-end_Ecov, ncol = data$n_Ecov, byrow=T))
      map.Ecov.obs.logsigma <- rbind(map.Ecov.obs.logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_year <- c(data$Ecov_year, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- dim(data$Ecov_obs)[1] # num years Ecov to model (padded)
    data$Ecov_use_re <- matrix(1, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
    
    #map of random effects
    #do this now in case of padding 
    map.Ecov.obs.logsigma.re = matrix(NA, data$n_years_Ecov, data$n_Ecov)
    #initial values of random effects
    par.Ecov.obs.logsigma.re = matrix(0, data$n_years_Ecov, data$n_Ecov)
    for(i in 1:data$n_Ecov) if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_re') {
      map.Ecov.obs.logsigma.re[,i] = max(0, map.Ecov.obs.logsigma.re, na.rm=T) + 1:data$n_years_Ecov
      par.Ecov.obs.logsigma.re[,i] <- par.Ecov.obs.logsigma[,i] # random effect initialize at values in matrix provided
    }
    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    data$ind_Ecov_out_start <- data$ind_Ecov_out_end <- matrix(NA, data$n_Ecov, n_effects)
    for(i in 1:data$n_Ecov) for(j in 1:n_effects) {
      data$ind_Ecov_out_start[i,j] <- which(data$Ecov_year==data$year1_model)-ecov$lag[[i]][j]-1 # -1 is for cpp indexing
      data$ind_Ecov_out_end[i,j] <- which(data$Ecov_year==end_model)-ecov$lag[[i]][j]-1 # -1 is for cpp indexing
    }

    if(!identical(length(ecov$lag), length(ecov$label), data$n_Ecov)) stop("Length of ecov$lag and ecov$label not equal to # Ecov")
    ecov$lag <- t(matrix(unlist(ecov$lag), n_effects, data$n_Ecov)) #just used on R side (e.g., plotting)
    
    data$Ecov_where = matrix(0, data$n_Ecov, n_effects)
    if(is.character(ecov$where)) ecov$where = as.list(ecov$where) #put in new allowed format so each ecov can have mulitple effects.
    for(i in 1:data$n_Ecov) {
      if(any(ecov$where[[i]] == "none")) if(!is.null(ecov$how)) if(ecov$how[i] != 0){
        warning(paste0("ecov$where[[", i, "]] is 'none', but ecov$how[", i, "] is not 0. Setting this Ecov to not have any effects"))
        ecov$how[i] = 0
      }
    }

    if(any(sapply(ecov$where, function(x) any(x == 'recruit'))) & data$n_NAA_sigma == 0){
      stop("Cannot estimate ecov effect on recruitment when
      recruitment in each year is estimated freely as fixed effect parameters.
      Either remove ecov-recruit effect or estimate recruitment
      (or all numbers-at-age) as random effects.")
    }
    if(is.null(ecov$where)) stop("ecov$where must be specified as 'none', 'recruit', 'M', and/or 'q' for each ecov.")
    if(any(sapply(ecov$where, function(x) any(!(x %in% c('recruit','M','q','none')))))){
      stop("Only ecov effects on recruitment, M, and catchability (q) currently implemented.
      Set ecov$where = 'none' or one or more of 'recruit', 'M', 'q'.")
    }

    if(!all(ecov$process_model %in% c(NA,"rw", "ar1"))){
      stop("ecov$process_model must be 'rw' (random walk), 'ar1', or NA (do not fit)")
    }
    for(i in 1:data$n_Ecov) {
      if(is.na(ecov$process_model[i]) & any(ecov$where[[i]] %in% c("recruit", "M", "q"))){ #ecov$how !=0){
      stop(paste0("ecov$process_model ", i, " is turned off (NA) but ecov$where includes 'M', 'recruit', or 'q'.
       Either 1) choose an ecov process model ('rw' or 'ar1'),
              2) turn off ecov (set ecov$where[i] = 'none' and ecov$process_model = NA),
           or 3) fit ecov but with no effect on population (ecov$where[i] = 'none', ecov$process_model[i] = 'rw' or 'ar1')."))
      }
    }
    data$Ecov_model <- sapply(ecov$process_model, match, c("rw", "ar1"))
    for(i in 1:data$n_Ecov){
      if(any(ecov$where[[i]] == "recruit")) data$Ecov_where[i,1] = 1 
      if(any(ecov$where[[i]] == "M")) data$Ecov_where[i,2] = 1 
      if(any(ecov$where[[i]] == "q")) {
        if(is.null(ecov$indices)) stop("ecov$indices must be specified if any ecov$where == 'q'") 
        if(!is.list(ecov$indices)) stop("ecov$indices must be a specified as a list (length = n_Ecov) of vectors of any indices each Ecov affects") 
        data$Ecov_where[i,2 + ecov$indices[[i]]] = 1 
      }
    }

    if(is.null(ecov$how) & ('recruit' %in% unlist(ecov$where))) stop("ecov$how must be specified when any ecov is affecting recruitment")
    if(length(ecov$how) != data$n_Ecov) stop("ecov$how must be a vector of length(n.ecov)")
    data$Ecov_how = rep(0, data$n_Ecov) 
    for(i in 1:data$n_Ecov){
      if(data$Ecov_where[i,2] == 1) if(!ecov$how[i] %in% c(0,1)){
        stop("Sorry, only the following ecov effects on M are currently implemented.
        Set ecov$how = 0 (no effect) or 1 (effect on mean M, shared across ages).")
      }
      if(data$Ecov_where[i,2] == 1 & ecov$how[i] == 0) data$Ecov_where[i,2] = 0 #change back to zero because data$Ecov_where now defines which Ecov_beta to estimate
      if(any(data$Ecov_where[i,index_effects] == 1)) if(!ecov$how[i] %in% c(0,1)){
        stop("Sorry, only the following ecov effects on q are currently implemented.
        Set ecov$how = 0 (no effect) or 1 (effect on q).")
      }
      if(any(data$Ecov_where[i,index_effects] == 1) & ecov$how[i] == 0) data$Ecov_where[i,index_effects] = 0 #change back to zero because data$Ecov_where now defines which Ecov_beta to estimate
      if(data$Ecov_where[i,1] == 1) if(!ecov$how[i] %in% c(0,1,2,4)){
        stop("Sorry, only the following ecov effects on recruitment are currently implemented.
        Set ecov$how = 0 (no effect), 1 (controlling), 2 (limiting, Bev-Holt only), or 4 (masking).")
      }
      if(data$Ecov_where[i,1] == 1 & ecov$how[i] == 0) data$Ecov_where[i,1] = 0 #change back to zero because data$Ecov_where now defines which Ecov_beta to estimate
      if(data$Ecov_where[i,1] == 1 & data$recruit_model == 1){
        stop("Random walk recruitment cannot have an ecov effect on recruitment.
        Either choose a different recruit_model (2, 3, or 4), or remove the Ecov effect.")
      }
      if(data$Ecov_where[i,1] == 1) if(data$recruit_model == 4 & ecov$how[i] == 2){
        stop("'Limiting' ecov effect on Ricker recruitment not implemented.
        Either set ecov$how = 0 (no effect), 1 (controlling), or 4 (masking)...
        Or set recruit_model = 3 (Bev-Holt).")
      }      
      #currently only need this if recruitment effects modeled.
      if(data$Ecov_where[i,1] == 1) data$Ecov_how[i] = ecov$how[i]
    }

    Ecov_poly <- matrix(1,data$n_Ecov, n_effects)
    if(!is.null(ecov$link_model)){
      if(!is.na(ecov$link_model)){
        if(!is.list(ecov$link_model)) ecov$link_model = lapply(ecov$link_model, rep, n_effects)
        for(i in 1:data$n_Ecov) for(j in 1:n_effects){
          ecov_str <- strsplit(ecov$link_model[[i]][j],"-")[[1]]
          if(!ecov_str[1] %in% c('linear','poly')) stop("Only 'linear' or 'poly-x' (x = 1, 2, ...) ecov link models implemented.")
          if(ecov_str[1]=='poly') Ecov_poly[i,j] <- as.numeric(ecov_str[2])
        }
      }
    }

    if(!is.null(ecov$ages)){
      if(length(ecov$ages) != data$n_Ecov) stop("ecov$ages must be a list of length n.ecov")
      for(i in 1:data$n_Ecov){
        if(!all(ecov$ages[[i]] %in% 1:data$n_ages)) stop("All ecov$ages must be in 1:n.ages")
      }
    } else {
      ecov$ages <- vector("list", data$n_Ecov)
      for(i in 1:data$n_Ecov) ecov$ages[[i]] <- 1:data$n_ages # default: ecov affects all ages
    }

    cat(paste0("Please check that the environmental covariates have been loaded and interpreted correctly.

      Model years: ", data$year1_model, " to ", end_model,"
      Ecov years: ", data$year1_Ecov, " to ", end_Ecov,"

    "))

    for(i in 1:data$n_Ecov){
      years <- data$Ecov_year[as.logical(data$Ecov_use_obs[,i])]
      lastyr <- tail(years,1)

      if(data$Ecov_where[i,1] == 1){ # recruitment
        cat(paste0("Ecov ",i,": ",ecov$label[i],"
          ",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how[i]+1]," (",ecov$link[[i]][1],") effect on: recruitment

          Model years:
        "))

        cat(years, fill=TRUE)

        cat(paste0("Lag: ",ecov$lag[i,1],"
        Ex: ",ecov$label[i]," in ",years[1]," affects recruitment in ",years[1+ecov$lag[i,1]],"
            ",ecov$label[i]," in ",lastyr," affects recruitment in ",lastyr+ecov$lag[i,1],"

        \n"))
      }

      if(data$Ecov_where[i,2] == 1){ # M
        cat(paste0("Ecov ",i,": ",ecov$label[i]," effect (", ecov$link[[i]][2], ") on: M 

          Model years:
        "))
        cat(years, fill=TRUE)

        cat(paste0("Lag: ",ecov$lag[i,2],"
        Ex: ",ecov$label[i]," in ",years[1]," affects M in ",years[1+ecov$lag[i,2]],"
            ",ecov$label[i]," in ",lastyr," affects M in ",lastyr+ecov$lag[i,2],"

        \n"))
      }

      for(j in index_effects) if(data$Ecov_where[i,j] == 1){ # q
        cat(paste0("Ecov ",i,": ",ecov$label[i]," effect (", ecov$link[[i]][j], ") on: q for index ", j + 1 - min(index_effects), " 

          Model years:
        "))
        cat(years, fill=TRUE)

        cat(paste0("Lag: ",ecov$lag[i,j],"
        Ex: ",ecov$label[i]," in ",years[1]," affects index ", j + 1 - min(index_effects), " in ",years[1+ecov$lag[i,j]],"
            ",ecov$label[i]," in ",lastyr," affects M index ", j + 1 - min(index_effects), " in ",lastyr+ecov$lag[i,j],"

        \n"))
      }
    }
    data$Ecov_label <- list(data$Ecov_label)
  } # end load Ecov

  # Ecov pars
  par$Ecov_re = matrix(rnorm(data$n_years_Ecov*data$n_Ecov), data$n_years_Ecov, data$n_Ecov)
  max.poly <- max(Ecov_poly)
  par$Ecov_beta = array(0, dim=c(n_effects, max.poly, data$n_Ecov, data$n_ages)) # beta_R in eqns 4-5, Miller et al. (2016)
  par$Ecov_process_pars = matrix(0, 3, data$n_Ecov) # nrows = RW: 2 par (Ecov1, log_sig), AR1: 3 par (mu, log_sig, phi); ncol = N_ecov
  #this row is for the mean not the sd of the process
  par$Ecov_process_pars[1,] = -1.3 # start sig_ecov at 0.27
  #changing the initial value for sig_ecov to the right place actually causes tests to not pass!
  #par$Ecov_process_pars[2,] = -1.3 # start sig_ecov at 0.27
  if(!is.null(ecov$process_mean_vals)){
    for(i in 1:data$n_Ecov) if(data$Ecov_model[i]==2) par$Ecov_process_pars[1,i] = ecov$process_mean_vals[i]
  }
  if(!is.null(ecov$process_sig_vals)){
    for(i in 1:data$n_Ecov) par$Ecov_process_pars[2,i] = log(ecov$process_sig_vals[i])
  }
  if(!is.null(ecov$process_cor_vals)){
    inv_trans_rho <- function(rho) log(rho+1) - log(1-rho) # not same transformation as other re on cpp side.
    for(i in 1:data$n_Ecov)  if(data$Ecov_model[i]==2) par$Ecov_process_pars[3,i] = inv_trans_rho(ecov$process_cor_vals[i])
  }

  par$Ecov_obs_logsigma <- par.Ecov.obs.logsigma
  par$Ecov_obs_sigma_par <- par.Ecov.obs.sigma.par
  par$Ecov_obs_logsigma_re = par.Ecov.obs.logsigma.re

  # turn off 3rd Ecov par if it's a RW
  tmp.pars <- par$Ecov_process_pars
  for(i in 1:data$n_Ecov) tmp.pars[3,i] <- ifelse(data$Ecov_model[i]==1, NA, 0)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)

  # turn off Ecov_beta to fit Ecov process model without effect on population
  tmp <- array(NA, dim=dim(par$Ecov_beta))
  ct = 1
  for(n in 1:n_effects) for(j in 1:data$n_Ecov){
    for(i in 1:max.poly){
      for(a in 1:data$n_ages){
        if(data$Ecov_where[j,n] == 1 & i <= Ecov_poly[j,n] & a %in% ecov$ages[[j]]) 
        {
          tmp[n,i,j,a] = ct # default share ecov_beta across ages
          if(!is.null(ecov$beta_vals[[j]][[n]])) {
             par$Ecov_beta[n,i,j,a] = ecov$beta_vals[[j]][[n]][i,a]
          }
        }
      }
      ct = ct+1
    }
  }
  map$Ecov_beta = factor(tmp)

  # turn off Ecov pars if no Ecov (re, process)
  # for any Ecov_model = NA, ecov$how must be 0 and beta is already turned off
  data$Ecov_model[is.na(data$Ecov_model)] = 0 # turn any NA into 0
  tmp.re <- matrix(1:length(par$Ecov_re), data$n_years_Ecov, data$n_Ecov, byrow=FALSE)
  for(i in 1:data$n_Ecov){
    tmp.pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else tmp.pars[,i]
    tmp.re[,i] <- if(data$Ecov_model[i]==0) rep(NA,data$n_years_Ecov) else tmp.re[,i]
    if(data$Ecov_model[i]==1) tmp.re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
  }
  ind.notNA <- which(!is.na(tmp.re))
  tmp.re[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_re = factor(tmp.re)
  ind.notNA <- which(!is.na(tmp.pars))
  tmp.pars[ind.notNA] <- 1:length(ind.notNA)
  map$Ecov_process_pars = factor(tmp.pars)
  map$Ecov_obs_logsigma <- factor(map.Ecov.obs.logsigma)
  map$Ecov_obs_sigma_par <- factor(map.Ecov.obs.sigma.par)
  map$Ecov_obs_logsigma_re = factor(map.Ecov.obs.logsigma.re)

  input$data = data
  input$par = par
  input$map = map
  
	# add vector of all observations for one step ahead residuals ==========================
	input = wham:::set_osa_obs(input)
	#print("osa_obs")

	# projection data will always be modified by 'prepare_projection'
	input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
	#print("proj")

	#set any parameters as random effects
	input$random = NULL
	input = wham:::set_random(input)
	#print("random")
	
  return(input)
}
