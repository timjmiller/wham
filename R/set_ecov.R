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
#'     \item{$use_obs}{T/F (or 1/0) vector/matrix of the same dimension as \code{$mean} and \code{$logsigma}.
#'     Use the observation? Can be used to ignore subsets of the ecov without changing data files.}
#'     \item{$process_model}{Process model for the ecov time-series. \code{"rw"} = random walk, \code{"ar1"} = 1st order autoregressive, 
#'        \code{NA} = do not fit}
#'     \item{$process_mean_vals}{vector of (initial) mean values for the ecov time-series.}
#'     \item{$process_sig_vals}{vector of (initial) standard deviation values for the ecov time-series.}
#'     \item{$process_cor_vals}{vector of (initial) correlation values for the ecov time-series.}
#'     \item{$recruitment_how}{character matrix (n_Ecov x n_stocks) indicating how each ecov affects recruitment for each stock. 
#'        Options are based on (see \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}) 
#'        combined with the order of orthogonal polynomial of the covariate and has the form "type-lag-order". "type" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "controlling"}{pre-recruit density-independent mortality.}
#'           \item{= "limiting"}{ maximum recruitment, e.g. ecov determines amount of suitable habitat)}
#'           \item{= "lethal"}{threshold, i.e. R --> 0 at some ecov value.}
#'           \item{= "masking"}{metabolic/growth, decreases dR/dS}
#'           \item{= "directive"}{e.g. behavioral}
#'          }
#'         for type other than "none", "lag" can be:
#'         \describe{
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects recruitment in year y+1.}
#'         }
#'         for "type" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed recruitment parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'         }
#'         so "limiting-lag-1-poly-2" would model the covariate affecting recruitment the next year (lag = 1) as a second order orthogonal 
#'          polynomial (\eqn{b_0 + b_1*ecov + b_2*ecov^2 + ...}) limiting effect.
#'      }
#'     \item{$M_how}{character array (n_Ecov x n_stocks x n_ages x n_regions) indicating how each ecov affects M by age,stock,region and 
#'        has the form "lag-order". "lag" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects M in year y+1.}
#'         }
#'         for "lag" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed M parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'          }
#'       }
#'     \item{$M_effect_map}{integer array (n_stocks x n_ages x n_regions x n_Ecov) indicating which estimated effects are common by age,stock,region.
#'       If not specified there the same effect is estimated for all M where $M_how is other than "none" for each covariate.}
#'     \item{$q_how}{character matrix (n_Ecov x n_indices) indicating whether each ecov affects catchability for each index. and has 
#'      the form "lag-order". "lag" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects catchability in year y+1.}
#'         }
#'         for "lag" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed catchability parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'          }
#'       }
#'     \item{$move_how}{character array (n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions - 1) indicating whether each ecov 
#'        affects movement from one region to the others by stock,age,season. and has the form "lag-order". "lag" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects a movement parameter in year y+1.}
#'         }
#'         for "lag" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed movement parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'          }
#'        }
#'     \item{$move_effect_map}{integer array (n_stocks x n_ages x n_seasons x n_regions x n_regions-1 x n_Ecov) indicating which estimated 
#'        effects are common by age,stock,region, season etc. If not specified the same effect is estimated for all movement parameters
#'        where $move_how is other than "none" for each covariate.}
#'     \item{$beta_R_vals}{n_stocks x n_ecov x max(n_poly_R) array of initial values for effects on recruitment.}
#'     \item{$beta_M_vals}{n_stocks x n_ages x n_regions x n_ecov x max(n_poly_M) array of initial values for effects on natural mortality.}
#'     \item{$beta_q_vals}{n_indices x n_ecov x max(n_poly_q) array of initial values for effects on catchability.}
#'     \item{$beta_mu_vals}{n_stocks x n_ages x n_seasons x n_regions x n_regions - 1 x n_ecov x max(n_poly_move) array of initial values for effects on movement parameters.}
#'   }
#'
#' @return a named list with same elements as the input provided with environmental covariate observations, effects, and model options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)
#' input <- prepare_wham_input(asap3, NAA_re = list(sigma = "rec"))
#' ecov <- list(
#'  label = "GSI",
#'  mean = as.matrix(env.dat$GSI),
#'  logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
#'  year = env.dat$year,
#'  use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
#'  process_model = 'ar1') # "rw" or "ar1"
#' input <- set_ecov(input, ecov = ecov) #GSI in the model without any effects
#' }
#'
#' @export
set_ecov = function(input, ecov) {
  data = input$data
  par = input$par
  map = input$map
  input$log$ecov <- list()
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  ecov_pars = c("Ecov_re", "Ecov_beta_R", "Ecov_beta_M", "Ecov_beta_mu", "Ecov_beta_q", 
    "Ecov_process_pars", "Ecov_obs_log_sigma", "Ecov_obs_logsigma_re", "Ecov_obs_sigma_par")
  map <- map[(!names(map) %in% ecov_pars)]
  
  #define new dimensions for all effects a given ecov can have on assessment model
  #currently the order is recruitment, M, index 1, ..., n_indices
  #n_effects = 2 + data$n_indices
  #index_effects = 2+1:data$n_indices

  # --------------------------------------------------------------------------------
  # Environmental covariate data
  #set up for ecov == NULL
  data$Ecov_obs <- matrix(1, nrow=1, ncol=1)
  data$Ecov_obs_sigma_opt <- 1
  data$n_Ecov <- 1
  data$Ecov_model <- rep(0, data$n_Ecov)
  input$years_Ecov <- input$years[1]
  data$n_years_Ecov <- 1
  data$years_use_Ecov <- 0
  data$Ecov_use_obs <- matrix(0, nrow=1, ncol=1)
  
  data$Ecov_how_R <- matrix(0, data$n_Ecov, data$n_stocks)
  data$Ecov_how_q <- matrix(0, data$n_Ecov, data$n_indices)
  data$Ecov_how_M <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages,data$n_regions))
  data$Ecov_how_mu <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
  data$ind_Ecov_out_start_R <- data$ind_Ecov_out_end_R <- matrix(0, data$n_Ecov, data$n_stocks)
  data$ind_Ecov_out_start_q <- data$ind_Ecov_out_end_q <- matrix(0, data$n_Ecov, data$n_indices)
  data$ind_Ecov_out_start_M <- data$ind_Ecov_out_end_M <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
  data$ind_Ecov_out_start_mu <- data$ind_Ecov_out_end_mu <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, 
    data$n_regions, data$n_regions-1))
  
  input$Ecov_names <- "none"
  data$Ecov_use_re <- rep(0, data$n_Ecov)

  data$n_poly_Ecov_R <- matrix(1,data$n_Ecov, data$n_stocks)
  data$n_poly_Ecov_M <- array(1,dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
  data$n_poly_Ecov_mu <- array(1,dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
  #print(data$n_Ecov)
  #print(dim(data$n_poly_Ecov_mu))
  #print(max(data$n_poly_Ecov_mu))
  #print(data$n_poly_Ecov_mu)
  data$n_poly_Ecov_q <- matrix(1,data$n_Ecov, data$n_indices)
  
  par$Ecov_obs_logsigma <- matrix(-2.3, nrow=1, ncol=1)
  par$Ecov_obs_logsigma_re <- matrix(-2.3, nrow=1, ncol=1)
  par$Ecov_obs_sigma_par <- matrix(0, nrow=1, ncol=1)
  par$Ecov_process_pars <- matrix(0, 3, data$n_Ecov)  
  par$Ecov_re <- matrix(0, data$n_years_Ecov, data$n_Ecov)
  par$Ecov_beta_R <- array(0, dim = c(data$n_stocks, data$n_Ecov, max(data$n_poly_Ecov_R)))
  par$Ecov_beta_q <- array(0, dim = c(data$n_indices, data$n_Ecov, max(data$n_poly_Ecov_q)))
  par$Ecov_beta_M <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov, max(data$n_poly_Ecov_M)))
  par$Ecov_beta_mu <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1, data$n_Ecov, max(data$n_poly_Ecov_mu,0)))
  map$Ecov_obs_logsigma <- rep(NA, length(par$Ecov_obs_logsigma))
  map$Ecov_obs_logsigma_re <- rep(NA, length(par$Ecov_obs_logsigma_re))
  map$Ecov_obs_sigma_par <- rep(NA, length(par$Ecov_obs_sigma_par))
  map$Ecov_process_pars <- rep(NA, length(par$Ecov_process_pars))
  map$Ecov_re <- rep(NA, length(par$Ecov_re))
  map$Ecov_beta_R <- rep(NA, length(par$Ecov_beta_R))
  map$Ecov_beta_q <- rep(NA, length(par$Ecov_beta_q))
  map$Ecov_beta_M <- rep(NA, length(par$Ecov_beta_M))
  map$Ecov_beta_mu <- rep(NA, length(par$Ecov_beta_mu))

  if(!is.null(ecov)){

    #FIRST: configure observation models, likelihoods, and process models, likelihoods for ecovs
    if(class(ecov$mean)[1] == "matrix") {data$Ecov_obs <- ecov$mean} else{
      input$log$ecov <- c(input$log$ecov, "NOTE: ecov$mean is not a matrix. Coercing to a matrix...")
      data$Ecov_obs <- as.matrix(ecov$mean)
    }
    if(class(ecov$use_obs)[1] == "matrix"){
      data$Ecov_use_obs <- ecov$use_obs
    } else{
      input$log$ecov <- c(input$log$ecov, "NOTE: ecov$use_obs is not a matrix with same dimensions as ecov$mean. Coercing to a matrix...")
      data$Ecov_use_obs <- as.matrix(as.integer(ecov$use_obs))
    }
    if(!identical(dim(data$Ecov_use_obs), dim(data$Ecov_obs))) stop("Dimensions of ecov$use_obs != dimensions of ecov$mean")

    # Handle Ecov sigma options
    data$n_Ecov <- dim(data$Ecov_obs)[2] # num Ecovs
    n_Ecov_obs <- dim(data$Ecov_obs)[1] # num Ecov obs
    data$Ecov_obs_sigma_opt = rep(1, data$n_Ecov) # Ecov sigma given, initialized at given values, not estimated by default

    if(length(ecov$year) != n_Ecov_obs) stop("ecov$year is not the same length as # rows in ecov$mean")
    #data$Ecov_year <- as.numeric(ecov$year)
    input$years_Ecov <- as.numeric(ecov$year)
    end_model <- tail(input$years,1)
    end_Ecov <- tail(ecov$year,1)
    
    #define labels and report info about effects of each
    if(length(ecov$label) == data$n_Ecov){
      input$Ecov_names <- ecov$label
    } else {
      input$log$ecov <- c(input$log$ecov, "NOTE: Number of Ecov labels not equal to number of Ecovs
              Setting Ecov labels = 'Ecov 1', 'Ecov 2', ...")
      input$Ecov_names = paste0("Ecov ",1:data$n_Ecov)
    }    
    
    
    #default NAs for parameter matrix of observation standard deviations
    par$Ecov_obs_logsigma = matrix(NA, n_Ecov_obs, data$n_Ecov)
    map$Ecov_obs_logsigma <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation

    logsig_more = list()
    if(is.list(ecov$logsigma)){
      n = length(ecov$logsigma)
      if(n>1) logsig_more = ecov$logsigma[[2]] #further elements than the second are ignored.
      ecov$logsigma = ecov$logsigma[[1]]
    }

    if(class(ecov$logsigma)[1] == "matrix"){
      par$Ecov_obs_logsigma <- ecov$logsigma
      if(!identical(dim(par$Ecov_obs_logsigma), dim(data$Ecov_obs))) stop("Dimensions of ecov$mean != dimensions of ecov$logsigma")
    }
    if(class(ecov$logsigma)[1] == 'numeric'){
      #data$Ecov_obs_sigma_opt[] = 1 #defined above
      input$log$ecov <- c(input$log$ecov, "ecov$logsigma is numeric. Coercing to a matrix... \n")
      if(length(ecov$logsigma) == data$n_Ecov) par$Ecov_obs_logsigma <- matrix(rep(ecov$logsigma, each=n_Ecov_obs), ncol=data$n_Ecov)
      if(length(ecov$logsigma) == n_Ecov_obs & data$n_Ecov == 1) par$Ecov_obs_logsigma <- matrix(ecov$logsigma, ncol=1)
      if(length(ecov$logsigma) != data$n_Ecov & length(ecov$logsigma) != n_Ecov_obs) stop("ecov$logsigma is numeric but length is not equal to # of ecovs or ecov observations")
    }

    #set up and check length of ecov$process_model
    if(length(ecov$process_model) == 1) ecov$process_model = rep(ecov$process_model, data$n_Ecov) #use the single value for all Ecovs
    if(length(ecov$process_model) != data$n_Ecov) stop("length of ecov$process_model must be either 1 or the number of Ecovs")

    #now over write ecov$logsigma with logsig_more if available because the fixed obs var matrices have been defined
    if(length(logsig_more)) ecov$logsigma = logsig_more


    if(class(ecov$logsigma)[1] == 'character'){
      #check that estimation options are right
      if(!all(ecov$logsigma %in% c("est_1", "est_re", NA))){
        stop("ecov$logsigma or ecov$logsigma[[2]] is character and must be NA (do not estimate), 'est_1' (single variance parameter), or 'est_re' (iid re annual variance parameters)")
      }
      if(length(ecov$logsigma) == 1) ecov$logsigma = rep(ecov$logsigma, data$n_Ecov) #use the single value for all Ecovs
      #check length of estimation options
      if(length(ecov$logsigma) != data$n_Ecov) stop("length of ecov$logsigma when character must be either 1 or the number of Ecovs")


      for(i in 1:data$n_Ecov) {
        if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_1'){ # estimate 1 Ecov obs sigma for each Ecov
          data$Ecov_obs_sigma_opt[i] = 2
          par$Ecov_obs_logsigma[,i] <- -1.3 
          map$Ecov_obs_logsigma[,i] <- i 
        }
        if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_re'){
          data$Ecov_obs_sigma_opt[i] = 4
          map$Ecov_obs_logsigma[,i] <- NA # turn off estimation of fixed effects
          par$Ecov_obs_sigma_par[,i] <- c(-1.3, -2.3) # random effect pars
          map$Ecov_obs_sigma_par[,i] <- max(0, map$Ecov_obs_sigma_par, na.rm =T) + 1:2 
        }
      }

      #default values for: 
      #map of observation error variance parameters
      map$Ecov_obs_logsigma_re <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
      #initial values of random effects
      par$Ecov_obs_logsigma_re = matrix(0, n_Ecov_obs, data$n_Ecov)

      for(i in 1:data$n_Ecov) if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_re') {
        map$Ecov_obs_logsigma_re[,i] = max(0, map$Ecov_obs_logsigma_re, na.rm=T) + 1:data$n_years_Ecov
        par$Ecov_obs_logsigma_re[,i] <- par$Ecov_obs_logsigma[,i] # random effect initialize at values in matrix provided
      }

    }

    
    if(!all(ecov$process_model %in% c(NA,"rw", "ar1"))){
      stop("ecov$process_model must be 'rw' (random walk), 'ar1', or NA (do not fit)")
    }
    data$Ecov_model <- sapply(ecov$process_model, match, c("rw", "ar1"))
    data$Ecov_model[is.na(data$Ecov_model)] <- 0 #don't fit for NA


    if(is.null(ecov$recruitment_how)) ecov$recruitment_how <- matrix("none", data$n_Ecov, data$n_stocks)
    if(is.null(ecov$q_how)) ecov$q_how <- matrix("none", data$n_Ecov, data$n_indices)
    if(is.null(ecov$M_how)) ecov$M_how <- array("none", dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
    if(is.null(ecov$M_effect_map)) ecov$M_effect_map <- array(NA, dim = c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov))
    if(is.null(ecov$move_how)) ecov$move_how <- array("none", dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    if(is.null(ecov$move_effect_map)) ecov$move_effect_map <- array(NA, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1, data$n_Ecov))

    for(i in 1:data$n_Ecov) {
      ecov_used <- any(c(ecov$recruitment_how[i,],ecov$q_how[i,],ecov$M_how[i,,,],ecov$move_how[i,,,,,]) != "none")
      if(data$Ecov_model[i] == 0 & ecov_used) {
      stop(paste0("ecov$process_model ", i, " is turned off (NA) but an effect is specified in ecov$recruitment_how, ecov$q_how, ecov$M_how, and/or ecov$move_how.
       Either 1) choose an ecov process model ('rw' or 'ar1'),
              2) turn off ecov (set the ecov$type_how = 'none' and ecov$process_model = NA),
           or 3) fit ecov but with no effect on population (set the ecov$type_how = 'none' and ecov$process_model[i] = 'rw' or 'ar1')."))
      }
    }

    #make R
    ecov$lag_R <- matrix(1, data$n_Ecov, data$n_stocks)
    data$n_poly_Ecov_R <- matrix(1, data$n_Ecov, data$n_stocks)
    data$Ecov_how_R <- matrix(0, data$n_Ecov, data$n_stocks)
    if(!is.null(ecov$recruitment_how)) {
      for(s in 1:data$n_stocks) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$recruitment_how[i,s], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,4:5))) stop(paste0("Form of ecov$recruitment_how[",i,",", s, "] for ecov, ", i, " and stock , ", s, 
          " should be 'none', 'type-lag-l-poly-p' or 'type-lag-l-linear'"))
        if(tmp[1] != "none") ecov$lag_R[i,s] <- as.integer(tmp[3])
        if(length(tmp) == 5) data$n_poly_Ecov_R[i,s] <- as.integer(tmp[5])
        if(tmp[1] != "none" & data$recruit_model[s] == 1) stop(paste0("Random walk recruitment for stock ", s, " cannot have an ecov effect on recruitment.
            Either choose a different recruit_model (2, 3, or 4), or remove the Ecov effect."))
        if(data$recruit_model[s] == 2 & !(tmp[1] %in% c("none", "controlling"))) stop(paste0("Random recruitment about mean for stock ", s, 
          "only allows effect type to be 'none' or 'controlling'."))
        if(data$recruit_model[s] == 4 & tmp[1] == "limiting") stop(paste0("'Limiting' ecov effect on Ricker recruitment for stock ", s, " not implemented.
          Either set ecov$how_R = 0 (no effect), 1 (controlling), or 4 (masking)...
          Or set recruit_model = 3 (Bev-Holt)."))
        if(tmp[1] != "none" & data$NAA_re_model[s] == 0) stop(paste0("Cannot estimate ecov effect on recruitment for stock ", s, " when
          recruitment in each year is estimated freely as fixed effect parameters.
          Either remove ecov-recruit effect or estimate recruitment
          (or all numbers-at-age) as random effects."))
        if(!(tmp[1] %in% c("none", "controlling", "limiting", "lethal", "masking", "directive"))) stop(paste0("The first component of each
          element of ecov$recruitment_how must be one of the following: 'none', 'controlling', 'limiting', 'lethal', 'masking', 'directive'"))
        if(tmp[1] == "controlling") data$Ecov_how_R[i,s] <- 1
        if(tmp[1] == "limiting") data$Ecov_how_R[i,s] <- 2
        if(tmp[1] == "lethal") data$Ecov_how_R[i,s] <- 3
        if(tmp[1] == "masking") data$Ecov_how_R[i,s] <- 4
        if(tmp[1] == "directive") data$Ecov_how_R[i,s] <- 5
        #it is already set to 1
        #if(length(tmp)== 4 & tmp[4] == "linear") ecov_poly_R[i,s] <- 1 
      }
    } else {
      input$log$ecov <- c(input$log$ecov, "no ecov$recruitment_how is provided so no effects on recruitment will be estimated. \n")
    }

    par$Ecov_beta_R <- array(0, dim = c(data$n_stocks, data$n_Ecov, max(data$n_poly_Ecov_R)))
    map$Ecov_beta_R <- array(NA, dim = dim(par$Ecov_beta_R))
    k <- 0
    for(s in 1:data$n_stocks) for(i in 1:data$n_Ecov) if(data$Ecov_how_R[i,s]>0){
      n <- data$n_poly_Ecov_R[i,s]
      map$Ecov_beta_R[s,i,1:n] <- k + 1:n
      k <- k + n
    }
    
    #make M
    ecov$lag_M <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
    data$n_poly_Ecov_M <- array(1, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_regions))
    data$Ecov_how_M <- array(0, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_regions))
    if(!is.null(ecov$M_how)) {
      for(s in 1:data$n_stocks) for(r in 1:data$n_regions) for(a in 1:data$n_ages) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$M_how[i,s,a,r], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,3:4))) stop(paste0("Form of ecov$M_how[",i,",", s, ",", a, ",", r, "] for ecov, ", i, ", stock , ", s, 
          ", age , ", a, ", region , ", r," should be 'none', 'lag-l-poly-p' or 'lag-l-linear'"))
        if(tmp[1] != "none"){
          data$Ecov_how_M[i,s,a,r] <- 1
          ecov$lag_M[i,s,a,r] <- as.integer(tmp[2])
          if(length(tmp) == 4) {
            #print(tmp)
            data$n_poly_Ecov_M[i,s,a,r] <- as.integer(tmp[4])
          }
        }
        #it is already set to 1
        #if(length(tmp)== 3 & tmp[4] == "linear") ecov_poly_R[i,s] <- 1 
      }
    } else {
      input$log$ecov <- c(input$log$ecov, "no ecov$M_how is provided so no effects on natural mortality will be estimated. \n")
    }
    #print(c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov, max(data$n_poly_Ecov_M)))
     # stop()
    par$Ecov_beta_M <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov, max(data$n_poly_Ecov_M)))
    map$Ecov_beta_M <- array(NA, dim = dim(par$Ecov_beta_M))
    
    if(all(is.na(ecov$M_effect_map))){ #wasn't provided or some non-NA values were provided
      k <- 0 #same beta for all effects of this covariate
      for(i in 1:data$n_Ecov) {
        for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) {
          if(data$Ecov_how_M[i,s,a,r] != 0) {
            map$Ecov_beta_M[s,a,r,i,] <- k + 1:data$n_poly_Ecov_M[i,s,a,r]
            #ecov$M_effect_map[s,a,r,i] <- i #same beta for all M for each Ecov
          }
        }
        k <- max(map$Ecov_beta_M, 0, na.rm = TRUE)
      }
    } else { #user provided 
      ecov$M_effect_map[] <- as.integer(factor(ecov$M_effect_map)) #reset to 1:n_unique_effects_M
      n_eff_tot <- max(ecov$M_effect_map, na.rm = TRUE)
      for(i in 1:data$n_Ecov) {
        for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) {
          if(data$Ecov_how_M[i,s,a,r] != 0) {
            n_eff <- data$n_poly_Ecov_M[i,s,a,r]
            k <- ecov$M_effect_map[s,a,r,i]
            if(n_eff>1) {
              k <- c(k, n_eff_tot+1:(n_eff-1))
              n_eff_tot <- max(k)
            }
            map$Ecov_beta_M[s,a,r,i,1:n_eff] <- k
          }
        }
      }
    }
    # unique_effects_M <- unique(ecov$M_effect_map[which(!is.na(ecov$M_effect_map))]) #will be at least 1 value
    # n_unique_effects_M <- length(unique_effects_M)
    # for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) for(i in 1:data$n_Ecov) if(data$Ecov_how_M[i,s,a,r]>0){
    #   n <- data$n_poly_Ecov_M[i,s,a,r]
    #   if(n>1) { #need unique parameters for higher order polynomials
    #     for(eff in 2:n) map$Ecov_beta_M[s,a,r,i,eff] <- map$Ecov_beta_M[s,a,r,i,1] + n_unique_effects_M * (eff-1)
    #   }
    # }

    #make q
    ecov$lag_q <- matrix(0, data$n_Ecov, data$n_indices)
    data$n_poly_Ecov_q <- matrix(1, data$n_Ecov, data$n_indices)
    data$Ecov_how_q <- matrix(0, data$n_Ecov, data$n_indices)
    if(!is.null(ecov$q_how)) {
      for(ind in 1:data$n_indices) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$q_how[i,ind], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,3:4))) stop(paste0("Form of ecov$q_how[",i,",", ind, "] for ecov, ", i, ", index , ", ind,
        " should be 'none', 'lag-l-poly-p' or 'lag-l-linear'"))
        if(tmp[1] != "none"){
          data$Ecov_how_q[i,ind] <- 1
          ecov$lag_q[i,ind] <- as.integer(tmp[2])
          if(length(tmp) == 4) data$n_poly_Ecov_q[i,ind] <- as.integer(tmp[4])
        }
      }
    } else {
      input$log$ecov <- c(input$log$ecov, "no ecov$q_how is provided so no effects on catchability will be estimated. \n")
    }

    par$Ecov_beta_q <- array(0, dim = c(data$n_indices, data$n_Ecov, max(data$n_poly_Ecov_q)))
    map$Ecov_beta_q <- array(NA, dim = dim(par$Ecov_beta_q))
    k <- 1
    for(ind in 1:data$n_indices) for(i in 1:data$n_Ecov) if(data$Ecov_how_q[i,ind]>0){
      n <- data$n_poly_Ecov_q[i,ind]
      map$Ecov_beta_q[ind,i,1:n] <- k + 1:n
      k <- k + n
    }
    
    #make mu
    ecov$lag_mu <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))

    data$n_poly_Ecov_mu <- array(1, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    data$Ecov_how_mu <- array(0, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    if(data$n_regions > 1 & !is.null(ecov$move_how)) {
      for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$move_how[i,s,a,t,r,rr], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,3:4))) stop(paste0("Form of ecov$move_how[" ,i, ",", s, ",", a, ",", t, ",", r, ",", rr, "] for ecov, ", i, 
          ", stock , ", s, ", age , ", a, ", season , ", t, ", from region , ", r, " to region , ", rr, " should be 'none', 'lag-l-poly-p' or 'lag-l-linear'"))
        if(tmp[1] != "none"){
          data$Ecov_how_mu[i,s,a,t,r,rr] <- 1
          ecov$lag_mu[i,s,a,t,r,rr] <- as.integer(tmp[2])
          if(length(tmp) == 4) data$n_poly_Ecov_mu[i,s,a,t,r,rr] <- as.integer(tmp[4])
        }
        #it is already set to 1
        #if(length(tmp)== 3 & tmp[4] == "linear") ecov_poly_R[i,s] <- 1 
      }
    } else {
      if(data$n_regions> 1) input$log$ecov <- c(input$log$ecov, "There is more than 1 region, but no ecov$move_how is provided so no effects on movement parameters will be estimated.\n")
    }
    par$Ecov_beta_mu <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1, data$n_Ecov, 
      max(data$n_poly_Ecov_mu,0)))
    map$Ecov_beta_mu <- array(NA, dim = dim(par$Ecov_beta_mu))

    if(data$n_regions>1) {
      if(all(is.na(ecov$move_effect_map))){ #wasn't provided or some non-NA values were provided
        k <- 0 #same beta for all effects of this covariate
        for(i in 1:data$n_Ecov) {
          for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)){
            if(data$Ecov_how_mu[i,s,a,t,r,rr] != 0) {
              #ecov$move_effect_map[s,a,t,r,rr,i] <- NA #no effect so map it so.
              map$Ecov_beta_mu[s,a,t,r,rr,i,] <- k + 1:data$n_poly_Ecov_mu[i,s,a,t,r,rr]
            }
          }
          k <- max(map$Ecov_beta_mu, 0, na.rm = TRUE)
        }
      } else { #user provided 
        ecov$move_effect_map[] <- as.integer(factor(ecov$move_effect_map))
        n_eff_tot <- max(ecov$move_effect_map, na.rm = TRUE)
        for(i in 1:data$n_Ecov) {
          for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)){
            if(data$Ecov_how_mu[i,s,a,t,r,rr] != 0) {
              n_eff <- data$n_poly_Ecov_mu[i,s,a,t,r,rr]
              k <- ecov$move_effect_map[s,a,t,r,rr,i]
              if(n_eff>1) {
                k <- c(k, n_eff_tot+1:(n_eff-1))
                n_eff_tot <- max(k)
              }
              map$Ecov_beta_mu[s,a,t,r,rr,i,] <- k
            }
          }
        }
      }
    }

    #   map$Ecov_beta_mu[,,,,,,1] <- ecov$move_effect_map# k + 1:n
    #   unique_effects_mu <- unique(ecov$move_effect_map[which(!is.na(ecov$move_effect_map))]) #will be at least 1 value
    #   n_unique_effects_mu <- length(unique_effects_mu)
    #   for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) for(i in 1:data$n_Ecov) if(data$Ecov_how_M[i,s,a,r]>0){
    #     n <- data$n_poly_Ecov_mu[i,s,a,t,r,rr]
    #     if(n>1) { #need unique parameters for higher order polynomials
    #       for(eff in 2:n) map$Ecov_beta_mu[s,a,t,r,rr,i,eff] <- map$Ecov_beta_mu[s,a,t,r,rr,i,1] + n_unique_effects_mu * (eff-1)
    #     }
    #   }
    # }
    
    # # check that Ecov year vector doesn't have missing gaps
    # pad Ecov if it starts after model year1 - max(lag)
    max.lag = max(c(ecov$lag_R,ecov$lag_M,ecov$lag_mu,ecov$lag_q))
    #if(is.null(ecov$lag)) stop("ecov$lag needs to be provided for each ecov")
    #if(!is.list(ecov$lag)) ecov$lag = lapply(ecov$lag, function(x) rep(x,n_effects))
    if(input$years_Ecov[1] > input$years[1] - max.lag){
      input$log$ecov <- c(input$log$ecov, "one or more ecov does not start by model year 1 - max(lag). Padding ecov... \n")
      data$Ecov_obs <- rbind(matrix(0, nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov), data$Ecov_obs)
      par$Ecov_obs_logsigma <- rbind(matrix(par$Ecov_obs_logsigma[1,], nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov, byrow=T), par$Ecov_obs_logsigma)
      map$Ecov_obs_logsigma <- rbind(matrix(NA, nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov), map$Ecov_obs_logsigma)
      data$Ecov_use_obs <- rbind(matrix(0, nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov), data$Ecov_use_obs)
      input$years_Ecov <- c(seq(input$years[1] - max.lag, input$years_Ecov[1]-1), input$years_Ecov)
    }

    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      input$log$ecov <- c(input$log$ecov, "Ecov last year is before model last year. Padding Ecov... \n")
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      par$Ecov_obs_logsigma <- rbind(par$Ecov_obs_logsigma, matrix(par$Ecov_obs_logsigma[NROW(par$Ecov_obs_logsigma),], nrow = end_model-end_Ecov, ncol = data$n_Ecov, byrow=T))
      map$Ecov_obs_logsigma <- rbind(map$Ecov_obs_logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      input$years_Ecov <- c(input$years_Ecov, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- dim(data$Ecov_obs)[1] # num years Ecov to model (padded)
    data$years_use_Ecov <- 1:data$n_years_Ecov - 1
    
    data$Ecov_use_re <- rep(1, data$n_Ecov)
    for(i in 1:data$n_Ecov){
      if(all(data$Ecov_use_obs[,i]==0)) {
        data$Ecov_use_re[i] <- data$Ecov_model[i] <- 0
        input$log$ecov <- c(input$log$ecov, paste0("No observations for Ecov", i, " so no latent process will be estimated.\n"))
      }
    }


    #set up Ecov_re with padded dimensions
    #par$Ecov_re = matrix(rnorm(data$n_years_Ecov*data$n_Ecov), data$n_years_Ecov, data$n_Ecov)
    par$Ecov_re = matrix(0, data$n_years_Ecov, data$n_Ecov)
    map$Ecov_re <- matrix(1:length(par$Ecov_re), data$n_years_Ecov, data$n_Ecov, byrow=FALSE)
    for(i in 1:data$n_Ecov){
      #tmp.pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else tmp.pars[,i]
      map$Ecov_re[,i] <- if(data$Ecov_model[i]==0) rep(NA,data$n_years_Ecov) else map$Ecov_re[,i]
      if(data$Ecov_model[i]==1) map$Ecov_re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
    }
    ind.notNA <- which(!is.na(map$Ecov_re))
    map$Ecov_re[ind.notNA] <- 1:length(ind.notNA)


    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    # print(end_model)
    # print(ecov$lag_R)
    #stop()
    data$ind_Ecov_out_start_R <- which(input$years_Ecov == input$years[1]) - ecov$lag_R - 1
    data$ind_Ecov_out_end_R <- which(input$years_Ecov==end_model)- ecov$lag_R - 1 # -1 is for cpp indexing
    data$ind_Ecov_out_start_M <- which(input$years_Ecov == input$years[1]) - ecov$lag_M - 1
    data$ind_Ecov_out_end_M <- which(input$years_Ecov==end_model)- ecov$lag_M - 1 # -1 is for cpp indexing
    data$ind_Ecov_out_start_q <- which(input$years_Ecov == input$years[1]) - ecov$lag_q - 1
    data$ind_Ecov_out_end_q <- which(input$years_Ecov==end_model)- ecov$lag_q - 1 # -1 is for cpp indexing
    data$ind_Ecov_out_start_mu <- which(input$years_Ecov == input$years[1]) - ecov$lag_mu - 1
    data$ind_Ecov_out_end_mu <- which(input$years_Ecov==end_model)- ecov$lag_mu - 1 # -1 is for cpp indexing
    
    input$log$ecov <- c(input$log$ecov, paste0("Please check that the environmental covariates have been loaded and interpreted correctly.

      Model years: ", input$years[1], " to ", end_model,"
      Ecov years: ", input$years_Ecov[1], " to ", end_Ecov,"

    "))

    #if(!identical(length(ecov$lag), length(ecov$label), data$n_Ecov)) stop("Length of ecov$lag and ecov$label not equal to # Ecov")
    for(i in 1:data$n_Ecov){
      years <- input$years_Ecov[as.logical(data$Ecov_use_obs[,i])]
      lastyr <- tail(years,1)

      # recruitment
      for(s in 1:data$n_stocks) if(data$Ecov_how_R[i,s] > 0){ 
        input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i],"
          ",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how_R[i,s]+1]," (",
            ifelse(data$n_poly_Ecov_R[i,s] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_R[i,s])),
            ") effect on: recruitment for stock ", s, "

          Model years:
        "))

        input$log$ecov <- c(input$log$ecov, paste0(years, collapse=', '))

        input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_R[i,s],"
        Ex: ",ecov$label[i]," in ",years[1]," affects recruitment in ",years[1+ecov$lag_R[i,s]],"
            ",ecov$label[i]," in ",lastyr," affects recruitment in ",lastyr+ecov$lag_R[i,s],"

        \n"))
      }

      #M
      for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) if(data$Ecov_how_M[i,s,a,r] == 1){
        input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i]," effect (", 
          ifelse(data$n_poly_Ecov_M[i,s,a,r] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_M[i,s,a,r])), 
          ") on: M for stock ", s, " at age ", a, " in region ", r,
          "

          Model years:
        "))
        input$log$ecov <- c(input$log$ecov, paste0(years, collapse= ", "))

        input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_M[i,s,a,r],"
        Ex: ",ecov$label[i]," in ",years[1]," affects M in ",years[1+ecov$lag_M[i,s,a,r]],"
            ",ecov$label[i]," in ",lastyr," affects M in ",lastyr+ecov$lag_M[i,s,a,r],"

        \n"))
      }

       # q
      for(j in 1:data$n_indices) if(data$Ecov_how_q[i,j] == 1){
        input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i]," effect (", 
          ifelse(data$n_poly_Ecov_q[i,j] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_q[i,j])), 
          ") on: q for index ", j, 
          " 

          Model years:
        "))
        input$log$ecov <- c(input$log$ecov, paste0(years, collapse = ", "))

        input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_q[i,j],"
        Ex: ",ecov$label[i]," in ",years[1]," affects index ", j, " in ",years[1+ecov$lag_q[i,j]],"
            ",ecov$label[i]," in ",lastyr," affects index ", j, " in ",lastyr+ecov$lag_q[i,j],"

        \n"))
      }

      #movement
      if(data$n_regions>1) for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) if(data$n_regions>1) {
        for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) if(data$Ecov_how_mu[i,s,a,t,r,rr] == 1){
          input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i]," effect (", 
            ifelse(data$n_poly_Ecov_mu[i,s,a,t,r,rr] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_mu[i,s,a,t,r,rr])), 
            ") on: movement for stock ", s, " at age ", a, "in season ", t, "from region ", r, "to ", rr, "of the other regions
            
            Model years:
          "))
          input$log$ecov <- c(input$log$ecov, paste0(years, collapse=", "))

          input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_mu[i,s,a,t,r,rr],"
          Ex: ",ecov$label[i]," in ",years[1]," affects movement in ",years[1+ecov$lag_mu[i,s,a,t,r,rr]],"
              ",ecov$label[i]," in ",lastyr," affects movement in ",lastyr+lag_mu[i,s,a,t,r,rr],"

          \n"))
        }
      }
    }
    #input$Ecov_names <- list(input$Ecov_names)

    # Ecov process pars
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
      inv_trans_rho <- function(rho) log(rho+1) - log(1-rho) 
      for(i in 1:data$n_Ecov)  if(data$Ecov_model[i]==2) par$Ecov_process_pars[3,i] = inv_trans_rho(ecov$process_cor_vals[i])
    }

    # turn off Ecov pars if no Ecov (re, process)
    # turn off 3rd Ecov par if it's a RW
    map$Ecov_process_pars <- par$Ecov_process_pars
    for(i in 1:data$n_Ecov) {
      map$Ecov_process_pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else map$Ecov_process_pars[,i]
      map$Ecov_process_pars[3,i] <- ifelse(data$Ecov_model[i]==1, NA, 0)
    }
    ind.notNA <- which(!is.na(map$Ecov_process_pars))
    map$Ecov_process_pars[ind.notNA] <- 1:length(ind.notNA)

    #fill in initial values for any covariate effects
    for(i in c("beta_q", "beta_R", "beta_M", "beta_mu")){
      inp <- paste0(i, "_vals")
      if(!is.null(ecov[[inp]])){
        if(!length(dim(ecov[[inp]])) == length(dim(par[[paste0("Ecov_",i)]])) | !all(dim(ecov[[inp]]) == dim(par[[paste0("Ecov_",i)]]))){
          stop(paste0("ecov$",i, "_vals is not an array of the correct dimensions."))
        }
        par[[paste0("Ecov_",i)]] <- ecov[[inp]]
      }
    }

  } # end load Ecov


  map$Ecov_obs_logsigma <- factor(map$Ecov_obs_logsigma)
  map$Ecov_obs_sigma_par <- factor(map$Ecov_obs_sigma_par)
  map$Ecov_obs_logsigma_re <- factor(map$Ecov_obs_logsigma_re)
  map$Ecov_process_pars = factor(map$Ecov_process_pars)
  map$Ecov_re <- factor(map$Ecov_re)
  map$Ecov_beta_R <- factor(map$Ecov_beta_R)
  map$Ecov_beta_M <- factor(map$Ecov_beta_M)
  map$Ecov_beta_q <- factor(map$Ecov_beta_q)
  map$Ecov_beta_mu <- factor(map$Ecov_beta_mu)

  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$ecov)) input$log$ecov <- c("Ecov: \n", input$log$ecov)
	# add vector of all observations for one step ahead residuals ==========================
  if(!is_internal_call()) { #check whether called by prepare_wham_input
    input <- set_osa_obs(input)
    cat(unlist(input$log$ecov, recursive=T))
  }
	#print("osa_obs")

	# projection data will always be modified by 'prepare_projection'
	input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
	#print("proj")

	#set any parameters as random effects
	input$random = NULL
	input = set_random(input)
	#print("random")
  input$options$ecov <- ecov
	
  return(input)
}
