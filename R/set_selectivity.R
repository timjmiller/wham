#' Specify model and parameter configuration for selectivity
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param selectivity (optional) list specifying options for selectivity blocks, models, initial parameter values, parameter fixing/mapping, and random effects (see details)
#'
#' \code{set_selectivity} specifies options for selectivity and allows you to overwrite existing options 
#' in the \code{input} list or as specified in the ASAP data file. If \code{selectivity = NULL}, selectivity options from 
#' \code{input} are used. 
#'
#' \code{\link{prepare_wham_input}(..., selectivity=selectivity)} calls \code{set_selectivity(..., selectivity=selectivity)}.
#' If you already have created \code{input} with \code{prepare_wham_input}, you can also use \code{set_selectivity(input, selectivity=selectivity)}
#' to modify the selectivity specification.
#'
#' \code{selectivity} is a list with the following entries:
#'   \describe{
#'     \item{$n_selblocks}{How many selectivity blocks. Optional. If unspecified and no asap3 object, then this is set to the number 
#'                of fleets + indices. If specified, ensure other components of \code{selectivity} are consistent.}
#'     \item{$model}{Selectivity model for each block. Vector with length = number of selectivity blocks. Each entry must be one of: "age-specific", "logistic", "double-logistic", or "decreasing-logistic".}
#'     \item{$re}{Time-varying (random effects) for each block. Vector with length = number of selectivity blocks.
#'                  If \code{NULL}, selectivity parameters in all blocks are constant over time and uncorrelated.
#'                  Each entry of \code{selectivity$re} must be one of the following options, where the selectivity parameters are:
#'                  \describe{
#'                    \item{"none"}{(default) are constant and uncorrelated}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by age/par (AR1), but not year}
#'                    \item{"ar1_y"}{correlated by year (AR1), but not age/par}
#'                    \item{"2dar1"}{correlated by year and age/par (2D AR1)}
#'                  }
#'                 }
#'     \item{$initial_pars}{Initial parameter values for each block. List of length = number of selectivity blocks. Each entry must be
#'                a vector of length # parameters in the block, i.e. \code{c(2,0.2)} for logistic (a50 and 1/slope) or \code{c(0.5,0.5,0.5,1,1,0.5)} for
#'                age-specific parameters when there are 6 ages. Default is to set at middle of parameter range. This is 0.5 for age-specific and n.ages/2 
#'                or logistic, double-logistic, and decreasing-logistic.}
#'     \item{$fix_pars}{Alternative to \code{$map_pars} for specifying which selectivity parameters (only fixed effects) to fix at initial values. 
#'                List of length = number of selectivity blocks. E.g. model with 3 age-specific blocks and 6 ages, 
#'                \code{list(4:5, 4, 2:4))} will fix ages 4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in block 3.
#'                Use NULL to not fix any parameters for a block, e.g. list(NULL, 4, 2) does not fix any pars in block 1.}
#'     \item{$par_min}{The lower bound for selectivity parameters and is used to populate \code{data$selpars_lower}.
#'                List of length = number of selectivity blocks, where each item is a 
#'                vector of length = number of selectivity parameters (age-specific: n.ages, logistic: 2, 
#'                double-logistic: 4).}
#'     \item{$par_max}{The upper bound for selectivity parameters and is used to populate \code{data$selpars_upper}.
#'                List of length = number of selectivity blocks, where each item is a 
#'                vector of length = number of selectivity parameters (age-specific: n.ages, logistic: 2, 
#'                double-logistic: 4).}
#'     \item{$map_pars}{Alternative to \code{$fix_pars} for specifying how to fix selectivity parameters (only fixed effects), corresponds 
#'                to \code{map$logit_selpars}. List of length = number of selectivity blocks, where each item is a 
#'                vector of length = number of selectivity parameters (age-specific: n.ages, logistic: 2, 
#'                double-logistic: 4). Use \code{NA} to fix a parameter and integers to estimate. Use the same integer
#'                for multiple ages or fleets/indices to estimate a shared parameter. E.g. for a model with 3 age-specific 
#'                blocks (1 fleet, 2 indices) and 6 ages, \code{$map_pars = list(c(1,2,3,NA,NA,4), c(5,6,7,NA,8,8), c(1,2,3,NA,NA,4))}
#'                will estimate ages 1-3 and 6 in block 1 (fleet), ages 1-3 and 4-5 (shared) in block 2 (index 1), and then set the
#'                index 2 (block 3) selectivity equal to the fleet.}
#'     \item{$sigma_vals}{Initial standard deviation values to use for the random effect deviations. Must be a vector with length = number of blocks. 
#'                Use natural (not log) scale, must be positive. \code{par$sel_repars[,1]} will be estimated on log-scale. Not used if \code{re = 'none'} for all blocks.}
#'     \item{$map_sigma}{Specify which SD parameters to fix for the random effect deviations. Must be a vector with length = number of blocks. 
#'                Use \code{NA} to fix a parameter and integers to estimate. Use the same integer for multiple blocks to estimate a shared SD parameter.
#'                Not used if \code{re = 'none'} for all blocks.}
#'     \item{$cor_vals}{Initial correlation values to use for the random effect deviations. Must be a n_selblocks x 2 integer matrix. Columns correspond to correlation by age and year,
#'                respectively. If \code{re = 'ar1'} or \code{re = 'ar1_y'} only the corresponding values are used.
#'                Values must be between -1 and 1, but parameters are estimated on a logit transformed scale internally. 
#'                Not used if \code{re = 'none'} or \code{re = 'iid'} for all blocks.}
#'     \item{$map_cor}{Specify which correlation parameters to fix for the random effect deviations. Must be a n_selblocks x 2 matrix. Columns correspond to correlation by age
#'                and year,respectively. Parameters can be shared by setting corresponding values of \code{map_cor} to the same integer. Use \code{NA} to fix a parameter.
#'                If \code{re = 'ar1'} or \code{re = 'ar1_y'}, only the column for the corresponding correlation are used.
#'                Not used if \code{re = 'none'} or \code{re = 'iid'} for all blocks.}
#'     \item{$map_re}{Specify whether and how to estimate selectivity random effects. Structure is the same as $map_pars: list of length = number of selectivity blocks, where each item is a 
#'                vector of length = number of selectivity parameters (age-specific: n.ages, logistic: 2, 
#'                double-logistic: 4). Use \code{NA} to not estimate random effects for a corresponding mean parameter and integers to estimate. Use the same integer
#'                for multiple ages to estimate shared parameters within a selectivity block. For \code{re = "ar1_y" or "2dar1"}, the annual RE are unique, but identical for parameters specified to be equal. E.g. for selectivity block 1 with age-specific 
#'                selectivity and annual random effect (and n_ages = 6), \code{$map_pars[[1]] = c(1,2,2,NA,NA,4)}
#'                will estimate 2 sets of unique annual RE for ages 1, 6, and another unique set applied to both ages 2 and 3.
#'                No RE will be estimated for any mean parameters that are at bounds (\code{par_min} or \code{par_max}). Lack of convergence should be expected when there is only 1 unique set of annual RE and correlation across parameters is estimated (i.e., \code{re = "2dar1"}) or when the ratio of unique RE and mean parameters is low without annual RE (i.e., \code{re = "ar1"}).}
#'   }
#'
#' @return a named list with same elements as the input provided with selectivity options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3, NAA_re = list(sigma = "rec"))
#' sel <- list(model=rep("logistic",input$data$n_selblocks),
#'    initial_pars=rep(list(c(3,3)),input$data$n_selblocks),
#'    fix_pars=rep(list(NULL),input$data$n_selblocks))
#' input <- set_selectivity(input, selectivity = sel) #logistic selectivity for all selectivity blocks
#' }
#'
#' @export
set_selectivity <- function(input, selectivity){
  data <- input[["data"]]
  if(is.null(data)) data <- list()
  par <- input$par
  if(is.null(par)) par <- list()
  map <- input$map
  if(is.null(map)) map <- list()
  
  par_index <- list(
    1:data$n_ages,
    data$n_ages + 1:2,
    data$n_ages + 3:6,
    data$n_ages + 1:2
  )
  selopts <- c("age-specific","logistic","double-logistic","decreasing-logistic")
  default_selpars <- list(rep(0.5,data$n_ages),rep(data$n_ages/2,2),rep(data$n_ages/2,4),rep(data$n_ages/2,2)) #For age-specfic, need to fix one par at 1.
  n_selpars <- c(data[["n_ages"]],2,4,2)
  asap3 <- input$asap3
  if(is.null(input$use_asap3)){
    input$use_asap3 <- FALSE
  }
  if(input$use_asap3 & is.null(asap3)) input$use_asap3 <- FALSE
  estimate_selpars <- selpars_ini <- NULL
  if(input$use_asap3) {
    x <- set_asap3_selectivity(asap3, input) 
    input <- x$input
    selpars_ini <- x$selpars_ini
    estimate_selpars <- x$estimate_selpars
  }

  input$log$selectivity <- list()

  if(is.null(selectivity$n_selblocks) & is.null(data$n_selblocks)) {
    #selblock pointers are defined upstream in set_indices and set_catch
    data$n_selblocks <- max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices)
    input$log$selectivity <- c(input$log$selectivity, paste0("selectivity$n_selblocks was not provided so number of selblocks: ", data$n_selblocks, 
      ", is being determined by max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices).\n"))
  }

  if(!is.null(selectivity$n_selblocks)){ #override asap structure
    data$n_selblocks <- selectivity$n_selblocks
  }
  if(is.null(data$selblock_models) & is.null(selectivity$model)) {
    data$selblock_models <- rep(2, data$n_selblocks)
    input$log$selectivity <- c(input$log$selectivity, paste0("selectivity$selblock_models was not provided so logistic selectivity is being set for mean models of all ", 
      data$n_selblocks, " selblocks.\n"))
  }
  if(!is.null(selectivity$model)){
    if(length(selectivity$model) != data$n_selblocks) stop(paste0("Length of selectivity$model: ", length(selectivity[["model"]])," must equal number of selectivity blocks: ", data[["n_selblocks"]]))
    if(!all(selectivity$model %in% selopts)) stop("Each model entry must be one of the following: 'age-specific','logistic','double-logistic','decreasing-logistic'")
    data$selblock_models <- match(selectivity$model, selopts)
  }
  input$log$selectivity <- c(input$log$selectivity, paste0("(Mean) selectivity block models are:\n",
    paste0("Block ", 1:data$n_selblocks, ": ", selopts[data$selblock_models], collapse = "\n"), "\n\n")
  )

  ##########################################################################################################################
  # Starting/Fixed values, and upper and lower bounds for mean selectivity parameters

  if(is.null(selpars_ini)) { #wasn't provided by set_asap3_selectivity()
    selpars_ini <- matrix(NA, data$n_selblocks, data$n_ages + 6) 
    for(b in 1:data$n_selblocks) selpars_ini[b,par_index[[data$selblock_models[b]]]] <- default_selpars[[data$selblock_models[b]]]
  }

  # initial values
  if(!is.null(selectivity$initial_pars)){ #override anything from asap3 file or defaults
    if(!is.list(selectivity$initial_pars)) stop("selectivity$initial_pars must be a list of length = number of selectivity blocks.")
    if(length(selectivity$initial_pars) != data$n_selblocks) stop(paste0("Length of selectivity$initial_pars: ", length(selectivity[["initial_pars"]])," must equal number of selectivity blocks: ", data[["n_selblocks"]]))
    for(b in 1:data$n_selblocks){
      if(length(selectivity$initial_pars[[b]]) != n_selpars[data$selblock_models[b]]) stop(paste0("Length of vector for selectivity block ",b," in the selectivity$initial_pars list is not equal to the number of selectivity parameters for block ",b,": ",n_selpars[data$selblock_models[b]]))
      selpars_ini[b,par_index[[data$selblock_models[b]]]] <- selectivity$initial_pars[[b]]
    }
  }

  # initial values on logit scale, par$logit_selpars
  selpars_lo <- selpars_hi <- matrix(0, data$n_selblocks, data$n_ages + 6)
  selpars_hi[,1:data$n_ages] <- 1
  selpars_hi[,data$n_ages + 1:6] <- data$n_ages
  
  #user-specified upper and lower bounds
  if(!is.null(selectivity$par_min)){
    for(b in 1:data$n_selblocks) {
      if(any(selectivity$par_min[[b]]<0)) stop(paste0("At least one of the lower bounds provided for selectivity block ", b,": ", , " is less than 0. Parameters must be positive.\n"))
      if(any(selectivity$par_min[[b]] > selpars_ini[b, par_index[[data$selblock_models[b]]]])){
        stop(paste0("At least one of the lower bounds provided for selectivity block ", b,": ", , " is greater than the corresonding initial value.\n",
        paste0("Initial values ", selpars_ini[b, par_index[[data$selblock_models[b]]]], ", and par_min: ", selectivity$par_min[[b]], collapse = "\n"), "\n\n"))
      }
      selpars_lo[b, par_index[[data$selblock_models[b]]]] <- selectivity$par_min[[b]]
    }
  }
  if(!is.null(selectivity$par_max)){
    for(b in 1:data$n_selblocks) {
      if(any(selectivity$par_max[[b]]<0)) stop(paste0("At least one of the upper bounds provided for selectivity block ", b,": ", , " is less than 0. Parameters must be positive.\n"))
      if(any(selectivity$par_max[[b]] < selpars_ini[b, par_index[[data$selblock_models[b]]]])){
        stop(paste0("At least one of the upper bounds provided for selectivity block ", b,": ", , " is less than the corresonding initial value.\n",
        paste0("Initial values ", selpars_ini[b, par_index[[data$selblock_models[b]]]], ", and par_max: ", selectivity$par_min[[b]], collapse = "\n"), "\n\n"))
      }
      selpars_hi[b, par_index[[data$selblock_models[b]]]] <- selectivity$par_max[[b]]
    }
  }

  #if selpars_ini are still outside of bounds, then set them at the closest bound
  selpars_ini[which(selpars_ini > selpars_hi)] <- selpars_hi[which(selpars_ini > selpars_hi)]
  selpars_ini[which(selpars_ini < selpars_lo)] <- selpars_lo[which(selpars_ini < selpars_lo)]
  data$selpars_lower <- selpars_lo #only need these for estimated parameters
  data$selpars_upper <- selpars_hi
  
  ##########################################################################################################################
  #MAPPING mean selectivity parameters, (but also checking for 0 observations to set age-specific pars to 0, so par$logit_selpars not set yet.)

  if(!is.null(selectivity$model) | is.null(estimate_selpars)){
    #redefining which parameters to estimate since selblock models may be changing or there was no asap3 input
    estimate_selpars <- matrix(0, data$n_selblocks, data$n_ages + 6)
    for(b in 1:data$n_selblocks) estimate_selpars[b,par_index[[data$selblock_models[b]]]] <- 1
  }

  map_logit_selpars <- matrix(NA, data$n_selblocks, data$n_ages + 6)
  map_logit_selpars[which(estimate_selpars > 0)] <- 1:sum(estimate_selpars>0)
  # which selpars to fix, either $fix_pars or $map_pars
  if(!is.null(selectivity$fix_pars) & !is.null(selectivity$map_pars)) {
    stop("Cannot specify $fix_pars and $map_pars (both set which pars to estimate). Choose one or the other.\n")
  }
  # $fix_pars
  if(!is.null(selectivity$fix_pars)){
    if(!is.list(selectivity$fix_pars)) stop("selectivity$fix_pars must be a list.\n")
    if(length(selectivity$fix_pars) != data$n_selblocks) stop(paste0("Length of selectivity$fix_pars: ",  length(selectivity$fix_pars), " must equal number of selectivity blocks: ", data$n_selblocks,".\n
      Use 'NULL' to not fix any parameters for a block, e.g. list(NULL,4,2) does not fix any pars in block 1.\n"))
    for(b in 1:data$n_selblocks){
      if(data$selblock_models[b] == 1) map_logit_selpars[b,selectivity$fix_pars[[b]]] <- NA
      if(data$selblock_models[b] %in% c(2,4)) map_logit_selpars[b,data$n_ages+selectivity$fix_pars[[b]]] <- NA
      if(data$selblock_models[b] == 3) map_logit_selpars[b,data$n_ages+2+selectivity$fix_pars[[b]]] <- NA
    }
  }
  # $map_pars 
  if(!is.null(selectivity$map_pars)){
    if(!is.list(selectivity$map_pars)) stop("selectivity$map_pars must be a list.\n")
    if(length(selectivity$map_pars) != data$n_selblocks) stop(paste0("Length of selectivity$map_pars", length(selectivity$map_pars), " must equal number of selectivity blocks ", data$n_selblocks,".\n"))    
    #start map from scratch
    map_logit_selpars <- matrix(NA, data$n_selblocks, data$n_ages + 6)
    all_map <- unlist(selectivity$map_pars)
    unique_selpars <- unique(all_map[!is.na(all_map)])
    for(b in 1:data$n_selblocks) {
      if(length(selectivity$map_pars[[b]] != length(par_index[[data$selblock_models[b]]]))) stop(paste0("Length of Length of selectivity$map_pars[[", b, 
        "]] is not equal to the number required by the specified selectivity model: ", selopts[data$selblock_models[b]], ".\n"))
      for(ind in unique_selpars) {
        ind_b <- par_index[[data$selblock_models[b]]][which(selectivity$map_pars[[b]] == ind)]
        if(length(ind_b)) map_logit_selpars[b,ind_b] <- ind
      }
    }
  }
  # For age-specific selectivity blocks, check for ages with ~zero catch and fix these at 0
  age_specific <- which(data$selblock_models==1)
  for(b in age_specific){
    ind <- list(fleets = which(apply(data$selblock_pointer_fleets == b,2,sum) > 0))
    ind$indices <- which(apply(data$selblock_pointer_indices == b,2,sum) > 0)
    paa <- matrix(nrow = 0, ncol = data$n_ages)
    if(length(ind$fleets)) for(f in ind$fleets) {
      y <- data$catch_paa[f,which(data$selblock_pointer_fleets[,f] == b & data$use_catch_paa[,f] == 1),]
      paa <- rbind(paa,y)
    }
    if(length(ind$indices)) for(i in ind$indices) {
      y <- data$index_paa[i,which(data$selblock_pointer_indices[,i] == b & data$use_index_paa[,i] == 1),]
      paa <- rbind(paa,y)
    }
    y <- apply(paa,2,sum)
    ind <- which(y < 1e-5 & !is.na(map_logit_selpars[b,1:data$n_ages])) #if phase is set to estimate par and 0s in all years for an age, then change input to fix those pars at 0.
    selpars_ini[b, ind] <- 0
    map_logit_selpars[b, ind] <- NA
    if(length(ind)){
      input$log$selectivity <- c(input$log$selectivity, paste0("Age-specific mean model set for selectivity block ", b, 
        ", but there are no observations > 0 for some ages corresponding to parameters specified to be estimated so they have been fixed at 0.\n"))
    }
  }
  par$logit_selpars <- log(selpars_ini-data$selpars_lower) - log(data$selpars_upper - selpars_ini)
  par$logit_selpars[!is.na(map_logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] <- -10
  par$logit_selpars[!is.na(map_logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] <- 10  

  ##########################################################################################################################
  #RANDOM EFFECTS:
  
  selblock_pointers <- cbind(data[["selblock_pointer_fleets"]], data[["selblock_pointer_indices"]])
  if(is.null(data[["selblock_years"]]) | !identical(dim(data[["selblock_years"]]), as.integer(c(data[["n_years_model"]], data[["n_selblocks"]])))) {
    data[["selblock_years"]] <- matrix(0L, data[["n_years_model"]], data[["n_selblocks"]])
    #for each selblock find which years it is used in any fleet or index
    for(b in 1:data[["n_selblocks"]]) data[["selblock_years"]][,b] <- apply(selblock_pointers, 1, function(x) b %in% x)
  }
  data[["n_years_selblocks"]] <- apply(data[["selblock_years"]], 2, sum)

  data[["selblock_models_re"]] <- rep(1, data[["n_selblocks"]]) # default: no RE on selectivity parameters
  if(!is.null(selectivity[["re"]])){
    if(length(selectivity[["re"]]) != data[["n_selblocks"]]) stop(paste0("Length of selectivity$re: ", length(selectivity[["re"]])," must equal number of selectivity blocks: ", data[["n_selblocks"]],".\n"))
    if(!all(selectivity[["re"]] %in% c("none","iid","ar1","ar1_y","2dar1"))) stop("Each selectivity$re entry must be one of the following: 'none','iid','ar1','ar1_y','2dar1'.\n")
    data$selblock_models_re <- match(selectivity[["re"]], c("none","iid","ar1","ar1_y","2dar1"))
  } 

  input$log$selectivity <- c(input$log$selectivity, paste0("Random effects options for each selectivity block are:\n",
    paste0("Block ", 1:data$n_selblocks, ": ", c("none","iid","ar1","ar1_y","2dar1")[data$selblock_models_re], collapse ="\n"), "\n\n")
  )
  
  for(b in 1:data$n_selblocks){ 
    if(data$selblock_models_re[b] == 3){ #ar1(age)
      if(data$selblock_models[b] != 1) stop(paste0("'ar1' (AR1(age)) random effects specified for selectivity block ",b,", but the mean model is not specified as 'age-specific'.\n"))
      pars_map <- map_logit_selpars[b, par_index[[data$selblock_models[b]]]]
      #only allow random effects for mean parameters that are not (fixed at) at bounds
      pars_with_re <- which(!is.infinite(par$logit_selpars[b,1:data[["n_ages"]]]))
      if(length(pars_with_re) < 3) {
        input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) RE with age-specific selectivity are specified for block ",b,", 
          but number of mean parameters that will have RE is <= 2, which may not be estimable if the RE variance/correlation parameters are estimated. \n",
          "RE are not allowed for mean parameters set at their bounds (0 or 1) so perhaps redefine those initial values.\n"))
      }
      unique_pars <- unique(pars_map[!is.na(pars_map)])
      if(length(unique_pars)>1){ #more than 1 mean parameter. multiple mean pars may make convergence difficult
        if(!is.null(selectivity$map_pars)) {
          input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) RE with age-specific selectivity are specified for block ",b,
          ", but there are multiple mean age-specific parameters being estimated which may not allow convergence.\n", 
          "Consider configuring $map_pars for this block to estimate just one mean parameter shared across ages that are not fixed at 1 or 0(there are at least 3 ages that are free here).\n"))
        } else { #map_pars not used, so do default behavior
          input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) RE with age-specific selectivity are specified for block ",b,
          ", $map_pars was not specified so default behavior is to estimate just one mean parameter shared across ages that are not fixed at 1 or 0.\n"))
          # but don't overwrite fixed pars (likely will be fixing one age at 1)
          pars_map[!is.na(pars_map)] <- min(pars_map, na.rm=TRUE)
          map_logit_selpars[b, par_index[[data$selblock_models[b]]]] <- pars_map
        }
      }
      # warning message if no mean sel pars (logit_selpars) are fixed
      # allow so user can fit model without fixing and then fix the age with highest sel at 1
      if(all(!is.na(pars_map))) {
        input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) RE with age-specific selectivity are specified for block ",b,
        ", but no age fixed at 1 which may not converge. Advised to fit the current model and then fix the age with highest selectivity at 1. 
        Can use selectivity$fix_pars.\n"))
      }
      # if there are few ages ages, estimating ar1(age) may be difficult
      if(data[["n_ages"]]<6) {
        input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) RE with age-specific selectivity are specified for block ",b,
        ", but the number of age classes for the population is less than 6 which may make converge difficult.\n"))
      }
    }
    if(data$selblock_models_re[b] == 4){ #ar1(year)
      if(data[["n_years_selblocks"]][b]<10){
        input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1_y' (AR1(year)) RE spedified for block ",b,
        ", but the number of years where this block is used is less than 10 which may make converge difficult.\n"))
      }
    }
  }
  #NOTE: n_selpars_re and selpare_re_index are used to configure RE on TMB side
  par$selpars_re <- array(0, c(data$n_selblocks, data$n_years_model, data$n_ages+6))
  map$selpars_re <- array(NA, c(data$n_selblocks, data$n_years_model, data$n_ages+6))
  data$selpars_re_index <- matrix(0, data$n_selblocks, max(data$n_ages,4)) #configured like M_re_index
  data$n_selpars_re <- rep(1,data[["n_selblocks"]]) #configured like n_M_re
  ct <- 0
  for(b in 1:data$n_selblocks){
    #if map_re not provided, just make RE for any parameter not fixed at bounds
    par_index_b <- par_index[[data$selblock_models[b]]]
    is_re <- !is.infinite(par$logit_selpars[b,par_index_b])
    if(data$selblock_models_re[b] > 1) {
      if(!sum(is_re)) stop(paste0("RE are specified for selectivity block ", b, ", but all mean parameters are set at upper or lower bounds."))
    }
    ind_y <- which(data[["selblock_years"]][,b]==1)
    if(data$selblock_models_re[b] == 4) { #ar1(year)
      data$n_selpars_re[b] <- 1
      data$selpars_re_index[b,which(is_re)] <- 1 #default is 1 common RE for all selpars
      for(s in which(is_re)) map$selpars_re[b,ind_y,par_index_b[s]] <- ct + 1:data[["n_years_selblocks"]][b]
      ct <- ct + data[["n_years_selblocks"]][b]
    }
    if(data$selblock_models_re[b] == 3) { #ar1(age)
      data$n_selpars_re[b] <- length(which(is_re))
      data$selpars_re_index[b,which(is_re)] <- 1:data$n_selpars_re[b]
      for(s in which(is_re)) {
        map$selpars_re[b,ind_y,par_index_b[s]] <- ct + 1
        ct <- ct + 1
      }
    }
    if(data$selblock_models_re[b] %in% c(2,5)) { #2d iid or 2d ar1
      data$n_selpars_re[b] <- length(which(is_re))
      data$selpars_re_index[b,which(is_re)] <- 1:data$n_selpars_re[b]
      map$selpars_re[b,ind_y,par_index_b[which(is_re)]] <- ct + 1:(data[["n_years_selblocks"]][b]*sum(is_re))
      ct <- ct + (data[["n_years_selblocks"]][b]*sum(is_re))
    }
  }

  ct <- 0
  if(!is.null(selectivity$map_re)){
    #this can override data$selblock_models_re[b] == 4
    if(!is.list(selectivity$map_re)) stop("selectivity$map_re must be a list.\n")
    if(length(selectivity[["map_re"]]) != data[["n_selblocks"]]) stop(paste0("Length of selectivity$map_re: ", length(selectivity[["map_re"]])," must equal number of selectivity blocks: ", data[["n_selblocks"]],".\n"))
    for(b in 1:data[["n_selblocks"]]) if(data$selblock_models_re[b]>1){
      if(length(selectivity$map_re[[b]] != length(par_index[[data$selblock_models[b]]]))) stop(paste0("Length of selectivity$map_re[[", b, 
        "]] is not equal to the number required by the specified selectivity model: ", selopts[data$selblock_models[b]], ".\n"))
      par_index_b <- par_index[[data$selblock_models[b]]]
      if(any(!is.na(selectivity$map_re[[b]]) & is.infinite(par[["logit_selpars"]][b,par_index_b]))){
        stop(paste0("selectivity$map_re for block ", b, " specifies RE (not NA), but there are corresponding mean parameters that are set at their bounds which is not allowed."))
      }
      if(all(is.na(selectivity$map_re[[b]]))) stop(paste0("Random effects have been specified for selectivity model: ", b, ", but selectivity$map_re[[b]] are all NA."))
      
      ind_re <- unclass(factor(selectivity$map_re[[b]]))
      
      unique_ind_re <- unique(ind_re[which(!is.na(ind_re))])
      data$selpars_re_index[b,] <- ind_re
      data$n_selpars_re[b] <- length(unique_ind_re)
      ind_y <- which(data[["selblock_years"]][,b]==1)
      if(data$selblock_models_re[b] == 3){ #ar1(age): constant across years
        for(y in ind_y){
          map$selpars_re[b,y,] <- ind_re
        }
        ct <- ct + data[["n_selpars_re"]]
      }
      if(data$selblock_models_re[b] %in% c(2,4,5)){ #ar1(year): constant across age, 2d(iid), 2d(ar1). Note ar1(year) and 2d(iid) with $map_re can configure equivalent models
        for(y in ind_y){
          map$selpars_re[b,y,] <- ind_re + ct
          ct <- ct + data[["n_selpars_re"]]
        }
      }
    }
  }

  map$selpars_re <- factor(map$selpars_re)
  
  # initial and map for parameters controlling selectivity RE
  # default initial values: sigma = 0.1, rho = 0
  par$sel_repars <- matrix(0, nrow=data$n_selblocks, ncol=3)
  par$sel_repars[,1] <- log(0.1)
  if(!is.null(selectivity$sigma_vals)){
    if(any(selectivity$sigma_vals < 0)) stop('Variance controlling selectivity random effects must be positive.') 
    par$sel_repars[,1] <- log(selectivity$sigma_vals) # log scale
  }
  if(!is.null(selectivity$cor_vals)){
    if(!is.matrix(selectivity$cor_vals)) stop("selectivity$cor_vals must be an n_selbocks x 2 matrix.")
    if(!all(dim(selectivity$cor_vals) == c(data$n_selblocks,2))) stop("selectivity$cor_vals must be an n_selbocks x 2 matrix.")
    if(any(abs(selectivity$cor_vals)>=1)) stop('|Correlation parameters|<1 is required.')
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) par$sel_repars[b,2] <- gen.logit(selectivity$cor_vals[b,1], -1, 1) # if ar1 over ages, use specified initial
      if(data$selblock_models_re[b] == 4) par$sel_repars[b,3] <- gen.logit(selectivity$cor_vals[b,2], -1, 1) # if ar1 over years, use specified initial
      if(data$selblock_models_re[b] == 5) par$sel_repars[b,2:3] <- gen.logit(selectivity$cor_vals[b,1:2], -1, 1) # if 2dar1 over years, use both
    }
  }  

  # map selectivity RE
  tmp.sel.repars <- matrix(NA, data$n_selblocks, 3)
  if(!is.null(selectivity$map_sigma)){
    if(length(selectivity$map_sigma) != data$n_selblocks) stop("selectivity$map_sigma must be a vector of length = number of selectivity blocks")
    for(b in 1:data$n_selblocks) if(data$selblock_models_re[b]>1) tmp.sel.repars[b,1] <- selectivity$map_sigma[b]
  } else{ #default mapping: unique values for each selblock
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] > 1) tmp.sel.repars[b,1] <- max(c(0,tmp.sel.repars), na.rm = T) + 1 # estimate sigma
    }
  }
  if(!is.null(selectivity$map_cor)){
    if(!is.matrix(selectivity$map_cor)) stop("selectivity$map_cor must be an n_selbocks x 2 matrix.")
    if(!all(dim(selectivity$map_cor) == c(data$n_selblocks,2))) stop("selectivity$map_cor must be an n_selbocks x 2 matrix.")
    selectivity$map_cor[] <- as.interger(selectivity$map_cor)
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,2] <- selectivity$map_cor[b,1]
      if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,3] <- selectivity$map_cor[b,2] 
      if(data$selblock_models_re[b] == 5) tmp.sel.repars[b,2:3] <- selectivity$map_cor[b,1:2]
    }
  } else { #default mapping: unique values for each selblock
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,2] <- max(c(0,tmp.sel.repars), na.rm = T) + 1 # estimate sigma, rho
      if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,3] <- max(c(0,tmp.sel.repars), na.rm = T) + 1 # estimate sigma, rho_y
      if(data$selblock_models_re[b] == 5) tmp.sel.repars[b,2:3] <- max(c(0,tmp.sel.repars), na.rm = T) + 1:2 # estimate sigma, rho, rho_y
    }
  }
  map$sel_repars <- factor(tmp.sel.repars)
  map$logit_selpars <- factor(map_logit_selpars)

  input$data <- data
  input$par <- par
  input$map <- map
	if(length(input$log$selectivity))	input$log$selectivity <- c(
    "--Selectivity------------------------------------------------------------------------------------------------------------------------",
    "\n", input$log$selectivity, 
    "-------------------------------------------------------------------------------------------------------------------------------------",
    "\n\n")
 
  if(!is.null(input$options$catch) & !is.null(input$data$selblock_pointer_fleets)) {
    if(any(!input$data$selblock_pointer_fleets %in% 1:input$data$n_selblocks)){
    input$log$catch <- c(input$log$catch, 
      paste0("NOTE: set_catch has previously been called, but some fleet selblock_pointers are outside of the number of selectivity blocks currently specified.\n",
      "\tMake sure to run set_catch and set_selectivity with appropriate options before using the input.\n"))
    }
  }
  if(!is.null(input$options$index) & !is.null(input$data$selblock_pointer_indices)) {
    if(any(!input$data$selblock_pointer_indices %in% 1:input$data$n_selblocks)){
    input$log$index <- c(input$log$index, 
      paste0("NOTE: set_indices has previously been called, but some index selblock_pointers are outside of the number of selectivity blocks currently specified.\n",
      "\tMake sure to run set_catch and set_selectivity with appropriate options before using the input.\n"))
    }
  }
 
  #set any parameters as random effects
  input$random <- NULL
  input <- set_random(input)
  input$options$selectivity <- selectivity
  if(is.null(input[["by_pwi"]])) { #check whether called by prepare_wham_input
    input <- set_age_comp(input, input[["options"]][["age_comp"]])
    input <- set_osa_obs(input) #check for ages to omit due to any selectivity = 0
    message(unlist(input[["log"]][["selectivity"]], recursive=T))
  }
  return(input)

}


set_asap3_selectivity <- function(asap3, input){
  data <- input$data
  data$selblock_models <- integer(0)
  data$n_selblocks <- 0
  par_index <- list(
    1:data$n_ages,
    data$n_ages + 1:2,
    data$n_ages + 3:6,
    data$n_ages + 1:2
  )
  for(i in 1:length(asap3)){ #have to do all fleets first, then indices
    data$n_selblocks <- data$n_selblocks + asap3[[i]]$n_fleet_sel_blocks
    data$selblock_models <- c(data$selblock_models, asap3[[i]]$sel_block_option)
  }
  for(i in 1:length(asap3)){
    which_indices <- which(asap3[[i]]$use_index ==1) # length <- data$n_indices
    asap3[[i]]$index_sel_option <- asap3[[i]]$index_sel_option[which_indices]
    asap3[[i]]$index_sel_ini <- asap3[[i]]$index_sel_ini[which_indices]
    data$n_selblocks <- data$n_selblocks + length(which_indices)
    data$selblock_models <- c(data$selblock_models, asap3[[i]]$index_sel_option)
  }
  
  selpars_ini <- matrix(NA, data$n_selblocks, data$n_ages + 6)
  estimate_selpars <- matrix(0, data$n_selblocks, data$n_ages + 6)
  j <- 1
  for(k in 1:length(asap3)) for(i in 1:asap3[[k]]$n_fleet_sel_blocks) {
    selpars_ini[j,] <- asap3[[k]]$sel_ini[[i]][,1]
    which_sel_model <- asap3[[k]]$sel_block_option[i]
    which_estimate <- which(asap3[[k]]$sel_ini[[i]][par_index[[which_sel_model]],2] > 0)
    estimate_selpars[j,par_index[[which_sel_model]][which_estimate]] <- 1
    j <- j + 1
  }
  for(k in 1:length(asap3)) for(i in 1:length(asap3[[k]]$index_sel_ini)){
    selpars_ini[j,] <- asap3[[k]]$index_sel_ini[[i]][,1]
    which_sel_model <- asap3[[k]]$index_sel_option[i]
    which_estimate <- which(asap3[[k]]$index_sel_ini[[i]][par_index[[which_sel_model]],2] > 0)
    estimate_selpars[j,par_index[[which_sel_model]][which_estimate]] <- 1
    j <- j + 1
  }
  input$data <- data
  return(list(input = input, selpars_ini = selpars_ini, estimate_selpars = estimate_selpars))
}

