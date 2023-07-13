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
#'                a vector of length # parameters in the block, i.e. \code{c(2,0.2)} for logistic or \code{c(0.5,0.5,0.5,1,1,0.5)} for
#'                age-specific with 6 ages. Default is to set at middle of parameter range. This is 0.5 for age-specific and n.ages/2 
#'                or logistic, double-logistic, and decreasing-logistic.}
#'     \item{$fix_pars}{Alternative to \code{$map_pars} for specifying which selectivity parameters (only fixed effects) to fix at initial values. 
#'                List of length = number of selectivity blocks. E.g. model with 3 age-specific blocks and 6 ages, 
#'                \code{list(4:5, 4, 2:4))} will fix ages 4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in block 3.
#'                Use NULL to not fix any parameters for a block, e.g. list(NULL, 4, 2) does not fix any pars in block 1.}
#'     \item{$map_pars}{Alternative to \code{$fix_pars} for specifying how to fix selectivity parameters (only fixed effects), corresponds 
#'                to \code{map$logit_selpars}. List of length = number of selectivity blocks, where each item is a 
#'                vector of length = number of selectivity parameters (age-specific: n.ages, logistic: 2, 
#'                double-logistic: 4). Use \code{NA} to fix a parameter and integers to estimate. Use the same integer
#'                for multiple ages or fleets/indices to estimate a shared parameter. E.g. for a model with 3 age-specific 
#'                blocks (1 fleet, 2 indices) and 6 ages, \code{$map_pars = list(c(1,2,3,NA,NA,4), c(5,6,7,NA,8,8), c(1,2,3,NA,NA,4))}
#'                will estimate ages 1-3 and 6 in block 1 (fleet), ages 1-3 and 4-5 (shared) in block 2 (index 1), and then set the
#'                index 2 (block 3) selectivity equal to the fleet.}
#'     \item{$initial_sigma}{Initial standard deviation values to use for the random effect deviations. Must be a vector with length = number of blocks. 
#'                Use natural (not log) scale, must be positive. \code{par$sel_repars[,1]} will be estimated on log-scale. Not used if \code{re = 'none'} for all blocks.}
#'     \item{$map_sigma}{Specify which SD parameters to fix for the random effect deviations. Must be a vector with length = number of blocks. 
#'                Use \code{NA} to fix a parameter and integers to estimate. Use the same integer for multiple blocks to estimate a shared SD parameter.
#'                Not used if \code{re = 'none'} for all blocks.}
#'     \item{$initial_cor}{Initial correlation values to use for the random effect deviations. Must be a list with length = number of blocks. If \code{re = 'ar1'} or \code{ar1_y'}
#'                for a block, list element must be a single value. If \code{re = '2dar1'} for a block, list element must be a vector of length 2 (first is for "age", second is for "year") for each block. 
#'                Use natural scale, must be between -1 and 1. \code{par$sel_repars[,2:3]} will be estimated on a transform scale,  (2 / (1 + exp(-2x))) - 1. 
#'                Not used if \code{re = 'none'} or \code{re = 'iid'} for all blocks.}
#'     \item{$map_cor}{Specify which correlation parameters to fix for the random effect deviations. Must be a list with length = number of blocks. 
#'                If \code{re = 'ar1'} or \code{ar1_y'}, each list element (block) must be a single value. If \code{re = '2dar1'}, 
#'                must be a vector of length 2 (first is for "age", second is for "year"). Use \code{NA} to fix a parameter and integers to estimate. 
#'                Use the same integer for multiple blocks to estimate a shared correlation parameter. Not used if \code{re = 'none'} or \code{re = 'iid'} for all blocks.}
#'     \item{$n_selblocks}{How many selectivity blocks. Optional. If unspecified and no asap3 object, then this is set to the number 
#'                of fleets + indices. If specified, ensure other components of \code{selectivity} are consistent.}
#'   }
#'
#' @export
set_selectivity = function(input, selectivity)
{
  data = input$data
  par = input$par
  map = input$map

  par_index = list(
    1:data$n_ages,
    data$n_ages + 1:2,
    data$n_ages + 3:6,
    data$n_ages + 1:2
  )
  asap3 = input$asap3
  input$log$selectivity <- list()
  if(is.null(asap3)) {
    if(is.null(selectivity$n_selblocks)) {
      data$n_selblocks = max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices)
      input$log$selectivity <- c(input$log$selectivity, paste0("number of selblocks, ", data$n_selblocks, 
        ", is being determined by input$data$selblock_pointer_fleets and input$data$selblock_pointer_indices. Those should be 1,...,max(pointers)."))
      data$n_selblocks = max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices)
      #data$n_selblocks = data$n_fleets + data$n_indices  #1 for fleet, 1 for index
    } else data$n_selblocks = selectivity$n_selblocks
  } else {
    data$n_selblocks = 0
    data$selblock_models = integer(0)
    for(i in 1:length(asap3)){ #have to do all fleets first, then indices
      data$n_selblocks <- data$n_selblocks + asap3[[i]]$n_fleet_sel_blocks
      data$selblock_models <- c(data$selblock_models, asap3[[i]]$sel_block_option)
    }
    for(i in 1:length(asap3)){
      which_indices <- which(asap3[[i]]$use_index ==1) # length = data$n_indices
      asap3[[i]]$index_sel_option <- asap3[[i]]$index_sel_option[which_indices]
      asap3[[i]]$index_sel_ini = asap3[[i]]$index_sel_ini[which_indices]
      data$n_selblocks <- data$n_selblocks + length(which_indices)
      data$selblock_models <- c(data$selblock_models, asap3[[i]]$index_sel_option)
    }
    orig_sel_models <- data$selblock_models
  }
  no_asap = is.null(asap3)
  selopts <- c("age-specific","logistic","double-logistic","decreasing-logistic")
  
  if(is.null(selectivity$model)) {
    if(no_asap) data$selblock_models = rep(2, data$n_selblocks)
  } 
  if(!is.null(selectivity$model)){
    if(length(selectivity$model) != data$n_selblocks) stop("Length of selectivity$model must equal number of selectivity blocks (e.g., asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$model %in% selopts)) stop("Each model entry must be one of the following: 'age-specific','logistic','double-logistic','decreasing-logistic'")
    data$selblock_models <- match(selectivity$model, selopts)
  }
  
  if(is.null(selectivity$re)) data$selblock_models_re <- rep(1, data$n_selblocks) # default: no RE on selectivity parameters
  if(!is.null(selectivity$re)){
    if(length(selectivity$re) != data$n_selblocks) stop("Length of selectivity$re must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$re %in% c("none","iid","ar1","ar1_y","2dar1"))) stop("Each selectivity$re entry must be one of the following: 'none','iid','ar1','ar1_y','2dar1'")
    data$selblock_models_re <- match(selectivity$re, c("none","iid","ar1","ar1_y","2dar1"))
  }
  
  selblock_pointers <- cbind(data$selblock_pointer_fleets, data$selblock_pointer_indices)
  data$selblock_years <- matrix(0, nrow=data$n_years_model, ncol=data$n_selblocks)
  for(b in 1:data$n_selblocks) data$selblock_years[,b] <- apply(selblock_pointers, 1, function(x) b %in% x)
  data$n_years_selblocks <- apply(data$selblock_years, 2, sum)
  
  data$n_selpars <- c(data$n_ages,2,4,2)[data$selblock_models] # num selpars per block
  # Prep selectivity initial values  
  selpars_ini = matrix(NA, data$n_selblocks, data$n_ages + 6)
  # Prep selectivity map
  phase_selpars = matrix(-1, data$n_selblocks, data$n_ages + 6)
  for(b in 1:data$n_selblocks){
    phase_selpars[b,par_index[[data$selblock_models[b]]]] = 1
  }

  # initial values
  if(is.null(selectivity$initial_pars)) {
    if(!no_asap) {
      j = 1
      for(k in 1:length(asap3)) for(i in 1:asap3[[k]]$n_fleet_sel_blocks) {
        selpars_ini[j,] = asap3[[k]]$sel_ini[[i]][,1]
        j = j + 1
      }
      for(k in 1:length(asap3)) for(i in 1:length(asap3[[k]]$index_sel_ini)){
        selpars_ini[j,] = asap3[[k]]$index_sel_ini[[i]][,1]
        j = j + 1
      }
    }
    default_selpars <- list()
    dpars = c(0.5,data$n_ages/2)
    orig_selpars <- list()

    for(b in 1:data$n_selblocks){
      if(data$selblock_models[b] == 1) 
      {
        default_selpars[[b]] <- rep(0.5, data$n_ages) # default to middle of par range
      }
      if(data$selblock_models[b] %in% c(2,4)) {
        default_selpars[[b]] <- rep(data$n_ages/2, 2) # default to middle of par range
      }
      if(data$selblock_models[b] == 3) {
        default_selpars[[b]] <- rep(data$n_ages/2, 4) # default to middle of par range
      }
      if(!no_asap){
        orig_selpars[[b]] <- selpars_ini[b,par_index[[data$selblock_models[b]]]]
      }
    }
    if(no_asap) for(b in 1:data$n_selblocks){
      selpars_ini[b,] <- c(rep(0.5,data$n_ages), rep(data$n_ages/2, 6))#default_selpars[[b]] # default to middle of par range
    }
    if(!no_asap) {
      #defined above
      #orig_sel_models <- c(asap3[[i]]$sel_block_option, asap3[[i]]$index_sel_option)
      sel_mod_diff_warn <- NULL
      #for(b in 1:data$n_selblocks){
      for(b in 1:length(orig_sel_models)){
        if(data$selblock_models[b] != orig_sel_models[b]){
          selpars_ini[b,par_index[[data$selblock_models[b]]]] <- default_selpars[[b]] # default to middle of par range
          sel_mod_diff_warn <- paste(sel_mod_diff_warn, paste0("For Selectivity Block ",b,":"), 
            #paste0("  ASAP .dat file: ",selopts[orig_sel_models[b]]),
            paste0("  selectivity$models: ",selopts[data$selblock_models[b]]),
            paste0("  Changing values from ASAP3 .dat file ",
            #paste0(orig_selpars[[b]], collapse = ", "),
            " to inits in middle of par range ",
            paste0(default_selpars[[b]], collapse = ", ")), sep="\n")
        }
      }
      if(!is.null(sel_mod_diff_warn)){
        sel_mod_diff_warn <- paste("Selectivity models differ from ASAP .dat file but initial","parameter values not specified. Please check initial values","and specify with selectivity$initial_pars if desired.",sel_mod_diff_warn,sep="\n")
        input$log$selectivity <- c(input$log$selectivity, sel_mod_diff_warn, sep="\n")
      }
    }
  } else {
    if(!is.list(selectivity$initial_pars)) stop("selectivity$initial_pars must be a list")
    if(length(selectivity$initial_pars) != data$n_selblocks) stop("Length of selectivity$initial_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    for(b in 1:data$n_selblocks){
      if(length(selectivity$initial_pars[[b]]) != data$n_selpars[b]) stop(paste0("Length of vector ",b," in the selectivity$initial_pars list is not equal to the number of selectivity parameters for block ",b,": ",data$n_selpars[b]))
      selpars_ini[b,par_index[[data$selblock_models[b]]]] = selectivity$initial_pars[[b]]
    }
  }
  
  # which selpars to fix, either $fix_pars or $map_pars
  if(is.null(selectivity$fix_pars) & is.null(selectivity$map_pars)){
    if(!no_asap){
      j = 1
      for(k in 1:length(asap3)) for(i in 1:asap3[[k]]$n_fleet_sel_blocks) {
        phase_selpars[j,par_index[[asap3[[k]]$sel_block_option[i]]]] = asap3[[k]]$sel_ini[[i]][par_index[[asap3[[k]]$sel_block_option[i]]],2]
        j <- j + 1
      }
      for(k in 1:length(asap3)) for(i in 1:length(asap3[[k]]$index_sel_ini)) {
        phase_selpars[j,par_index[[asap3[[k]]$index_sel_option[i]]]] = asap3[[k]]$index_sel_ini[[i]][par_index[[asap3[[k]]$index_sel_option[i]]],2]
        j <- j + 1
      }
    } 
  } else {
    if(!is.null(selectivity$fix_pars) & !is.null(selectivity$map_pars)) stop("Cannot specify $fix_pars and $map_pars (both set which pars to estimate). Choose one or the other.")
    # fix_pars
    if(!is.null(selectivity$fix_pars)){
      if(!is.list(selectivity$fix_pars)) stop("selectivity$fix_pars must be a list")
      if(length(selectivity$fix_pars) != data$n_selblocks) stop("Length of selectivity$fix_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices).
        Use 'NULL' to not fix any parameters for a block, e.g. list(NULL,4,2) does not fix any pars in block 1")
      for(b in 1:data$n_selblocks){
        if(data$selblock_models[b] == 1) phase_selpars[b,selectivity$fix_pars[[b]]] = -1
        if(data$selblock_models[b] %in% c(2,4)) phase_selpars[b,data$n_ages+selectivity$fix_pars[[b]]] = -1
        if(data$selblock_models[b] == 3) phase_selpars[b,data$n_ages+2+selectivity$fix_pars[[b]]] = -1
      }
    }
    # map_pars
    if(!is.null(selectivity$map_pars)){
      if(!is.list(selectivity$map_pars)) stop("selectivity$map_pars must be a list")
      if(length(selectivity$map_pars) != data$n_selblocks) stop("Length of selectivity$map_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
      for(b in 1:data$n_selblocks){
        naind <- which(is.na(selectivity$map_pars[[b]]))
        if(data$selblock_models[b] == 1) phase_selpars[b, naind] = -1
        if(data$selblock_models[b] %in% c(2,4)) phase_selpars[b, data$n_ages + naind] = -1
        if(data$selblock_models[b] == 3) phase_selpars[b,data$n_ages + 2 + naind] = -1
      }
    }
  }
  # For age-specific selectivity blocks, check for ages with ~zero catch and fix these at 0
  age_specific <- which(data$selblock_models==1)
  for(b in age_specific){
    if(all(phase_selpars[b,] < 0)){ # if no selpars estimated, keep fixed at specified initial values
     phase_selpars[b,] = -1
    } else {
      ind = list(fleets = which(apply(data$selblock_pointer_fleets == b,2,sum) > 0))
      ind$indices = which(apply(data$selblock_pointer_indices == b,2,sum) > 0)
      paa = matrix(nrow = 0, ncol = data$n_ages)
      if(length(ind$fleets)) for(f in ind$fleets) {
        y = data$catch_paa[f,which(data$selblock_pointer_fleets[,f] == b & data$use_catch_paa[,f] == 1),]
        paa = rbind(paa,y)
      }
      if(length(ind$indices)) for(i in ind$indices) {
        y = data$index_paa[i,which(data$selblock_pointer_indices[,i] == b & data$use_index_paa[,i] == 1),]
        paa = rbind(paa,y)
      }
      y = apply(paa,2,sum)
      ind <- which(y < 1e-5 & phase_selpars[b,1:data$n_ages] > 0) #if phase is set to -1 and 0s in all years for an age and selpars_ini >0, then keep selpars as is.
      selpars_ini[b, ind] = 0
      phase_selpars[b, ind] = -1
    }
  }
  temp <- matrix(NA, data$n_selblocks, data$n_ages + 6)
  # if(!is.null(selectivity$fix_pars)){ # use fix_pars
    temp[which(phase_selpars > 0)] = 1:sum(phase_selpars>0)
  # }
  if(!is.null(selectivity$map_pars)){ # use map_pars directly
    for(b in 1:data$n_selblocks) temp[b, par_index[[data$selblock_models[b]]]] = selectivity$map_pars[[b]]
  }
  for(b in 1:data$n_selblocks){ 
    if(data$selblock_models_re[b] == 3){
      bl <- temp[b, par_index[[data$selblock_models[b]]]]
      # print(b)
      # print(data$selblock_models_re[b])
      # print(data$selblock_models[b])
      # print(bl)
      if(sum(!is.na(bl)) < 3) {
        #data$selblock_models_re[b] <- 1 #no RE for this block
        stop(paste0("'ar1' (AR1(age)) selectivity random effects specified for block ",b," but no number of free mean parameters is <= 2, 
        which will not be identifiable. Use age-specific selectivity with more free mean selectivity parameters (par$logit_selpars) or use 
        a different random effects specification."))
      } else {
      # if re='ar1' (by age) with age-specific selectivity, only estimate one mean shared across ages
      # but don't overwrite fixed pars (likely will be fixing one age at 1)
        bl[!is.na(bl)] = min(bl, na.rm=TRUE)
        temp[b, par_index[[data$selblock_models[b]]]] = bl
      }
      # warning message if no mean sel pars (logit_selpars) are fixed
      # allow so user can fit model without fixing and then fix the age with highest sel at 1
      if(all(!is.na(bl))  & data$selblock_models[b] == 1) input$log$selectivity <- c(input$log$selectivity, paste0("NOTE: 'ar1' (AR1(age)) with age-specfici selectivity for block ",b,
        " but no age fixed at 1. Advised to fit the current model and then fix the age with highest selectivity at 1. 
        Can use selectivity$fix_pars."))
      # could add warning message if most ages are fixed, leaving less than xx ages to estimate with the AR1 re
    }
  }
  #data$selpars_est <- phase_selpars
  #data$selpars_est[data$selpars_est == -1] = 0
  data$selpars_est <- matrix(0, data$n_selblocks, data$n_ages + 6)
  data$selpars_est[which(!is.na(temp))] = 1
  # print(data$selpars_est)
  data$n_selpars_est <- apply(data$selpars_est > 0, 1, sum)
  map$logit_selpars = factor(temp)
  # print(temp)
  # print(data$n_selpars_est)


  # initial values on logit scale, par$logit_selpars
  selpars_lo = selpars_hi = matrix(0, data$n_selblocks, data$n_ages + 6)
  selpars_hi[,1:data$n_ages] = 1
  selpars_hi[,data$n_ages + 1:6] = data$n_ages
  data$selpars_lower = selpars_lo #only need these for estimated parameters
  data$selpars_upper = selpars_hi

  selpars_ini[which(selpars_ini > selpars_hi)] <- selpars_hi[which(selpars_ini > selpars_hi)]
  selpars_ini[which(selpars_ini < selpars_lo)] <- selpars_lo[which(selpars_ini < selpars_lo)]
  par$logit_selpars = log(selpars_ini-selpars_lo) - log(selpars_hi - selpars_ini)
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] = -10
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] = 10
  
  # random effects, selpars_re
  # number of estimated selpars per block * number of years per block (only if that block has re)
  if(any(data$selblock_models_re > 1)){
    par$selpars_re <- rep(0, sum((data$selblock_models_re > 1)*data$n_selpars_est*data$n_years_selblocks))
    tmp_vec <- c()
    ct <- 0
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] > 1){
        tmp <- matrix(0, nrow=data$n_years_selblocks[b], ncol=data$n_selpars_est[b])
        if(data$selblock_models_re[b] %in% c(2,5)){ # 2d ar1
          tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) + ct # all y,a estimated
        }
        if(data$selblock_models_re[b] == 3){ # ar1_a (devs by age, constant by year)
          for(i in 1:dim(tmp)[2]) tmp[,i] = (i + ct)
        }
        if(data$selblock_models_re[b] == 4){ # ar1_y (devs by year, constant by age)
          for(i in 1:dim(tmp)[1]) tmp[i,] = (i + ct)
        }
        ct = max(tmp)
        tmp_vec = c(tmp_vec, as.vector(tmp))
      }
    }
    map$selpars_re <- factor(tmp_vec)
  } else {
    par$selpars_re <- matrix(0)
    map$selpars_re <- factor(NA)
  }
  
  # initial and map for parameters controlling selectivity RE
  trans <- function(x) return((2/(1 + exp(-2*x))) - 1) # transform for correlation par
  # default initial values: sigma = 0.1, rho = 0
  par$sel_repars <- matrix(0, nrow=data$n_selblocks, ncol=3)
  par$sel_repars[,1] <- log(0.1)
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 3) par$sel_repars[b,3] <- 0 # if ar1 over ages only, fix rho_y = 0
    if(data$selblock_models_re[b] == 4) par$sel_repars[b,2] <- 0 # if ar1 over years only, fix rho = 0
    # check if only 1 estimated sel par (e.g. because all but 1 age is fixed), can't estimate rho
    if(data$n_selpars_est[b] < 2) par$sel_repars[b,2] <- 0
  }
  if(!is.null(selectivity$initial_sigma)){
    if(any(selectivity$initial_sigma < 0)) stop('Variance controlling selectivity random effects must be positive.') 
    par$sel_repars[,1] = log(selectivity$initial_sigma) # log scale
  }
  if(!is.null(selectivity$initial_cor)){
    if(any(sapply(selectivity$initial_cor, function(x) any(x < -1)))) stop('Correlation parameters controlling selectivity random effects must be between -1 and 1.')
    if(any(sapply(selectivity$initial_cor, function(x) any(x > 1)))) stop('Correlation parameters controlling selectivity random effects must be between -1 and 1.')
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) par$sel_repars[b,2] <- trans(selectivity$initial_cor[[b]]) # if ar1 over ages, use specified initial
      if(data$selblock_models_re[b] == 4) par$sel_repars[b,3] <- trans(selectivity$initial_cor[[b]]) # if ar1 over years, use specified initial
      if(data$selblock_models_re[b] == 5) par$sel_repars[b,2:3] <- trans(selectivity$initial_cor[[b]]) # if 2dar1 over years, use both
    }
  }  

  # map selectivity RE
  tmp.sel.repars <- par$sel_repars
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 1) tmp.sel.repars[b,] <- rep(NA,3) # no RE pars to estimate
    if(data$selblock_models_re[b] == 2) tmp.sel.repars[b,2:3] <- rep(NA,2) # estimate sigma
    if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,3] <- NA # estimate sigma, rho
    if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,2] <- NA # estimate sigma, rho_y
    if(data$n_selpars_est[b] < 2) tmp.sel.repars[b,2] <- NA # can't estimate rho if only 1 selpar estimated
  }
  if(!is.null(selectivity$map_sigma)){
    if(length(selectivity$map_sigma) != data$n_selblocks) stop("selectivity$map_sigma must be a vector of length = number of selectivity blocks")
    tmp.sel.repars[is.na(selectivity$map_sigma),1] = NA
  }
  if(!is.null(selectivity$map_cor)){
    if(length(selectivity$map_cor) != data$n_selblocks) stop("selectivity$map_cor must be a list of length = number of selectivity blocks")
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) tmp.sel.repars[is.na(selectivity$map_cor[[b]]),2] = NA
      if(data$selblock_models_re[b] == 4) tmp.sel.repars[is.na(selectivity$map_cor[[b]]),3] = NA
      if(data$selblock_models_re[b] == 5) tmp.sel.repars[is.na(selectivity$map_cor[[b]]),2:3] = NA
    }
  }
  ind.notNA <- which(!is.na(tmp.sel.repars))
  tmp.sel.repars[ind.notNA] <- 1:length(ind.notNA)
  if(!is.null(selectivity$map_sigma)){
    tmp.sel.repars[,1] = selectivity$map_sigma
  }
  if(!is.null(selectivity$map_cor)){
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,2] = selectivity$map_cor[[b]]
      if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,3] = selectivity$map_cor[[b]]
      if(data$selblock_models_re[b] == 5) tmp.sel.repars[b,2:3] = selectivity$map_cor[[b]]
    }
    st <- ifelse(!all(is.na(tmp.sel.repars[,1])), max(tmp.sel.repars[,1], na.rm=T), 0)
    tmp.sel.repars[,2:3] = tmp.sel.repars[,2:3] + st
  }
  map$sel_repars = factor(tmp.sel.repars)

  input$data = data
  input$par = par
  input$map = map
	if(length(input$log$selectivity))	input$log$selectivity <- c("Selectivity: \n", input$log$selectivity)
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  return(input)

}