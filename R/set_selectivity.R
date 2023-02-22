set_selectivity = function(input, selectivity)
{
  data = input$data
  par = input$par
  map = input$map

  par_index = list(
    1:data$n_ages,
    data$n_ages + 1:2,
    data$n_ages + 3:6,
    data$n_ages + 1:2,
	  data$n_ages + 7:12,
    data$n_ages + 13:14,
    data$n_ages + 13:14,
    data$n_ages + 15:20
  )

  if(is.null(input$asap3)) {
    asap3 = NULL
    if(is.null(selectivity$n_selblocks)) {
      data$n_selblocks = max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices)
      print(paste0("number of selblocks, ", data$n_selblocks, 
        ", is being determined by input$data$selblock_pointer_fleets and input$data$selblock_pointer_indices. Those should be 1,...,max(pointers)."))
      data$n_selblocks = max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices)
      #data$n_selblocks = data$n_fleets + data$n_indices  #1 for fleet, 1 for index
    }
    else data$n_selblocks = selectivity$n_selblocks
  }
  else {
    asap3 = input$asap3
    which_indices <- which(asap3$use_index ==1) # length = data$n_indices
    asap3$index_sel_option <- asap3$index_sel_option[which_indices]
    asap3$index_sel_ini = asap3$index_sel_ini[which_indices]
    data$n_selblocks <- asap3$n_fleet_sel_blocks + data$n_indices
    data$selblock_models <- c(asap3$sel_block_option, asap3$index_sel_option)  
  }
  no_asap = is.null(asap3)
  selopts <- c("age-specific","logistic","double-logistic","decreasing-logistic","double-normal","len-logistic","len-decreasing-logistic","len-double-normal")
  # if(!no_asap) data$n_selblocks <- asap3$n_fleet_sel_blocks + asap3$n_indices
  # if(no_asap) data$n_selblocks = data$n_fleets + data$n_indices
  
  if(is.null(selectivity$model)) {
    #if(!no_asap) data$selblock_models <- c(asap3$sel_block_option, asap3$index_sel_option)
    #else 
    if(no_asap) data$selblock_models = rep(2, data$n_selblocks)
  } 
  if(!is.null(selectivity$model)){
    if(length(selectivity$model) != data$n_selblocks) stop("Length of selectivity$model must equal number of selectivity blocks (e.g., asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$model %in% selopts)) stop("Each model entry must be one of the following: 'age-specific','logistic','double-logistic','decreasing-logistic','double-normal','len-logistic','len-decreasing-logistic','len-double-normal'")
    data$selblock_models <- match(selectivity$model, selopts)
  }
  
  if(is.null(selectivity$re)) data$selblock_models_re <- rep(1, data$n_selblocks) # default: no RE on selectivity parameters
  if(!is.null(selectivity$re)){
    if(length(selectivity$re) != data$n_selblocks) stop("Length of selectivity$re must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$re %in% c("none","iid","ar1","ar1_y","2dar1"))) stop("Each selectivity$re entry must be one of the following: 'none','iid','ar1','ar1_y','2dar1'")
    data$selblock_models_re <- match(selectivity$re, c("none","iid","ar1","ar1_y","2dar1"))
  }

  # if(!no_asap)
  # {
  #   data$selblock_pointer_fleets = cbind(sapply(asap3$sel_block_assign, function(x) return(x)))
  #   data$selblock_pointer_indices = matrix(rep(asap3$n_fleet_sel_blocks + 1:data$n_indices, each = data$n_years_model), data$n_years_model, data$n_indices)
  # }
  # else{
  #   data$selblock_pointer_fleets = matrix(rep(1:data$n_fleets, each = data$n_years_model), data$n_years_model, data$n_fleets)
  #   data$selblock_pointer_indices = matrix(rep(1:data$n_indices, each = data$n_years_model), data$n_years_model, data$n_indices) + data$n_fleets
  # }
  
  selblock_pointers <- cbind(data$selblock_pointer_fleets, data$selblock_pointer_indices)
  data$selblock_years <- matrix(0, nrow=data$n_years_model, ncol=data$n_selblocks)
  for(b in 1:data$n_selblocks) data$selblock_years[,b] <- apply(selblock_pointers, 1, function(x) b %in% x)
  data$n_years_selblocks <- apply(data$selblock_years, 2, sum)
  
  data$n_selpars <- c(data$n_ages,2,4,2,6,2,2,6)[data$selblock_models] # num selpars per block
  # Prep selectivity initial values  
  selpars_ini = matrix(NA, data$n_selblocks, data$n_ages + 20)
  # Prep selectivity map
  phase_selpars = matrix(-1, data$n_selblocks, data$n_ages + 20)
  for(b in 1:data$n_selblocks){
    phase_selpars[b,par_index[[data$selblock_models[b]]]] = 1
  }
  if(is.null(selectivity$initial_pars)) {
    if(!no_asap) {
      for(i in 1:asap3$n_fleet_sel_blocks) selpars_ini[i,] = asap3$sel_ini[[i]][,1]
      for(i in 1:data$n_indices) selpars_ini[i+asap3$n_fleet_sel_blocks,] = asap3$index_sel_ini[[i]][,1]
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
      if(data$selblock_models[b] == 5) {
        default_selpars[[b]] <- c(data$n_ages/2, -2, 0.5, 0.5, -6, -6) # default to middle of par range
      }
      if(data$selblock_models[b] %in% c(6,7)) {
        default_selpars[[b]] <- rep(data$n_lengths/2, 2) # default to middle of par range
      }
      if(data$selblock_models[b] == 8) {
        default_selpars[[b]] <- c(data$n_lengths/2, -2, 0.5, 0.5, -6, -6) # default to middle of par range
      }
      if(!no_asap){
        orig_selpars[[b]] <- selpars_ini[b,par_index[[data$selblock_models[b]]]]
      }
    }
    if(no_asap) for(b in 1:data$n_selblocks){
      selpars_ini[b,] <- c(rep(0.5,data$n_ages), rep(data$n_ages/2, 6), c(data$n_ages/2, -2, 0.5, 0.5, -6, -6), 
                          rep(data$n_lengths/2, 2), 
                          c(data$n_lengths/2, -2, 0.5, 0.5, -6, -6))#default_selpars[[b]] # default to middle of par range
    }
    if(!no_asap) {
      orig_sel_models <- c(asap3$sel_block_option, asap3$index_sel_option)
      sel_mod_diff_warn <- NULL
      for(b in 1:data$n_selblocks){
        if(data$selblock_models[b] != orig_sel_models[b]){
          selpars_ini[b,par_index[[data$selblock_models[b]]]] <- default_selpars[[b]] # default to middle of par range
          sel_mod_diff_warn <- paste(sel_mod_diff_warn, paste0("Block ",b,":"), paste0("  ASAP .dat file: ",selopts[orig_sel_models[b]]),paste0("  selectivity$models: ",selopts[data$selblock_models[b]]),paste0("  Changing .dat file inits ",paste0(orig_selpars[[b]], collapse = ", ")," to inits in middle of par range ",paste0(default_selpars[[b]], collapse = ", ")), sep="\n")
        }
      }
      if(!is.null(sel_mod_diff_warn)){
        sel_mod_diff_warn <- paste("Selectivity models differ from ASAP .dat file but initial","parameter values not specified. Please check initial values","and specify with selectivity$initial_pars if desired.",sel_mod_diff_warn,sep="\n")
        cat(sel_mod_diff_warn, sep="\n")
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
  
  if(is.null(selectivity$fix_pars)){
    if(!no_asap){
      for(i in 1:asap3$n_fleet_sel_blocks) phase_selpars[i,par_index[[asap3$sel_block_option[i]]]] = asap3$sel_ini[[i]][par_index[[asap3$sel_block_option[i]]],2]
      for(i in (1:data$n_indices)) phase_selpars[i+asap3$n_fleet_sel_blocks,par_index[[asap3$index_sel_option[i]]]] = asap3$index_sel_ini[[i]][par_index[[asap3$index_sel_option[i]]],2]
    } 
  } else {
    if(!is.list(selectivity$fix_pars)) stop("selectivity$fix_pars must be a list")
    if(length(selectivity$fix_pars) != data$n_selblocks) stop("Length of selectivity$fix_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices).
      Use 'NULL' to not fix any parameters for a block, e.g. list(NULL,4,2) does not fix any pars in block 1")
    for(b in 1:data$n_selblocks){
      if(data$selblock_models[b] == 1) phase_selpars[b,selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] %in% c(2,4)) phase_selpars[b,data$n_ages+selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] == 3) phase_selpars[b,data$n_ages+2+selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] == 5) phase_selpars[b,data$n_ages+6+selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] %in% c(6,7)) phase_selpars[b,data$n_ages+12+selectivity$fix_pars[[b]]] = -1
      if(data$selblock_models[b] == 8) phase_selpars[b,data$n_ages+14+selectivity$fix_pars[[b]]] = -1
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
      if(length(ind$fleets)) for(f in ind$fleets)
      {
        y = data$catch_paa[f,which(data$selblock_pointer_fleets[,f] == b & data$use_catch_paa[,f] == 1),]
        paa = rbind(paa,y)
      }
      if(length(ind$indices)) for(i in ind$indices)
      {
        y = data$index_paa[i,which(data$selblock_pointer_indices[,i] == b & data$use_index_paa[,i] == 1),]
        paa = rbind(paa,y)
      }
      y = apply(paa,2,sum)
      selpars_ini[b, which(y < 1e-5)] = 0
      phase_selpars[b, which(y < 1e-5)] = -1
    }
  }
  data$selpars_est <- phase_selpars
  data$selpars_est[data$selpars_est == -1] = 0
  data$n_selpars_est <- apply(data$selpars_est > 0, 1, sum)
  selpars_lo = selpars_hi = matrix(0, data$n_selblocks, data$n_ages + 20)
  selpars_lo[,data$n_ages + 7] = 1 # par1 age double-normal
  selpars_lo[,data$n_ages + 8:12] = -20 
  selpars_lo[,data$n_ages + 13:15] = min(data$lengths) 
  selpars_lo[,data$n_ages + 16:20] = -20 
  selpars_hi[,1:data$n_ages] = 1
  selpars_hi[,data$n_ages + 1:7] = data$n_ages
  selpars_hi[,data$n_ages + 8:12] = 11 
  selpars_hi[,data$n_ages + 13:15] = max(data$lengths) 
  selpars_hi[,data$n_ages + 16:20] = 11 

  temp = matrix(NA, data$n_selblocks, data$n_ages + 20)
  temp[which(phase_selpars>0)] = 1:sum(phase_selpars>0)
  map$logit_selpars = factor(temp)
  
  data$selpars_lower = selpars_lo #only need these for estimated parameters
  data$selpars_upper = selpars_hi
  #map = list(logit_selpars = factor(temp))

  # selectivity pars
  selpars_ini[which(selpars_ini > selpars_hi)] <- selpars_hi[which(selpars_ini > selpars_hi)]
  selpars_ini[which(selpars_ini < selpars_lo)] <- selpars_lo[which(selpars_ini < selpars_lo)]
  par$logit_selpars = log(selpars_ini-selpars_lo) - log(selpars_hi - selpars_ini)
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] = -10
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] = 10
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
  par$sel_repars <- matrix(0, nrow=data$n_selblocks, ncol=3)
  par$sel_repars[,1] <- log(0.1) # start sigma at 0.1, rho at 0
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 3) par$sel_repars[b,3] <- 0 # if ar1 over ages only, fix rho_y = 0
    if(data$selblock_models_re[b] == 4) par$sel_repars[b,2] <- 0 # if ar1 over years only, fix rho = 0
    # check if only 1 estimated sel par (e.g. because all but 1 age is fixed), can't estimate rho
    if(data$n_selpars_est[b] < 2) par$sel_repars[b,2] <- 0
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
  ind.notNA <- which(!is.na(tmp.sel.repars))
  tmp.sel.repars[ind.notNA] <- 1:length(ind.notNA)
  map$sel_repars = factor(tmp.sel.repars)

  input$data = data
  input$par = par
  input$map = map
  return(input)

}