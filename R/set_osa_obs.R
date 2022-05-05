set_osa_obs = function(input)
{
  data = input$data

  obs.colnames <- c("year","fleet","bin","type","val") # change age -> bin
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 5 components: Ecov, fleet catch (log), index catch (log), paa catch, paa index
  # do all Ecov first
  # 1. Ecov
  if(!all(data$Ecov_use_obs==0)){
    x <- as.data.frame(data$Ecov_obs)
    x[data$Ecov_use_obs==0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("Ecov_", 1:data$n_Ecov)
    x$year <- seq(from=data$year1_Ecov-data$year1_model+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
    tmp <- tmp[complete.cases(tmp),]
    tmp$bin <- NA
    tmp$type <- "Ecov"
    obs <- rbind(obs, tmp[, obs.colnames])
  }


  #x <- data.frame(NA,nrow = data$n_years_model, ncol = data$n_indices)
  #x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s

  #step through time adding data
  for(y in 1:data$n_years_model){
    # 2. log index catch
    x <- as.data.frame(data$agg_indices)
    x[data$use_indices == 0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("index_", 1:data$n_indices)
    for(i in 1:data$n_indices) {
      if(!is.na(x[y,i])) {
        tmp = data.frame(year = y, fleet = colnames(x)[i], bin = NA, type = 'logindex', val = log(x[y,i]))
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }
    #tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
      #tmp <- tmp[complete.cases(tmp),]
      #tmp$val <- log(tmp$val) # all obs of 0 catch should have use_indices==0, turned to NA, and already removed
      #tmp$age <- NA
      #tmp$type <- "logindex"

    # 3. log fleet catch
    x <- as.data.frame(data$agg_catch)
    x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
    colnames(x) <- paste0("fleet_", 1:data$n_fleets)
    for(i in 1:data$n_fleets) {
      if(!is.na(x[y,i])) {
        tmp = data.frame(year = y, fleet = colnames(x)[i], bin = NA, type = 'logcatch', val = log(x[y,i]))
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }
    #x$year <- 1:data$n_years_catch
    #tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
    #tmp <- tmp[complete.cases(tmp),]  
    #tmp$val <- log(tmp$val) # all obs of 0 catch should have use_agg_catch==0, turned to NA, and removed
    #tmp$age <- NA
    #tmp$type <- "logcatch"
    #obs <- rbind(obs, tmp[, obs.colnames])
  
    # 4. paa index
    for(i in 1:data$n_indices){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$index_paa[i,,]
      x[which(data$use_index_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
      indices = paste0("index_", 1:data$n_indices)
      if(data$use_index_paa[y,i]!=0) 
      {
        tmp = data.frame(year = y, fleet = indices[i], bin = 1:data$n_ages, type = 'indexpaa', val = x[y,])
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }

    # 4.1 pal index
    for(i in 1:data$n_indices){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$index_pal[i,,]
      x[which(data$use_index_pal[,i]==0),] <- NA # only include catch data to fit in obsvec
      indices = paste0("index_", 1:data$n_indices)
      if(data$use_index_pal[y,i]!=0) 
      {
        tmp = data.frame(year = y, fleet = indices[i], bin = data$lengths, type = 'indexpal', val = x[y,])
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }


    # 5. paa catch
    #dimnames(data$catch_paa) <- list(
    #  fleet=paste0("fleet_", 1:data$n_fleets),
    #  year=1:data$n_years_catch,
    #  age=1:data$n_ages)  
    for(i in 1:data$n_fleets){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$catch_paa[i,,]
      x[which(data$use_catch_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
      #x = as.data.frame(x)
      fleets = paste0("fleet_", 1:data$n_fleets)
      if(data$use_catch_paa[y,i]!=0) {
        tmp = data.frame(year = y, fleet = fleets[i], bin = 1:data$n_ages, type = 'catchpaa', val = x[y,])
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }
    #colnames(x) <- paste0(1:data$n_ages)
    #x$year <- 1:data$n_years_catch # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
    #tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="age")
    #tmp$fleet = dimnames(data$catch_paa)[[1]][i]
    #tmp$type <- "catchpaa"
    #tmp <- tmp[complete.cases(tmp),]
    #obs = rbind(obs, tmp[,obs.colnames])

  # 5.1 pal catch
    for(i in 1:data$n_fleets){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$catch_pal[i,,]
      x[which(data$use_catch_pal[,i]==0),] <- NA # only include catch data to fit in obsvec
      #x = as.data.frame(x)
      fleets = paste0("fleet_", 1:data$n_fleets)
      if(data$use_catch_pal[y,i]!=0) {
        tmp = data.frame(year = y, fleet = fleets[i], bin = data$lengths, type = 'catchpal', val = x[y,])
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }
  # # 5. paa index
  #dimnames(data$index_paa) <- list(
  #  fleet=paste0("index_", 1:data$n_indices),
  #  year=1:data$n_years_indices,
  #  age=1:data$n_ages)
  #for(i in 1:data$n_indices){
  #  x <- data$index_paa[i,,]
  #  x[which(data$use_index_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
  #  x = as.data.frame(x)
  #  colnames(x) <- paste0(1:data$n_ages)
  #  x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
  #  tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="age")
  #  tmp$fleet = dimnames(data$index_paa)[[1]][i]
  #  tmp$type <- "indexpaa"
  #  tmp <- tmp[complete.cases(tmp),]
  #  obs = rbind(obs, tmp[,obs.colnames])
  }
  obs_levels <- paste0("fleet_", 1:data$n_fleets)
  obs_levels <- c(obs_levels, paste0("index_", 1:data$n_indices))
  obs_levels <- c(obs_levels, paste0("Ecov_", 1:data$n_Ecov))
  # order by year, fleet, age, type
  #obs$fleet <- factor(obs$fleet, levels=obs_levels)
  #o <- order(obs$fleet, obs$type, as.numeric(obs$year), as.numeric(obs$age))
  #obs <- obs[o,]

  # calculate obsvec indices in keep arrays
  obs$ind <- 1:dim(obs)[1]


  data$keep_E <- matrix(NA, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
  for(i in 1:data$n_Ecov){
    ind <- which(data$Ecov_use_obs[,i]==1) 
    data$keep_E[ind,i] <- subset(obs, type=='Ecov' & fleet == paste0("Ecov_",i))$ind
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_E <- data$keep_E - 1

  data$keep_C <- matrix(NA, nrow=data$n_years_catch, ncol=data$n_fleets)
  for(y in 1:data$n_years_model) for(i in 1:data$n_fleets){
    if(data$use_agg_catch[y,i]==1){ 
      data$keep_C[y,i] <- subset(obs, year == y & type=='logcatch' & fleet == paste0("fleet_",i))$ind
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_C <- data$keep_C - 1

  data$keep_I <- matrix(NA, nrow=data$n_years_indices, ncol=data$n_indices)
  for(y in 1:data$n_years_model) for(i in 1:data$n_indices){
    if(data$use_indices[y,i]==1){ 
      data$keep_I[y,i] <- subset(obs, year == y & type=='logindex' & fleet == paste0("index_",i))$ind
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_I <- data$keep_I - 1

  data$condition_no_osa = NULL #to condition on age comps for likelihoods we don't have osa residuals set up for.
  data$subset_discrete_osa = NULL #age comp obs for likelihoods we need to specify as discrete obs.

  # keep_Cpaa:
  data$keep_Cpaa <- array(NA, dim=c(data$n_fleets, data$n_years_model, data$n_ages))
  for(i in 1:data$n_fleets) {
    for(y in 1:data$n_years_model) if(data$use_catch_paa[y,i]==1){
      data$keep_Cpaa[i,y,] <- subset(obs, year == y & type=='catchpaa' & fleet==paste0("fleet_",i))$ind
      #subset for oneStepPredict can't include these
      if(data$age_comp_model_fleets[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, data$keep_Cpaa[i,y,])
      if(data$age_comp_model_fleets[i] %in% 8:9) data$condition_no_osa = c(data$condition_no_osa, data$keep_Cpaa[i,y,])
    }
  }

  # subtract 1 bc TMB indexes from 0
  data$keep_Cpaa <- data$keep_Cpaa - 1
  
  # keep_Cpal:
  data$keep_Cpal <- array(NA, dim=c(data$n_fleets, data$n_years_model, data$n_lengths))
  for(i in 1:data$n_fleets) {
    for(y in 1:data$n_years_model) if(data$use_catch_pal[y,i]==1){
      data$keep_Cpal[i,y,] <- subset(obs, year == y & type=='catchpal' & fleet==paste0("fleet_",i))$ind
      #subset for oneStepPredict can't include these
      if(data$len_comp_model_fleets[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, data$keep_Cpal[i,y,])
      if(data$len_comp_model_fleets[i] %in% 8:9) data$condition_no_osa = c(data$condition_no_osa, data$keep_Cpal[i,y,])
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Cpal <- data$keep_Cpal - 1

  # keep_Ipaa:
  data$keep_Ipaa <- array(NA, dim=c(data$n_indices, data$n_years_model, data$n_ages))
  for(i in 1:data$n_indices) {
    for(y in 1:data$n_years_model) if(data$use_index_paa[y,i]==1){
      data$keep_Ipaa[i,y,] <- subset(obs, year == y & type=='indexpaa' & fleet==paste0("index_",i))$ind 
      #subset for oneStepPredict can't include these
      if(data$age_comp_model_indices[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, data$keep_Ipaa[i,y,])
      if(data$age_comp_model_indices[i] %in% 8:9) data$condition_no_osa = c(data$condition_no_osa, data$keep_Ipaa[i,y,])
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Ipaa <- data$keep_Ipaa - 1
  
  # keep_Ipal:
  data$keep_Ipal <- array(NA, dim=c(data$n_indices, data$n_years_model, data$n_lengths))
  for(i in 1:data$n_indices) {
    for(y in 1:data$n_years_model) if(data$use_index_pal[y,i]==1){
      data$keep_Ipal[i,y,] <- subset(obs, year == y & type=='indexpal' & fleet==paste0("index_",i))$ind 
      #subset for oneStepPredict can't include these
      if(data$len_comp_model_indices[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, data$keep_Ipal[i,y,])
      if(data$len_comp_model_indices[i] %in% 8:9) data$condition_no_osa = c(data$condition_no_osa, data$keep_Ipal[i,y,])
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Ipal <- data$keep_Ipal - 1

  #these are the models in order
  #all_models <- c("multinomial","dir-mult","dirichlet-miss0","dirichlet-pool0","logistic-normal-miss0","logistic-normal-ar1-miss0",
  #  "logistic-normal-pool0","logistic-normal-01-infl","logistic-normal-01-infl-2par")

  obs$cohort = NA
  find_paa = grep(pattern = 'paa', x = obs$type) # only for paa
  obs$cohort[find_paa] = as.numeric(obs$year[find_paa]) - as.numeric(obs$bin[find_paa])   #make cohort column. could be useful for analyzing age comp OSA residuals
  data$obs <- obs
  data$obsvec <- obs$val
  data$do_osa = 0 #this will be changed when TMB::oneStepPredict is called by fit_wham
  data$do_post_samp = rep(0,5) #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated


  input$data = data
  return(input)
}