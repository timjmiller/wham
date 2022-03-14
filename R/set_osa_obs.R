set_osa_obs = function(input)
{
  data = input$data

  # 5 components: fleet catch (log), index catch (log), Ecov, paa catch, paa index
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 1. log fleet catch
  x <- as.data.frame(data$agg_catch)
  x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
  colnames(x) <- paste0("fleet_", 1:data$n_fleets)
  obs_levels <- paste0("fleet_", 1:data$n_fleets)
  x$year <- 1:data$n_years_catch
  tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
  tmp <- tmp[complete.cases(tmp),]  
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_agg_catch==0, turned to NA, and removed
  tmp$age <- NA
  tmp$type <- "logcatch"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 2. log index catch
  x <- as.data.frame(data$agg_indices)
  x[data$use_indices==0] <- NA # only include index data to fit in obsvec
  colnames(x) <- paste0("index_", 1:data$n_indices)
  obs_levels <- c(obs_levels, paste0("index_", 1:data$n_indices))
  x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
  tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
  tmp <- tmp[complete.cases(tmp),]
  tmp$val <- log(tmp$val) # all obs of 0 catch should have use_indices==0, turned to NA, and already removed
  tmp$age <- NA
  tmp$type <- "logindex"
  obs <- rbind(obs, tmp[, obs.colnames])

  # 3. Ecov
  if(!all(data$Ecov_use_obs==0)){
    x <- as.data.frame(data$Ecov_obs)
    x[data$Ecov_use_obs==0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("Ecov_", 1:data$n_Ecov)
    obs_levels <- c(obs_levels, paste0("Ecov_", 1:data$n_Ecov))
    x$year <- seq(from=data$year1_Ecov-data$year1_model+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
    tmp <- tmp[complete.cases(tmp),]
    tmp$age <- NA
    tmp$type <- "Ecov"
    obs <- rbind(obs, tmp[, obs.colnames])
  }

  # # 4. paa catch
  dimnames(data$catch_paa) <- list(
    fleet=paste0("fleet_", 1:data$n_fleets),
    year=1:data$n_years_catch,
    age=1:data$n_ages)  
  for(i in 1:data$n_fleets){
    #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
    x <- data$catch_paa[i,,]
    x[which(data$use_catch_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
    x = as.data.frame(x)
    colnames(x) <- paste0(1:data$n_ages)
    x$year <- 1:data$n_years_catch # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="age")
    tmp$fleet = dimnames(data$catch_paa)[[1]][i]
    tmp$type <- "catchpaa"
    tmp <- tmp[complete.cases(tmp),]
    obs = rbind(obs, tmp[,obs.colnames])
  }

  # # 5. paa index
  dimnames(data$index_paa) <- list(
    fleet=paste0("index_", 1:data$n_indices),
    year=1:data$n_years_indices,
    age=1:data$n_ages)
  for(i in 1:data$n_indices){
    #obs_levels <- c(obs_levels, paste0("index_",i, "_paa"))
    x <- data$index_paa[i,,]
    x[which(data$use_index_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
    x = as.data.frame(x)
    colnames(x) <- paste0(1:data$n_ages)
    x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="age")
    tmp$fleet = dimnames(data$index_paa)[[1]][i]
    tmp$type <- "indexpaa"
    tmp <- tmp[complete.cases(tmp),]
    obs = rbind(obs, tmp[,obs.colnames])
  }

  # order by year, fleet, age, type
  obs$fleet <- factor(obs$fleet, levels=obs_levels)
  o <- order(obs$fleet, obs$type, as.numeric(obs$year), as.numeric(obs$age))
  obs <- obs[o,]

  # calculate obsvec indices in keep arrays
  obs$ind <- 1:dim(obs)[1]
  data$keep_C <- matrix(NA, nrow=data$n_years_catch, ncol=data$n_fleets)
  for(i in 1:data$n_fleets){
    ind = which(data$use_agg_catch[,i]==1) 
    data$keep_C[ind,i] <- subset(obs, type=='logcatch' & fleet == paste0("fleet_",i))$ind
  }
  #xl <- lapply(seq_len(nrow(data$use_agg_catch)), function(r) which(data$use_agg_catch[r,]==1))
  #Col <- unlist(xl)
  #Row <- rep(1:data$n_years_catch, times=sapply(xl, length))
  #data$keep_C[cbind(Row,Col)] <- subset(obs, type=='logcatch')$ind
  # subtract 1 bc TMB indexes from 0
  data$keep_C <- data$keep_C - 1

  data$keep_I <- matrix(NA, nrow=data$n_years_indices, ncol=data$n_indices)
  for(i in 1:data$n_indices){
    ind = which(data$use_indices[,i]==1) 
    data$keep_I[ind,i] <- subset(obs, type=='logindex' & fleet == paste0("index_",i))$ind
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_I <- data$keep_I - 1

  data$keep_E <- matrix(NA, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
  for(i in 1:data$n_Ecov){
    ind = which(data$Ecov_use_obs[,i]==1) 
    data$keep_E[ind,i] <- subset(obs, type=='Ecov' & fleet == paste0("Ecov_",i))$ind
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_E <- data$keep_E - 1

  data$condition_no_osa = NULL #to condition on age comps for likelihoods we don't have osa residuals set up for.
  data$subset_discrete_osa = NULL #age comp obs for likelihoods we need to specify as discrete obs.

  data$keep_Cpaa <- array(NA, dim=c(data$n_fleets, data$n_years_catch, data$n_ages))
  for(i in 1:data$n_fleets) {
    ind = which(data$use_catch_paa[,i]==1)
    data$keep_Cpaa[i,ind,] <- matrix(subset(obs, type=='catchpaa' & fleet==paste0("fleet_",i))$ind, 
      nrow=length(ind), ncol=data$n_ages, byrow=TRUE)
    # subtract 1 bc TMB indexes from 0
    data$keep_Cpaa[i,,] <- data$keep_Cpaa[i,,] - 1
    #subset for oneStepPredict can't include these
    if(data$age_comp_model_fleets[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, data$keep_Ipaa[i,,]+1)
    if(data$age_comp_model_fleets[i] %in% 8:9) data$condition_no_osa = c(data$condition_no_osa, data$keep_Ipaa[i,,]+1)
  }
  data$keep_Ipaa <- array(NA, dim=c(data$n_indices, data$n_years_indices, data$n_ages))
  for(i in 1:data$n_indices) {
    ind = which(data$use_index_paa[,i]==1)
    data$keep_Ipaa[i,ind,] <- matrix(subset(obs, type=='indexpaa' & fleet==paste0("index_",i))$ind, 
      nrow=length(ind), ncol=data$n_ages, byrow=TRUE)
    # subtract 1 bc TMB indexes from 0
    data$keep_Ipaa[i,,] <- data$keep_Ipaa[i,,] - 1
    #subset for oneStepPredict can't include these
    if(data$age_comp_model_indices[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, data$keep_Ipaa[i,,]+1)
    if(data$age_comp_model_indices[i] %in% 8:9) data$condition_no_osa = c(data$condition_no_osa, data$keep_Ipaa[i,,]+1)
  }
  #these are the models in order
  #all_models <- c("multinomial","dir-mult","dirichlet-miss0","dirichlet-pool0","logistic-normal-miss0","logistic-normal-ar1-miss0",
  #  "logistic-normal-pool0","logistic-normal-01-infl","logistic-normal-01-infl-2par")

  obs$cohort = as.numeric(obs$year) - as.numeric(obs$age)   #make cohort column. could be useful for analyzing age comp OSA residuals
  data$obs <- obs
  data$obsvec <- obs$val
  data$do_osa = 0 #this will be changed when TMB::oneStepPredict is caled by fit_wham
  data$do_post_samp = rep(0,5) #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated


  input$data = data
  return(input)
}