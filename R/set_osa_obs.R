#' Set up observation vector that is used by the model for likelihood calculations and one-step-ahead residuals.
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#'
#' @return the same input list as provided, but with $obs and $obsvec configured.
#' This is run after any changes have been made to the data
#' @export
set_osa_obs = function(input)
{
  data = input$data
  par = input$par
  map = input$map
  input$log$osa_obs <- list()
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames

  # 5 components: Ecov, fleet catch (log), index catch (log), paa catch, paa index
  # do all Ecov first
  # 1. Ecov
  if(!all(data$Ecov_use_obs==0)){
    x <- as.data.frame(data$Ecov_obs)
    x[data$Ecov_use_obs==0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("Ecov_", 1:data$n_Ecov)
    x$year <- seq(from=input$years_Ecov[1]-input$years[1]+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
    tmp <- tmp[complete.cases(tmp),]
    tmp$age <- NA
    tmp$type <- "Ecov"
    obs <- rbind(obs, tmp[, obs.colnames])
  }


  #x <- data.frame(NA,nrow = data$n_years_model, ncol = data$n_indices)
  #x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s

  #find out if there are any age-specific selectivity parameters set to 0 and not estimated
  #if the selblock is used, then we will reduce the observation vector to omit these ages. 
  #an OSA residual for ages with selectivity assumed to be zero will be useless and causes estimation problems with conditional construction of 
  #age comp likelihoods.
  ind <- which(data$selblock_models == 1) #age-specific
  ages_omit <- lapply(1:data$n_selblocks, function(x) integer(0))
  if(length(ind)){ #some
    if(is.null(map$logit_selpars)) {
      input$log$osa_obs <- c(input$log$osa_obs, "NOTE: input$map$logit_selapars has not been specified. This will only work for simulation.")
    } else {
      lsp_map <- matrix(as.integer(map$logit_selpars),nrow = data$n_selblocks, ncol = data$n_ages+6)
      for(i in ind){
        age_pars <- par$logit_selpars[i,1:data$n_ages]  
        age_map <- lsp_map[i,1:data$n_ages]
        ages_omit[[i]] <- which(age_pars == -Inf & is.na(age_map)) #if length then observations for these ages will be omitted.
        if(length(ages_omit[[i]])){
          if(i %in% data$selblock_pointer_indices) for(j in 1:data$n_indices){ #some indices affected
            bad_year = which(data$selblock_pointer_indices[,j] == i & data$use_index_paa[,j] == 1)
            if(length(bad_year)){
              input$log$osa_obs <- c(input$log$osa_obs, paste0(
                "Selectivity block ", i, " has parameters for ages ", paste0(input$ages.lab[ages_omit[[i]]], collapse = ", "), 
                " fixed at 0 and is used by index ", j, "\n",
                " in years ", paste0(input$years[bad_year], collapse = ", "), ".\n", 
                "Observations at these ages will be omitted and remaining observations will be rescaled.\n\n"))
            }
          }
          if(i %in% data$selblock_pointer_fleets) for(j in 1:data$n_fleets){ #some fleets affected
            bad_year = which(data$selblock_pointer_fleets[,j] == i & data$use_catch_paa[,j] == 1)
            if(length(bad_year)){
              input$log$osa_obs <- c(input$log$osa_obs, paste0(
                "Selectivity block ", i, " has parameters for ages ", paste0(input$ages.lab[ages_omit[[i]]], collapse = ", "), 
                " fixed at 0 and is used by fleet ", j, "\n",
                " in years ", paste0(input$years[bad_year], collapse = ", "), ".\n", 
                "Observations at these ages will be omitted and remaining observations will be rescaled.\n\n"))
            }
          }
        }
      }
    }
  }

  #step through time adding data
  for(y in 1:data$n_years_model){
    
    # 2. log index catch
    x <- as.data.frame(data$agg_indices)
    x[data$use_indices == 0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("index_", 1:data$n_indices)
    for(i in 1:data$n_indices) {
      if(!is.na(x[y,i])) {
        tmp = data.frame(year = y, fleet = colnames(x)[i], age = NA, type = 'logindex', val = log(x[y,i]))
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }

    # 3. log fleet catch
    x <- as.data.frame(data$agg_catch)
    x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
    colnames(x) <- paste0("fleet_", 1:data$n_fleets)
    for(i in 1:data$n_fleets) {
      if(!is.na(x[y,i])) {
        tmp = data.frame(year = y, fleet = colnames(x)[i], age = NA, type = 'logcatch', val = log(x[y,i]))
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }

    # 4. paa index
    for(i in 1:data$n_indices){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$index_paa[i,,]
      x[which(data$use_index_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
      indices = paste0("index_", 1:data$n_indices)
      if(data$use_index_paa[y,i] == 1)
      {
        obs_y = x[y,]
        tmp <- ages_omit[[data$selblock_pointer_indices[y,i]]]
        res = transform_paa_obs(obs_y, data$age_comp_model_indices[i], ages_omit = tmp)
        obs_y = res[[1]]
        ind = res[[2]] #now the ages to use is specified for all likelihods by transform_paa_obs
        #multinom, D-M, mvtweedie
        if(data$age_comp_model_indices[i] %in% c(1:2,10,11)) obs_y = obs_y * data$index_Neff[y,i]

        #if(data$age_comp_model_indices[i] %in% 3:7) {
        #  ind = res[[2]]
        #} else ind = 1:data$n_ages
        if(length(ind)) {
          tmp = data.frame(year = y, fleet = indices[i], age = (1:data$n_ages)[ind], type = 'indexpaa', val = obs_y[ind])
          #tmp = data.frame(year = y, fleet = indices[i], age = (1:data$n_ages), type = 'indexpaa', val = obs_y)
          obs <- rbind(obs, tmp[, obs.colnames])
        } else {
          data$use_index_paa[y,i] <- 0 #set to not use because there are not enough positive values 
          input$log$osa_obs <- c(input$log$osa_obs, paste0("Setting data$use_index_paa to 0 for index ", i, " and year ", y, "because not enough positive values. \n"))
        }
      }
    }


    # 5. paa catch
    for(i in 1:data$n_fleets){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$catch_paa[i,,]
      x[which(data$use_catch_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
      #x = as.data.frame(x)
      fleets = paste0("fleet_", 1:data$n_fleets)
      if(data$use_catch_paa[y,i] == 1) {
        obs_y = x[y,]
        tmp <- ages_omit[[data$selblock_pointer_fleets[y,i]]]
        #multinom, D-M, mvtweedie
        res = transform_paa_obs(obs_y, data$age_comp_model_fleets[i], ages_omit = tmp)
        obs_y = res[[1]]
        ind = res[[2]] #now the ages to use is specified for all likelihods by transform_paa_obs
        if(data$age_comp_model_fleets[i] %in% c(1:2,10,11)) obs_y = obs_y * data$catch_Neff[y,i]
        #if(data$age_comp_model_fleets[i] %in% 3:7) {
        #  ind = res[[2]]
        #} else ind = 1:data$n_ages

        if(length(ind)) {
          tmp = data.frame(year = y, fleet = fleets[i], age = (1:data$n_ages)[ind], type = 'catchpaa', val = obs_y[ind])
          obs <- rbind(obs, tmp[, obs.colnames])
        } else {
          data$use_catch_paa[y,i] <- 0 #set to not use because there are not enough positive values
          input$log$osa_obs <- c(input$log$osa_obs, paste0("Setting data$use_catch_paa to 0 for fleet ", i, " and year ", y, "because not enough positive values. \n"))
        }
      }
    }

  }

  # calculate obsvec indices in keep arrays
  obs$ind <- 1:dim(obs)[1]


  data$keep_E <- matrix(NA, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
  for(i in 1:data$n_Ecov){
    ind <- which(data$Ecov_use_obs[,i]==1) 
    data$keep_E[ind,i] <- subset(obs, type=='Ecov' & fleet == paste0("Ecov_",i))$ind
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_E <- data$keep_E - 1

  data$keep_C <- matrix(NA, nrow=data$n_years_model, ncol=data$n_fleets)
  for(y in 1:data$n_years_model) for(i in 1:data$n_fleets){
    if(data$use_agg_catch[y,i]==1){ 
      data$keep_C[y,i] <- subset(obs, year == y & type=='logcatch' & fleet == paste0("fleet_",i))$ind
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_C <- data$keep_C - 1

  data$keep_I <- matrix(NA, nrow=data$n_years_model, ncol=data$n_indices)
  for(y in 1:data$n_years_model) for(i in 1:data$n_indices){
    if(data$use_indices[y,i]==1){ 
      data$keep_I[y,i] <- subset(obs, year == y & type=='logindex' & fleet == paste0("index_",i))$ind
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_I <- data$keep_I - 1

  data$condition_no_osa = NULL #to condition on age comps for likelihoods we don't have osa residuals set up for.
  data$subset_discrete_osa = NULL #age comp obs for likelihoods we need to specify as discrete obs.

  data$keep_Cpaa <- array(NA, dim=c(data$n_fleets, data$n_years_model, 2))
  for(i in 1:data$n_fleets) {
    for(y in 1:data$n_years_model) if(data$use_catch_paa[y,i]==1){
      tmp = subset(obs, year == y & type=='catchpaa' & fleet==paste0("fleet_",i))
      if(length(tmp$ind)) #should always be TRUE because use_paa changed above
      {
        data$keep_Cpaa[i,y,1:2] <- c(tmp$ind[1], length(tmp$ind))
        #if(data$age_comp_model_fleets[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, tmp$ind)
        #subset for oneStepPredict can't include these
        if(data$age_comp_model_fleets[i] %in% 8:10) data$condition_no_osa = c(data$condition_no_osa, tmp$ind)
      }
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Cpaa[,,1] <- data$keep_Cpaa[,,1] - 1
  
  data$keep_Ipaa <- array(NA, dim=c(data$n_indices, data$n_years_model, 2))
  for(i in 1:data$n_indices) {
    for(y in 1:data$n_years_model) if(data$use_index_paa[y,i]==1){
      tmp = subset(obs, year == y & type=='indexpaa' & fleet==paste0("index_",i))
      if(length(tmp$ind)) #should always be TRUE because use_paa changed above
      {
        data$keep_Ipaa[i,y,1:2] <- c(tmp$ind[1], length(tmp$ind))
        #if(data$age_comp_model_indices[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, tmp$ind)
        #subset for oneStepPredict can't include these
        if(data$age_comp_model_indices[i] %in% 8:10) data$condition_no_osa = c(data$condition_no_osa, tmp$ind)
      }
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Ipaa[,,1] <- data$keep_Ipaa[,,1] - 1
  
  obs$cohort = as.numeric(obs$year) - as.numeric(obs$age)   #make cohort column. could be useful for analyzing age comp OSA residuals
  data$obs <- obs
  data$obsvec <- obs$val
  data$agesvec <- obs$age #potentially needed for AR1 sigma correlation of logistic-normal paa obs. 
  data$do_osa = 0 #this will be changed when TMB::oneStepPredict is called by fit_wham
  attr(data, "check.passed") <- NULL #if data have been generated from obj$simulate(complete=T), this can be problematic
  input$data = data

  #data$do_post_samp = rep(0,5) #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
  if(length(input$log$osa_obs)) input$log$osa_obs <- c("OSA obs: \n", input$log$osa_obs)
  if(!is_internal_call()) cat(unlist(input$log$osa_obs, recursive=T))
  return(input)
}

transform_paa_obs = function(x, model, zero.criteria = 1e-15, do_mult = FALSE, ages_omit = integer(0)){
    #transforms paa obs for obsvec and appropriate for OSA residuals once model is fit.
    #remove zeros for dirichlet and logistic-normal
    #remove ages/classes where predicted probability = 0 for all likelihoods (i.e., surveys that only observe a subset of ages).
    #transform logistic-normal obs to MVN obs (analogous to log-catch and log-indices)
    all_models <- c("multinomial","dir-mult","dirichlet-miss0","dirichlet-pool0",
    "logistic-normal-miss0", "logistic-normal-ar1-miss0", "logistic-normal-pool0",
    "logistic-normal-01-infl","logistic-normal-01-infl-2par", "mvtweedie", "dir-mult-linear")
  # if model %in% 1:2 do nothing for multinomial and D-m
  is_pred_pos = !(1:length(x) %in% ages_omit)
  x[which(!is_pred_pos)] = 0 #if obs in omitted ages are zero this does nothing
  x =  x/sum(x) #rescale to a reduced set of ages, if obs in omitted ages are zero this does nothing
  if(model>2 & model<10){ #not multinom, D-M, mvtweedie or DM-linear
    is_pos = x> zero.criteria # looking for zeros here will also omit the ages with predicted probability = 0
    pos_ind = which(is_pos)
    npos = sum(is_pos)
    zero_ind = which(!is_pos)
    if(npos>1){ #need at least 2 categories
      if(length(zero_ind)){
        x_pos = x[pos_ind] # 0s missing or pooled
        # if(model %in% c(3,5:6)) x_pos = x[pos_ind] # 0s missing
        # if(model %in% c(4,7)) { # 0s pooled
        #   index = 0
        #   x_pos = rep(NA,npos)
        #   for(a in 1:length(x))
        #   {
        #     x_pos[index] = x_pos[index] + x[a]; #just adding zeros. Not necessary? Only necessary for expected proportions
        #     if((x[a] > zero.criteria) & index < npos) index = index+ 1
        #   }
        # }
        x[zero_ind] <- NA
      } else x_pos = x

      #transform to multivariate normal obs for logistic-normal.
      #Necessary because TMB::oneStepPredict needs observation transformations on R side.
      #Conditional distribution format for logistic-normal would be needed to avoid this. 
      if(model %in% 5:7){
        y = log(x_pos[-length(x_pos)])
        if(do_mult){ #multiplicative
          for(i in 1:length(y)) {
            y[i] = y[i] - log(1-sum(x_pos[1:i]))
          }
        } else { #additive
          y = y - log(x_pos[length(x_pos)])
        }
        y = c(y,NA) #NA for last category which = 1 - sum of the other categories
        x_pos = y
      }
      x[pos_ind] = x_pos
    } else { #only 1 positive category
      x[] = NA
    }
  } else { #multinom, D-M, mvtweedie
    x[which(!is_pred_pos)] = NA
    pos_ind = which(is_pred_pos)
  }
  return(list(x, pos_ind))
}
