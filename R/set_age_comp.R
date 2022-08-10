set_age_comp = function(input, age_comp)
{
	data = input$data
	par = input$par
  map = input$map
  all_models <- c( "multinomial",
                   "dir-mult",
                   "dirichlet-miss0",
                   "dirichlet-pool0",
                   "logistic-normal-miss0",
                   "logistic-normal-ar1-miss0",
                   "logistic-normal-pool0",
                   "logistic-normal-01-infl",
                   "logistic-normal-01-infl-2par",
                   "mvtweedie" )
  n_pars <- c(0,1,1,1,1,2,1,3,2,2)
  if(is.null(age_comp)){
    data$age_comp_model_fleets = rep(1, data$n_fleets) # multinomial by default
    data$age_comp_model_indices = rep(1, data$n_indices) # multinomial by default
  } else {
    if(is.character(age_comp) & length(age_comp)==1){ # all use the same model
      name_change = age_comp == "dirichlet"
      if(any(name_change)){
        cat("'dirichlet' is no longer an option and the old option is equivalent to 'dirichlet-pool0' so using that.\n")
        age_comp[which(name_change)] = "dirichlet-pool0"
      }
      themod <- match(age_comp, all_models)
      if(is.na(themod)) stop("age_comp option not recognized. See ?prepare_wham_input.")
      data$age_comp_model_fleets = rep(themod, data$n_fleets)
      data$age_comp_model_indices = rep(themod, data$n_indices)
    } else {
      if(all(names(age_comp) %in% c("fleets","indices"))){
        name_change = list(fleets = age_comp$fleets == "dirichlet", indices = age_comp$indices == "dirichlet")
        if(any(unlist(name_change))){
          cat("'dirichlet' is no longer an option and the old option is equivalent to 'dirichlet-pool0' so using that.\n")
          age_comp$fleets[which(name_change$fleets)] = "dirichlet-pool0"
          age_comp$indices[which(name_change$indices)] = "dirichlet-pool0"
        }
        themods <- match(age_comp$fleets, all_models)
        if(any(is.na(themods))) stop("age_comp$fleets option not recognized. See ?prepare_wham_input for available options.")
        if(length(themods) != data$n_fleets) stop("age_comp$fleets must have length = the number of fleets")
        data$age_comp_model_fleets = themods

        themods <- match(age_comp$indices, all_models)
        if(any(is.na(themods))) stop("age_comp$indices option not recognized. See ?prepare_wham_input for available options.")
        if(length(themods) != data$n_indices) stop("age_comp$indices must have length = the number of indices")
        data$age_comp_model_indices = themods
      } else {
        stop("age_comp must either be a character or a named ('fleets','indices') list. See ?prepare_wham_input.")
      }
    }
  }
  #data$n_age_comp_pars_fleets = c(0,1,1,3,1,2,1)[data$age_comp_model_fleets]
  #data$n_age_comp_pars_indices = c(0,1,1,3,1,2,1)[data$age_comp_model_indices]

  # age comp pars
  #n_catch_acomp_pars = c(0,1,1,3,1,2,1)[data$age_comp_model_fleets[which(apply(data$use_catch_paa,2,sum)>0)]]
  #n_index_acomp_pars = c(0,1,1,3,1,2,1)[data$age_comp_model_indices[which(apply(data$use_index_paa,2,sum)>0)]]
  par$catch_paa_pars = matrix(0,data$n_fleets, 3) #rep(0, sum(n_catch_acomp_pars))
  par$index_paa_pars = matrix(0,data$n_indices, 3) #rep(0, sum(n_catch_acomp_pars))
  #par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  #par$index_paa_pars = rep(0, sum(n_index_acomp_pars))
  neff <- data$catch_Neff
  neff[neff <= 0] <- 1
  catch_neff <- apply(neff,2,mean, na.rm=TRUE)#[which(apply(data$use_catch_paa,2,sum)>0)]
  ind = which(data$age_comp_model_fleets %in% 5:7)
  par$catch_paa_pars[ind,1] <- 0.5*log(catch_neff[ind])
  neff <- data$index_Neff
  neff[neff <= 0] <- 1
  index_neff <- apply(neff,2,mean, na.rm=TRUE)#[which(apply(data$use_index_paa,2,sum)>0)]
  ind = which(data$age_comp_model_indices %in% 5:7)
  par$index_paa_pars[ind,1] <- 0.5*log(index_neff[ind])

  map$index_paa_pars = matrix(NA,data$n_indices, 3)
  for(i in 1:data$n_indices) if(sum(data$use_index_paa[,i])){
    if(data$age_comp_model_indices[i] %in% c(2:5,7))  map$index_paa_pars[i,1] = 1
    if(data$age_comp_model_indices[i] %in% c(6,9,10))  map$index_paa_pars[i,1:2] = 1
    if(data$age_comp_model_indices[i] %in% 8)  map$index_paa_pars[i,1:3] = 1
  }
  map$catch_paa_pars = matrix(NA,data$n_fleets, 3)
  for(i in 1:data$n_fleets) if(sum(data$use_catch_paa[,i])){
    if(data$age_comp_model_fleets[i] %in% c(2:5,7))  map$catch_paa_pars[i,1] = 1
    if(data$age_comp_model_fleets[i] %in% c(6,9,10))  map$catch_paa_pars[i,1:2] = 1
    if(data$age_comp_model_fleets[i] %in% 8)  map$catch_paa_pars[i,1:3] = 1
  }
  nest = sum(map$index_paa_pars,na.rm=TRUE)
  if(nest) map$index_paa_pars[which(!is.na(map$index_paa_pars))] = 1:nest
  map$index_paa_pars = factor(map$index_paa_pars)
  nest = sum(map$catch_paa_pars,na.rm=TRUE)
  if(nest) map$catch_paa_pars[which(!is.na(map$catch_paa_pars))] = 1:nest
  map$catch_paa_pars = factor(map$catch_paa_pars)

#  if(all(data$age_comp_model_fleets %in% c(5,7))){ # start tau/neff at 0
#    neff <- data$catch_Neff
#    neff[neff <= 0] <- NA
#    par$catch_paa_pars = 0.5*log(neff) # exp(age_comp_pars(0)-0.5*log(Neff))
#  }  
#  if(all(data$age_comp_model_indices %in% c(5,7))){ # start tau/neff at 0
#    neff <- data$index_Neff
#    neff[neff <= 0] <- NA
#    neff <- apply(neff,2,mean, na.rm=TRUE)[which(apply(data$use_index_paa,2,sum)>0)]
#    par$index_paa_pars = 0.5*log(neff) # exp(age_comp_pars(0)-0.5*log(Neff))
#  }

	input$data = data
	input$par = par
  input$map = map
	return(input)
}
