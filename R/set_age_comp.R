set_age_comp = function(input, age_comp)
{
	data = input$data
	par = input$par
  if(is.null(age_comp)){
    data$age_comp_model_fleets = rep(1, data$n_fleets) # multinomial by default
    data$age_comp_model_indices = rep(1, data$n_indices) # multinomial by default
  } else {
    if(is.character(age_comp)){ # all use the same model
      themod <- match(age_comp, c("multinomial","dir-mult","dirichlet","logistic-normal-01-infl","logistic-normal-pool0","logistic-normal-01-infl-2par","logistic-normal-miss0"))
      if(is.na(themod)) stop("age_comp option not recognized. See ?prepare_wham_input.")
      data$age_comp_model_fleets = rep(themod, data$n_fleets)
      data$age_comp_model_indices = rep(themod, data$n_indices)
    } else {
      if(all(names(age_comp) == c("fleets","indices"))){
        themods <- match(age_comp$fleets, c("multinomial","dir-mult","dirichlet","logistic-normal-01-infl","logistic-normal-pool0","logistic-normal-01-infl-2par","logistic-normal-miss0"))
        if(any(is.na(themods))) stop("age_comp$fleets option not recognized. See ?prepare_wham_input for available options.")
        if(length(themods) != data$n_fleets) stop("age_comp$fleets must have length = the number of fleets")
        data$age_comp_model_fleets = themods

        themods <- match(age_comp$indices, c("multinomial","dir-mult","dirichlet","logistic-normal-01-infl","logistic-normal-pool0","logistic-normal-01-infl-2par","logistic-normal-miss0"))
        if(any(is.na(themods))) stop("age_comp$indices option not recognized. See ?prepare_wham_input for available options.")
        if(length(themods) != data$n_indices) stop("age_comp$indices must have length = the number of indices")
        data$age_comp_model_indices = themods
      } else {
        stop("age_comp must either be a character or a named list. See ?prepare_wham_input.")
      }
    }
  }
  data$n_age_comp_pars_fleets = c(0,1,1,3,1,2,1)[data$age_comp_model_fleets]
  data$n_age_comp_pars_indices = c(0,1,1,3,1,2,1)[data$age_comp_model_indices]

  # age comp pars
  n_catch_acomp_pars = c(0,1,1,3,1,2,1)[data$age_comp_model_fleets[which(apply(data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2,1)[data$age_comp_model_indices[which(apply(data$use_index_paa,2,sum)>0)]]
  par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  par$index_paa_pars = rep(0, sum(n_index_acomp_pars))  
  if(all(data$age_comp_model_fleets %in% c(5,7))){ # start tau/neff at 0
    neff <- data$catch_Neff
    neff[neff <= 0] <- NA
    neff <- apply(neff,2,mean, na.rm=TRUE)[which(apply(data$use_catch_paa,2,sum)>0)]
    par$catch_paa_pars = 0.5*log(neff) # exp(age_comp_pars(0)-0.5*log(Neff))
  }  
  if(all(data$age_comp_model_indices %in% c(5,7))){ # start tau/neff at 0
    neff <- data$index_Neff
    neff[neff <= 0] <- NA
    neff <- apply(neff,2,mean, na.rm=TRUE)[which(apply(data$use_index_paa,2,sum)>0)]
    par$index_paa_pars = 0.5*log(neff) # exp(age_comp_pars(0)-0.5*log(Neff))
  }

	input$data = data
	input$par = par
	return(input)
}
