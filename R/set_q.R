set_q = function(input, q_opts = NULL){
	
  data = input$data
  par = input$par
  map = input$map
  
	if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3
  
  data$q_lower <- rep(0,data$n_indices)
  data$q_upper <- rep(1000,data$n_indices)
  if(!is.null(q_opts$q_lower)) {
  	ind = which(!is.null(q_opts$q_lower))
  	data$q_lower[ind] = q_opts$q_lower[ind]
  }
  if(!is.null(q_opts$q_upper)) {
  	ind = which(!is.null(q_opts$q_upper))
  	data$q_upper[ind] = q_opts$q_upper[ind]
  }
	par$logit_q = rep(gen.logit(0.3, data$q_lower, data$q_upper), data$n_indices)
  if(!is.null(asap3)) par$logit_q = gen.logit(asap3$q_ini[which(asap3$use_index ==1)], data$q_lower, data$q_upper) # use q_ini values from asap3 file
	if(!is.null(q_opts$initial_q)) {
  	ind = which(!is.null(q_opts$initial_q))
		par$logit_q[ind] = gen.logit(q_opts$initial_q[ind], data$q_lower[ind], data$q_upper[ind])
	}
	data$q_re_model = rep(0, data$n_indices)
	par$q_re = matrix(0, data$n_years_model + data$n_years_proj, data$n_indices)
	map$q_re = matrix(NA, data$n_years_model + data$n_years_proj, data$n_indices)
	par$q_re_pars = matrix(0, data$n_indices, 2)
	map$q_re_pars = matrix(NA, data$n_indices, 2)
	if(!is.null(q_opts$model)){
  	ind = which(q_opts$model %in% c("iid", "ar1")
		map$q_re[,ind] = 1:(length(ind) * (data$n_years_model + data$n_years_proj)) #turn on appropriate columns of q_re
		data$q_re_model[ind] = 1
		map$q_re_pars[ind,1] = 1:length(ind)
  	ind = which(q_opts$model == "ar1")
		map$q_re_pars[ind,2] = sum(!is.na(map$q_re_pars[,1])) + 1:length(ind)
	}
	map$q_re_pars = factor(map$q_re_pars)

	data$use_q_prior = rep(0, data$n_indices)
	par$q_prior_re = rep(0, data$n_indices)
	map$q_prior_re = rep(NA, data$n_indices)
	data$logit_q_prior_sigma = rep(1, data$n_indices)
	if(!is.null(q_opts$prior)){
  	ind = which(!is.null(q_opts$prior))
		data$use_q_prior[ind] = 1
		data$logit_q_prior_sigma[ind] = q_opts$prior[ind]
		map$q_prior_re[ind] = 1:length(ind)
	}
	map$q_prior_re = factor(map$q_prior_re)

  input$data = data
  input$par = par
  input$map = map

	return(input)
}