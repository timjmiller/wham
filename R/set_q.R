set_q = function(input, q_opts = NULL){
	
  data = input$data
  par = input$par
  map = input$map
  
	if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3
  
  data$q_lower <- rep(0,data$n_indices)
  data$q_upper <- rep(1000,data$n_indices)
  if(!is.null(q_opts$q_lower)) {
		if(length(q_opts$q_lower) != data$n_indices) stop("the length of q_opts$q_lower is not equal to the number of indices")
  	data$q_lower = q_opts$q_lower
  }
  if(!is.null(q_opts$q_upper)) {
		if(length(q_opts$q_upper) != data$n_indices) stop("the length of q_opts$q_upper is not equal to the number of indices")
  	data$q_upper = q_opts$q_upper
  }
	par$logit_q = gen.logit(0.3, data$q_lower, data$q_upper)
  if(!is.null(asap3)) par$logit_q = gen.logit(asap3$q_ini[which(asap3$use_index ==1)], data$q_lower, data$q_upper) # use q_ini values from asap3 file
	if(!is.null(q_opts$initial_q)) {
		if(length(q_opts$initial_q) != data$n_indices) stop("the length of q_opts$initial_q is not equal to the number of indices")
		par$logit_q = gen.logit(q_opts$initial_q, data$q_lower, data$q_upper)
	}
	data$use_q_re = rep(0, data$n_indices)
	par$q_re = matrix(0, data$n_years_model, data$n_indices)
	map$q_re = matrix(NA, data$n_years_model, data$n_indices)
	par$q_repars = matrix(0, data$n_indices, 2)
	map$q_repars = matrix(NA, data$n_indices, 2)
	if(!is.null(q_opts$re)){
  	ind = which(q_opts$re %in% c("iid", "ar1"))
		map$q_re[,ind] = 1:(length(ind) * (data$n_years_model)) #turn on appropriate columns of q_re
		data$use_q_re[ind] = 1
		map$q_repars[ind,1] = 1:length(ind)
  	ind = which(q_opts$re == "ar1")
		map$q_repars[ind,2] = sum(!is.na(map$q_repars[,1])) + 1:length(ind)
	}
	map$q_repars = factor(map$q_repars)
  map$q_re = factor(map$q_re)

	data$use_q_prior = rep(0, data$n_indices)
	par$q_prior_re = rep(0, data$n_indices)
	map$q_prior_re = rep(NA, data$n_indices)
	data$logit_q_prior_sigma = rep(1, data$n_indices)
	map$logit_q = 1:data$n_indices
	if(!is.null(q_opts$prior)){
  	ind = which(!is.na(q_opts$prior_sd))
		data$use_q_prior[ind] = 1
		data$logit_q_prior_sigma[ind] = q_opts$prior_sd[ind]
		par$q_prior_re[ind] = par$logit_q[ind]
		map$q_prior_re[ind] = 1:length(ind)
		map$logit_q[ind] = NA #turn off logit q (mean of the prior) when random effects are used
	}
	map$q_prior_re = factor(map$q_prior_re)
	map$logit_q = factor(map$logit_q)

  input$data = data
  input$par = par
  input$map = map

	return(input)
}