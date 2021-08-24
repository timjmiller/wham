set_q = function(input, q_opts = NULL){
	
	if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3
  
  input$data$q_lower <- rep(0,input$data$n_indices)
  input$data$q_upper <- rep(1000,input$data$n_indices)
  if(!is.null(q_opts$q_lower)) input$data$q_lower = q_opts$q_lower
  if(!is.null(q_opts$q_upper)) input$data$q_upper = q_opts$q_upper
	if(is.null(q_opts$q)) {
  	if(!is.null(asap3)) input$par$logit_q = gen.logit(asap3$q_ini[which(asap3$use_index ==1)], input$data$q_lower, input$data$q_upper) # use q_ini values from asap3 file
		else input$par$logit_q = rep(gen.logit(0.3, input$data$q_lower, input$data$q_upper), input$data$n_indices)
	}
	else input$par$logit_q = gen.logit(q_opts$q, input$data$q_lower, input$data$q_upper)
	return(input)
}