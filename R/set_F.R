set_F = function(input, F_opts = NULL)
{
  asap3 = input$asap3
  input$par$log_F <- matrix(NA,input$data$n_years_model, input$data$n_fleets)
  F_devs = matrix(0, input$data$n_years_model-1, input$data$n_fleets)
	if(!is.null(asap3)) {
    k = 1
    for(i in 1:length(asap3))for(j in 1:length(asap3[[i]]$F1_ini)){
      input$par$log_F[,k] = log(asap3[[i]]$F1_ini[j]) # use F1_ini values from asap3 file  
      k = k + 1
    }
	} else {
  	input$par$log_F[] = log(0.2) # old
  }  
  if(!is.null(F_opts$F)) {
    input$par$log_F[] = log(F_opts$F)
    #input$par$F_devs = cbind(apply(log(cbind(F_opts$F)),2,diff))
  }
  return(input)
}