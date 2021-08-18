set_F = function(input, F_opts = NULL)
{
  asap3 = if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3

  input$par$F_devs = matrix(0, input$data$n_years_model-1, input$data$n_fleets)
	if(!is.null(asap3))
	{
	  input$par$log_F1 = log(asap3$F1_ini) # use F1_ini values from asap3 file  
	}
  else
  {
  	input$par$log_F1 = rep(log(0.2), input$data$n_fleets) # old
    if(!is.null(F_opts$F)) {
      input$par$log_F1 = log(F_opts$F[1,])
      input$par$F_devs = cbind(apply(log(cbind(F_opts$F)),2,diff))
    }
  }  
  return(input)
}