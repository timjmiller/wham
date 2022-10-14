set_random <- function(input){
	random = input$random
  #if(input$data$Ecov_obs_sigma_opt == 4) random = "Ecov_obs_logsigma"
  if(any(input$data$Ecov_obs_sigma_opt==4)) random = c(random, "Ecov_obs_logsigma_re")
  if(any(input$data$selblock_models_re > 1)) random = c(random, "selpars_re")
  if(input$data$M_re_model > 1) random = c(random, "M_re")
  if(sum(input$data$Ecov_model) > 0) random = c(random, "Ecov_re")
  if(any(input$data$NAA_re_mode > 0) random = c(random, "log_NAA")
  if(sum(input$data$use_q_prior)) random = c(random, "q_prior_re")
  if(sum(input$data$use_q_re)) random = c(random, "q_re")
  input$random = random	
  return(input)
}