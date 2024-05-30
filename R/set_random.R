set_random <- function(input){
	random = NULL
  #if(input$data$Ecov_obs_sigma_opt == 4) random = "Ecov_obs_logsigma"
  if(any(input$data$Ecov_obs_sigma_opt==4)) random = c(random, "Ecov_obs_logsigma_re")
  if(any(input$data$selblock_models_re > 1)) random = c(random, "selpars_re")
  if(any(input$data$M_re_model > 1)) random = c(random, "M_re")
  if(any(input$data$use_b_prior>0)) random <- c(random, "log_b")
  if(any(input$data$Ecov_model > 0)) random = c(random, "Ecov_re")
  if(any(input$data$NAA_re_model > 0)) random = c(random, "log_NAA")
  if(any(input$data$N1_model == 2)) random = c(random, "log_N1")
  if(any(input$data$use_mu_prior > 0)) random = c(random, "mu_prior_re")
  #print("here")
  #print(input$data$use_mu_prior)
  #print(input$data$mu_model)
  if(any(input$data$mu_model %in% c(2:4,6:8,10:12,14:16))) random = c(random, "mu_re")
  if(sum(input$data$use_q_prior)) random = c(random, "q_prior_re")
  if(sum(input$data$use_q_re)) random = c(random, "q_re")
  input$random = random	
  return(input)
}