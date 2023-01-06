#adapted from stockassessment::procres
get_post_samp <- function(fit){
  
  input = fit$input
  input$par = fit$parList
  re_names = c("log_NAA", "M_re", "selpars_re", "Ecov_re", "q_re", "growth_re", "LAA_re", "LW_re", "WAA_re")
  input$data$do_post_samp[] = as.integer(re_names %in% input$random)
  res = list()
  if(sum(input$data$do_post_samp)){
    if(fit$is_sdrep) if(!fit$na_sdrep){
      tfit = fit_wham(input, do.fit = FALSE)
      sdrep = sdreport(tfit)
      for(i in re_names){
        if(re_names[i] %in% input$random) {
          idx <- which(names(sdrep$value)==re_names[i])
          if(res_names[i] %in% c("log_NAA", "M_re", "q_re", "LAA_re", "WAA_re")) {
            res[[re_names[i]]] <- rmvnorm(1,mu=sdrep$value[idx], Sigma=sdrep$cov[idx,idx]) # does this work when M_Re = ar1_a?
            res[[re_names[i]]] <- matrix(res[[re_names[i]]], nrow=length(fit$years_full))
          }
          if(res_names[i] %in% c("growth_re")) {
            mod_re = input$data$growth_re_model
            n_par = input$data$n_growth_par
          }
          if(res_names[i] %in% c("LW_re")) {
            mod_re = input$data$LW_re_model
            n_par = input$data$n_LW_par
          }
          if(res_names[i] == "Ecov_re") { 
            res[[re_names[i]]] <- rmvnorm(1,mu=sdrep$value[idx], Sigma=sdrep$cov[idx,idx])
            res[[re_names[i]]] <- matrix(res[[re_names[i]]], nrow=input$data$n_years_Ecov + input$data$n_years_proj_Ecov)
          }
        }
      }
      return(res)
    } else {
      warning("sdreport did not complete successfully. Therefore sample of posterior for random effects not possible.")
    }
  } else {
    cat("No process errors specified in model. Therefore no posterior to sample.\n")
  }
}
