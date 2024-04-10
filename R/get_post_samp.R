#adapted from stockassessment::procres
get_post_samp <- function(fit){
  
  input = fit$input
  input$par = fit$parList
  re_names = c("log_NAA", "M_re", "mu_re", "selpars_re", "Ecov_re", "q_re")
  post_samp_names <- paste0("do_post_samp_", c("N", "M","mu", "sel", "Ecov", "q"))
  ind <- which(re_names %in% input$random)
  for(i in ind) input$data[[post_samp_names[i]]] = 1
  res = list()
  if(length(ind)){
    if(fit$is_sdrep) if(!fit$na_sdrep){
      # tfit = fit_wham(input, do.fit = FALSE)
      # sdrep = TMB::sdreport(tfit)
      for(i in re_names[ind]){
        if(i %in% input$random) {
          idx <- which(names(sdrep$value)==i)
          Sig <- sdrep$cov[idx,idx]
          L <- t(chol(Sig))
          Z <- rnorm(length(idx))
          res[i] <- L * Z + sdrep$value[idx]
          #res[i] <- rmvnorm(1,mu=sdrep$value[idx], Sigma=)
          # if(res_names[i] %in% c("log_NAA", "M_re", "q_re")) {
          #   res[[re_names[i]]] <- matrix(res[[re_names[i]]], nrow=length(fit$years_full))
          # }
          # if(res_names[i] == "Ecov_re") {
          #   res[[re_names[i]]] <- matrix(res[[re_names[i]]], nrow=input$data$n_years_Ecov + input$data$n_years_proj_Ecov)
          # }
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
