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
          res[[re_names[i]]] <- rmvnorm(1,mu=sdrep$value[idx], Sigma=sdrep$cov[idx,idx]) # does this work when M_Re = ar1_a?
          if(res_names[i] %in% c("log_NAA", "M_re", "q_re", "LAA_re", "WAA_re")) {
            res[[re_names[i]]] <- matrix(res[[re_names[i]]], nrow=length(fit$years_full))
          }
          if(res_names[i] %in% c("growth_re")) {
            in_res = res[[re_names[i]]]
            i_min = 1
            i_max = 0
            mod_re = input$data$growth_re_model
            n_par = input$data$n_growth_par
            out_array = array(0, dim = c(length(fit$years_full), input$data$n_ages, n_par))
            for(k in 1:n_par) {
              if(mod_re[k] == 1) next
              if(mod_re[k] %in% c(2,4)) {# year effect
                n_re = length(fit$years_full) 
                i_max = i_max + n_re
                in_res = in_res[i_min:i_max]
                out_array[,,k] = matrix(in_res[i_min:i_max], nrow=fit$years_full, ncol = input$data$n_ages)
                i_min = i_max + 1
              }
              if(mod_re[k] %in% c(3,5)) {# cohort effect
                n_re = input$data$n_ages + length(fit$years_full) - 1 
                i_max = i_max + n_re
                in_res = in_res[i_min:i_max]
                out_array[,,k] = sort_cohort_matrix(in_res, fit$years_full, input$data$n_ages)
                i_min = i_max + 1
              }
            }
            res[[re_names[i]]] = out_array
          }
          if(res_names[i] %in% c("LW_re")) {
            in_res = res[[re_names[i]]]
            i_min = 1
            i_max = 0
            mod_re = input$data$LW_re_model
            n_par = input$data$n_LW_par
            out_array = array(0, dim = c(length(fit$years_full), input$data$n_ages, n_par))
            for(k in 1:n_par) {
              if(mod_re[k] == 1) next
              if(mod_re[k] %in% c(2,4)) {# year effect
                n_re = length(fit$years_full) 
                i_max = i_max + n_re
                in_res = in_res[i_min:i_max]
                out_array[,,k] = matrix(in_res[i_min:i_max], nrow=fit$years_full, ncol = input$data$n_ages)
                i_min = i_max + 1
              }
              if(mod_re[k] %in% c(3,5)) {# cohort effect
                n_re = input$data$n_ages + length(fit$years_full) - 1 
                i_max = i_max + n_re
                in_res = in_res[i_min:i_max]
                out_array[,,k] = sort_cohort_matrix(in_res, fit$years_full, input$data$n_ages)
                i_min = i_max + 1
              }
            }
            res[[re_names[i]]] = out_array
          }
          if(res_names[i] == "Ecov_re") { 
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

#auxiliary function for cohort effects (growth and LW)
sort_cohort_matrix = function(vector, nrows, ncols) {

      tmp1 = matrix(NA, ncol = ncols, nrow = nrows)
      loop_row = rep(0, times = ncols + nrows - 1)
      loop_row[1:ncols] = (ncols - 1):0
      loop_col = rep(ncols - 1, times = ncols + nrows - 1)
      loop_col[(length(loop_col) - ncols + 1):length(loop_col)] = (ncols - 1):0

      for(j in seq_along(loop_col)) {
        tmp1[(loop_row[j]:loop_col[j])*(nrows + 1) + (j - ncols + 1)] <- vector[j]
      }

      return(tmp1)
}