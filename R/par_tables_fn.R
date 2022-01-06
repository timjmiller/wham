par_tables_fn = function(mod, do.tex=FALSE, do.html=FALSE, od)
{
  library(rmarkdown)

  ci = function(par,se, p=0.975, lo = 0, hi = 1, type = "I"){
    ci = par + c(-1,1) * qnorm(0.975) * se

    if(type == "I") {
      return(c(se,ci))
    }
    if(type == "exp") {
      return(c(exp(par)*se, exp(ci)))
    }
    if(type == "expit") { #Delta-method: V(lo + (hi-lo)/(1 + exp(-x))) ~ ((hi-lo) * p * (1-p))^2 * V(x)
      p = 1/(1 + exp(-par))
      dm.se = abs(hi-lo)*p*(1-p)*se
      return(c(dm.se, lo + (hi-lo)/(1+ exp(-ci))))
    }
  }
  data = mod$env$data
  if(mod$is_sdrep) {
    sdrep = mod$sdrep
    pars = as.list(sdrep, "Est")
    sd = as.list(sdrep, "Std")
  } else {
    pars = mod$parList
    sd = lapply(pars, function(x) x[] = NA)
  }

  fe.names = character()
  fe.vals = numeric()
  fe.cis = matrix(nrow = 0,ncol = 3)

  if(data$n_NAA_sigma>0) {
    if(data$recruit_model == 2) { #random about mean
      tvar = length(unique(mod$rep$pred_NAA[-1,1])) != 1 #see if anything is causing mean recruitment to vary over time
      if(!tvar) {
        fe.names = c(fe.names, "Mean Recruitment")
        fe.vals = c(fe.vals, exp(pars$mean_rec_pars[1]))
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1], type = "exp"))
      } else{
        fe.names = c(fe.names, "mean log(R) intercept")
        fe.vals = c(fe.vals, pars$mean_rec_pars[1])
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1]))
      }
    }
    if(data$recruit_model == 3){
      if(data$use_steepness == 0){
        tvar_a = length(unique(mod$rep$log_SR_a)) != 1 #see if anything is causing mean recruitment to vary over time
        tvar_b = length(unique(mod$rep$log_SR_b)) != 1 #see if anything is causing mean recruitment to vary over time
        if(!tvar_a){
          fe.names = c(fe.names, "B-H a")
          fe.vals = c(fe.vals, exp(pars$mean_rec_pars[1]))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1], type = "exp"))
        } else{
          fe.names = c(fe.names, "mean log(B-H a) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[1])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1]))          
        }
        if(!tvar_b){
          fe.names = c(fe.names, "B-H b")
          fe.vals = c(fe.vals, exp(pars$mean_rec_pars[2]))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2], type = "exp"))
        } else{
          fe.names = c(fe.names, "mean log(B-H b) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[2])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2]))          
        }
      }
      if(data$use_steepness == 1){
        tvar_h = length(unique(mod$rep$SR_h_tf)) != 1 #see if anything is causing mean recruitment to vary over time
        tvar_RO = length(unique(mod$rep$log_SR_R0)) != 1 #see if anything is causing mean recruitment to vary over time
        if(!tvar_h){
          fe.names = c(fe.names, "B-H h")
          fe.vals = c(fe.vals, 0.2 + 0.8/(1+exp(-pars$mean_rec_pars[1])))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1],lo=0.2, type = "expit")) 
        } else{
          fe.names = c(fe.names, "mean logit(B-H h) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[1])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1]))          
        }
        if(!tvar_R0){
          fe.names = c(fe.names, "B-H $R_0$")
          fe.vals = c(fe.vals, exp(pars$mean_rec_pars[2]))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2], type = "exp"))
        } else{
          fe.names = c(fe.names, "mean log(B-H $R_0$) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[2])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2]))          
        }
      }
    }
    if(data$recruit_model == 4) { #Ricker
      if(data$use_steepness == 0){
        tvar_a = length(unique(mod$rep$log_SR_a)) != 1 #see if anything is causing mean recruitment to vary over time
        tvar_b = length(unique(mod$rep$log_SR_b)) != 1 #see if anything is causing mean recruitment to vary over time
        if(!tvar_a){
          fe.names = c(fe.names, "Ricker a")
          fe.vals = c(fe.vals, exp(pars$mean_rec_pars[1]))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1], type = "exp"))
        } else{
          fe.names = c(fe.names, "mean log(Ricker a) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[1])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1]))          
        }
        if(!tvar_b){
          fe.names = c(fe.names, "Ricker b")
          fe.vals = c(fe.vals, exp(pars$mean_rec_pars[2]))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2], type = "exp"))
        } else{
          fe.names = c(fe.names, "mean log(Ricker b) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[2])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2]))          
        }
      }
      if(data$use_steepness == 1){
        tvar_h = length(unique(mod$rep$SR_h_tf)) != 1 #see if anything is causing mean recruitment to vary over time
        tvar_RO = length(unique(mod$rep$log_SR_R0)) != 1 #see if anything is causing mean recruitment to vary over time
        if(!tvar_h){
          fe.names = c(fe.names, "Ricker h")
          fe.vals = c(fe.vals, 0.2 + exp(pars$mean_rec_pars[1]))
          fe.cis = rbind(fe.cis, 0.2 + ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1], type = "exp")) 
        } else{
          fe.names = c(fe.names, "mean log(Ricker h - 0.2) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[1])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[1],sd$mean_rec_pars[1]))          
        }
        if(!tvar_R0){
          fe.names = c(fe.names, "Ricker $R_0$")
          fe.vals = c(fe.vals, exp(pars$mean_rec_pars[2]))
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2], type = "exp"))
        } else{
          fe.names = c(fe.names, "mean log(Ricker R_0) intercept")
          fe.vals = c(fe.vals, pars$mean_rec_pars[2])
          fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[2],sd$mean_rec_pars[2]))          
        }
      }
    }

    ages = mod$ages
    ind = unique(data$NAA_sigma_pointer)
    npar = length(ind)
    al = ah = integer()
    if(data$n_NAA_sigma>1){
      for(i in 1:npar) {
        al = c(al, ages[min(which(data$NAA_sigma_pointer == ind[i]))])
        ah = c(ah, ages[max(which(data$NAA_sigma_pointer == ind[i]))])
        more = ah[i] != al[i]
        fe.names = c(fe.names, paste0("NAA $\\sigma$ (age ", al[i], ifelse(more, paste0("-", ah[i], ")"), ")")))
        fe.vals = c(fe.vals, exp(pars$log_NAA_sigma[i]))
        fe.cis = rbind(fe.cis, ci(pars$log_NAA_sigma[i], sd$log_NAA_sigma[i], type = "exp"))
      }
      fe.names = c(fe.names, paste("NAA residual AR1 $\\rho$", c("age", "year")))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(-pars$trans_NAA_rho)))
      fe.cis = rbind(fe.cis, 
        ci(pars$trans_NAA_rho[1], sd$trans_NAA_rho[1], lo = -1, hi = 1, type = "expit"),
        ci(pars$trans_NAA_rho[2], sd$trans_NAA_rho[2], lo = -1, hi = 1, type = "expit"))
    }
    else {
      fe.names = c(fe.names, paste0("NAA $\\sigma$ (age ", ages[1], ")"))
      fe.vals = c(fe.vals, exp(pars$log_NAA_sigma[1]))
      fe.cis = rbind(fe.cis, ci(pars$log_NAA_sigma[1], sd$log_NAA_sigma[1], type = "exp"))
      fe.names = c(fe.names, paste("NAA residual AR1 $\\rho$", "year"))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(-pars$trans_NAA_rho[2])))
      fe.cis = rbind(fe.cis, 
        ci(pars$trans_NAA_rho[2], sd$trans_NAA_rho[2], lo = -1, hi = 1, type = "expit"))
    }
  }

  for(i in 1:data$n_indices) {
    tvar_q = length(unique(mod$rep$logit_q_mat[,i])) != 1  #see if anything is causing q to vary over time
    parname = ifelse(data$use_q_prior[i]== 0, "logit_q", "q_prior_re")
    if(!tvar_q){
      fe.names = c(fe.names, paste("Index" , i, "fully selected q"))
      fe.vals = c(fe.vals, data$q_lower[i] + (data$q_upper[i]-data$q_lower[i])/(1 + exp( - pars[[parname]][i])))
      fe.cis = rbind(fe.cis, ci(pars[[parname]][i],sd[[parname]][i], lo = data$q_lower[i],hi = data$q_upper[i], type = "expit"))
    } else {
      fe.names = c(fe.names, paste("Index" , i, "logit(q) intercept"))
      fe.vals = c(fe.vals, pars[[parname]][i])
      fe.cis = rbind(fe.cis, ci(pars[[parname]][i],sd[[parname]][i]))
    }
  }
  for(i in 1:data$n_selblocks){
    if(all(apply(mod$rep$selAA[[i]], 2, function(x) length(unique(x))) == 1)){
      extra.name = ""
    } else extra.name = "Mean "

    if(data$selblock_models[i] == 1) {
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.name, "Selectivity for age ", mod$ages))
      ind = 1:data$n_ages
    }
    if(data$selblock_models[i] == 2){ #increasing logistic
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.name, c("$a_{50}$", "1/slope (increasing)")))
      ind = data$n_ages + 1:2
    }
    if(data$selblock_models[i] == 3){ #double logistic
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.name, c("$a_{50}$ (1)", "1/slope (1)","$a_{50}$ (2)", "1/slope (2)")))
      ind = data$n_ages + 3:6
    }
    if(data$selblock_models[i] == 4){ #increasing logistic
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.name, c("$a_{50}$", "-1/slope (decreasing)")))
      ind = data$n_ages + 1:2
    }
    fe.vals = c(fe.vals, ((data$selpars_lower + data$selpars_upper-data$selpars_lower)/(1 + exp(-pars$logit_selpars)))[i,ind])
    for(a in ind) {
      fe.cis = rbind(fe.cis, ci(pars$logit_selpars[i,a], sd$logit_selpars[i,a], lo = data$selpars_lower[i,a], hi = data$selpars_upper[i,a], type = "expit"))
    }
  }
  for(i in 1:data$n_selblocks) if(data$selblock_models_re[i]>1){
    fe.names = c(fe.names, paste0("Block ", i , ": Selectivity RE $\\sigma$"))
    fe.vals = c(fe.vals, exp(pars$sel_repars[i,1]))
    fe.cis = rbind(fe.cis, ci(pars$sel_repars[i,1], sd$sel_repars[i,1], type = "exp"))
    if(data$selblock_models_re[i] %in% c(3,5)) {
      modify = ""
      if(data$selblock_models[i] == 1) modify = " AR1 $\\rho$ (age)"
      if(data$selblock_models[i] %in% c(2,4)) modify = " $\\rho$ for $a_{50}$ and 1/slope" 
      if(data$selblock_models[i] == 3) modify = " AR1 $\\rho$ for double-logistic pars"
      fe.names = c(fe.names, paste0("Block ", i , ": Selectivity RE", modify))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(-pars$sel_repars[i,2])))
      fe.cis = rbind(fe.cis, ci(pars$sel_repars[i,2], sd$sel_repars[i,2], lo = -1, hi = 1, type = "expit"))
    }
    if(data$selblock_models_re[i] %in% c(4,5)) {
      fe.names = c(fe.names, paste0("Block ", i , ": Selectivity RE AR1 $\\rho$ (year)"))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(-pars$sel_repars[i,3])))
      fe.cis = rbind(fe.cis, ci(pars$sel_repars[i,3], sd$sel_repars[i,3], lo = -1, hi = 1, type = "expit"))
    }
  }
  acomp_par_count = 0
  for(i in 1:data$n_fleets){
    if(sum(data$use_catch_paa[,i]) > 0){
      if(data$age_comp_model_fleets[i] %in% c(2:3,5:7)){
        if(data$age_comp_model_fleets[i] == 2){
          fe.names = c(fe.names, paste0("Fleet ", i , " age comp, Dirichlet-multinomial: dispersion ($\\phi$)"))
          ind = acomp_par_count+1
        }
        if(data$age_comp_model_fleets[i] == 3){
          fe.names = c(fe.names, paste0("Fleet ", i , " age comp, Dirichlet: dispersion ($\\phi$)"))
          ind = acomp_par_count+1
        }
        if(data$age_comp_model_fleets[i] == 5){
          fe.names = c(fe.names, paste0("Fleet ", i , " age comp, logistic-normal (0-pooled): $\\sigma$"))
          ind = acomp_par_count+1
        }
        if(data$age_comp_model_fleets[i] == 6){
          fe.names = c(fe.names, paste0("Fleet ", i , " age comp, 0/1-inflated logistic-normal: ",
            c("Binomial N parameter probablity of 0",
              "logistic-normal $\\sigma$")))
          ind = acomp_par_count+1:2
        }
        if(data$age_comp_model_fleets[i] == 7){
          fe.names = c(fe.names, paste0("Fleet ", i , " age comp, logistic-normal (0-missing): $\\sigma$"))
          ind = acomp_par_count+1
        }
        fe.vals = c(fe.vals, exp(pars$catch_paa_pars[ind]))
        fe.cis = rbind(fe.cis, ci(pars$catch_paa_pars[ind], sd$catch_paa_pars[ind], type = "exp"))
      }
      if(data$age_comp_model_fleets[i] == 4){
        fe.names = c(fe.names, paste0("Fleet ", i , " age comp, 0/1-inflated logistic-normal: ",
          c("Declining probablity of 0 parameter 1",
            "Declining probablity of 0 parameter 2",
            "logistic-normal $\\sigma$")))
        ind = acomp_par_count+1:2
        fe.vals = c(fe.vals, pars$catch_paa_pars[ind])
        fe.cis = rbind(fe.cis, ci(pars$catch_paa_pars[ind], sd$catch_paa_pars[ind]))
        ind = acomp_par_count+3
        fe.vals = c(fe.vals, exp(pars$catch_paa_pars[ind]))
        fe.cis = rbind(fe.cis, ci(pars$catch_paa_pars[ind], sd$catch_paa_pars[ind], type = "exp"))
      }
      acomp_par_count = acomp_par_count + data$n_age_comp_pars_fleets[i]
    }
  }
  acomp_par_count = 0
  for(i in 1:data$n_indices){
    if(sum(data$use_index_paa[,i]) > 0){
      if(data$age_comp_model_indices[i] %in% c(2:3,5:7)){
        if(data$age_comp_model_indices[i] == 2){
          fe.names = c(fe.names, paste0("Index ", i , " age comp, Dirichlet-multinomial: dispersion ($\\phi$)"))
          ind = acomp_par_count+1
        }
        if(data$age_comp_model_indices[i] == 3){
          fe.names = c(fe.names, paste0("Index ", i , " age comp, Dirichlet: dispersion ($\\phi$)"))
          ind = acomp_par_count+1
        }
        if(data$age_comp_model_indices[i] == 5){
          fe.names = c(fe.names, paste0("Index ", i , " age comp, logistic-normal (0-pooled): $\\sigma$"))
          ind = acomp_par_count+1
        }
        if(data$age_comp_model_indices[i] == 6){
          fe.names = c(fe.names, paste0("Index ", i , " age comp, 0/1-inflated logistic-normal: ",
            c("Binomial N parameter probablity of 0",
              "logistic-normal $\\sigma$")))
          ind = acomp_par_count+1:2
        }
        if(data$age_comp_model_indices[i] == 7){
          fe.names = c(fe.names, paste0("Index ", i , " age comp, logistic-normal (0-missing): $\\sigma$"))
          ind = acomp_par_count+1
        }
        fe.vals = c(fe.vals, exp(pars$index_paa_pars[ind]))
        fe.cis = rbind(fe.cis, ci(pars$index_paa_pars[ind], sd$index_paa_pars[ind], type = "exp"))
      }
      if(data$age_comp_model_indices[i] == 4){
        fe.names = c(fe.names, paste0("Index ", i , " age comp, 0/1-inflated logistic-normal: ",
          c("Declining probablity of 0 parameter 1",
            "Declining probablity of 0 parameter 2",
            "logistic-normal $\\sigma$")))
        ind = acomp_par_count+1:2
        fe.vals = c(fe.vals, pars$index_paa_pars[ind])
        fe.cis = rbind(fe.cis, ci(pars$index_paa_pars[ind], sd$index_paa_pars[ind]))
        ind = acomp_par_count+3
        fe.vals = c(fe.vals, exp(pars$index_paa_pars[ind]))
        fe.cis = rbind(fe.cis, ci(pars$index_paa_pars[ind], sd$index_paa_pars[ind], type = "exp"))
      }
      acomp_par_count = acomp_par_count + data$n_age_comp_pars_indices[i]
    }
  }

  if(sum(!is.na(mod$input$map$M_a))){ #any M_a estimated?
    if(data$M_re_model == 1 & data$Ecov_where[2] == 0 & data$M_model %in% 1:2){ #no random effects, ecov or WAA effects on M
      modify = "M for ages("
    } else {
      if(data$M_model != 3) modify = "mean log(M) for ages ("
      if(data$M_model == 3 | data$Ecov_where[2] == 0) modify = "mean log(M) intercept for log(WAA) effects"
      if(data$M_model != 3 | data$Ecov_where[2] == 1) modify = "mean log(M) intercept for ages ("
    }
    age.list = M_a_point = list()
    M_map = as.integer(as.character(input$map$M_a))
    ind = unique(M_map[which(!is.na(M_map))])
    if(data$M_model == 1) {
      M_a_point[[1]] = 1
      ages.list = list(mod$ages)
    }
    if(data$M_model == 2){
      npar = length(ind)
      for(i in 1:npar) {
        M_a_point[[i]] = which(M_map == ind[i])[1]
        age.list[[i]] = mod$ages[which(M_map == ind[i])]
      }
    }
    if(length(age.list)){
      for(i in 1:length(age.list)){
        fe.names = c(fe.names, paste0(modify, paste0(age.list[[i]], collapse = ", "),")"))
        if(data$M_re_model == 1 & data$Ecov_where[2] == 0 & data$M_model %in% 1:2){
          fe.vals = c(fe.vals, exp(pars$M_a[M_a_point]))
          fe.cis = rbind(fe.cis, ci(pars$M_a[M_a_point], sd$M_a[M_a_point], type = "exp"))
        } else {
          fe.vals = c(fe.vals, pars$M_a[M_a_point])
          fe.cis = rbind(fe.cis, ci(pars$M_a[M_a_point], sd$M_a[M_a_point]))
        }
      }
    }
  }
  if(data$M_model == 3){
    fe.names = c(fe.names, "mean log(M) log(WAA) effect")
    fe.vals = c(fe.vals, exp(pars$log_b))
    fe.cis = rbind(fe.cis, ci(pars$log_b, sd$log_b, type = "exp"))
  }
  if(data$M_re_model>1){
    fe.names = c(fe.names, "M RE $\\sigma$")
    fe.vals = c(fe.vals, exp(pars$M_repars[1]))
    fe.cis = rbind(fe.cis, ci(pars$M_repars[1], sd$M_repars[1], type = "exp"))
    if(data$M_re_model %in% c(3,5)){
      fe.names = c(fe.names, "M RE AR1 $\\rho$ (age)")
      fe.vals = c(fe.vals, exp(pars$M_repars[2]))
      fe.cis = rbind(fe.cis, ci(pars$M_repars[2], sd$M_repars[2], lo = -1, hi = 1, type = "expit"))
    }
    if(data$M_re_model %in% c(4,5)) {
      fe.names = c(fe.names, "M RE AR1 $\\rho$ (year)")
      fe.vals = c(fe.vals, exp(pars$M_repars[3]))
      fe.cis = rbind(fe.cis, ci(pars$M_repars[3], sd$M_repars[3], lo = -1, hi = 1, type = "expit"))
    }
  }
  if(sum(!is.na(mod$input$map$log_catch_sig_scale))){ #any agg catch obs var estimated?
    sig_map = as.integer(as.character(mod$input$map$log_catch_sig_scale))
    for(i in 1:data$n_fleets){
      if(!is.na(sig_map[i])){
        fe.names = c(fe.names, paste0("Fleet ", i, ": log-catch observation SD scalar"))
        fe.vals = c(fe.vals, exp(pars$log_catch_sig_scale[i]))
        fe.cis = rbind(fe.cis, ci(pars$log_catch_sig_scale[i], sd$log_catch_sig_scale[i], type = "exp"))
      }
    }
  }
  if(sum(!is.na(mod$input$map$log_index_sig_scale))){ #any agg index obs var estimated?
    sig_map = as.integer(as.character(mod$input$map$log_index_sig_scale))
    for(i in 1:data$n_fleets){
      if(!is.na(sig_map[i])){
        fe.names = c(fe.names, paste0("Index ", i, ": log-index observation SD scalar"))
        fe.vals = c(fe.vals, exp(pars$log_index_sig_scale[i]))
        fe.cis = rbind(fe.cis, ci(pars$log_index_sig_scale[i], sd$log_index_sig_scale[i], type = "exp"))
      }
    }
  }

  for(i in 1:data$n_Ecov){
    if(data$Ecov_model[i] > 0){
      fe.vals = c(fe.vals, pars$Ecov_process_pars[1,i])
      fe.cis = rbind(fe.cis, ci(pars$Ecov_process_pars[1,i], sd$Ecov_process_pars[1,i]))
      if(data$Ecov_model[i] == 1){
        fe.names = c(fe.names, paste0("Ecov ", data$Ecov_label[[1]][i], ": ", c("Ecov$_1$", "RW $\\sigma$")))
      }
      if(data$Ecov_model[i] == 2){
        fe.names = c(fe.names, paste0("Ecov ", data$Ecov_label[[1]][i], ": ", c("AR1 $\\mu$", "AR1 $\\rho$", "AR1 $\\sigma$")))
        fe.vals = c(fe.vals, -1 + 2/(1 + exp(-pars$Ecov_process_pars[3,i])))
        fe.cis = rbind(fe.cis, ci(pars$Ecov_process_pars[3,i], sd$Ecov_process_pars[3,i], lo = -1, hi = 1, type = "expit"))
      }
      fe.vals = c(fe.vals, exp(pars$Ecov_process_pars[2,i]))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_process_pars[2,i], sd$Ecov_process_pars[2,i], type = "exp"))

    }
  }

  ecov_beta_map = array(dim = dim(pars$Ecov_beta)) #(2 + n_indices) x n_poly x n_ecov x n_ages
  ecov_beta_map[] = as.integer(as.character(mod$input$map$Ecov_beta))
  for(i in 1:data$n_Ecov){
    if(data$Ecov_where[i,1] == 1){ #Recruitment
      if(any(!is.na(ecov_beta_map[1,,i,1]))){
        npoly = sum(!is.na(ecov_beta_map[1,,i,1]))
        fe.names = c(fe.names, paste0("Recruitment Ecov: ", data$Ecov_label[[1]][i], " $\\beta_", 1:npoly, "$"))  
        for(p in 1:npoly) {
          fe.vals = c(fe.vals, pars$Ecov_beta[1,p,i,1])
          fe.cis = rbind(fe.cis, ci(pars$Ecov_beta[1,p,i,1], sd$Ecov_beta[1,p,i,1]))
        }
      }
    }
    if(data$Ecov_where[i,2] == 1){ #M
      if(any(!is.na(ecov_beta_map[2,,i,]))){
        npoly = dim(ecov_beta_map)[2]
        ecov_beta_map_i = matrix(ecov_beta_map[2,,i,],nrow = npoly)
        ecov_beta_par_i = matrix(pars$Ecov_beta[2,,i,],nrow = npoly)
        npoly = max(apply(ecov_beta_map_i,2, function(x) sum(!is.na(x))))
        beta_ind = unique(ecov_beta_map_i[!is.na(ecov_beta_map_i)])
        ages.list = list()
        for(k in 1:length(beta_ind)){
          #find the ages which the beta is being applied to
          ages.list[[k]] = which(apply(ecov_beta_map_i,2, function(x) any(x == k)))
          #find the order of the polynomial. Ecov_beta should never be mapped to apply to more than one order of the orthogonal polynomial
          poly.ind = which(apply(ecov_beta_map_i,1, function(x) any(x == k)))[1]
          if(length(ages.list[[k]])){
            modify = ifelse(length(ages.list[[k]])>1, "M at ages (", "M at age (")
            fe.names = c(fe.names, paste0(modify, paste0(ages.list[[k]], collapse = ", "), " Ecov: ", data$Ecov_label[[1]][i], " $\\beta_", poly.ind, "$"))
            fe.vals = c(fe.vals, pars$Ecov_beta[1,poly.ind,i,ages.list[[k]][1]])
            fe.cis = rbind(fe.cis, ci(pars$Ecov_beta[1,poly.ind,i,ages.list[[k]][1]], sd$Ecov_beta[1,poly.ind,i,ages.list[[k]][1]]))
          }
        }
      }
    }
    for(j in 1:data$n_indices){
      if(data$Ecov_where[i,2+j] == 1){ #catchability. No reason Ecov_beta parameters would be shared across indices 
        if(any(!is.na(ecov_beta_map[2+j,,i,1]))){
          npoly = sum(!is.na(ecov_beta_map[2+j,,i,1]))
          fe.names = c(fe.names, paste0("Catchability Ecov: ", data$Ecov_label[[1]][i], " $\\beta_", 1:npoly, "$"))  
          for(p in 1:npoly) {
            fe.vals = c(fe.vals, pars$Ecov_beta[2+j,p,i,1])
            fe.cis = rbind(fe.cis, ci(pars$Ecov_beta[2+j,p,i,1], sd$Ecov_beta[2+j,p,i,1]))
          }
        }
      }
    }
  }

  for(i in 1:data$n_Ecov){
    if(data$Ecov_obs_sigma_opt[i] == 2){ #single ecov obs sd estimated
      fe.names = c(fe.names, paste0("Ecov: ", data$Ecov_label[[1]][i], " obs. sd."))
      ind = which(!is.na(matrix(input$map$Ecov_obs_logsigma, NROW(input$par$Ecov_obs_logsigma))[,i]))[1]
      fe.vals = c(fe.vals, exp(pars$Ecov_obs_logsigma[ind,i]))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_obs_logsigma[ind,i], sd$Ecov_obs_logsigma[ind,i], type = "exp"))
    }
    if(data$Ecov_obs_sigma_opt[i] == 4){
      fe.names = c(fe.names, paste0("Ecov: ", data$Ecov_label[[1]][i], " obs. log(sd.) RE ", c("$\\mu$", "$\\sigma$")))
      fe.vals = c(fe.vals, pars$Ecov_obs_sigma_pars[1,i])
      fe.cis = rbind(fe.cis, ci(pars$Ecov_obs_sigma_pars[1,i], sd$Ecov_obs_sigma_pars[1,i]))
      fe.vals = c(fe.vals, exp(pars$Ecov_obs_sigma_par[2,i]))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_obs_sigma_par[2,i], sd$Ecov_obs_sigma_par[2,i], type = "exp"))
    }
  }

  fe = cbind(fe.vals, fe.cis)
  rownames(fe) = fe.names
  colnames(fe) = c("Estimate", "Std. Error", "95\\% CI lower", "95\\% CI upper")
  saveRDS(fe, file = file.path(od,"parameter_estimates_table.RDS"))
  
  #numbers at age
  NAA = NAA.cv = mod$rep$NAA
  if(!mod$na_sdrep) NAA.cv[] = mod$sdrep$sd["log_NAA_rep"]
  NAA.sd = NAA * NAA.cv
  NAA.lo = exp(log(NAA) - qnorm(0.975) * NAA.cv)
  NAA.hi = exp(log(NAA) + qnorm(0.975) * NAA.cv)
  rownames(NAA) = rownames(NAA.cv) = rownames(NAA.lo) = rownames(NAA.hi) = mod$years_full
  colnames(NAA) = colnames(NAA.cv) = colnames(NAA.lo) = colnames(NAA.hi) = mod$ages
  saveRDS(NAA, file = file.path(od,"NAA_table.RDS"))
  saveRDS(NAA.sd, file = file.path(od,"NAA_sd_table.RDS"))
  saveRDS(NAA.lo, file = file.path(od,"NAA_lo_table.RDS"))
  saveRDS(NAA.hi, file = file.path(od,"NAA_hi_table.RDS"))
    
  #Total F at age
  FAA_tot = mod$rep$FAA_tot
  FAA_tot.cv = matrix(NA, NROW(FAA_tot), NCOL(FAA_tot))
  if(!mod$na_sdrep) FAA_tot.cv[] = mod$sdrep$sd["log_FAA_tot"]
  FAA_tot.sd = FAA_tot * FAA_tot.cv
  FAA_tot.lo = exp(mod$rep$FAA_tot - qnorm(0.975) * FAA_tot.cv)
  FAA_tot.hi = exp(mod$rep$FAA_tot + qnorm(0.975) * FAA_tot.cv)
  rownames(FAA_tot) = rownames(FAA_tot.cv) = rownames(FAA_tot.lo) = rownames(FAA_tot.hi) = mod$years_full
  colnames(FAA_tot) = colnames(FAA_tot.cv) = colnames(FAA_tot.lo) = colnames(FAA_tot.hi) = mod$ages
  saveRDS(FAA_tot, file = file.path(od,"FAA_tot_table.RDS"))
  saveRDS(FAA_tot.sd, file = file.path(od,"FAA_tot_sd_table.RDS"))
  saveRDS(FAA_tot.lo, file = file.path(od,"FAA_tot_lo_table.RDS"))
  saveRDS(FAA_tot.hi, file = file.path(od,"FAA_tot_hi_table.RDS"))

  #F at age
  FAA = mod$rep$FAA
  FAA.cv = array(NA, dim = dim(FAA))
  if(!mod$na_sdrep) FAA.cv[] = mod$sdrep$sd["log_FAA"]
  FAA.sd = FAA * FAA.cv
  FAA.lo = exp(mod$rep$FAA - qnorm(0.975) * FAA.cv)
  FAA.hi = exp(mod$rep$FAA + qnorm(0.975) * FAA.cv)
  dimnames(FAA) = dimnames(FAA.cv) = dimnames(FAA.lo) = dimnames(FAA.hi) = list(mod$years_full, paste0("Fleet ", 1:data$n_fleets), mod$ages)
  saveRDS(FAA, file = file.path(od,"FAA_table.RDS"))
  saveRDS(FAA.sd, file = file.path(od,"FAA_sd_table.RDS"))
  saveRDS(FAA.lo, file = file.path(od,"FAA_lo_table.RDS"))
  saveRDS(FAA.hi, file = file.path(od,"FAA_hi_table.RDS"))
  
  wham.dir <- find.package("wham")
  pt = list.files(find.package("wham"), pattern = "par_tables.Rmd", recursive = T, full.names = T)[1]
  file.copy(from=pt, to=od, overwrite=TRUE)
  if(do.html) rmarkdown::render(file.path(od,"par_tables.Rmd"), output_format = "html_document", output_file = file.path(od, "wham_par_tables.html"), quiet = T)
  if(do.tex) rmarkdown::render(file.path(od,"par_tables.Rmd"), output_format = "pdf_document", output_file = file.path(od,"wham_par_tables.pdf"), quiet = T)

  #delete par_tables.Rmd from od
  #file.remove(file.path(od,"par_tables.Rmd"))

}
