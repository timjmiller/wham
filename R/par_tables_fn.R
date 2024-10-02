par_tables_fn = function(mod, do.tex=FALSE, do.html=FALSE, od = NULL)
{
  #require(rmarkdown)

  ci = function(par,se, p=0.975, lo = 0, hi = 1, type = "I", k = 1){
    ci = par + c(-1,1) * qnorm(0.975) * se

    if(type == "I") {
      return(c(se,ci))
    }
    if(type == "exp") {
      return(c(exp(par)*se, exp(ci)))
    }
    if(type == "expit") { #Delta-method: V(lo + (hi-lo)/(1 + exp(-x))) ~ ((hi-lo) * p * (1-p))^2 * V(x)
      p = 1/(1 + exp(- k * par))
      dm.se = abs(k) * abs(hi-lo)*p*(1-p)*se
      return(c(dm.se, lo + (hi-lo)/(1+ exp(-k * ci))))
    }
  }
  data = mod$env$data
  if(mod$is_sdrep) {
    sdrep = mod$sdrep
    pars = TMB:::as.list.sdreport(sdrep, "Est")
    sd = TMB:::as.list.sdreport(sdrep, "Std")
  } else {
    pars = mod$parList
    sd = lapply(pars, function(x) {
      x[] <- NA
      return(x)})
  }
  stock.names.tab <- gsub("_", " ", mod$input$stock_names, fixed = TRUE)
  region.names.tab <- gsub("_", " ", mod$input$region_names, fixed = TRUE)
  index.names.tab <- gsub("_", " ", mod$input$index_names, fixed = TRUE)
  fleet.names.tab <- gsub("_", " ", mod$input$fleet_names, fixed = TRUE)
  ecov.names.tab <- gsub("_", " ", mod$input$Ecov_names, fixed = TRUE)
  stock.names.f <- gsub(" ", "_", mod$input$stock_names, fixed = TRUE)
  region.names.f <- gsub(" ", "_", mod$input$region_names, fixed = TRUE)
  index.names.f <- gsub(" ", "_", mod$input$index_names, fixed = TRUE)
  fleet.names.f <- gsub(" ", "_", mod$input$fleet_names, fixed = TRUE)
  ecov.names.f <- gsub(" ", "_", mod$input$Ecov_names, fixed = TRUE)

  fe.names = character()
  fe.vals = numeric()
  fe.cis = matrix(nrow = 0,ncol = 3)

  ages = mod$ages.lab
  map_sig <- array(mod$input$map$log_NAA_sigma, dim = dim(pars$log_NAA_sigma))
  map_rho <- array(mod$input$map$trans_NAA_rho, dim = dim(pars$trans_NAA_rho))
  for(s in 1:data$n_stocks) if(data$NAA_re_model[s]>0) {
    if(data$recruit_model[s] == 2) { #random about mean
      tvar = sum(data$Ecov_how_R[,s]) == 1
      #tvar = length(unique(mod$rep$pred_NAA[s,data$spawn_regions[s],-1,1])) != 1 #see if anything is causing mean recruitment to vary over time
      if(!tvar) {
        fe.names = c(fe.names, paste0(stock.names.tab[s]," Mean Recruitment"))
        fe.vals = c(fe.vals, exp(pars$mean_rec_pars[s,1]))
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,1],sd$mean_rec_pars[s,1], type = "exp"))
      } else{
        fe.names = c(fe.names, paste0(stock.names.tab[s]," mean log(R) intercept"))
        fe.vals = c(fe.vals, pars$mean_rec_pars[s,1])
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,1],sd$mean_rec_pars[s,1]))
      }
    }
    if(data$recruit_model[s] == 3){
      tvar_a = length(unique(mod$rep$log_SR_a[,s])) != 1 #see if anything is causing mean recruitment to vary over time
      tvar_b = length(unique(mod$rep$log_SR_b[,s])) != 1 #see if anything is causing mean recruitment to vary over time
      if(!tvar_a){
        fe.names = c(fe.names, paste0(stock.names.tab[s]," B-H a"))
        fe.vals = c(fe.vals, exp(pars$mean_rec_pars[s,1]))
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,1],sd$mean_rec_pars[s,1], type = "exp"))
      } else{
        fe.names = c(fe.names, paste0(stock.names.tab[s]," mean log(B-H a) intercept"))
        fe.vals = c(fe.vals, pars$mean_rec_pars[s,1])
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,1],sd$mean_rec_pars[s,1]))          
      }
      if(!tvar_b){
        fe.names = c(fe.names, paste0(stock.names.tab[s]," B-H b"))
        fe.vals = c(fe.vals, exp(pars$mean_rec_pars[s,2]))
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,2],sd$mean_rec_pars[s,2], type = "exp"))
      } else{
        fe.names = c(fe.names, paste0(stock.names.tab[s]," mean log(B-H b) intercept"))
        fe.vals = c(fe.vals, pars$mean_rec_pars[s,2])
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,2],sd$mean_rec_pars[s,2]))          
      }
    }
    if(data$recruit_model[s] == 4) { #Ricker
      tvar_a = length(unique(mod$rep$log_SR_a[,s])) != 1 #see if anything is causing mean recruitment to vary over time
      tvar_b = length(unique(mod$rep$log_SR_b[,s])) != 1 #see if anything is causing mean recruitment to vary over time
      if(!tvar_a){
        fe.names = c(fe.names, paste0(stock.names.tab[s]," Ricker a"))
        fe.vals = c(fe.vals, exp(pars$mean_rec_pars[s,1]))
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,1],sd$mean_rec_pars[s,1], type = "exp"))
      } else{
        fe.names = c(fe.names, paste0(stock.names.tab[s]," mean log(Ricker a) intercept"))
        fe.vals = c(fe.vals, pars$mean_rec_pars[s,1])
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,1],sd$mean_rec_pars[s,1]))          
      }
      if(!tvar_b){
        fe.names = c(fe.names, paste0(stock.names.tab[s]," Ricker b"))
        fe.vals = c(fe.vals, exp(pars$mean_rec_pars[s,2]))
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,2],sd$mean_rec_pars[s,2], type = "exp"))
      } else{
        fe.names = c(fe.names, paste0(stock.names.tab[s]," mean log(Ricker b) intercept"))
        fe.vals = c(fe.vals, pars$mean_rec_pars[2])
        fe.cis = rbind(fe.cis, ci(pars$mean_rec_pars[s,2],sd$mean_rec_pars[s,2]))          
      }
    }

    for(r in 1:data$n_regions){
      ind = unique(map_sig[s,r,])
      use.reg.name <- ""
      if(data$n_regions>1) use.reg.name <- paste0(" in ", region.names.tab[r])
      #npar = length(ind)
      #al = ah = integer()
      for(i in ind) {
        #al = c(al, ages[min(which(ind == ind[i]))])
        #ah = c(ah, ages[max(which(ind == ind[i]))])
        as <- ages[which(map_sig[s,r,] == i)]
        diff.ages <- diff(which(map_sig[s,r,] == i))
        if(length(as)){
          if(length(as)==1) fe.names = c(fe.names, paste0(stock.names.tab[s], use.reg.name, " NAA $\\sigma$ (age ", as[1], ")"))
          else{
            if(all(diff.ages==1)) fe.names = c(fe.names, paste0(stock.names.tab[s], use.reg.name, " NAA $\\sigma$ (ages ", as[1], "-", as[length(as)], ")"))
            else fe.names = c(fe.names, paste0(stock.names.tab[s], use.reg.name, " NAA $\\sigma$ (ages ", paste0(as, collapse = ","), ")"))
          }
          #more = ah[i] != al[i]
          #fe.names = c(fe.names, paste0("NAA $\\sigma$ (age ", al[i], ifelse(more, paste0("-", ah[i], ")"), ")")))
          fe.vals = c(fe.vals, exp(pars$log_NAA_sigma[s,r,which(ind==i)[1]]))
          fe.cis = rbind(fe.cis, ci(pars$log_NAA_sigma[s,r,which(ind==i)[1]], sd$log_NAA_sigma[s,r,which(ind==i)[1]], type = "exp"))
        }
      }
    }
    # if(is.null(mod$input$map$trans_NAA_rho)) mod$input$map$trans_NAA_rho <- 1:length(mod$input$par$trans_NAA_rho)
    for(r in 1:data$n_regions){
      use.reg.name <- ""
      if(data$n_regions>1) use.reg.name <- paste0(" in ", region.names.tab[r])
      rho.names <- paste0(" ", c("NAA", "NAA","Recruitment"), " AR1 $\\rho$ ", c("age", "year", "year"))
      for(i in 1:3) if(!is.na(map_rho[s,r,i])){
        # fe.name <- paste(stock.names.tab[s], use.reg.name, rho.names[i])
        fe.names = c(fe.names, paste(stock.names.tab[s], use.reg.name, rho.names[i]))
        fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$trans_NAA_rho[s,r,i])))
        fe.cis = rbind(fe.cis, 
          ci(pars$trans_NAA_rho[s,r,i], sd$trans_NAA_rho[s,r,i], lo = -1, hi = 1, type = "expit") #see trans_rho in helper.cpp
        )
      }
    }
  }

  for(i in 1:data$n_indices) {
    tvar_q = length(unique(mod$rep$logit_q_mat[,i])) != 1  #see if anything is causing q to vary over time
    parname = ifelse(data$use_q_prior[i]== 0, "logit_q", "q_prior_re")
    if(!tvar_q){
      fe.names = c(fe.names, paste(index.names.tab[i], "fully selected q"))
      fe.vals = c(fe.vals, data$q_lower[i] + (data$q_upper[i]-data$q_lower[i])/(1 + exp( - pars[[parname]][i])))
      fe.cis = rbind(fe.cis, ci(pars[[parname]][i],sd[[parname]][i], lo = data$q_lower[i],hi = data$q_upper[i], type = "expit"))
    } else {
      fe.names = c(fe.names, paste(index.names.tab[i], "logit(q) intercept"))
      fe.vals = c(fe.vals, pars[[parname]][i])
      fe.cis = rbind(fe.cis, ci(pars[[parname]][i],sd[[parname]][i]))
      if(mod$input$data$use_q_re[i]){
        fe.names = c(fe.names, paste(index.names.tab[i], "q RE $\\sigma$"))
        fe.vals = c(fe.vals, exp(pars$q_repars[i,1]))
        fe.cis = rbind(fe.cis, ci(pars$q_repars[i,1], sd$q_repars[i,1], type = "exp"))
        if(mod$input$options$q$re[i] == "ar1"){
          fe.names = c(fe.names, paste0(index.names.tab[i], "q RE AR1 $\\rho$ (year)"))
          fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$q_repars[i,2])))
          fe.cis = rbind(fe.cis, ci(pars$q_repars[i,2], sd$q_repars[i,2], lo = -1, hi = 1, type = "expit"))
        }
      }
    }
  }
  block.fleets.indices <- lapply(1:data$n_selblocks, function(x){
    y <- data$selblock_pointer_fleets
    z <- matrix(as.integer(y == x), NROW(y), NCOL(y))
    fleet_ind <- apply(z,2,any)
    out <- fleet.names.tab[which(fleet_ind)]
    y <- data$selblock_pointer_indices
    z <- matrix(as.integer(y == x), NROW(y), NCOL(y))
    index_ind <- apply(z,2,any)
    out <- c(out, index.names.tab[which(index_ind)])
  })
  include.selblock <- sapply(block.fleets.indices, length) > 0
  extra.names = rep("", data$n_selblocks)
  for(i in 1:data$n_selblocks) if(include.selblock[i]){
    extra.names[i] <- paste0(extra.names[i], paste(block.fleets.indices[[i]], collapse = ", "), " ")
  }
  extra.names.mean <- extra.names
  for(i in 1:data$n_selblocks) if(include.selblock[i]){
    if(!all(apply(mod$rep$selAA[[i]], 2, function(x) length(unique(x))) == 1)){
      extra.names.mean[i] <- paste0(extra.names.mean[i], "Mean ")
    }

    if(data$selblock_models[i] == 1) {
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.names.mean[i], "Selectivity for age ", mod$ages.lab))
      ind = 1:data$n_ages
    }
    if(data$selblock_models[i] == 2){ #increasing logistic
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.names.mean[i], c("$a_{50}$", "1/slope (increasing)")))
      ind = data$n_ages + 1:2
    }
    if(data$selblock_models[i] == 3){ #double logistic
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.names.mean[i], c("$a_{50}$ (1)", "1/slope (1)","$a_{50}$ (2)", "1/slope (2)")))
      ind = data$n_ages + 3:6
    }
    if(data$selblock_models[i] == 4){ #increasing logistic
      fe.names = c(fe.names, paste0("Block ", i, ": ", extra.names.mean[i], c("$a_{50}$", "-1/slope (decreasing)")))
      ind = data$n_ages + 1:2
    }
    fe.vals = c(fe.vals, ((data$selpars_lower + data$selpars_upper-data$selpars_lower)/(1 + exp(-pars$logit_selpars)))[i,ind])
    for(a in ind) {
      fe.cis = rbind(fe.cis, ci(pars$logit_selpars[i,a], sd$logit_selpars[i,a], lo = data$selpars_lower[i,a], hi = data$selpars_upper[i,a], type = "expit"))
    }
  }
  for(i in 1:data$n_selblocks) if(data$selblock_models_re[i]>1){
    fe.names = c(fe.names, paste0("Block ", i , ": ", extra.names[i], "Selectivity RE $\\sigma$"))
    fe.vals = c(fe.vals, exp(pars$sel_repars[i,1]))
    fe.cis = rbind(fe.cis, ci(pars$sel_repars[i,1], sd$sel_repars[i,1], type = "exp"))
    if(data$selblock_models_re[i] %in% c(3,5)) {
      modify = ""
      if(data$selblock_models[i] == 1) modify = " AR1 $\\rho$ (age)"
      if(data$selblock_models[i] %in% c(2,4)) modify = " $\\rho$ for $a_{50}$ and 1/slope" 
      if(data$selblock_models[i] == 3) modify = " AR1 $\\rho$ for double-logistic pars"
      fe.names = c(fe.names, paste0("Block ", i ,": ", extra.names[i], "Selectivity RE", modify))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$sel_repars[i,2])))
      fe.cis = rbind(fe.cis, ci(pars$sel_repars[i,2], sd$sel_repars[i,2], lo = -1, hi = 1, type = "expit"))
    }
    if(data$selblock_models_re[i] %in% c(4,5)) {
      fe.names = c(fe.names, paste0("Block ", i ,": ", extra.names[i], "Selectivity RE AR1 $\\rho$ (year)"))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$sel_repars[i,3])))
      fe.cis = rbind(fe.cis, ci(pars$sel_repars[i,3], sd$sel_repars[i,3], lo = -1, hi = 1, type = "expit"))
    }
  }
  #acomp_par_count = 0
  add_age_comp_pars = function(age_comp_models, use_paa, pars, pars_sd, is_fleet = TRUE, fe.names, fe.vals, fe.cis) {
    n_mods = length(age_comp_models)
    if(is_fleet) startnames <- paste0(fleet.names.tab, " in ", region.names.tab[mod$input$data$fleet_regions])
    else startnames <- paste0(index.names.tab, " in ", region.names.tab[mod$input$data$index_regions])
    for(i in 1:n_mods){
      if(sum(use_paa[,i]) > 0){
        if(age_comp_models[i] %in% c(2:5,7,11)){
          if(age_comp_models[i] == 2){
            fe.names = c(fe.names, paste0(startnames[i] , " age comp, Dirichlet-multinomial: dispersion ($\\phi$)"))
            #ind = acomp_par_count+1
          }
          if(age_comp_models[i] %in% 3:4){
            fe.names = c(fe.names, paste0(startnames[i] , " age comp, Dirichlet: dispersion ($\\phi$)"))
            #ind = acomp_par_count+1
          }
          if(age_comp_models[i] %in% c(5,7)){
            fe.names = c(fe.names, paste0(startnames[i] , " age comp, logistic-normal: $\\sigma$"))
            #ind = acomp_par_count+1
          }
          if(age_comp_models[i] == 11){
            fe.names = c(fe.names, paste0(startnames[i] , " age comp, Dirichlet-multinomial (linearized): dispersion ($\\phi$)"))
            #ind = acomp_par_count+1
          }
          fe.vals = c(fe.vals, exp(pars[i,1]))
          fe.cis = rbind(fe.cis, ci(pars[i,1], pars_sd[i,1], type = "exp"))
        }
        if(age_comp_models[i] == 6){
          fe.names = c(fe.names, paste0(startnames[i] , " age comp, logistic-normal: $", c("\\sigma", "\\rho"), "$"))
          fe.vals = c(fe.vals, exp(pars[i,1]), 1/(1+exp(-pars[i,2])))
          fe.cis = rbind(fe.cis, 
            ci(pars[i,1], pars_sd[i,1], type = "exp"),
            ci(pars[i,2], pars_sd[i,2], lo = 0, hi = 1, type = "expit"))
          #ind = acomp_par_count+1
        }
        if(age_comp_models[i] == 8){
          fe.names = c(fe.names, paste0(startnames[i], " age comp, 0/1-inflated logistic-normal: ",
            c("Declining probablity of 0 parameter 1",
              "Declining probablity of 0 parameter 2",
              "logistic-normal $\\sigma$")))
          fe.vals = c(fe.vals, pars[i,1:2], exp(pars[i,3]))
          fe.cis = rbind(fe.cis, 
            ci(pars[i,1], pars_sd[i,1]),
            ci(pars[i,2], pars_sd[i,2]),
            ci(pars[i,3], pars_sd[i,3], type = "exp"))
        }
        if(age_comp_models[i] == 9){
          fe.names = c(fe.names, paste0(startnames[i], " age comp, 0/1-inflated logistic-normal: ",
            c("Binomial N parameter probablity of 0",
              "logistic-normal $\\sigma$")))
          fe.vals = c(fe.vals, exp(pars[i,1:2]))
          fe.cis = rbind(fe.cis, 
            ci(pars[i,1], pars_sd[i,1], type = "exp"),
            ci(pars[i,2], pars_sd[i,2], type = "exp"))
        }
        if(age_comp_models[i] == 10){
          fe.names = c(fe.names, paste0(startnames[i], " age comp, MV Tweedie: ",
            c("$\\phi$", "power")))
          fe.vals = c(fe.vals, exp(pars[i,1:2]))
          fe.cis = rbind(fe.cis, 
            ci(pars[i,1], pars_sd[i,1], type = "exp"),
            1 + ci(pars[i,2], pars_sd[i,2], lo = 0, hi = 1, type = "expit")) # 1 < x < 2
        }
      }
    }
    return(list(fe.names, fe.vals, fe.cis))
  }
  temp = add_age_comp_pars(data$age_comp_model_fleets, data$use_catch_paa, pars$catch_paa_pars, sd$catch_paa_pars, 
      is_fleet = TRUE, fe.names, fe.vals, fe.cis)
  fe.names = temp[[1]]
  fe.vals = temp[[2]]
  fe.cis = temp[[3]]
  temp = add_age_comp_pars(data$age_comp_model_indices, data$use_index_paa, pars$index_paa_pars, sd$index_paa_pars, 
      is_fleet = FALSE, fe.names, fe.vals, fe.cis)
  fe.names = temp[[1]]
  fe.vals = temp[[2]]
  fe.cis = temp[[3]]


  if(sum(!is.na(mod$input$map$Mpars))){ #any (mean) M estimated?
    for(s in 1:data$n_stocks) for(r in 1:data$n_regions) for(a in 1:data$n_ages){
      if(!is.na(sd$Mpars[s,r,a])){
        if(data$M_re_model[s,r] == 1) {# & data$M_model == 1){ #no random effects, ecov or WAA effects on M
          if(data$M_model == 1) {
            if(sum(data$Ecov_how_M[,s,a,r]) == 0) modify = "M for age "
            if(sum(data$Ecov_how_M[,s,a,r]) >0) modify = "mean log(M) for age "
          }
          if(data$M_model == 2) modify = "mean log(M) intercept for log(WAA) effects" #whether Ecov or not
        }
        if(data$M_re_model[s,r] >1){
          if(data$M_model == 1) modify = "mean log(M) intercept for age " #whether Ecov or not
        }
        modify <- paste(stock.names.tab[s], region.names.tab[r], modify)
        if(data$M_model == 1) modify <- paste0(modify, mod$ages.lab[a])
        fe.names = c(fe.names, modify)
        if(data$M_re_model[s,r]>1 | data$M_model == 2 | sum(data$Ecov_how_M[,s,a,r]) >0){
          fe.vals = c(fe.vals, pars$Mpars[s,r,a])
          fe.cis = rbind(fe.cis, ci(pars$Mpars[s,r,a], sd$Mpars[s,r,a]))
        } else {
          fe.vals = c(fe.vals, exp(pars$Mpars[s,r,a]))
          fe.cis = rbind(fe.cis, ci(pars$Mpars[s,r,a], sd$Mpars[s,r,a], type = "exp"))
        }
      }
    }
  }
  for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
    modify <- paste(stock.names.tab[s], region.names.tab[r])
    if(!is.na(sd$log_b[s,r])){
      fe.names = c(fe.names, paste0(modify, "log(M) slope for log(WAA) effect"))
      fe.vals = c(fe.vals, exp(pars$log_b[s,r]))
      fe.cis = rbind(fe.cis, ci(pars$log_b[s,r], sd$log_b[s,r], type = "exp"))
    }
  }
  for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
    modify <- paste(stock.names.tab[s], region.names.tab[r])
    if(!is.na(sd$M_repars[s,r,1])) {
      fe.names = c(fe.names, paste0(modify, "M RE $\\sigma$"))
      fe.vals = c(fe.vals, exp(pars$M_repars[s,r,1]))
      fe.cis = rbind(fe.cis, ci(pars$M_repars[s,r,1], sd$M_repars[s,r,1], type = "exp"))
    }
    if(!is.na(sd$M_repars[s,r,2])) {
    #if(data$M_re_model %in% c(3,5)){
      fe.names = c(fe.names, paste0(modify, "M RE AR1 $\\rho$ (age)"))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$M_repars[s,r,2])))
      fe.cis = rbind(fe.cis, ci(pars$M_repars[s,r,2], sd$M_repars[s,r,2], lo = -1, hi = 1, type = "expit"))
    }
    if(!is.na(sd$M_repars[s,r,3])) {
    #if(data$M_re_model %in% c(4,5)) {
      fe.names = c(fe.names, paste0(modify, "M RE AR1 $\\rho$ (year)"))
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$M_repars[s,r,3])))
      fe.cis = rbind(fe.cis, ci(pars$M_repars[s,r,3], sd$M_repars[s,r,3], lo = -1, hi = 1, type = "expit"))
    }
  }

  #movement fixed effects or prior-based RE
  if(data$n_regions>1) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) {
    k <- rr
    if(rr>=r) k <- k + 1
    if(data$mu_model[r,rr] %in% 1:4) modify <- paste("$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
    if(data$mu_model[r,rr] %in% 5:8) modify <- paste("stock", stock.names.tab, "$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
    if(data$mu_model[r,rr] %in% 9:12) modify <- paste("season", 1:data$n_seasons, "$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
    if(data$mu_model[r,rr] %in% 13:16) modify <- paste(rep(stock.names.tab, each = data$n_seasons), "season", 1:data$n_seasons, 
      "$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
    if(data$mu_model[r,rr] %in% c(1:4,9:12)) ns <- 1
    else ns <- data$n_stocks
    if(data$mu_model[r,rr] %in% c(1:4,5:8)) nt <- 1
    else nt <- data$n_seasons
    modify <- matrix(modify, nt, ns)
    for(s in 1:ns) for(t in 1:nt) {
      if(data$use_mu_prior[s,t,r,rr]) parname <- "mu_prior_re"
      else parname <- "trans_mu"
      
      if(!is.na(sd[[parname]][s,t,r,rr])) {
        fe.names = c(fe.names, paste0(modify[t,s], " (intercept)"))
        if(data$mig_type[s] == 1){ #instantaneous: log transform
          fe.vals = c(fe.vals, exp(pars[[parname]][s,t,r,rr]))
          fe.cis = rbind(fe.cis, ci(pars[[parname]][s,t,r,rr], sd[[parname]][s,t,r,rr], type = "exp"))
        } else { #sequential: logit transform
          fe.vals = c(fe.vals, 1/(1 + exp(- pars[[parname]][s,t,r,rr])))
          fe.cis = rbind(fe.cis, ci(pars[[parname]][s,t,r,rr], sd[[parname]][s,t,r,rr], lo = 0, hi = 1, type = "expit"))
        }
      }
    }
  }

  #movement random effects variance, correlation
  if(data$n_regions>1) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) {
    k <- rr
    if(rr>=r) k <- k + 1
    if(data$mu_model[r,rr] %in% c(2:4,6:8,10:12,14:16)){ #some movement random effects
      if(data$mu_model[r,rr] %in% 2:4) modify <- paste("$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
      if(data$mu_model[r,rr] %in% 6:8) modify <- paste("stock", stock.names.tab, "$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
      if(data$mu_model[r,rr] %in% 10:12) modify <- paste("season", 1:data$n_seasons, "$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
      if(data$mu_model[r,rr] %in% 14:16) modify <- paste(rep(stock.names.tab, each = data$n_seasons), "season", 1:data$n_seasons, 
        "$\\mu$ from",  region.names.tab[r], "to", region.names.tab[k])
      if(data$mu_model[r,rr] %in% c(2:4,10:12)) ns <- 1
      else ns <- data$n_stocks
      if(data$mu_model[r,rr] %in% c(2:4,6:8)) nt <- 1
      else nt <- data$n_seasons
      modify <- matrix(modify, nt, ns)
      for(s in 1:ns) for(t in 1:nt) {
        if(!is.na(sd$mu_repars[s,t,r,rr,1])) {
          fe.names = c(fe.names, paste0(modify[t,s], " RE $\\sigma$"))
          fe.vals = c(fe.vals, exp(pars$mu_repars[s,t,r,rr,1]))
          fe.cis = rbind(fe.cis, ci(pars$mu_repars[s,t,r,rr,1], sd$mu_repars[s,t,r,rr,1], type = "exp"))
        }
        if(!is.na(sd$mu_repars[s,t,r,rr,2])) {
          fe.names = c(fe.names, paste0(modify[t,s], " RE AR1 $\\rho$ (age)"))
          fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$mu_repars[s,t,r,rr,2])))
          fe.cis = rbind(fe.cis, ci(pars$mu_repars[s,t,r,rr,2], sd$mu_repars[s,t,r,rr,2], lo = -1, hi = 1, type = "expit"))
        }
        if(!is.na(sd$mu_repars[s,t,r,rr,3])) {
          fe.names = c(fe.names, paste0(modify[t,s], " RE AR1 $\\rho$ (year)"))
          fe.vals = c(fe.vals, -1 + 2/(1 + exp(- pars$mu_repars[s,t,r,rr,3])))
          fe.cis = rbind(fe.cis, ci(pars$mu_repars[s,t,r,rr,3], sd$mu_repars[s,t,r,rr,3], lo = -1, hi = 1, type = "expit"))
        }
      }
    }
  }
  for(i in 1:data$n_fleets){
    if(!is.na(sd$log_catch_sig_scale[i])){
      fe.names = c(fe.names, paste0(fleet.names.tab[i], " log-catch observation SD scalar"))
      fe.vals = c(fe.vals, exp(pars$log_catch_sig_scale[i]))
      fe.cis = rbind(fe.cis, ci(pars$log_catch_sig_scale[i], sd$log_catch_sig_scale[i], type = "exp"))
    }
  }
  for(i in 1:data$n_indices){
    if(!is.na(sd$log_index_sig_scale[i])){
      fe.names = c(fe.names, paste0(index.names.tab[i], " log-index observation SD scalar"))
      fe.vals = c(fe.vals, exp(pars$log_index_sig_scale[i]))
      fe.cis = rbind(fe.cis, ci(pars$log_index_sig_scale[i], sd$log_index_sig_scale[i], type = "exp"))
    }
  }

  ecov.labels <- paste0("Ecov ", ecov.names.tab, ": ")
  for(i in 1:data$n_Ecov){
    #if any parameters sds are not NA, then Ecov_model != 0
    if(!is.na(sd$Ecov_process_pars[1,i])){
      fe.vals = c(fe.vals, pars$Ecov_process_pars[1,i])
      fe.cis = rbind(fe.cis, ci(pars$Ecov_process_pars[1,i], sd$Ecov_process_pars[1,i]))
      if(data$Ecov_model[i] == 1) fe.names = c(fe.names, paste0(ecov.labels[i], "RW Ecov$_1$"))
      if(data$Ecov_model[i] == 2) fe.names = c(fe.names, paste0(ecov.labels[i], "$\\mu$"))
    }
    if(!is.na(sd$Ecov_process_pars[2,i])){
      fe.vals = c(fe.vals, exp(pars$Ecov_process_pars[2,i]))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_process_pars[2,i], sd$Ecov_process_pars[2,i], type = "exp"))
      if(data$Ecov_model[i] == 1) fe.names = c(fe.names, paste0(ecov.labels[i], "RW $\\sigma$"))
      if(data$Ecov_model[i] == 2) fe.names = c(fe.names, paste0(ecov.labels[i], "$\\sigma$"))
    }
    if(!is.na(sd$Ecov_process_pars[3,i])){
      fe.vals = c(fe.vals, -1 + 2/(1 + exp(-pars$Ecov_process_pars[3,i])))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_process_pars[3,i], sd$Ecov_process_pars[3,i], lo = -1, hi = 1, type = "expit")) #doesn't use rho_trans
      fe.names = c(fe.names, paste0(ecov.labels[i], "AR1 $\\rho$"))
    }
  }
  ecov_beta_names = paste0("Ecov_beta_", c("R", "M", "mu","q"))
  for(et in 1:4){
  #   ecov_beta_map = array(mod$input$map[[ecov_beta_names[i]]], dim = dim(pars[[ecov_beta_names[i]]])) 
  #   # R: n_stocks x n_ecov x n_poly
  #   # M: n_stocks x n_ages x n_regions x n_ecov x n_poly
  #   #mu: n_stocks x n_ages x n_seasons x n_regions x n_regions-1 x n_ecov x n_poly
  #   #q : n_indices x n_ecov x n_poly
  # }
    for(i in 1:data$n_Ecov){
      if(any(!is.na(sd[[ecov_beta_names[et]]]))){
        if(et<4) for(s in 1:data$n_stocks){
          if(et==1) if(any(!is.na(sd[[ecov_beta_names[et]]][s,i,]))){ #Recruitment
            ind <- which(!is.na(sd[[ecov_beta_names[et]]][s,i,]))
            for(k in ind) {
              fe.vals = c(fe.vals, pars[[ecov_beta_names[et]]][s,i,k])
              fe.cis = rbind(fe.cis, ci(pars[[ecov_beta_names[et]]][s,i,k], sd[[ecov_beta_names[et]]][s,i,k]))
              fe.names = c(fe.names, paste0(stock.names.tab[s], " Recruitment Ecov: ", ecov.names.tab[i], " $\\beta_", k, "$"))
            }
          }
          if(et==2) for(a in 1:data$n_ages) for(r in 1:data$n_regions) if(any(!is.na(sd[[ecov_beta_names[et]]][s,a,r,i,]))){ #M
            ind <- which(!is.na(sd[[ecov_beta_names[et]]][s,a,r,i,]))
            modify <- paste(stock.names.tab[s], region.names.tab[r], "M at age", mod$ages.lab[a])
            for(k in ind) {
              fe.vals = c(fe.vals, pars[[ecov_beta_names[et]]][s,a,r,i,k])
              fe.cis = rbind(fe.cis, ci(pars[[ecov_beta_names[et]]][s,a,r,i,k], sd[[ecov_beta_names[et]]][s,a,r,i,k]))
              fe.names = c(fe.names, paste0(modify, " Ecov: ", ecov.names.tab[i], " $\\beta_", k, "$"))
            }
          }
          if(et==3) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)){#movement
            if(rr != r) if(any(!is.na(sd[[ecov_beta_names[et]]][s,a,t,r,,rr,i,]))){ 
              ind <- which(!is.na(sd[[ecov_beta_names[et]]][s,a,t,r,,rr,i,]))
              from <- region.names.tab[r]
              to <- rr
              if(rr >= r) to <- rr+1
              to <- region.names.tab[to]
              modify <- paste(stock.names.tab[s], " age ", mod$ages.lab[a], "season", t, "move from", from, "to", to)
              for(k in ind) {
                fe.vals = c(fe.vals, pars[[ecov_beta_names[et]]][s,a,t,r,,rr,i,k])
                fe.cis = rbind(fe.cis, ci(pars[[ecov_beta_names[et]]][s,a,t,r,,rr,i,k], sd[[ecov_beta_names[et]]][s,a,t,r,,rr,i,k]))
                fe.names = c(fe.names, paste0(modify, " Ecov: ", ecov.names.tab[i], " $\\beta_",k, "$"))
              }
            }
          }
        }
        if(et == 4) for(j in 1:data$n_indices){ #q
          if(any(!is.na(sd[[ecov_beta_names[et]]][j,i,]))){
            ind <- which(!is.na(sd[[ecov_beta_names[et]]][j,i,]))
            for(k in ind) {
              fe.vals = c(fe.vals, pars[[ecov_beta_names[et]]][j,i,k])
              fe.cis = rbind(fe.cis, ci(pars[[ecov_beta_names[et]]][j,i,k], sd[[ecov_beta_names[et]]][j,i,k]))
              fe.names = c(fe.names, paste0(index.names.tab[j], " Catchability Ecov: ", ecov.names.tab[i], " $\\beta_", k, "$"))
            }
          }
        }
      }
    }
  }

  for(i in 1:data$n_Ecov){
    if(data$Ecov_obs_sigma_opt[i] == 2){ #single ecov obs sd estimated
      fe.names = c(fe.names, paste0("Ecov: ", ecov.names.tab[i], " obs. sd."))
      ind = which(!is.na(matrix(mod$input$map$Ecov_obs_logsigma, NROW(mod$input$par$Ecov_obs_logsigma))[,i]))[1]
      fe.vals = c(fe.vals, exp(pars$Ecov_obs_logsigma[ind,i]))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_obs_logsigma[ind,i], sd$Ecov_obs_logsigma[ind,i], type = "exp"))
    }
    if(data$Ecov_obs_sigma_opt[i] == 4){
      fe.names = c(fe.names, paste0("Ecov: ", ecov.names.tab[i], " obs. log(sd.) RE ", c("$\\mu$", "$\\sigma$")))
      fe.vals = c(fe.vals, pars$Ecov_obs_sigma_pars[1,i])
      fe.cis = rbind(fe.cis, ci(pars$Ecov_obs_sigma_pars[1,i], sd$Ecov_obs_sigma_pars[1,i]))
      fe.vals = c(fe.vals, exp(pars$Ecov_obs_sigma_par[2,i]))
      fe.cis = rbind(fe.cis, ci(pars$Ecov_obs_sigma_par[2,i], sd$Ecov_obs_sigma_par[2,i], type = "exp"))
    }
  }

  fe = cbind(fe.vals, fe.cis)
  rownames(fe) = fe.names
  colnames(fe) = c("Estimate", "Std. Error", "95\\% CI lower", "95\\% CI upper")
  if(!is.null(od)) {
    saveRDS(fe, file = file.path(od,"parameter_estimates_table.RDS"))
    saveRDS(mod$input, file = file.path(od,"fit_input.RDS"))
  }
  
  if(mod$is_sdrep) {
    sdrep = mod$sdrep
    pars = TMB:::as.list.sdreport(sdrep, "Est", report = TRUE)
    sd = TMB:::as.list.sdreport(sdrep, "Std", report = TRUE)
  }
  #Jan 1 numbers at age, by stock and region
  for(s in 1:data$n_stocks) for(r in 1:data$n_regions){
    NAA = NAA.cv = mod$rep$NAA[s,r,,]
    rownames(NAA) =mod$years_full
    colnames(NAA) = mod$ages.lab
    if(!is.null(od)) saveRDS(NAA, file = file.path(od,paste0(stock.names.f[s],"_", region.names.f[r], "_NAA_table.RDS")))
    NAA.cv[] <- NA
    if(!is.null(mod$opt)) if(!is.na(mod$na_sdrep)) if(mod$is_sdrep) {
      NAA.cv[] = sd[["log_NAA_rep"]][s,r,,]
      NAA.sd = NAA * NAA.cv
      NAA.lo = exp(log(NAA) - qnorm(0.975) * NAA.cv)
      NAA.hi = exp(log(NAA) + qnorm(0.975) * NAA.cv)
      rownames(NAA.sd) = rownames(NAA.cv) = rownames(NAA.lo) = rownames(NAA.hi) = mod$years_full
      colnames(NAA.sd) = colnames(NAA.cv) = colnames(NAA.lo) = colnames(NAA.hi) = mod$ages.lab
      if(!is.null(od)){
        saveRDS(NAA.sd, file = file.path(od,paste0(stock.names.f[s],"_", region.names.f[r], "_NAA_sd_table.RDS")))
        saveRDS(NAA.lo, file = file.path(od,paste0(stock.names.f[s],"_", region.names.f[r], "_NAA_lo_table.RDS")))
        saveRDS(NAA.hi, file = file.path(od,paste0(stock.names.f[s],"_", region.names.f[r], "_NAA_hi_table.RDS")))
      }
    }
  }
      
  #Total F at age by region
  for(r in 1:data$n_regions){
    FAA_r <- FAA_r.cv <- mod$rep$FAA_by_region[r,,]
    rownames(FAA_r) = mod$years_full
    colnames(FAA_r) = mod$ages.lab
    if(!is.null(od)) saveRDS(FAA_r, file = file.path(od,paste0(region.names.f[r], "_FAA_tot_table.RDS")))
    FAA_r.cv[] <- NA
    if(!is.null(mod$opt)) if(!is.na(mod$na_sdrep)) if(mod$is_sdrep) {
      FAA_r.cv = sd[["log_FAA_by_region"]][r,,]
      FAA_r.sd = FAA_r * FAA_r.cv
      FAA_r.lo = FAA_r * exp(- qnorm(0.975) * FAA_r.cv)
      FAA_r.hi = FAA_r * exp( qnorm(0.975) * FAA_r.cv)
      rownames(FAA_r.sd) = rownames(FAA_r.cv) = rownames(FAA_r.lo) = rownames(FAA_r.hi) = mod$years_full
      colnames(FAA_r.sd) = colnames(FAA_r.cv) = colnames(FAA_r.lo) = colnames(FAA_r.hi) = mod$ages.lab
      if(!is.null(od)){ 
        saveRDS(FAA_r.sd, file = file.path(od,paste0(region.names.f[r], "_FAA_tot_sd_table.RDS")))
        saveRDS(FAA_r.lo, file = file.path(od,paste0(region.names.f[r], "_FAA_tot_lo_table.RDS")))
        saveRDS(FAA_r.hi, file = file.path(od,paste0(region.names.f[r], "_FAA_tot_hi_table.RDS")))
      }
    }
  }

  #F at age by fleet
  for(f in 1:data$n_fleets){
    FAA = FAA.cv = mod$rep$FAA[f,,]
    dimnames(FAA) = list(mod$years_full, mod$ages.lab)
    if(!is.null(od)) saveRDS(FAA, file = file.path(od,paste0(fleet.names.f[f], "_FAA_table.RDS")))
    FAA.cv[] <- NA
    if(!is.null(mod$opt)) if(!is.na(mod$na_sdrep)) if(mod$is_sdrep) {
      FAA.cv[] = sd[["log_FAA"]][f,,]
      FAA.sd = FAA * FAA.cv
      FAA.lo = exp(mod$rep$FAA[f,,] - qnorm(0.975) * FAA.cv)
      FAA.hi = exp(mod$rep$FAA[f,,] + qnorm(0.975) * FAA.cv)
      dimnames(FAA.sd) = dimnames(FAA.cv) = dimnames(FAA.lo) = dimnames(FAA.hi) = list(mod$years_full, mod$ages.lab)
      if(!is.null(od)){ 
        saveRDS(FAA.sd, file = file.path(od,paste0(fleet.names.f[f], "_FAA_sd_table.RDS")))
        saveRDS(FAA.lo, file = file.path(od,paste0(fleet.names.f[f], "_FAA_lo_table.RDS")))
        saveRDS(FAA.hi, file = file.path(od,paste0(fleet.names.f[f], "_FAA_hi_table.RDS")))
      }
    }
  }
  
  wham.dir <- find.package("wham")
  pt = list.files(find.package("wham"), pattern = "par_tables.Rmd", recursive = T, full.names = T)[1]
  if(!is.null(od)){ 
    od <- normalizePath(od) #some latex code doesn't like tildes in abreviated directories
    file.copy(from=pt, to=od, overwrite=TRUE)
    #print(dir(od))
    
    if(do.html) rmarkdown::render(file.path(od,"par_tables.Rmd"), output_format = "html_document", output_file = file.path(od, "wham_par_tables.html"), 
      quiet = T, envir = new.env())
    #if(do.tex) rmarkdown::render(file.path(od,"par_tables.Rmd"), output_format = "pdf_document", output_file = file.path(od,"wham_par_tables.pdf"), quiet = T)
    if(do.tex) { #for some reason on windows working outside of the temp directory was causing issues for tinytex::latexmf.
      origdir = getwd()
      setwd(od)
      rmarkdown::render("par_tables.Rmd", output_format = "pdf_document", output_file = file.path(od, "wham_par_tables.pdf"), quiet = T, envir = new.env())
      setwd(origdir)
    }
  }

  #delete par_tables.Rmd from od
  #file.remove(file.path(od,"par_tables.Rmd"))

}
