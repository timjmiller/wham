
fit.wham.fn = function(input, version="wham_v0", n.newton = 3, do.sdrep = TRUE, do.retro = TRUE, n.peels = 7)
{
  mod <- TMB::MakeADFun(input$data,input$par, DLL=version, random = input$random, map = input$map)
  mod = fit.tmb.fn(mod, n.newton = n.newton, do.sdrep = do.sdrep)
  if(do.retro) mod$peels = retro.fn(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels)
  mod$years = input$years
  mod$ages.lab = input$ages.lab
  mod$model_name = input$model_name
  return(mod)
}

fit.tmb.fn = function(model, n.newton=3, do.sdrep = TRUE)
{
  model$opt <- nlminb(model$par, model$fn,model$gr, control = list(iter.max = 1000, eval.max = 1000))
  if(n.newton) for(i in 1:n.newton) { # Take a few extra newton steps
    g <- as.numeric(model$gr(model$opt$par))
    h <- optimHess(model$opt$par, model$fn, model$gr)
    model$opt$par <- model$opt$par - solve(h, g)
    model$opt$objective <- model$fn(model$opt$par)
  }
  model$date = Sys.time()
  model$dir = getwd()
  model$rep <- model$report()
  model$TMB_version = packageVersion("TMB")
  model$parList = model$env$parList()
  model$final_gradient = model$gr()
  if(do.sdrep)
  {
    model$sdrep <- try(sdreport(model))
    model$is_sdrep = !is.character(model$sdrep)
    if(model$is_sdrep) model$na_sdrep = any(is.na(summary(model$sdrep)[,2]))
    else model$na_sdrep = NA
  }
  return(model)
}

peel.fit.fn = function(peel, model, do.sdrep = FALSE, version="wham_v0", n.newton = 3)
{
  out = list()
  print(peel)
  temp = model
  n_years = temp$dat$n_years_catch = temp$dat$n_years_indices = temp$dat$n_years_model = temp$dat$n_years_model - peel
  log_NAA_na_ind = rbind(matrix(1:(temp$dat$n_ages*(n_years-1)), n_years-1), matrix(rep(NA, peel*temp$dat$n_ages), peel))
  F_devs_na_ind = rbind(matrix(1:(temp$dat$n_fleets * (n_years-1)), n_years-1), matrix(rep(NA, peel * temp$dat$n_fleets), peel))
  log_R_na_ind = c(1:(n_years-1), rep(NA, peel))
  if("log_R" %in% model$random) temp$map$log_R = factor(log_R_na_ind)
  if("log_NAA" %in% model$random) temp$map$log_NAA = factor(log_NAA_na_ind)
  temp$map$F_devs = factor(F_devs_na_ind)
  temp.mod <- TMB::MakeADFun(temp$dat,temp$par,DLL=version, random = temp$random, map = temp$map)
  out = fit.tmb.fn(temp.mod, do.sdrep = do.sdrep, n.newton = n.newton)
  return(out)
}

retro.fn = function(model, n.peels = 7, ran = "log_NAA", do.sdrep = FALSE, n.newton = 0)
{
  temp = list(dat = model$env$data, par = model$parList, map = model$env$map, random = ran)
  if(n.peels>0) peels = list(peel.fit.fn(1, model = temp, do.sdrep = do.sdrep, n.newton = n.newton))
  if(n.peels>1) for(i in 2:n.peels) peels[[i]] = peel.fit.fn(i, model = temp, do.sdrep = do.sdrep, n.newton = n.newton)
  return(peels)
}

rho.fn = function(model)
{
  npeels = length(model$peels)
  ny = model$env$data$n_years_model
  na = model$env$data$n_ages
  if(npeels)
  {
    rho = c(
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[ny-x]/model$rep$SSB[ny-x] - 1)),
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$Fbar[ny-x]/model$rep$Fbar[ny-x] - 1)))#,
      #mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[ny-x,1]/model$rep$NAA[ny-x,1] - 1)))
    names(rho) = c("SSB","Fbar")#,"R")
    rho.naa = sapply(1:na, function(y)
    {
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[ny-x,y]/model$rep$NAA[ny-x,y] - 1))
    })
    names(rho.naa) = c("R", paste0("N", model$ages.lab[2:na]))
    rho = c(rho, rho.naa)
    return(rho)
  }
  else stop("There are no peels in this model")
}

retro.res.fn = function(model) #get time series for retro plots
{
  npeels = length(model$peels)
  ny = model$env$data$n_years_model
  if(npeels)
  {
    SSB = lapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[1:(ny-x)]/model$rep$SSB[1:(ny-x)] - 1)
    Fbar = lapply(1:npeels, function(x) model$peels[[x]]$rep$Fbar[1:(ny-x)]/model$rep$Fbar[1:(ny-x)] - 1)
    na = model$env$data$n_ages
    NAA = lapply(1:na, function(a)
    {
      lapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[1:(ny-x),a]/model$rep$NAA[1:(ny-x),a] - 1)
    })
    res = list(SSB = SSB, Fbar = Fbar, NAA = NAA)
    return(res)
  }
  else stop("There are no peels in this model")
}
