plot.ecov <- function(mod, plot.pad = FALSE, do.tex=FALSE, do.png=FALSE, fontfam="", res=72, od){
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  ecov.pred = mod$rep$Ecov_x
  ecov.pred.low <- ecov.pred.high <- ecov.pred.se <- matrix(NA, NROW(ecov.pred), NCOL(ecov.pred))
  ecov.obs = dat$Ecov_obs[1:dat$n_years_Ecov,,drop=F]
  years <- seq(from=mod$input$years_Ecov[1], by=1, length.out=NROW(ecov.obs))
  years_full <- seq(from=mod$input$years_Ecov[1], by=1, length.out=NROW(ecov.pred))#dat$n_years_Ecov+dat$n_years_proj_Ecov)

  # ecov.obs.sig = mod$rep$Ecov_obs_sigma # Ecov_obs_sigma now a derived quantity in sdrep
  # if(class(mod$sdrep)[1] == "sdreport"){
  #   sdrep = summary(mod$sdrep)
  # } else {
  #   sdrep = mod$sdrep
  # }
  if(class(mod$sdrep)[1] == "sdreport"){
    temp <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = T),
      TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = T))
    ecov.pred.se[] <- temp[[2]]$Ecov_x #TMB:::as.list.sdreport(mod$sdrep, what = "Std.", report=TRUE)$Ecov_x
    ecov.pred.low[] <- ecov.pred - 1.96 * ecov.pred.se
    ecov.pred.high[] <- ecov.pred + 1.96 * ecov.pred.se
  }

  ecov.obs.sig = mod$rep$Ecov_obs_sigma # Ecov_obs_sigma is filled with fixed, or estimated values (fe or re) for each covariate depending on the respective options
  # if("Ecov_obs_logsigma" %in% names(mod$env$par)){
  #   ecov.obs.sig = matrix(exp(sdrep[rownames(sdrep) %in% "Ecov_obs_logsigma",1]), ncol=dat$n_Ecov) # all the same bc obs_sig_var --> 0
  #   if(dim(ecov.obs.sig)[1] == 1) ecov.obs.sig = matrix(rep(ecov.obs.sig, dim(ecov.obs)[1]), ncol=dat$n_Ecov, byrow=T)
  # } else {
  #   ecov.obs.sig = exp(mod$input$par$Ecov_obs_logsigma)
  # }
  ecov.use = dat$Ecov_use_obs[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig = ecov.obs.sig[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig[ecov.use == 0] <- NA
  #ecov.pred.se = matrix(sdrep[rownames(sdrep) %in% "Ecov_x",2], ncol=dat$n_Ecov)

  # default: don't plot the padded entries that weren't used in ecov likelihood
  if(!plot.pad) ecov.obs[ecov.use == 0] <- NA

  ecov.res = (ecov.obs - ecov.pred[1:dat$n_years_Ecov,]) / ecov.obs.sig # standard residual (obs - pred)

  ecovs <- 1:dat$n_Ecov
  plot.colors = mypalette(dat$n_Ecov)
  for (i in ecovs)
  {
    if(do.tex) cairo_pdf(file.path(od, paste0("Ecov_",i,"_",mod$input$Ecov_names[i],".pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Ecov_",i, "_", mod$input$Ecov_names[i],'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)

    # ecov.pred.low <- ecov.pred[,i] - 1.96 * ecov.pred.se[,i]
    # ecov.pred.high <- ecov.pred[,i] + 1.96 * ecov.pred.se[,i]
    ecov.low <- ecov.obs[,i] - 1.96 * ecov.obs.sig[,i]
    ecov.high <- ecov.obs[,i] + 1.96 * ecov.obs.sig[,i]
    y.min <- ifelse(min(ecov.low,na.rm=T) < 0, 1.1*min(ecov.low,na.rm=T), 0.9*min(ecov.low,na.rm=T))
    y.max <- ifelse(max(ecov.high,na.rm=T) < 0, 0.9*max(ecov.high,na.rm=T), 1.1*max(ecov.high,na.rm=T))
    if(max(ecov.pred[,i],na.rm=T) > y.max) y.max <- max(ecov.pred[,i],na.rm=T)
    if(min(ecov.pred[,i],na.rm=T) < y.min) y.min <- min(ecov.pred[,i],na.rm=T)
    plot(years_full, ecov.pred[,i], type='n', xlab="Year", ylab=mod$input$Ecov_names[i],
         ylim=c(y.min, y.max))
    polygon(c(years_full,rev(years_full)), c(ecov.pred.low[,i], rev(ecov.pred.high[,i])), col=adjustcolor(plot.colors[i], alpha.f=0.4), border = "transparent")
    arrows(years, ecov.low, years, ecov.high, length=0)
    points(years, ecov.obs[,i], pch=19)
    lines(years_full, ecov.pred[,i], col=plot.colors[i], lwd=3)
    title (paste0("Ecov ",i, ": ",mod$input$Ecov_names[i]), outer=T, line=-1)
    if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2)

    if(do.tex | do.png) dev.off() else par(origpar)
  }
  #plot standard deviations if estimated as random effects
  ecov.obs.sigma.cv = as.list(mod$sdrep, "Std")$Ecov_obs_logsigma_re
  for (i in ecovs) if(dat$Ecov_obs_sigma_opt[i]==4)
  {
    if(do.tex) cairo_pdf(file.path(od, paste0("Ecov_",i, "_", mod$input$Ecov_names[i],"_sig_re.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Ecov_",i, "_", mod$input$Ecov_names[i], '_sig_re.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    sigmas = ecov.obs.sig[,i]
    sigmas.cv = ecov.obs.sigma.cv[,i]
    sigmas.low <- exp(log(sigmas) - 1.96 * sigmas.cv)
    sigmas.high <- exp(log(sigmas) + 1.96 * sigmas.cv)
    y.min = min(sigmas.low, na.rm =T)
    y.max = max(sigmas.high, na.rm =T)
    plot(years, sigmas, type='n', xlab="Year", ylab=mod$input$Ecov_names[i],
         ylim=c(y.min, y.max))
    polygon(c(years,rev(years)), c(sigmas.low, rev(sigmas.high)), col=adjustcolor(plot.colors[i], alpha.f=0.4), border = "transparent")
    points(years, sigmas, pch=19)
    lines(years, sigmas, col=plot.colors[i], lwd=3)
    title (paste0("Ecov ",i, ": ",mod$input$Ecov_names[i], " SD"), outer=T, line=-1)
    if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2)

    if(do.tex | do.png) dev.off() else par(origpar)
  }
}

#revised

plot.osa.residuals <- function(mod, do.tex=FALSE, do.png=FALSE, fontfam="", res=72, od){
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  if("logcatch" %in% mod$osa$type){
    dat <- subset(mod$osa, type == "logcatch")
    dat$fleet <- sapply(dat$fleet, function(x) as.integer(strsplit(x, "_")[[1]][2]))
    dat$fleet <- factor(dat$fleet) #have to make sure integers are in order
    n.fleets <- length(levels(dat$fleet))
    plot.colors = mypalette(n.fleets)
    ffns <- chartr(" ", "_", mod$input$fleet_names)
    for(f in 1:n.fleets){
      tmp <- subset(dat, fleet==as.character(f))
      if(do.tex) cairo_pdf(file.path(od, paste0("OSA_resid_catch_4panel_", ffns[f],".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSA_resid_catch_4panel_", ffns[f],'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = fontfam)
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(dat$residual, na.rm=TRUE)))
      if(is.infinite(ylim.max)) {
        cat("Infinite osa residuals for aggregate catch in fleet ", f, ", so using +/-10 for range of y axis \n")
        ylim.max = 10
      }
      ylims <- c(-ylim.max, ylim.max)

      # 1. trend vs. year
      plot(years, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 2. trend vs. fitted val
      plot(mod$rep$pred_log_catch[1:length(mod$years),f], tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Log(Predicted Catch)", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 3. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(tmp$residual, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)

      # 4. QQ plot modified from car:::qqPlot.default
      notNA <- tmp$residual[!is.na(tmp$residual)]
      ord.x <- notNA[order(notNA)]
      n <- length(ord.x)
      P <- ppoints(n)
      z <- qnorm(P, mean=0, sd=1)
      plot(z, ord.x, xlab="Std Normal Quantiles", ylab="OSA Residual Quantiles", main="", type = "n")
      grid(lty = 1, equilogs = FALSE)
      box()
      points(z, ord.x, col=plot.colors[f], pch=19)
      abline(0,1, col=plot.colors[f])
      conf = 0.95
      zz <- qnorm(1 - (1 - conf)/2)
      SE <- (1/dnorm(z)) * sqrt(P * (1 - P)/n)
      upper <- z + zz * SE
      lower <- z - zz * SE
      lines(z, upper, lty=2, col=plot.colors[f])
      lines(z, lower, lty=2, col=plot.colors[f])

      title (paste0("OSA residual diagnostics: ", mod$input$fleet_names[f]), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }
  if("catchpaa" %in% mod$osa$type) {
    dat <- subset(mod$osa, type == "catchpaa")
    dat$fleet <- sapply(dat$fleet, function(x) as.integer(strsplit(x, "_")[[1]][2]))
    dat$fleet <- factor(dat$fleet) #have to make sure integers are in order
    n.fleets <- length(levels(dat$fleet))
    plot.colors = mypalette(n.fleets)
    ffns <- chartr(" ", "_", mod$input$fleet_names)
    for(f in 1:n.fleets) if(any(mod$input$data$use_catch_paa[,f]==1)){
      if(do.tex) cairo_pdf(file.path(od, paste0("OSA_resid_paa_6panel_", ffns[f],".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSA_resid_paa_6panel_", ffns[f],'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = fontfam)
      
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(3,2))
      yind = which(mod$input$data$use_catch_paa[,f] ==1)
      resids = vals = ages = pyears = cohorts = matrix(NA, nrow = mod$input$data$n_years_model, 
        ncol = mod$input$data$n_ages)
      for(j in yind){
        tmp = subset(dat, year == j & fleet == as.character(f))
        resids[j,tmp$age] = tmp$residual
        vals[j,tmp$age] = tmp$val
        ages[j,tmp$age] = tmp$age
        pyears[j,tmp$age] = mod$years[tmp$year]
        cohorts[j,tmp$age] = tmp$cohort
      }

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(resids, na.rm=TRUE)))
      if(is.infinite(ylim.max)) {
        cat("Infinite osa residuals for catch proportions at age in fleet ", mod$input$fleet_names[f], ", so using +/-10 for range \n")
        ylim.max = 10
      }
      ylims <- c(-ylim.max, ylim.max)

      #plot_years = rep(years, NCOL(resids))#years[as.integer(tmp$year)]
      # 1. trend vs. year
      plot(as.integer(pyears), resids, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[f], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      # 2. trend vs. age
      plot(as.integer(ages), resids, type='p', col=plot.colors[f], pch=19, xlab="Age", ylab="OSA Residuals",
           ylim=ylims, xaxt ="n")
      axis(1, at = 1:mod$input$data$n_ages, labels = mod$ages.lab)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[f], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      # 3. trend vs. cohort
      plot(min(years) + as.integer(cohorts), resids, type='p', col=plot.colors[f], pch=19, xlab="Cohort", ylab="OSA Residuals",
           ylim=ylims)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[f], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      # 4. trend vs. observed val
      plot(vals, resids, type='p', col=plot.colors[f], pch=19, xlab="Observed", ylab="OSA Residuals",
           ylim=ylims)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[f], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[f], lwd=2, lty = 2)
      # 5. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(resids, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)
      # 6. QQ plot modified from car:::qqPlot.default
      notNA <- resids[!is.na(resids)]
      ord.x <- notNA[order(notNA)]
      n <- length(ord.x)
      P <- ppoints(n)
      z <- qnorm(P, mean=0, sd=1)
      plot(z, ord.x, xlab="Std Normal Quantiles", ylab="OSA Residual Quantiles", main="", type = "n")
      grid(lty = 1, equilogs = FALSE)
      box()
      points(z, ord.x, col=plot.colors[f], pch=19)
      abline(0,1, col=plot.colors[f])
      conf = 0.95
      zz <- qnorm(1 - (1 - conf)/2)
      SE <- (1/dnorm(z)) * sqrt(P * (1 - P)/n)
      upper <- z + zz * SE
      lower <- z - zz * SE
      lines(z, upper, lty=2, col=plot.colors[f])
      lines(z, lower, lty=2, col=plot.colors[f])

      title (paste0("age composition OSA residual diagnostics: ", mod$input$fleet_names[f]), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }

  if("logindex" %in% mod$osa$type){
    dat <- subset(mod$osa, type == "logindex")
    dat$fleet <- sapply(dat$fleet, function(x) as.integer(strsplit(x, "_")[[1]][2]))
    dat$fleet <- factor(dat$fleet) #have to make sure integers are in order
    n.fleets <- length(levels(dat$fleet))
    plot.colors = mypalette(n.fleets)
    ifns <- chartr(" ", "_", mod$input$index_names)
    for(f in 1:n.fleets){
      tmp <- subset(dat, fleet==as.character(f))
      if(do.tex) cairo_pdf(file.path(od, paste0("OSA_resid_catch_4panel_", ifns[f],".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSA_resid_catch_4panel_", ifns[f],'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = fontfam)
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(dat$residual, na.rm=TRUE)))
      if(is.infinite(ylim.max)) {
        cat("Infinite osa residuals for aggregate indices in  index ", mod$input$index_names[f], ", so using +/-10 for range of y axis \n")
        ylim.max = 10
      }
      ylims <- c(-ylim.max, ylim.max)

      # make robust to missing years
      tmp$year <- years[tmp$year]
      tmp <- rbind(tmp[c("year","residual")], data.frame(year=years[!(years %in% tmp$year)], residual=rep(NA, length(years[!(years %in% tmp$year)]))))
      tmp <- tmp[order(tmp$year),]
      # 1. trend vs. year
      plot(tmp$year, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 2. trend vs. fitted val
      plot(mod$rep$pred_log_indices[1:length(mod$years),f], tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Log(Predicted Index)", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 3. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(tmp$residual, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)

      # 4. QQ plot modified from car:::qqPlot.default
      notNA <- tmp$residual[!is.na(tmp$residual)]
      ord.x <- notNA[order(notNA)]
      n <- length(ord.x)
      P <- ppoints(n)
      z <- qnorm(P, mean=0, sd=1)
      plot(z, ord.x, xlab="Std Normal Quantiles", ylab="OSA Residual Quantiles", main="", type = "n")
      grid(lty = 1, equilogs = FALSE)
      box()
      points(z, ord.x, col=plot.colors[f], pch=19)
      abline(0,1, col=plot.colors[f])
      conf = 0.95
      zz <- qnorm(1 - (1 - conf)/2)
      SE <- (1/dnorm(z)) * sqrt(P * (1 - P)/n)
      upper <- z + zz * SE
      lower <- z - zz * SE
      lines(z, upper, lty=2, col=plot.colors[f])
      lines(z, lower, lty=2, col=plot.colors[f])

      title (paste0("OSA residual diagnostics: ", mod$input$index_names[f]), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }
  if("indexpaa" %in% mod$osa$type) {
    dat <- subset(mod$osa, type == "indexpaa")
    dat$fleet <- sapply(dat$fleet, function(x) as.integer(strsplit(x, "_")[[1]][2]))
    dat$fleet <- factor(dat$fleet) #have to make sure integers are in order
    n.fleets <- length(levels(dat$fleet))
    dat$residual[which(is.infinite(dat$residual))] = NA #some happen for zeros or last age class
    #dat$residual[which(as.integer(dat$age) == mod$input$data$n_ages)] = NA #remove last age class
    
    n.indices <- mod$input$data$n_indices
    plot.colors = mypalette(n.indices)
    ifns <- chartr(" ", "_", mod$input$index_names)
    for(i in 1:n.indices) if(any(mod$input$data$use_index_paa[,i]==1)){
      if(do.tex) cairo_pdf(file.path(od, paste0("OSA_resid_paa_6panel_", ifns[i],".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSA_resid_paa_6panel_", ifns[i],'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = fontfam)

      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(3,2))
      yind = which(mod$input$data$use_index_paa[,i] ==1)
      resids = vals = ages = pyears = cohorts = matrix(NA, nrow = mod$input$data$n_years_model, 
        ncol = mod$input$data$n_ages)
      for(j in yind){
        tmp = subset(dat, year == j & fleet == as.character(i))
        resids[j,tmp$age] = tmp$residual
        vals[j,tmp$age] = tmp$val
        ages[j,tmp$age] = tmp$age
        pyears[j,tmp$age] = mod$years[tmp$year]
        cohorts[j,tmp$age] = tmp$cohort
      }

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(resids, na.rm=TRUE)))
      if(is.infinite(ylim.max)) {
        cat("Infinite osa residuals for proportions at age in  index ", mod$input$index_names[i], ", so using +/-10 for range \n")
        ylim.max = 10
      }
      ylims <- c(-ylim.max, ylim.max)

      #plot_years = rep(years, NCOL(resids))#years[as.integer(tmp$year)]
      # 1. trend vs. year
      plot(as.integer(pyears), resids, type='p', col=plot.colors[i], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[i], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      # 2. trend vs. age
      plot(as.integer(ages), resids, type='p', col=plot.colors[i], pch=19, xlab="Age", ylab="OSA Residuals",
           ylim=ylims, xaxt ="n")
      axis(1, at = 1:mod$input$data$n_ages, labels = mod$ages.lab)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[i], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      # 3. trend vs. cohort
      plot(min(years) + as.integer(cohorts), resids, type='p', col=plot.colors[i], pch=19, xlab="Cohort", ylab="OSA Residuals",
           ylim=ylims)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[i], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      # 4. trend vs. observed val
      plot(vals, resids, type='p', col=plot.colors[i], pch=19, xlab="Observed", ylab="OSA Residuals",
           ylim=ylims)
      grid(lty = 2, col = "grey")
      abline(h=0, col=plot.colors[i], lwd=2)
      abline(h=-qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      abline(h=qnorm(0.975), col=plot.colors[i], lwd=2, lty = 2)
      # 5. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(resids, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[i], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)
      # 6. QQ plot modified from car:::qqPlot.default
      notNA <- resids[!is.na(resids)]
      ord.x <- notNA[order(notNA)]
      n <- length(ord.x)
      P <- ppoints(n)
      z <- qnorm(P, mean=0, sd=1)
      plot(z, ord.x, xlab="Std Normal Quantiles", ylab="OSA Residual Quantiles", main="", type = "n")
      grid(lty = 1, equilogs = FALSE)
      box()
      points(z, ord.x, col=plot.colors[i], pch=19)
      abline(0,1, col=plot.colors[i])
      conf = 0.95
      zz <- qnorm(1 - (1 - conf)/2)
      SE <- (1/dnorm(z)) * sqrt(P * (1 - P)/n)
      upper <- z + zz * SE
      lower <- z - zz * SE
      lines(z, upper, lty=2, col=plot.colors[i])
      lines(z, lower, lty=2, col=plot.colors[i])

      title (paste0("age composition OSA residual diagnostics: ", mod$input$index_names[i]), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }

  if(!all(mod$env$data$Ecov_model == 0)){
    dat <- subset(mod$osa, type == "Ecov")
    dat$fleet <- sapply(dat$fleet, function(x) as.integer(strsplit(x, "_")[[1]][2]))
    dat$fleet <- factor(dat$fleet) #have to make sure integers are in order
    n.fleets <- length(levels(dat$fleet))
    plot.colors = mypalette(n.fleets)
    for(f in 1:n.fleets){
      tmp <- subset(dat, fleet==names(table(dat$fleet))[f])
      obs_ind <- which(mod$input$data$Ecov_use_obs[,f]==1)
      ecov_years <- mod$input$years_Ecov[obs_ind] 
      #ecov_years = mod$input$years_Ecov[1] + 0:(mod$env$data$n_years_Ecov-1)
      tmp$year <- ecov_years#[tmp$year+1] # year in osa is MODEL year, not Ecov year
      tmp$pred <- mod$rep$Ecov_x[obs_ind,f]
      tmp <- subset(tmp, !is.nan(tmp$residual))
      if(do.tex) cairo_pdf(file.path(od, paste0("OSA_resid_ecov_4panel_", chartr(" ", "_", mod$input$Ecov_names[f]),".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSA_resid_ecov_4panel_", chartr(" ", "_", mod$input$Ecov_names[f]),'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = fontfam)
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(dat$residual, na.rm=TRUE)))
      if(is.infinite(ylim.max)) {
        cat("Infinite osa residuals for Environmental observations in series ", mod$input$Ecov_names[f], ", so using +/-10 for range of y axis \n")
        ylim.max = 10
      }
      ylims <- c(-ylim.max, ylim.max)

      # 1. trend vs. year
      plot(tmp$year, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 2. trend vs. observed
      plot(tmp$val, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab=paste0("Observed ", mod$input$Ecov_names[f]), ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 3. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(tmp$residual, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)

      # 4. QQ plot modified from car:::qqPlot.default
      notNA <- tmp$residual[!is.na(tmp$residual)]
      ord.x <- notNA[order(notNA)]
      n <- length(ord.x)
      P <- ppoints(n)
      z <- qnorm(P, mean=0, sd=1)
      plot(z, ord.x, xlab="Std Normal Quantiles", ylab="OSA Residual Quantiles", main="", type = "n")
      grid(lty = 1, equilogs = FALSE)
      box()
      points(z, ord.x, col=plot.colors[f], pch=19)
      abline(0,1, col=plot.colors[f])
      conf = 0.95
      zz <- qnorm(1 - (1 - conf)/2)
      SE <- (1/dnorm(z)) * sqrt(P * (1 - P)/n)
      upper <- z + zz * SE
      lower <- z - zz * SE
      lines(z, upper, lty=2, col=plot.colors[f])
      lines(z, lower, lty=2, col=plot.colors[f])

      title (paste0("OSA residual diagnostics: ", mod$input$Ecov_names[f]), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }
}

#revised

mypalette = function(n){
  palette.fn <- colorRampPalette(c("dodgerblue","green","red"), space = "Lab")
  palette.fn(n)
}

fit.summary.text.plot.fn <- function(mod){
  acm = c("Multinomial", "Dirichlet-multinomial", "Dirichlet (miss0)", "Dirichlet (pool0)","Logistic normal (miss0)",
    "Logistic normal AR1 corr (miss0)", "Logistic normal (pool0)", "ZI-logistic normal(1)","ZI-logistic normal(2)", "MV Tweedie", 
    "Dirichlet-multinomial (linearized)")
  selmods = c("Age-specific", "Logistic(+)", "Double-Logistic", "Logistic(-)")
  recs <- c("Random walk","Random about mean","Bev-Holt","Ricker")
  env.mod <- c("RW", "AR1")
  # env.where <- c('Recruitment','Growth','Mortality')
  env.where <- c('Recruitment','Mortality',paste0("q for index ", 1:mod$env$data$n_indices))  
  env.how <- c("controlling", "limiting", "lethal", "masking", "directive")
  fleet_selblocks = lapply(1:mod$env$data$n_fleets, function(x) unique(mod$env$data$selblock_pointer_fleets[,x]))
  index_selblocks = lapply(1:mod$env$data$n_indices, function(x) unique(mod$env$data$selblock_pointer_indices[,x]))
  plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
  nl = 10
  text(5,nl <- nl-0.5,mod$model_name)
  text(5,nl <- nl-0.5,paste0("Model years: ", min(mod$years), "-", max(mod$years)))
  if(mod$env$data$n_years_proj>0){
    proj.yrs <- tail(mod$years_full, mod$env$data$n_years_proj)
    proj.txt <- ifelse(mod$env$data$n_years_proj>0, "none", paste0(min(proj.yrs), "-", max(proj.yrs)))
    text(5,nl <- nl-0.5,paste0("Projection years: ", proj.txt))
  }
  text(5,nl <- nl-0.5,paste0("Number of stocks: ", mod$env$data$n_stocks))
  text(5,nl <- nl-0.5,paste0("Number of regions: ", mod$env$data$n_regions))
  text(5,nl <- nl-0.5,paste0("Number of fleets: ", mod$env$data$n_fleets))
  text(5,nl <- nl-0.5, paste0("Fleet Age Comp Models: ", paste(acm[mod$env$data$age_comp_model_fleets], collapse = ", ")))
  text(5,nl <- nl-0.5,paste0("Number of indices: ", mod$env$data$n_indices))
  text(5,nl <- nl-0.5, paste0("Index Age Comp Models: ", paste(acm[mod$env$data$age_comp_model_indices], collapse = ", ")))
  text(5,nl <- nl-0.5,paste0("Recruitment model for each stock: ", recs[mod$env$data$recruit_model]))
  dat <- mod$env$data
  if(!all(mod$env$data$Ecov_model == 0)){
    for(ec in 1:mod$env$data$n_Ecov) {
      ec.mod = env.mod[mod$env$data$Ecov_model[ec]]
      ec.label = mod$input$Ecov_names[ec]
      if(length(ec.mod)) {
        out = paste0("Environmental covariate ", ec.label, " modeled as ",ec.mod,".\n")
        for(i in 1:dat$n_stocks){
          if(dat$Ecov_how_R[ec,i]>0) out = paste0(out, paste0(" Effects on recruitment for stock ", i, " assumed. \n"))
          for(a in 1:dat$n_ages) for(r in 1:dat$n_regions) if(dat$Ecov_how_M[ec,i,a,r]>0){
            out = paste0(out, paste0(" Effects on M for stock ", i, "at age ", a, "in region ", r, " assumed. \n"))
          }
          if(dat$n_regions>1) for(a in 1:dat$n_ages) for(s in 1:dat$n_seasons) for(r in 1:dat$n_regions) for(rr in 1:(dat$n_regions-1)){
            if(dat$Ecov_how_mu[ec,i,a,s,r,rr]>0) {
              out = paste0(out, paste0(" Effects on movement for stock ", i, "at age ", a, "in season ", s, " from region ", r, "to region ",
              ifelse(rr>=r, rr+1, rr), "  assumed. \n"))
            }
          }
        }
        for(i in 1:dat$n_indices){
          if(dat$Ecov_how_q[ec,i]>0) out = paste0(out, paste0(" Effects on catchability for index ", i, " assumed. \n"))
        }
        if(sum(dat$Ecov_how_R,dat$Ecov_how_M, dat$Ecov_how_mu,dat$Ecov_how_q)==0) out = paste(out, " No effects estimated.")
        text(5,nl <- nl-0.5, out)
      }
    }
  } else {
    text(5,nl <- nl-0.5, "No Environmental covariates.")
  }
  text(5,nl <- nl-0.5,paste0("Number of Selectivity blocks: ", mod$env$data$n_selblocks))
  text(5,nl <- nl-0.5, paste0("Selectivity Block Types: ", paste(selmods[mod$env$data$selblock_models], collapse = ", ")))
  for(i in 1:length(fleet_selblocks)) text(5,nl <- nl-0.5, paste0("Fleet ", i, " Selectivity Blocks: ", paste(fleet_selblocks[i], collapse = ", ")))
  for(i in 1:length(index_selblocks)) text(5,nl <- nl-0.5, paste0("Index ", i, " Selectivity Blocks: ", paste(index_selblocks[i], collapse = ", ")))
  text(5,nl <- nl-0.5,paste0("Run date: ",format(mod$date, usetz = TRUE)))
  text(5,nl <- nl-0.5,paste0("Run directory: ",mod$dir))
  text(5,nl <- nl-0.5,paste0("WHAM version: ",mod$wham_version)) 
  text(5,nl <- nl-0.5,paste0("TMB version: ",mod$TMB_version)) 
  if(is.null(mod$opt)){
      text(5,nl <- nl-0.5,paste0("Model has NOT been fit."))
  } else {
    if(mod$is_sdrep){
      text(5,nl <- nl-0.5,paste0("sdreport() performed",
        ifelse(mod$na_sdrep, ", but with NAs for some variance estimates.", " with all variance estimates provided.")))
    }
    else {
      text(5,nl <- nl-0.5,"Warning: run did not provide pos-def Hessian or sdreport() not performed", col="red")
    }
    mgind = which(abs(mod$final_gradient) == max(abs(mod$final_gradient)))
    text(5,nl <- nl-0.5,paste0("Maximum absolute gradient: ", names(mod$opt$par)[mgind]," ", format(mod$final_gradient[mgind],digits =4)))
    text(5,nl <- nl-0.5,paste0("Number of fixed effects = ",length(mod$opt$par),", Number of random effects = ", length(mod$env$random)))
  }

  return()
}

#revised

plot.ll.table.fn <- function(mod,plot.colors){
  par(mfrow=c(1,1))

  npar <- length(mod$opt$par)
  lls <- mod$rep[grep("nll",names(mod$rep))]
  ll.names <- names(lls)
  #n.like <- length(lls)
  n_fleets <- mod$env$data$n_fleets
  n_indices <- mod$env$data$n_indices
  obs.ll.names <-c("nll_agg_catch", "nll_catch_acomp", "nll_agg_indices", "nll_index_acomp", "nll_Ecov_obs")
  obs.lls <- lls[names(lls) %in% obs.ll.names]
  obs.lls <- lapply(obs.lls, function(x) apply(x,2,sum))
  fleet.names <-  paste0(mod$input$fleet_names, " ", mod$input$region_names[mod$input$data$fleet_regions])
  index.names <- paste0(mod$input$index_names, " ", mod$input$region_names[mod$input$data$index_regions])
  names(obs.lls$nll_agg_catch) <- paste0(fleet.names, " Catch")
  names(obs.lls$nll_catch_acomp) <- paste0(fleet.names, " Age Comp")
  names(obs.lls$nll_agg_indices) <- paste0(index.names, " Catch")
  names(obs.lls$nll_index_acomp) <- paste0(index.names, " Age Comp")
  #if(sum(mod$env$data$Ecov_use_obs)>0) {
    obs.lls$nll_Ecov_obs <- apply(lls$nll_Ecov_obs,2,sum)
    names(obs.lls$nll_Ecov_obs) <- paste0("Ecov: ", mod$input$Ecov_names, " observations")
  #}
  #print(obs.lls)
  #stop()
  #names(obs.lls) = NULL
  #obs.lls = unlist(obs.lls)
  #print(obs.lls)
  #n.obs.ll = length(obs.lls)
  #obs.dists = character(n.obs.ll)
  obs.dists <- obs.lls
  obs.dists$nll_agg_catch[] <- "log(x) ~ Gaussian"
  obs.dists$nll_agg_indices[] <- "log(x) ~ Gaussian"
  acm = c("Multinomial", "Dirichlet-multinomial", "Dirichlet (miss0)", "Dirichlet (pool0)","Logistic normal (miss0)",
    "Logistic normal AR1 corr (miss0)", "Logistic normal (pool0)", "ZI-logistic normal(1)","ZI-logistic normal(2)", "MV Tweedie", 
    "Dirichlet-multinomial (linearized)")
  obs.dists$nll_catch_acomp[] <- paste0("x ~ ", acm[mod$env$data$age_comp_model_fleets])
  obs.dists$nll_index_acomp[] <- paste0("x ~ ", acm[mod$env$data$age_comp_model_indices])
  if(!is.null(obs.dists$nll_Ecov_obs)) obs.dists$nll_Ecov_obs[] <- "x ~ Gaussian"

  proc.rep.names <- c("nll_NAA", "nll_N1", "nll_mu_re", "nll_mu_prior", "nll_M", "nll_log_b", "nll_L", "nll_sel", 
    "nll_q_re", "nll_q_prior", "nll_Ecov", "nll_Ecov_obs_sig")
  proc.names <- c("NAA", "N1", "move", "move_prior", "M", "M_b_prior", "Extra_M", "Selex", 
    "q", "q_prior", "Ecov", "Ecov_obs_sig")
  these.proc.names <- names(lls)[which(names(lls) %in% proc.rep.names)]
  proc.lls = lls[these.proc.names]
  names(proc.lls) = proc.names[match(names(proc.lls),proc.rep.names)]
  proc.lls = unlist(lapply(proc.lls, sum, na.rm = TRUE))
  n.proc.ll = length(proc.lls)
  proc.dists = character(n.proc.ll)
  proc.dists <- character()
  proc.dists[which(names(proc.lls) %in% c("NAA","N1", "M", "M_b_prior", "Extra_M", "Ecov_obs_sig"))] <- "log(x) ~ MVN"
  proc.dists[which(names(proc.lls) %in% c("move","move_prior"))] <- "add_logit(x) ~ MVN"
  proc.dists[which(names(proc.lls) %in% c("Selex","q", "q_prior"))] <- "logit(x) ~ MVN"
  proc.dists[which(names(proc.lls) %in% c("Ecov"))] <- "x ~ MVN"
  names(proc.dists) = names(proc.lls)
  likes = -c(obs.lls[[1]],obs.lls[[2]],obs.lls[[3]],obs.lls[[4]], obs.lls[[5]], proc.lls)
  dists = c(obs.dists[[1]],obs.dists[[2]],obs.dists[[3]],obs.dists[[4]],obs.dists[[5]], proc.dists)
  n.likes = length(likes)
  my.range <- c(min(likes)-50,max(likes)+50)#1.2*range(likes)#c(min(likes), 1.2*max(likes))
  par(mar=c(5,20,1,1), oma=c(1,0,0,0))
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=n.likes) #mypalette(n.likes)
  barplot(horiz=TRUE, likes, beside=FALSE, col=plot.colors, xlab="Joint log-likelihood components",  axisnames=FALSE,  axes=FALSE,  space=0,
    xlim=my.range)
  axis(side=1, at=pretty(seq(my.range[1],my.range[2]), n=10), labels=pretty(seq(my.range[1],my.range[2]), n=10), cex=.75 )
  axis(side=2, at=seq(0.5,(n.likes-0.5)), labels= paste0(names(likes), "\n", dists), las=2)
  #axis(side=2, at=seq(0.5,(n.likes-0.5))-0.2, labels = dists, las = 2, tick = FALSE)
  text(x= likes, y=seq(0.5,(n.likes-0.5)), labels=round(likes,0), cex=0.8, pos=ifelse(likes>0, 4, 2))
  box()
  #title(paste0("Components of Obj. Function (", round(as.numeric(asap$like[1]),0), "), npar=", npar), cex=0.9 )
  title(sub=paste0("Model: ", mod$model_name, "     ", mod$date))
}

#revised

residual.bubbleplot.fn = function(x, x.cats, y.cats, scale.bubble = 0.25*25, ylab = "Pearson Residuals", xlab = "Age"){
  resids <- x  # NOTE obs-pred
  n_y = NCOL(x)
  n_x = NROW(x)
  pos.col = "#ffffffaa"
  neg.col = "#8c8c8caa"
  range.resids<-range(abs((as.vector(x))), na.rm=T)
  #scale.resid.bubble <- 25

  z <- x * scale.bubble #* scale.resid.bubble
  resid.col = ifelse(z > 0.0, pos.col, neg.col)

  plot(1:n_x, 1:n_y,  xlim = c(1,n_x), ylim = c(1,n_y), xlab = x_lab, ylab = ylab, type = "n", axes=F)
  axis(1, lab = x.cats, lwd = 2)
  axis(2, lab = y.cats, cex.axis=0.75, las=1, lwd = 2)
  box(lwd = 2)
  abline(h = 1:n_y, col="lightgray")
  segments(x0 = 1:n_x, y0 = rep(1,n_x), x1= 1:n_x, y1 = rep(n_y,n_x), col = "lightgray", lty = 1)
  for (j in 1:n_y)
  {
    points(1:n_x, rep(j, n_x), cex = abs(z[j,]), col="black", bg = resid.col[j,],  pch = 21)
  }

  bubble.legend1 = c(0.1, 1, 2)
  #bubble.legend1 <- c(0.01,0.05,0.1)
  bubble.legend2 <- bubble.legend1 * scale.bubble#scale.resid.bubble*scale.bubble
  legend("topright", xpd = TRUE, legend = bubble.legend1, pch = rep(1, 3), pt.cex = bubble.legend2, horiz = TRUE, col = 'black')
  legend("topleft", xpd = TRUE, legend = c("Neg.", "Pos."), pch = rep(21, 2), pt.cex = 3, horiz = TRUE,
    pt.bg = c(neg.resid.col, pos.resid.col), col = "black")
  text(x = trunc(n_x/2), y = 0, cex = 0.8, label=paste0("Max(abs(resid)) = ",round(range.resids[2],2)))
}

hi.cor.fn <- function(mod, out.dir = "", do.tex = FALSE, do.csv = FALSE, cor.limit = 0.99)
{
  cor = mod$sdrep$cov.fixed
  cor = cor/sqrt(diag(cor) %*% t(diag(cor)))
  npar = NROW(cor)
	which.hi <- matrix(nrow = 0, ncol = 2)
	for(i in 2:npar) for(j in 1:(i-1))  if(cor[i,j]>cor.limit) which.hi = rbind(which.hi, c(i, j))

  hi.cor = cor[which.hi[,1],which.hi[,2]]
  if(length(which.hi))
  {
    out = cbind.data.frame(
      "Parameter 1" = rownames(cor)[which.hi[,1]], "Parameter 1 value" = round(mod$opt$par[which.hi[,1]],3),
      "Parameter 2" = colnames(cor)[which.hi[,2]], "Parameter 2 value" = round(mod$opt$par[which.hi[,2]],3),
      "Correlation" = round(cor[which.hi],3))

    if(do.tex) x = latex(out, file = paste0(out.dir,"hi_cor.tex"), rowlabel = '', table.env = FALSE)
    if(do.csv) write.csv(out, file = paste0(out.dir,"hi_cor.csv"))
    return(out)
  }
  else cat(paste0("No fixed effects had correlations estimated greater than ", cor.limit, " \n"))
}

get.RMSEs.fn <- function(model)
{
  out <- list(RMSE=list(), RMSE_n = list())
  if(class(mod$sdrep)[1] == "sdreport"){
    #sdrep = summary(model$sdrep)
    # catch_stdresid <- as.list(model$sdrep, "Est", report=TRUE)$log_catch_resid/as.list(model$sdrep, "Std", report=TRUE)$log_catch_resid
    # index_stdresid <- as.list(model$sdrep, "Est", report=TRUE)$log_index_resid/as.list(model$sdrep, "Std", report=TRUE)$log_index_resid
    catch_stdresid <- model$rep$log_catch_resid/model$input$data$agg_catch_sigma
    index_stdresid <- model$rep$log_index_resid/model$input$data$agg_index_sigma
  # } else {
  #   sdrep = model$sdrep
  # }  
    # temp = sdrep[rownames(sdrep) %in% "log_catch_resid",]
    # catch_stdresid <- matrix(temp[,1]/temp[,2], model$env$data$n_years_model, model$env$data$n_fleets)
    # temp = sdrep[rownames(sdrep) %in% "log_index_resid",]
    # index_stdresid <- matrix(temp[,1]/temp[,2], model$env$data$n_years_model, model$env$data$n_indices)
    # temp = model$env$data$catch_paa - aperm(model$rep$pred_catch_paa[1:model$env$data$n_years_model,,,drop=FALSE],c(2,1,3))
    #temp = model$env$data$catch_paa - aperm(model$rep$pred_catch_paa,c(2,1,3))
    if(length(catch_stdresid)){
      out$RMSE$catch <- sqrt(apply(catch_stdresid^2,2,mean, na.rm = TRUE))
      out$RMSE_n$catch <- apply(catch_stdresid^2,2,function(x) sum(!is.na(x)))
      out$RMSE$catch_total <- sqrt(mean(catch_stdresid^2, na.rm = TRUE))
      out$RMSE_n$catch_total <- sum(!is.na(catch_stdresid^2))
    }
    if(length(index_stdresid)){
      out$RMSE$index <- sqrt(apply(index_stdresid^2,2,mean, na.rm = TRUE))
      out$RMSE_n$index <- apply(index_stdresid^2,2,function(x) sum(!is.na(x)))
      out$RMSE$index_total <- sqrt(mean(index_stdresid^2, na.rm = TRUE))
      out$RMSE_n$index_total <- sum(!is.na(index_stdresid^2))
    }
    if(length(out$RMSE)) return(out)
  }
}
RMSE.table.fn <- function(model)
{
  origpar <- par(no.readonly = TRUE)
  RMSEs <- get.RMSEs.fn(model)
  if(length(RMSEs)){
    par(mfrow=c(1,1), oma=rep(2,4), mar=c(0,0,1,0))
    max.txt <- 16
    rmses <- which(unlist(RMSEs$RMSE)>0)
    n.rmses <- length(rmses)
    plot(seq(1,15), seq(1,15), type='n', axes=F, xlab="", ylab="", xlim=c(1,max.txt+2+ 6+2+ 8+2), ylim=c(n.rmses+4, 1))
    text(rep(1, n.rmses), seq(3,n.rmses+2), labels=names(unlist(RMSEs$RMSE)[rmses]), pos=4)
    text(rep(max.txt+2, n.rmses), seq(3,n.rmses+2), labels=unlist(RMSEs$RMSE_n)[rmses], pos=4)
    text(rep(max.txt+2+ 6+2, n.rmses), seq(3,n.rmses+2), labels=signif(as.numeric(unlist(RMSEs$RMSE)[rmses]),3), pos=4)
    text(c(1, max.txt+2, max.txt+2+ 6+2), rep( 2, 3), labels=c("Component","# resids","RMSE"), font=2, pos=4)
    title(main="Root Mean Square Error computed from Standardized Residuals", outer=T, cex=0.85)
    par(origpar)
  }
	#if (save.plots) savePlot(paste(od, "RMSE_Comp_Table.",plotf, sep=""), type=plotf)
}

#revised

get.wham.results.fn <- function(mod, out.dir, do.tex = FALSE, do.png = FALSE)
{
  years <- mod$years_full #not used in cpp
  ages <- mod$ages.lab #not used in cpp, assumes last age is a plus group
  ny <- length(years)
  na <- length(ages)
  ni <- mod$env$data$n_indices
  nf <- mod$env$data$n_fleets
  ns <- mod$env$data$n_stocks
  nr <- mod$env$data$n_regions

  res <- list()
  res$ll <- -mod$opt$obj
  res$np <- length(mod$opt$par)
  res$aic <- 2*(mod$opt$obj + res$np)
  #rho <- mohns_rho(mod) #if no peels, then this will be NULL
  # tcol <- col2rgb('black')
  # black.poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  # tcol <- col2rgb('red')
  # red.poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  # tcol <- col2rgb('blue')
  # blue.poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')

  origpar <- par(no.readonly = TRUE)
  par(mar = c(5,5,1,1), oma = c(1,1,1,1), mfrow = c(3,1))
  use_outer <- FALSE
  x_line <- -1
  y_line <- 3
  fn <- paste0(out.dir, '/SSB')
  if(do.tex | do.png) {
    if(do.tex) cairo_pdf(file.path(out.dir,'SSB.pdf'), family = fontfam, height = 10, width = 10)
    else png(filename = file.path(out.dir,, 'SSB.png'), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mar = c(0,0,0,0), oma = c(4,4,1,1), mfrow = c(1,1))
    use_outer <- TRUE
    x_line <- y_line <- 2.5
  }
  SSB <- SSB.lo <- SSB.hi <- R <- R.lo <- R.hi <- matrix(NA, nrow = ny, ncol = ns)
  Fbar <- Fbar.lo <- Fbar.hi <- matrix(NA, nrow = ny, ncol = nr)
  if(class(mod$sdrep)[1] == "sdreport"){
    temp <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = T),
      TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = T))
    SSB <- exp(temp[[1]]$log_SSB)
    SSB.lo <- SSB * exp( - qnorm(0.975)*temp[[2]]$log_SSB)
    SSB.hi <- SSB * exp( qnorm(0.975)*temp[[2]]$log_SSB)
    for(i in 1:ns) {
      R[,i] <- exp(temp[[1]]$log_NAA_rep[i,mod$env$data$spawn_regions[i],,1])
      R.lo[,i] <- R[,i] * exp(- qnorm(0.975)*temp[[2]]$log_NAA_rep[i,mod$env$data$spawn_regions[i],,1])
      R.hi[,i] <- R[,i] * exp( qnorm(0.975)*temp[[2]]$log_NAA_rep[i,mod$env$data$spawn_regions[i],,1])
    }
    Fbar <- exp(temp[[1]]$log_Fbar)
    Fbar.lo <- Fbar * exp(- qnorm(0.975)*temp[[2]]$log_Fbar)
    Fbar.hi <- Fbar * exp( qnorm(0.975)*temp[[2]]$log_Fbar)
  } else {
    SSB[] <- mod$rep$SSB
    for(i in 1:ns) R[,i] <- mod$rep$NAA[i,mod$env$data$spawn_regions[i],,1]
    Fbar[] <- mod$rep$Fbar
  }

  max.y <- max(SSB.hi)
  na.se <- is.na(max.y)
  if(na.se) max.y <- max(SSB)
  pal = viridisLite::viridis(n=ns)
  plot(years,SSB[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  grid(col = gray(0.7))
  for(i in 1:ns) {
    lines(years,SSB[,i], lwd = 2, col = pal[i])
    if(!na.se) polygon(c(years,rev(years)), c(SSB.lo[,3],rev(SSB.hi[,4])), col = adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
  }
  if(ns>1) legend("topright", col = pal, lty = 1, lwd = 2, legend = mod$stock.names)
  box(lwd = 2)
  mtext(side = 1, "Year", cex = 2, outer = TRUE, line = x_line)
  mtext(side = 2, "SSB (kmt)", cex = 2, outer = use_outer, line = y_line)
  if(do.tex | do.png) dev.off() else par(origpar)

  if(do.tex | do.png) {
    if(do.tex) cairo_pdf(file.path(out.dir,'recruits.pdf'), family = fontfam, height = 10, width = 10)
    else png(filename = file.path(out.dir, 'recruits.png'), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mar = c(0,0,0,0), oma = c(4,4,1,1), mfrow = c(1,1))
  }
  max.y <- max(R.hi)
  na.se <- is.na(max.y)
  if(na.se) max.y <- max(R)
  plot(years,R[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  grid(col = gray(0.7))
  for(i in 1:ns){
    lines(years,R[,i], lwd = 2, col = pal[i])
    if(!na.se) polygon(c(years,rev(years)), c(R.lo[,i],rev(R.hi[,i])), col = adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
  }
  if(ns>1) legend("topright", col = pal, lty = 1, lwd = 2, legend = mod$stock.names)
  box(lwd = 2)
  if(use_outer) mtext(side = 1, "Year", cex = 2, outer = TRUE, line = x_line)
  mtext(side = 2, "Recruits (1000s)", cex = 2, outer = use_outer, line = y_line)
  if(do.tex | do.png) dev.off() else par(origpar)

  if(do.tex | do.png) {
    if(do.tex) cairo_pdf(file.path(out.dir,'Fbar.pdf'), family = fontfam, height = 10, width = 10)
    else png(filename = file.path(out.dir, 'Fbar.png'), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mar = c(0,0,0,0), oma = c(4,4,1,1), mfrow = c(1,1))
  }
  max.y <- max(Fbar.hi[,4])
  na.se <- is.na(max.y)
  if(na.se) max.y <- max(Fbar)
  plot(years,Fbar[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  grid(col = gray(0.7))
  for(i in 1:nr){
    lines(years,Fbar[,i], lwd = 2, col = pal[i])
    if(!na.se) polygon(c(years,rev(years)), c(Fbar.lo[,3],rev(Fbar.hi[,4])), col = adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
  }
  box(lwd = 2)
  if(use_outer) mtext(side = 1, "Year", cex = 2, outer = TRUE,line = x_line)
  ar <- c(min(mod$env$data$Fbar_ages), max(mod$env$data$Fbar_ages))
  mtext(side = 2, paste0("Average F (",ages[ar[1]],"-",ages[ar[2]],")"), cex = 2, outer = use_outer, line = y_line)
  if(do.tex | do.png) dev.off() else par(origpar)
  # par(origpar)
}

#revised

plot.all.stdresids.fn = function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
  years = mod$years
  # load Ecov residuals
  xe <- NULL
  if(!all(mod$env$data$Ecov_model == 0)){
    ny = mod$env$data$n_years_Ecov
    ni = mod$env$data$n_Ecov
    Ecov_resid <- Ecov_resid.lo <- Ecov_resid.hi <- matrix(NA, nrow = ny, ncol = ni)
    if(class(mod$sdrep)[1] == "sdreport"){
      temp <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = T),
        TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = T))
      Ecov_resid[] <- temp[[1]]$Ecov_resid
      Ecov_resid.lo[] <- Ecov_resid - qnorm(0.975)*temp[[2]]$Ecov_resid
      Ecov_resid.hi[] <- Ecov_resid + qnorm(0.975)*temp[[2]]$Ecov_resid
    } else {
      Ecov_resid <- mod$rep$Ecov_resid
    }
    # ind = rownames(temp) == "Ecov_resid"
    # templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
    # temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
    # temp = matrix(temp[ind,1], ny, ni)
    xe = data.frame(Label = integer(),
      Year = numeric(),
      Stdres = numeric(),
      lo = numeric(),
      hi = numeric())
    for(i in 1:ni)
    {
      ind = which(mod$env$data$Ecov_use_obs[,i] == 1)
      td = data.frame(Label = rep(i,length(ind)),
        Year = mod$input$years_Ecov[ind],
        Stdres = Ecov_resid[ind,i],
        lo = Ecov_resid.lo[ind,i],
        hi = Ecov_resid.hi[ind,i])
      xe <- rbind(xe, td)
    }
    xe$row = xe$Label
    xe$Label = factor(xe$Label)
    levels(xe$Label) = mod$input$Ecov_names
    xe$type = "Ecov"
  }

  # load Index residuals
  ny = mod$env$data$n_years_model
  ni = mod$env$data$n_indices
  nf = mod$env$data$n_fleets
  log_index_resid <- log_index_resid.lo <- log_index_resid.hi <- matrix(NA, ny, ni)
  log_catch_resid <- log_catch_resid.lo <- log_catch_resid.hi <- matrix(NA, ny, nf)
  log_catch_resid[] <- mod$rep$log_catch_resid
  log_index_resid[] <- mod$rep$log_index_resid
  log_catch_resid.lo[] <- log_catch_resid - qnorm(0.975)*mod$input$data$agg_catch_sigma
  log_catch_resid.hi[] <- log_catch_resid + qnorm(0.975)*mod$input$data$agg_catch_sigma
  log_index_resid.lo[] <- log_index_resid - qnorm(0.975)*mod$input$data$agg_index_sigma
  log_index_resid.hi[] <- log_index_resid + qnorm(0.975)*mod$input$data$agg_index_sigma

  if(class(mod$sdrep)[1] == "sdreport"){
    temp <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = T),
      TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = T))
      if(!is.null(temp[[1]]$log_index_resid)){
        log_index_resid[] <- temp[[1]]$log_index_resid
        log_index_resid.lo[] <- log_index_resid - qnorm(0.975)*temp[[2]]$log_index_resid
        log_index_resid.hi[] <- log_index_resid + qnorm(0.975)*temp[[2]]$log_index_resid
        log_catch_resid[] <- temp[[1]]$log_catch_resid
        log_catch_resid.lo[] <- log_catch_resid - qnorm(0.975)*temp[[2]]$log_catch_resid
        log_catch_resid.hi[] <- log_catch_resid + qnorm(0.975)*temp[[2]]$log_catch_resid
      }
    } else {
  }

  xi = data.frame(Label = integer(),
    Year = numeric(),
    Stdres = numeric(),
    lo = numeric(),
    hi = numeric())
  for(i in 1:ni)
  {
    ind = which(mod$env$data$use_indices[,i] == 1)
    td = data.frame(Label = rep(i,length(ind)),
      Year = years[ind],
      Stdres = log_index_resid[ind,i],
      lo = log_index_resid.lo[ind,i],
      hi = log_index_resid.hi[ind,i])
    xi <- rbind(xi, td)
  }
  xi$row = xi$Label
  xi$Label = factor(xi$Label)
  levels(xi$Label) = mod$input$index_names #paste0("Index ",1:length(table(xi$Label)))
  xi$type = "Index"
  # if(!is.null(index.names)) levels(x$Index) = index.names

  # load catch data (fleet)
  ni = mod$env$data$n_fleets
  xc = data.frame(Label = integer(),
    Year = numeric(),
    Stdres = numeric(),
    lo = numeric(),
    hi = numeric())
  for(i in 1:ni)
  {
    td = data.frame(Label = rep(i,ny),
      Year = years,
      Stdres = log_catch_resid[,i],
      lo = log_catch_resid.lo[,i],
      hi = log_catch_resid.hi[,i])
    xc <- rbind(xc, td)
  }
  xc$row = xc$Label
  xc$Label = factor(xc$Label)
  levels(xc$Label) = mod$input$fleet_names #paste0("Fleet ",1:length(table(xc$Label)))
  xc$type = "Catch"
  # if(!is.null(fleet.names)) levels(xc$Fleet) = fleet.names

  x <- rbind(xe, xi, xc)
  x$row = factor(x$row)
  x$type = factor(x$type)

  ggp = ggplot2::ggplot(x, ggplot2::aes(x=Year, y = Stdres, color=type, fill=type)) +
    # ggplot2::geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=type), alpha=0.3, linetype = 0) +
    # ggplot2::geom_line(size=1.1) +
    #ggplot2::geom_smooth(method = "lm", alpha=0.2) +
    ggplot2::geom_smooth(formula = y ~ x, method = "lm", alpha=0.2) +
    ggplot2::geom_point(size=0.8) +
    ggplot2::ylab("Standardized residual") +
    # expand_limits(y=0) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    # ggplot2::scale_color_manual(values=plot.colors) +
    # ggplot2::scale_fill_manual(values=plot.colors) +
    ggplot2::facet_wrap(~Label)
    # ggplot2::facet_grid(type ~ row)
    # ggplot2::facet_grid(row ~ type)
  if(do.tex) cairo_pdf(file.path(od, paste0("Residuals_time.pdf")), family = fontfam, height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("Residuals_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
  print(ggp)
  if(do.tex | do.png) dev.off()
  # return(ggp)
}

#revised

plot.catch.4.panel <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  dat = mod$env$data
  years_full = mod$years_full
  pred_log_catch = mod$rep$pred_log_catch
  pred_catch = mod$rep$pred_catch
  sigma = dat$agg_catch_sigma %*% diag(exp(mod$parList$log_catch_sig_scale), nrow = length(mod$parList$log_catch_sig_scale)) # dims: [ny,nf] x [nf]
  catch = dat$agg_catch
  log_stdres = (log(catch) - pred_log_catch[1:length(years),])/sigma # cpp already bias-corrects if bias_correct_oe = 1
  if(!missing(use.i)) fleets <- use.i
  else fleets <- 1:dat$n_fleets
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=dat$n_fleets) #mypalette(dat$n_fleets)
  ffns <- chartr(" ", "_", mod$input$fleet_names)
  rfns <- chartr(" ", "_", mod$input$region_names)
	for (i in fleets)
	{
    fleet.name.fn <- paste0(ffns[i], "_", rfns[mod$input$data$fleet_regions[i]])
    fleet.name.plt <- paste0(mod$input$fleet_names[i], " in ", mod$input$region_names[mod$input$data$fleet_regions[i]])
		if(do.tex) cairo_pdf(file.path(od, paste0("Catch_4panel_", fleet.name.fn,".pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_4panel_fleet_",fleet.name.fn,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
		plot(years_full, pred_catch[,i], col=plot.colors[i], lwd=2, type='l', xlab="Year", ylab="Total Catch",
			ylim=c(0, 1.1*max(c(catch[,i],pred_catch[,i]))))
    points(years, catch[,i], col=plot.colors[i], pch=1)
    if(mod$env$data$n_fleets == 1){
      if(length(years_full) > length(years)){
        abline(v=tail(years,1), lty=2, lwd=1)
      }
    }
		log.ob.min <- log(catch[,i])-1.96*sigma[,i]
		log.ob.max <- log(catch[,i])+1.96*sigma[,i]
		plot(years_full, log(pred_catch[,i]), col=plot.colors[i], lwd=2, type='l', xlab="Year", ylab="Ln(Total Catch)",
			ylim=c(min(log.ob.min,log(pred_catch[,i]), na.rm=T), 1.1*max(log.ob.max,log(pred_catch[,i]), na.rm=T)))
		points(years, log(catch[,i]), pch=1, col=plot.colors[i])
    if(mod$env$data$n_fleets == 1) abline(v=tail(years,1), lty=2, lwd=1)
		arrows(years, log.ob.min, years, log.ob.max, length=0)
		#title (paste0("Fleet ",i, " Catch"), outer=T, line=-1)
		title(fleet.name.plt, outer=T, line=-1)
		plot(years, log_stdres[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Log-scale Std. Residual")
		abline(h=0)
		hist(log_stdres[,i], plot=T, xlab="Std. Residual", ylab="Probability Density", freq=F, main=NULL)
		if(do.tex | do.png) dev.off() else par(origpar)
	}
}

#revised

plot.index.4.panel <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  dat = mod$env$data
  pred_index = exp(mod$rep$pred_log_indices[1:length(years),,drop=F])
  index = dat$agg_indices
  # index[index < 0] = NA # robustify to missing values entered as negative
  index[dat$use_indices == 0] = NA # don't plot unused values
  sigma = dat$agg_index_sigma %*% diag(exp(mod$parList$log_index_sig_scale), nrow = length(mod$parList$log_index_sig_scale)) # dims: [ny,ni] x [ni]
  log_stdres = (log(index)-log(pred_index))/sigma
  if(!missing(use.i)) indices <- use.i
  else indices <- 1:dat$n_indices
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=dat$n_indices) #mypalette(dat$n_indices)
	for (i in indices)
	{
    index.name.fn <- paste0(chartr(" ", "_", mod$input$index_names[i]), "_", chartr(" ", "_", mod$input$region_names[mod$input$data$index_regions[i]]))
    index.name.plt <- paste0(mod$input$index_names[i], " in ", mod$input$region_names[mod$input$data$index_regions[i]])
		if(do.tex) cairo_pdf(file.path(od, paste0("Index_4panel_",index.name.fn,".pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Index_4panel_",index.name.fn,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
		plot(years, index[,i], type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Index",
			ylim=c(0, 1.1*max(index[,i], na.rm=T)))
		lines(years, pred_index[,i], col=plot.colors[i], lwd=2)
		log.ob.min <- log(index[,i])-1.96*sigma[,i]
		log.ob.max <- log(index[,i])+1.96*sigma[,i]
    y.min <- min(c(log.ob.min,log(pred_index[,i]))[is.finite(c(log.ob.min,log(pred_index[,i])))], na.rm=T)
    y.max <- 1.1*max(c(log.ob.max,log(pred_index[,i]))[is.finite(c(log.ob.max,log(pred_index[,i])))], na.rm=T)
		plot(years, log(index[,i]), type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Ln(Index)", ylim=c(y.min, y.max))
		lines(years, log(pred_index[,i]), col=plot.colors[i], lwd=2)
		arrows(years, log.ob.min, years, log.ob.max, length=0)
		#title (paste0("Index ",i), outer=T, line=-1)
		title (index.name.plt, outer=T, line=-1)
		plot(years, log_stdres[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Log-scale Std. Residual")
		abline(h=0)
		hist(log_stdres[,i], plot=T, xlab="Std. Residual", ylab="Probability Density", freq=F, main=NULL)
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

#revised

plot.NAA.4.panel <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, use.s, use.r, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
  years <- mod$years
  years_full = mod$years_full
  dat = mod$env$data
  pred_NAA = mod$rep$pred_NAA
  NAA = mod$rep$NAA
  sigma_all = exp(mod$parList$log_NAA_sigma) #n_stocks x n_ages
  stocks <- 1:dat$n_stocks
  if(!missing(use.s)) stocks <- use.s
  regions <- 1:dat$n_regions
  if(!missing(use.r)) regions <- use.r
  ages <- 1:dat$n_ages
  if(!missing(use.i)) ages <- use.i
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=length(ages)) #mypalette(dat$n_ages)
  sfns <- chartr(" ", "_", mod$input$stock_names)
  rfns <- chartr(" ", "_", mod$input$region_names)
	for(s in stocks) for(r in regions) for (i in ages) if(dat$NAA_where[s,r,i])
	{
    sigma = matrix(sigma_all[s,r,], length(years_full), dat$n_ages, byrow = TRUE)
    log_stdres = (log(NAA[s,r,,])-log(pred_NAA[s,r,,]))/sigma
		if(do.tex) cairo_pdf(file.path(od, paste0("NAA_4panel_", sfns[s], "_", rfns[r], "_age_",i,".pdf")), 
      family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("NAA_4panel_", sfns[s], "_", rfns[r], "_age_",i,'.png')), 
      width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
    y.max <- max(NAA[s,r,,i], na.rm = TRUE)
    if(!is.na(y.max)) {
      plot(years_full, NAA[s,r,,i], type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Abundance (1000s)",
        ylim=c(0, 1.1*y.max))
      lines(years_full, pred_NAA[s,r,,i], col=plot.colors[i], lwd=2)
      if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
    } else plot(years_full, NAA[s,r,,i], type='n', col=plot.colors[i], pch=1, xlab="Year", ylab="Abundance (1000s)")
    log.ob.min <- log(NAA[s,r,,i])-1.96*sigma[,i]
		log.ob.max <- log(NAA[s,r,,i])+1.96*sigma[,i]
    y.max <- max(log.ob.max,log(pred_NAA[s,r,,i]), na.rm = TRUE)
    y.min <- min(log.ob.min,log(pred_NAA[s,r,,i]), na.rm = TRUE)
    if(!is.na(y.max) & !is.na(y.min)) {
      plot(years_full, log(NAA[s,r,,i]), type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Ln(Abundance)",
        ylim=c(y.min,y.max))
      lines(years_full, log(pred_NAA[s,r,,i]), col=plot.colors[i], lwd=2)
      arrows(years_full, log.ob.min, years_full, log.ob.max, length=0)
      if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
      plot(years_full, log_stdres[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Log-scale (Conditional) Std. Residual")
      abline(h=0)
      hist(log_stdres[,i], plot=T, xlab="(Conditional) Std. Residual", ylab="Probability Density", freq=F, main=NULL)
    } else {
      plot(years_full, log(NAA[s,r,,i]), type='n', col=plot.colors[i], pch=1, xlab="Year", ylab="Ln(Abundance)")
    }
  	title (paste0("Conditional Expected and Posterior Estimates of Age ",i, " Abundance "), outer=T, line=-1)
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

#revised

plot.NAA.res <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  #ages = mod$ages.lab
  dat <- mod$env$data
  ages <- 1:dat$n_ages
	n.ages <- length(ages)
  stocks <- 1:dat$n_stocks
  regions <- 1:dat$n_regions

  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=length(ages)) #mypalette(n.ages)
  years = mod$years_full
	n.yrs <- length(years)
  pred_NAA = mod$rep$pred_NAA
  NAA = mod$rep$NAA
  sigma_all = exp(mod$parList$log_NAA_sigma) #n_stocks x n_regions x n_ages
  x.at <- seq(1,n.yrs,5)
  x.lab <- years[x.at]
  sfns <- chartr(" ", "_", mod$input$stock_names)
  rfns <- chartr(" ", "_", mod$input$region_names)
  for(s in stocks) for(r in regions){
    sigma = matrix(sigma_all[s,r,], length(years), length(ages), byrow = TRUE)
    log_stdres = (log(NAA[s,r,,])-log(pred_NAA[s,r,,]))/sigma
    ymin <- min(apply(log_stdres,1, function(x) sum(x[which(x<0)])))
    ymax <- max(apply(log_stdres,1, function(x) sum(x[which(x>0)])))
    dat <- as.data.frame.table(log_stdres)
    if(do.tex) cairo_pdf(file.path(od, paste0("NAA_res_barplot_stacked_", sfns[s], "_", rfns[r],".pdf")), 
      family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("NAA_res_barplot_stacked_", sfns[s], "_", rfns[r],".png")), 
      width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mfrow=c(1,1), mar=c(5,5,1,1), oma = c(0,0,0,0))
    naa_res_plot <- lattice::barchart(Freq ~ Var1, data = dat, groups = Var2, stack = TRUE, col = plot.colors, xlab = "Year", ylab = "Std. Abundance Residuals", box.ratio = 10, reference = TRUE,
      scales = list(x=list(at = x.at, labels = x.lab), alternating = FALSE),
      key = list(text = list(lab = mod$ages.lab), rectangles = list(col = plot.colors), columns = n.ages, title = "Age"))
    print(naa_res_plot)
    if(do.tex | do.png) dev.off() else par(origpar)
  }
  # par(origpar)
}

#revised

plot.catch.age.comp <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years = mod$years
  ages = 1:mod$env$data$n_ages
  ages.lab = mod$ages.lab
  fleets <- 1:mod$env$data$n_fleets
  if(!missing(use.i)) fleets <- use.i
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=length(fleets)) #mypalette(mod$env$data$n_fleets)
  ffns <- chartr(" ", "_", mod$input$fleet_names)
  rfns <- chartr(" ", "_", mod$input$region_names)
	for (i in fleets)
	{
    acomp.obs = mod$env$data$catch_paa[i,,]
    acomp.pred = mod$rep$pred_catch_paa[i,1:mod$env$data$n_years_model,]
    fleet.name.fn <- paste0(ffns[i], "_", rfns[mod$input$data$fleet_regions[i]])
    fleet.name.plt <- paste0(mod$input$fleet_names[i], " in ", mod$input$region_names[mod$input$data$fleet_regions[i]])
    if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_", fleet.name.fn, ".pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_fleet_", fleet.name.fn, ".png")), width = 10*144, height = 10*144, 
      res = 144, pointsize = 12, family = fontfam)
    par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
    my.title <- fleet.name.plt #paste0("Fleet ", i)
    for (j in 1:length(years))
    {
      plot(1:length(ages), acomp.obs[j,], type='p', col=plot.colors[which(fleets == i)], pch=1, xlab="", ylab="",
        ylim=c(0, 1), axes = FALSE)
      if(j %% 15 == 1)
      {
        title(my.title, outer=TRUE, line=0)
        mtext(side = 2, "Proportion", line = 2, outer = TRUE)
        mtext(side = 1, "Age", line = 2, outer = TRUE)
      }
      if(j %% 15 %in% 1:5) axis(2)
      else axis(2, labels = FALSE)
      if(j %in% seq(5,length(years)+1,5)) axis(1, at = ages, labels = ages.lab)
      else axis(1, labels = FALSE)
      grid()
      box()
      lines(1:length(ages), acomp.pred[j,], col=plot.colors[which(fleets == i)],  lwd=2)
      title(paste("Year = ", years[j], sep=""), outer = FALSE, line = 1)

      # if 5x3 multipanel is full, save png and open new one
      if((j %% 15 == 0) & (do.tex | do.png) & (j < length(years))){
        dev.off()
        if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_", fleet.name.fn,"_",letters[j/15],".pdf")), family = fontfam, 
          height = 10, width = 10)
        if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_", fleet.name.fn,"_",letters[j/15],".png")), width = 10*144, 
          height = 10*144, res = 144, pointsize = 12, family = fontfam)
        par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
      }
    }  #end loop on n_years
    if(length(years) %% 15 != 0) frame()
    if(do.tex | do.png) dev.off() else par(origpar)
	}  #end loop on n_fleets

  # par(origpar)
}

#revised

plot.index.age.comp <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years = mod$years
  ages = 1:mod$env$data$n_ages
  ages.lab = mod$ages.lab
  if(!missing(use.i)) indices <- use.i
  else indices <- 1:mod$env$data$n_indices
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=length(indices)) #mypalette(mod$env$data$n_indices)

	for (i in indices)
	{
    acomp.obs = mod$env$data$index_paa[i,,]
    acomp.pred = mod$rep$pred_IAA[i,1:length(years),] #biomass is accounted for on the cpp side
    acomp.pred = acomp.pred/apply(acomp.pred,1,sum)
    index.name.fn <- paste0(chartr(" ", "_", mod$input$index_names[i]), "_", chartr(" ", "_", mod$input$region_names[mod$input$data$index_regions[i]]))
    index.name.plt <- paste0(mod$input$index_names[i], " in ", mod$input$region_names[mod$input$data$index_regions[i]])
    if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_", index.name.fn,".pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_", index.name.fn,".png")), width = 10*144, height = 10*144, 
      res = 144, pointsize = 12, family = fontfam)
    par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
    my.title <- index.name.plt #paste0("Index ", i)
    for (j in 1:length(years))
    {
      plot(1:length(ages), acomp.obs[j,], type='p', col=plot.colors[which(indices == i)], pch=1, xlab="", ylab="",
        ylim=c(0, 1), axes = FALSE)
      if(j %% 15 == 1)
      {
        title(my.title, outer=TRUE, line=0)
        mtext(side = 2, "Proportion", line = 2, outer = TRUE)
        mtext(side = 1, "Age", line = 2, outer = TRUE)
      }
      if(j %% 15 %in% 1:5) axis(2)
      else axis(2, labels = FALSE)
      if(j %in% seq(5,length(years)+1,5)) axis(1, at = ages, labels = ages.lab)
      else axis(1, labels = FALSE)
      grid()
      box()
      lines(1:length(ages), acomp.pred[j,], col=plot.colors[which(indices == i)],  lwd=2)
      title(paste("Year = ", years[j], sep=""), outer = FALSE, line = 1)

      # if 5x3 multipanel is full, save png and open new one
      if((j %% 15 == 0) & (do.tex | do.png) & (j < length(years))){
        dev.off()
        if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_", index.name.fn,"_",letters[j/15],".pdf")), family = fontfam, 
          height = 10, width = 10)
        if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_", index.name.fn,"_",letters[j/15],".png")), width = 10*144, 
          height = 10*144, res = 144, pointsize = 12, family = fontfam)
        par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
      }
    }  #end loop on n_years
    if(length(years) %% 15 != 0) frame()
		if(do.tex | do.png) dev.off() else par(origpar)
	}  #end loop on n_indices
  # par(origpar)
}

#revised

pearson.fn = function(mod, index = NULL, fleet = NULL, age_comp_mod= 1, sims = NULL) {
  dat <- mod$env$data
  rep <- mod$rep
  if(!is.null(index)) {
    age_comp_mod <- dat$age_comp_model_indices[index]
    paa <- dat$index_paa[index,,,drop = FALSE]
    paa[which(dat$use_index_paa[,index] == 0),] = NA
    paahat <- rep$pred_index_paa[index,1:dat$n_years_model,,drop = FALSE]
    paa_pars <- mod$parList$index_paa_pars[index,]
    N <- dat$index_Neff[,index]
  } else if(!is.null(fleet)){
    age_comp_mod <- dat$age_comp_model_fleets[fleet]
    paa <- dat$catch_paa[fleet,,,drop = FALSE]
    paa[which(dat$use_catch_paa[,fleet] == 0),] = NA
    paahat <- rep$pred_catch_paa[fleet,1:dat$n_years_model,,drop = FALSE]
    paa_pars <- mod$parList$catch_paa_pars[fleet,]
    N <- dat$catch_Neff[,fleet]
  }
  if(age_comp_mod %in% 3:7) paa[which(paa<1e-10)] <- NA #no zeros
  r <- paa - paahat
  if(age_comp_mod == 1) var <- paahat*(1-paahat)/N #multinomial
  if(age_comp_mod == 2) var <- paahat*(1-paahat) * (N + exp(paa_pars[1]))/((1 + exp(paa_pars[1])) * N) #D-M regular
  if(age_comp_mod %in% 3:4) var <- paahat*(1-paahat) * exp(paa_pars[1])/(1 + exp(paa_pars[1])) #Dirichlet
  if(age_comp_mod == 11) var <- paahat*(1-paahat) * (1 + exp(paa_pars[1]))/(1 + N * exp(paa_pars[1])) #D-M linear
  #for age_comp_mod 5:7 (Logistic normal) do simulations to determine variance (and mean?)
  if(age_comp_mod %in% 5:7) if(!is.null(sims)){
    if(!is.null(index)) var <- t(sapply(1:length(input$years), \(z) apply(sapply(sims, \(x) x$indx_paa[index,z,]),1,var)))
    else if(!is.null(fleet)) var <- t(sapply(1:length(input$years), \(z) apply(sapply(sims, \(x) x$catch_paa[fleet,z,]),1,var)))
  }
  pearson = r/sqrt(var)
  return(pearson)
}
#mean(x < 0, na.rm = TRUE)
#revised

plot.catch.age.comp.resids <- function(mod, ages, ages.lab, scale.catch.bubble2 = 2, pos.resid.col = "#ffffffaa", neg.resid.col = "#8c8c8caa",
  do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, osa = FALSE, use.i, od) {

  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- dat$n_ages
  years = mod$years
	nyrs <- length(years)
	if(!missing(use.i)) fleets <- use.i
	else fleets <- 1:mod$env$data$n_fleets
  tylab <- "Year"
  ffns <- chartr(" ", "_", mod$input$fleet_names)

	for (i in fleets) {
    yind = which(dat$use_catch_paa[,i] ==1)
    if(length(yind)){
      if(osa & "catchpaa" %in% mod$osa$type){
        df = subset(mod$osa, type == "catchpaa")
        my.title <- "Age Comp OSA Quantile Residuals for "
        fname = paste0("Catch_age_comp_osa_resids_", ffns[i])
        resids = matrix(NA, nrow = dat$n_years_model, ncol = dat$n_ages)
        vals = resids
        for(j in yind){
          tmp = subset(df, year == j & fleet == paste0("fleet_",i))
          resids[j,tmp$age] = tmp$residual
          vals[j,tmp$age] = tmp$val
          if(dat$age_comp_model_fleets[i] %in% c(1:2,10,11)) vals[j,tmp$age]/sum(vals[j,tmp$age]) #obs are numbers not proportions
        }

        scale.resid.bubble.catch <- 2
      } else{
        acomp.obs = dat$catch_paa[i,,]
        acomp.pred = mod$rep$pred_catch_paa[i,1:length(years),]
        #acomp.pred = aperm(mod$rep$pred_catch_paa[1:length(years),,,drop=FALSE], c(2,1,3))[i,,] #biomass is accounted for on the cpp side
        #acomp.pred = acomp.pred/apply(acomp.pred,1,sum)
        my.title <- "Age Comp Residuals (Observed-Predicted) for "
        resids <- acomp.obs - acomp.pred  # NOTE obs-pred
        resids[dat$use_catch_paa[,i]==0,] = NA # don't plot residuals for catch paa not fit in model
        fname = paste0("Catch_age_comp_resids_", ffns[i])
        scale.resid.bubble.catch <- 25
      }
      if(do.tex) cairo_pdf(file.path(od, paste0(fname,".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0(fname,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
      par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
      z1 <- resids
      if(any(!is.na(resids))) range.resids<-range(abs((as.vector(z1))), na.rm=T)
      else range.resids = c(0,0)

      z3 <- z1 * scale.resid.bubble.catch * scale.catch.bubble2
      resid.col <- matrix(NA, nrow=nyrs, ncol=n_ages)   # set color for residual bubbles
      resid.col <- ifelse(z3 > 0.0, pos.resid.col, neg.resid.col)
      plot(ages, rev(ages),  xlim = c(1, n_ages+1), ylim = c(years[nyrs],(years[1]-2)), xlab = "Age", ylab = tylab,
        type = "n", axes=F)
      axis(1, at= ages, lab=ages.lab)
      axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
      box()
      abline(h=years, col="lightgray")
      segments(x0=ages, y0=rep(years[1],n_ages), x1=ages, y1=rep(years[nyrs],n_ages), col = "lightgray", lty = 1)

      for (j in 1:nyrs) points(ages, rep(years[j], n_ages), cex=abs(z3[j,]), col="black", bg = resid.col[j,],  pch = 21)

      bubble.legend1 <- round(quantile(abs(resids), probs = c(0.05,0.5,0.95), na.rm = TRUE),3)
      bubble.legend2 <- bubble.legend1 * scale.resid.bubble.catch*scale.catch.bubble2
      legend("topright", xpd=T, legend=bubble.legend1, pch=rep(1, 3), pt.cex=bubble.legend2, horiz=T , col='black')
      legend("topleft", xpd=T, legend=c("Neg.", "Pos."), pch=rep(21, 2), pt.cex=3, horiz=T, pt.bg=c(neg.resid.col, pos.resid.col), col="black")
      legend("top", xpd = TRUE, legend = paste("Max(resid)=",round(range.resids[2],2), sep=""), horiz = TRUE)
      title (paste0(my.title,mod$input$fleet_names[i]), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    } #end if any age comp for fleet
	}   #end loop n_fleets
  # par(origpar)
}

#revised

plot.index.age.comp.resids <- function(mod, ages, ages.lab, scale.catch.bubble2 = 2, pos.resid.col = "#ffffffaa", neg.resid.col = "#8c8c8caa",
  do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, osa=FALSE, use.i, od) {

  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- dat$n_ages
  years = mod$years
	nyrs <- length(years)
	if(!missing(use.i)) indices <- use.i
	else indices <- 1:dat$n_indices
  ifns <- chartr(" ", "_", mod$input$index_names)

	for (i in indices) {
    yind = which(dat$use_index_paa[,i] ==1)
    if(length(yind)){
      if(osa & "indexpaa" %in% mod$osa$type){
        df = subset(mod$osa, type == "indexpaa")
        my.title <- "Age Comp OSA Quantile Residuals for "
        fname = paste0("Catch_age_comp_osa_resids_",ifns[i])
        resids = matrix(NA, nrow = dat$n_years_model, ncol = dat$n_ages)
        vals = resids
        for(j in yind){
          tmp = subset(df, year == j & fleet == paste0("index_",i))
          resids[j,tmp$age] = tmp$residual
          vals[j,tmp$age] = tmp$val
          if(dat$age_comp_model_indices[i] %in% c(1:2,10,11)) vals[j,tmp$age]/sum(vals[j,tmp$age]) #obs are numbers not proportions
        }
        scale.resid.bubble.catch <- 2
      } else {
        acomp.obs = dat$index_paa[i,,]
        acomp.pred = mod$rep$pred_index_paa[i,1:length(years),]
        #acomp.pred = aperm(mod$rep$pred_IAA[1:length(years),,,drop=FALSE], c(2,1,3))[i,,] #biomass is accounted for on the cpp side
        #acomp.pred = acomp.pred/apply(acomp.pred,1,sum)
        my.title <- "Age Comp Residuals (Observed-Predicted) for "
        resids <- acomp.obs - acomp.pred  # NOTE obs-pred
        resids[dat$use_index_paa[,i]==0,] = NA # don't plot residuals for index paa not fit in model
        fname = paste0("Catch_age_comp_resids_",ifns[i])
        scale.resid.bubble.catch <- 25
      }
      tylab <- "Year"
      
      if(do.tex) cairo_pdf(file.path(od, paste0(fname,".pdf")), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0(fname,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
      par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
      z1 <- resids
      if(any(!is.na(resids))) range.resids<-range(abs((as.vector(z1))), na.rm=T)
      else range.resids <- c(0,0)

      z3 <- z1 * scale.resid.bubble.catch * scale.catch.bubble2
      resid.col <- matrix(NA, nrow=nyrs, ncol=n_ages)   # set color for residual bubbles
      resid.col <- ifelse(z3 > 0.0, pos.resid.col, neg.resid.col)
      plot(ages, rev(ages),  xlim = c(1, n_ages+1), ylim = c(years[nyrs],(years[1]-2)), xlab = "Age", ylab = tylab,
        type = "n", axes=F)
      axis(1, at= ages, lab=ages.lab)
      axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
      box()
      abline(h=years, col="lightgray")
      segments(x0=ages, y0=rep(years[1],n_ages), x1=ages, y1=rep(years[nyrs],n_ages), col = "lightgray", lty = 1)

      for (j in 1:nyrs) if(dat$use_index_paa[j,i] == 1)
        points(ages, rep(years[j], n_ages), cex=abs(z3[j,]), col="black", bg = resid.col[j,],  pch = 21)

      bubble.legend1 <- round(quantile(abs(resids), probs = c(0.05,0.5,0.95), na.rm = TRUE),3)
      bubble.legend2 <- bubble.legend1 * scale.resid.bubble.catch*scale.catch.bubble2
      legend("topright", xpd=T, legend=bubble.legend1, pch=rep(1, 3), pt.cex=bubble.legend2, horiz=T , col='black')
      legend("topleft", xpd=T, legend=c("Neg.", "Pos."), pch=rep(21, 2), pt.cex=3, horiz=T, pt.bg=c(neg.resid.col, pos.resid.col), col="black")
      legend("top", xpd = TRUE, legend = paste("Max(resid)=",round(range.resids[2],2), sep=""), horiz = TRUE)
      title (paste0(paste0(my.title,mod$input$index_names[i])), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    } #end if any age comp for index
	}   #end loop n_fleets
  # par(origpar)
}

#revised

plot.sel.blocks <- function(mod, ages, ages.lab, plot.colors, indices = FALSE, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, od)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
	cc<-0
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = ages
	years <- mod$years
	if(!missing(use.i)) fleets <- use.i
	else {
    fleets <- 1:mod$env$data$n_fleets
    if(indices) fleets <- 1:mod$env$data$n_indices
  }
  if(!indices) {
    fleet_names <- mod$input$fleet_names
    fleet_regions <- mod$input$region_names[mod$input$data$fleet_regions]
    sb_p = dat$selblock_pointer_fleets #selblock pointer by year and fleet
  } else {
    fleet_names <- mod$input$index_names
    fleet_regions <- mod$input$region_names[mod$input$data$index_regions]
    sb_p <- dat$selblock_pointer_indices
  }
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=length(unique(sb_p))) #mypalette(length(unique(sb_p)))
	for (i in fleets) {
    fn <- paste0("Selectivity_", chartr(" ", "_", fleet_names[i]), "_", chartr(" ", "_", fleet_regions[i]))
    pn <- paste0(fleet_names[i], " in ", fleet_regions[i])
	  if(do.tex) cairo_pdf(file.path(od, paste0(fn,".pdf")), family = fontfam, height = 10, width = 10)
	  if(do.png) png(filename = file.path(od, paste0(fn,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
	  blocks = unique(sb_p[,i])
		n.blocks <- length(blocks)
    # sel = rbind(mod$rep$selblocks[blocks,])
    sel = do.call(rbind, lapply(mod$rep$selAA, function(x) apply(x,2,mean)))[blocks,,drop=FALSE]
		minyr <- rep(NA, n.blocks)
		maxyr <- rep(NA, n.blocks)
		my.col <- rep(NA, n.blocks)
		for (j in 1:n.blocks) {
			cc<-cc+1
			my.col[j] <- plot.colors[cc]
			minyr[j] <- min(years[sb_p[,i] == blocks[j]])
			maxyr[j] <- max(years[sb_p[,i] == blocks[j]])
			if (j==1)
			{
				plot(ages, sel[j,], type='l', col=my.col[j], xlim=c(min(ages),max(ages)), ylim=c(0,1.1), xlab="Age", ylab="Selectivity",
					lwd=2, axes = FALSE)
				grid(col = gray(0.7), lwd = 2)
				axis(1, at = ages, labels = ages.lab, lwd = 2)
				axis(2, lwd = 2)
				box(lwd=2)
			}
			if (j>1) lines(ages, sel[j,], type='l', col=my.col[j], lwd=2)
		}
		title(pn, line = 1)
		legend("topright", col=my.col, legend=paste0(minyr, " - ", maxyr), lwd=2, bg = "white")
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

#revised

#not used
plot.fleet.sel.blocks <- function(mod, ages, ages.lab, plot.colors, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, od)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
	cc<-0
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = ages
  sb_p = dat$selblock_pointer_fleets #selblock pointer by year and fleet
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=length(unique(sb_p))) #mypalette(length(unique(sb_p)))
	years <- mod$years
	if(!missing(use.i)) fleets <- use.i
	else fleets <- 1:mod$env$data$n_fleets
  ffns <- chartr(" ", "_", mod$input$fleet_names)

	for (i in fleets)
	{
    fn <- paste0("Selectivity_",ffns[i])
	  if(do.tex) cairo_pdf(file.path(od, paste0(fn,".pdf")), family = fontfam, height = 10, width = 10)
	  if(do.png) png(filename = file.path(od, paste0(fn,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
	  blocks = unique(sb_p[,i])
		n.blocks <- length(blocks)
    # sel = rbind(mod$rep$selblocks[blocks,])
    sel = do.call(rbind, lapply(mod$rep$selAA, function(x) apply(x,2,mean)))[blocks,,drop=FALSE]
		minyr <- rep(NA, n.blocks)
		maxyr <- rep(NA, n.blocks)
		my.col <- rep(NA, n.blocks)
		for (j in 1:n.blocks)
		{
			cc<-cc+1
			my.col[j] <- plot.colors[cc]
			minyr[j] <- min(years[sb_p[,i] == blocks[j]])
			maxyr[j] <- max(years[sb_p[,i] == blocks[j]])
			if (j==1)
			{
				plot(ages, sel[j,], type='l', col=my.col[j], xlim=c(min(ages),max(ages)), ylim=c(0,1.1), xlab="Age", ylab="Selectivity",
					lwd=2, axes = FALSE)
				grid(col = gray(0.7), lwd = 2)
				axis(1, at = ages, labels = ages.lab, lwd = 2)
				axis(2, lwd = 2)
				box(lwd=2)
			}
			if (j>1) lines(ages, sel[j,], type='l', col=my.col[j], lwd=2)
		}
		title(paste0(mod$input$fleet_names[i]), line = 1)
		legend("topright", col=my.col, legend=paste0(minyr, " - ", maxyr), lwd=2, bg = "white")
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

#not used
plot.index.sel.blocks <- function(mod, ages, ages.lab, plot.colors, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, use.i, od)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
	cc<-0
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = ages
  sb_p = dat$selblock_pointer_indices #selblock pointer by year and index
  if(missing(plot.colors)) plot.colors = mypalette(length(unique(sb_p)))
	years <- mod$years
	if(!missing(use.i)) indices <- use.i
	else indices <- 1:mod$env$data$n_indices

	for (i in indices)
	{
	  if(do.tex) cairo_pdf(file.path(od, paste0("Selectivity_index",i,".pdf")), family = fontfam, height = 10, width = 10)
	  if(do.png) png(filename = file.path(od, paste0("Selectivity_index",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
	  blocks = unique(sb_p[,i])
		n.blocks <- length(blocks)
    # sel = rbind(mod$rep$selblocks[blocks,])
    sel = do.call(rbind, lapply(mod$rep$selAA, function(x) apply(x,2,mean)))[blocks,,drop=FALSE]
		minyr <- rep(NA, n.blocks)
		maxyr <- rep(NA, n.blocks)
		my.col <- rep(NA, n.blocks)
		for (j in 1:n.blocks)
		{
			cc<-cc+1
			my.col[j] <- plot.colors[cc]
			minyr[j] <- min(years[sb_p[,i] == blocks[j]])
			maxyr[j] <- max(years[sb_p[,i] == blocks[j]])
			if (j==1)
			{
				plot(ages, sel[j,], type='l', col=my.col[j], xlim=c(min(ages),max(ages)), ylim=c(0,1.1), xlab="Age", ylab="Selectivity",
					lwd=2, axes = FALSE)
				grid(col = gray(0.7), lwd = 2)
				axis(1, at = ages, labels = ages.lab, lwd = 2)
				axis(2, lwd = 2)
				box(lwd=2)
			}
			if (j>1) lines(ages, sel[j,], type='l', col=my.col[j], lwd=2)
		}
		title(paste0("Index ",i), line = 1)
		legend("topright", col=my.col, legend=paste0(minyr, " - ", maxyr), lwd=2, bg = "white")
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

#revised

plot.SSB.F.trend<-function(mod, alpha = 0.05)
{
  origpar <- par(no.readonly = TRUE)
  years_full <- mod$years_full
  years <- mod$years
  n_ages = mod$env$data$n_ages
  n_regions <- mod$env$data$n_regions

  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  ssb <- ssb.lo <- ssb.hi <- matrix(NA, nrow = length(years_full), ncol = mod$env$data$n_stocks)
  faa <- faa.lo <- faa.hi <- array(NA, dim = dim(mod$rep$log_FAA_by_region))
  full.f <- full.f.lo <- full.f.hi <- matrix(NA, length(years_full), n_regions)
  if(class(mod$sdrep)[1] == "sdreport"){
    std <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = TRUE), TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = TRUE))
    ssb[] <- exp(std[[1]]$log_SSB)
    ssb.lo[] <- ssb*exp(qnorm(alpha/2)*std[[2]]$log_SSB)
    ssb.hi[] <- ssb*exp( -qnorm(alpha/2)*std[[2]]$log_SSB)
    faa[] <- exp(std[[1]]$log_FAA_by_region)
    faa.lo[] <- faa*exp(qnorm(alpha/2)*std[[2]]$log_FAA_by_region)
    faa.hi[] <- faa*exp(-qnorm(alpha/2)*std[[2]]$log_FAA_by_region)
    for(r in 1:n_regions) {
      full.f.lo[,r] <- apply(faa.lo[r,,],1, function(x) x[max(which(x == max(x)))])
      full.f.hi[,r] <- apply(faa.hi[r,,],1, function(x) x[max(which(x == max(x)))])
    }
  } else {
    ssb[] <- mod$rep$SSB
    faa[] <- mod$rep$FAA_tot
    #std = mod$sdrep
  }
	for(r in 1:n_regions) full.f[,r] <- apply(faa[r,,],1, function(x) x[max(which(x == max(x)))])
	
  par(mfrow=c(2,1), mar=c(1,1,1,1), oma = c(4,4,0,0))
  max.y <- max(ssb.hi)
  na.se <- is.na(max.y)
  if(na.se) max.y <- max(ssb)
  pal <- viridisLite::viridis(n=mod$env$data$n_stocks)
  plot(years_full, ssb[,1], type='n', lwd=2, xlab="", ylab="", ylim=c(0,max.y), axes = FALSE)
  axis(1, labels = FALSE)
  axis(2)
  box()
  mtext(side = 2, "SSB (mt)", outer = FALSE, line = 3)
  grid(col = gray(0.7))
  for(s in 1:mod$env$data$n_stocks){
    lines(years_full, ssb[,s], col = pal[s], lwd = 2)
    if(!na.se) polygon(c(years_full,rev(years_full)), c(ssb.lo[,s],rev(ssb.hi[,s])), col = adjustcolor(pal[s], alpha.f=0.4), border = "transparent")
  }
  if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
  if(mod$env$data$n_stocks>1) legend("topright", lty = 1, col = pal, lwd = 2, legend = mod$input$stock_names)
  
  # F trend
  max.y <- max(full.f.hi)
  na.se <- is.na(max.y)
  if(na.se) max.y <- max(full.f)
  pal <- viridisLite::viridis(n_regions)
  #if(!no.f.ci){ # have CI
  plot(years_full, full.f[,1], type='n', lwd=2, xlab="", ylab="", ylim=c(0,max.y), axes = FALSE)
  axis(1)
  axis(2)
  box()
  mtext(side = 1, "Year", outer = FALSE, line = 3)
  mtext(side = 2, "Fully-selected F", outer = FALSE, line = 3)
  grid(col = gray(0.7))
  for(i in 1:n_regions){
    lines(years_full, full.f[,i], col = pal[i], lwd = 2)
  	polygon(c(years_full,rev(years_full)), c(full.f.lo[,i],rev(full.f.hi[,i])), col = adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
  }
  if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
  if(mod$env$data$n_regions>1) legend("topright", lty = 1, col = pal, lwd = 2, legend = mod$input$region_names)
  par(origpar)
}  #end function

#revised

plot.SSB.AA <- function(mod, ages, ages.lab, plot.colors, prop=FALSE, stock = 1)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n.ages <- length(ages)
  if(missing(plot.colors)) plot.colors = mypalette(n.ages)
	years<-  mod$years
  years_full <-  mod$years_full
	n.yrs <- length(years_full)
  ssbfrac = dat$fracyr_SSB
  ssb.aa <- mod$rep$NAA_spawn[stock,,] * mod$rep$waa_ssb[stock,,] * mod$rep$mature_all[stock,,]/1000
	#ssb.aa <- (mod$rep$NAA * exp(-ssbfrac * (mod$rep$FAA_tot + mod$rep$MAA)) * dat$waa[dat$waa_pointer_ssb,,] * dat$mature)/1000
	ssb.max <- max(apply(ssb.aa,1,sum))

	par(mfrow=c(1,1), mar=c(5,5,3,1), oma = c(0,0,1,0))
	if(!prop){ # plot SSB at age
  	res <- barplot(t(ssb.aa), beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year', ylab =paste('SSB at age (', "kmt", ')', sep = ''),
  		ylim = c(0,1.15*ssb.max), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
    if(length(years_full) > length(years)) abline(v=length(years), lty=2, lwd=2)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8)
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years_full[seq(5,n.yrs,5)])
    box()
	}
	if(prop){ # plot *proportion* SSB at age
  	barplot(t(ssb.aa/apply(ssb.aa,1,sum)), beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year', ylab ='Proportion SSB at age',
  		ylim = c(0,1.1), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
    if(length(years_full) > length(years)) abline(v=length(years), lty=2, lwd=2)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8)
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years_full[seq(5,n.yrs,5)])
    box()
	}
  title(mod$input$stock_names[stock], line = 1)

  par(origpar)
}  #end funciton
#plot.SSB.AA(ssm)

#------------------------------------
plot.NAA <- function(mod, ages, ages.lab, plot.colors, scale = 1000, units = expression(10^6), prop=FALSE, stock = 1, region = 1)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
  dat = mod$env$data
	## stacked barplot of NAA
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n.ages <- length(ages)
  if(missing(plot.colors)) plot.colors = mypalette(n.ages)
  years<-  mod$years
  years_full <-  mod$years_full
  n.yrs <- length(years_full)
	NAA <- mod$rep$NAA[stock,region,,]
	N.max=max(apply(NAA,1,sum))/scale
	par(mfrow=c(1,1), mar=c(5,5,3,1), oma = c(0,0,1,0))
	if(!prop){ # plot numbers at age
  	barplot(t(NAA)/scale, beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year',
  		ylab =as.expression(substitute(paste("January 1 numbers at age (", units, ")", sep = ''), list(units = units[[1]]))),
  		ylim = c(0,1.15*N.max), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
    if(length(years_full) > length(years)) abline(v=length(years), lty=2, lwd=2)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8)
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years_full[seq(5,n.yrs,5)])
    box()
	}
	if(prop){ # plot *proportion* of numbers at age
  	barplot(t(NAA/apply(NAA,1,sum)), beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year',
  		ylab ='January 1 proportions at age', ylim = c(0,1.1), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
    if(length(years_full) > length(years)) abline(v=length(years), lty=2, lwd=2)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8 )
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years_full[seq(5,n.yrs,5)])
    box()
	}
  title(paste(mod$input$stock_names[stock],mod$input$region_names[region]), line = 1)
  par(origpar)
} # end function
#plot.NAA(ssm)
#------------------------------------
#revised

#scatter plot of SSB, R with 2-digit year as symbol (lag by 1 year)
plot.recr.ssb.yr <- function(mod, ssb.units = "kmt", recruits.units = expression(10^6), alpha = 0.05,
  scale.ssb = 1000, scale.recruits = 1000, age.recruit = 1, plot.colors, loglog=FALSE, stock = 1)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1), mar = c(4,5,3,1), oma = c(1,1,1,1))

  # std <- list(TMB:::as.list.sdreport(mod$sdrep, what = "Est", report = T),
  #     TMB:::as.list.sdreport(mod$sdrep, what = "Std", report = T))

  std <- summary(mod$sdrep, "report")
  cov <- mod$sdrep$cov
  dat = mod$env$data
  years = mod$years
  nyrs <- length(years)
  nages <- dat$n_ages
  #nstates <- model$dimensions$n_states
  #nprojyrs <- model$dimensions$n_years_proj
	ssb.ind <- matrix(which(rownames(std) == "log_SSB"), length(mod$years_full), dat$n_stocks)[,stock]
	log.ssb <- std[ssb.ind,1]
	ssb.cv <- std[ssb.ind,2]
  log.ssb.lo <- log.ssb + qnorm(alpha/2)*ssb.cv
  log.ssb.hi <- log.ssb + qnorm(1-alpha/2)*ssb.cv
  R.ind = array(which(rownames(std) == "log_NAA_rep"),dim = dim(mod$rep$NAA))[stock, dat$spawn_regions[stock],,1]
  #NAA.ind = which(rownames(std) == "log_NAA_rep")
  log.R = std[R.ind,1]
  R.cv = std[R.ind,2]
	#ssb.R.cov <- cov[c(ssb.ind,R.ind),c(ssb.ind,R.ind)]

  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=nyrs-age.recruit)# mypalette(nyrs-age.recruit)
	SR <- matrix(NA, (nyrs-age.recruit), 3)
	SR[,1] <- years[1:(nyrs-age.recruit)]
	SR[,2] <- exp(log.ssb[1:(nyrs-age.recruit)])/scale.ssb
	SR[,3] <- exp(log.R[age.recruit +1:(nyrs-age.recruit)])/scale.recruits
	yr.text <- substr(years[1:(nyrs-age.recruit)],3,4)
	npts <- nyrs-age.recruit
	ci.regs <- lapply(1:(nyrs-age.recruit), function(x)
	{
	  tcov <- cov[c(ssb.ind[x],R.ind[x+age.recruit]),c(ssb.ind[x],R.ind[x+age.recruit])]
	  tsd <- std[c(ssb.ind[x],R.ind[x+age.recruit]),2]
    tcor = tcov/(tsd %*% t(tsd))
	  el <- exp(ellipse::ellipse(tcor, level = 1-alpha, scale = tsd, centre = c(log(SR[x,2]),log(SR[x,3]))))
	  return(el)
	})
  # if(missing(max.y)) max.y <- max(sapply(ci.regs, function(x) max(x[,2])))
  # if(missing(max.x)) max.x <- max(sapply(ci.regs, function(x) max(x[,1])))
  max.y <- max(sapply(ci.regs, function(x) max(x[,2])))
  max.x <- max(sapply(ci.regs, function(x) max(x[,1])))
  na.lims <- any(is.na(c(max.y,max.x))) # check for na lims

	if(!loglog){ # untransformed SSB and Rec
    if(!na.lims){
    	plot(SR[,2], SR[,3], type='n', col='black',
    		xlab=as.expression(substitute(paste("SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))),
    		ylab= as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
    			list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), ylim=c(0, max.y), xlim=c(0,max.x), axes = FALSE)
    	grid(col = gray(0.7), lty = 2)
    	axis(1)
    	axis(2)
    	box()
    	lines(SR[,2], SR[,3], col = gray(0.7), lwd =2)
    	for(i in 1:length(ci.regs))
      {
        tcol <- col2rgb(plot.colors[i])
        poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
        polygon(ci.regs[[i]][,1],ci.regs[[i]][,2], border = poly)
      }
      points(SR[npts,2], SR[npts,3], pch=19, col="#ffaa22", cex=2.5)
      text(SR[,2], SR[,3], yr.text, cex=0.9, col=plot.colors)
    } else { # NA lims
      max.y <- max(SR[,3])
      max.x <- max(SR[,2])
      plot(SR[,2], SR[,3], type='n', col='black',
        xlab=as.expression(substitute(paste("SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))),
        ylab= as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
          list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), ylim=c(0, max.y), xlim=c(0,max.x), axes = FALSE)
      grid(col = gray(0.7), lty = 2)
      axis(1)
      axis(2)
      box()
      lines(SR[,2], SR[,3], col = gray(0.7), lwd =2)
      points(SR[npts,2], SR[npts,3], pch=19, col="#ffaa22", cex=2.5)
      text(SR[,2], SR[,3], yr.text, cex=0.9, col=plot.colors)
    }
	}

	if(loglog){ # log(SSB) and log(Rec)
    if(!na.lims){
    	plot(log(SR[,2]), log(SR[,3]), type='n', col='black',
    		xlab=as.expression(substitute(paste("Log-SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))),
    		ylab= as.expression(substitute(paste("Age-", age.recruit, " Log-Recruits (", units, ")", sep = ''),
    			list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), ylim=c(log(min(SR[,3])), log(max.y)), xlim=c(log(min(SR[,2])),log(max.x)), axes = FALSE)
    	grid(col = gray(0.7), lty = 2)
    	axis(1)
    	axis(2)
    	box()
    	lines(log(SR[,2]), log(SR[,3]), col = gray(0.7), lwd =2)
    	for(i in 1:length(ci.regs))
      {
        tcol <- col2rgb(plot.colors[i])
        poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
        polygon(log(ci.regs[[i]][,1]),log(ci.regs[[i]][,2]), border = poly)
      }
      points(log(SR[npts,2]), log(SR[npts,3]), pch=19, col="#ffaa22", cex=2.5)
      text(log(SR[,2]), log(SR[,3]), yr.text, cex=0.9, col=plot.colors)
    } else { # na lims but correct max.x and max.y already calculated for untransformed plot SSB-Rec
      max.y <- max(SR[,3])
      max.x <- max(SR[,2])
      plot(log(SR[,2]), log(SR[,3]), type='n', col='black',
        xlab=as.expression(substitute(paste("Log-SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))),
        ylab= as.expression(substitute(paste("Age-", age.recruit, " Log-Recruits (", units, ")", sep = ''),
          list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), ylim=c(log(min(SR[,3])), log(max.y)), xlim=c(log(min(SR[,2])),log(max.x)), axes = FALSE)
      grid(col = gray(0.7), lty = 2)
      axis(1)
      axis(2)
      box()
      lines(log(SR[,2]), log(SR[,3]), col = gray(0.7), lwd =2)
      points(log(SR[npts,2]), log(SR[npts,3]), pch=19, col="#ffaa22", cex=2.5)
      text(log(SR[,2]), log(SR[,3]), yr.text, cex=0.9, col=plot.colors)
    }
	}
  title(paste0(mod$input$stock_names[stock], " recruitment in region ", mod$input$region_names[dat$spawn_regions[stock]]), line = 1)
  par(origpar)
}  #end function

# revised

#------------------------------------
plot.SARC.R.SSB <- function(mod, scale.ssb=1, scale.recruits=1, age.recruit = 1, ssb.units = 'mt', recruits.units = expression(10^3), stock = NULL)
{
  origpar <- par(no.readonly = TRUE)
  par(mar = c(5,5,3,5), oma = c(0,0,0,1), family='serif')
  years = mod$years
  years_full = mod$years_full
  nyrs <- length(years_full)
  dat <- mod$input$data
  if(is.null(stock) & dat$n_stocks> 1) {
    ssb <- apply(mod$rep$SSB,1, sum)
    R <- mod$rep$NAA[1,dat$spawn_regions[1],,1]
    for(i in 2:dat$n_stocks) R <- R + mod$rep$NAA[i,dat$spawn_regions[i],,1]
  } else {
    ssb <- mod$rep$SSB[,stock]
    R <- mod$rep$NAA[stock, dat$spawn_regions[stock],,1]
  }
  ssb.plot <- ssb[1:(nyrs-age.recruit)]/scale.ssb
  recr.plot <- R[age.recruit + 1:(nyrs-age.recruit)]/scale.recruits
  yr.text <- substr(years_full,3,4)
  plot.colors <- viridisLite::viridis(n=2)
  max.r <- max(recr.plot)
  max.ssb <- max(ssb.plot)
  scale.r <- max(ssb.plot)/max(recr.plot)
  ylimr <- c(0,1.1*max(recr.plot))
  barplot(recr.plot/scale.recruits, axisnames=FALSE, width=1, space=rep(0,nyrs-age.recruit), offset=rep(-0.5,nyrs-age.recruit), axes=FALSE, xpd=FALSE,
    xlab = '', ylab ='', ylim = ylimr, xlim=c(0.5,nyrs-age.recruit - 0.5), col=plot.colors[1])
  xr <-pretty(c(0,recr.plot/scale.recruits))
  axis(2, at = xr, lab = xr )
  axis(side=1, las=2, at=seq(0.5,nyrs-age.recruit-0.5, by=2),
  labels=as.character(seq(years_full[1],years_full[nyrs-age.recruit], by=2)), cex=0.75, las=2)

  y.ssb <- (ssb.plot)*max.r/max.ssb
  lines(seq(0.5,nyrs-age.recruit-0.5, by=1), y.ssb, lwd=2, col = plot.colors[2])
  x <- pretty(c(0,ssb.plot))
  axis(4, at = c(0,x*max.r/max.ssb), lab = c(0,x))#, col=plot.colors[2], col.axis=plot.colors[2])
  box()
  mtext(side = 1, 'Year', line = 3)
  mtext(side = 4, as.expression(substitute(paste("SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))), line = 3)#, col=plot.colors[2])
  mtext(side = 2, as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
    list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), line = 3)
  if(length(years_full) > length(years)) abline(v=length(years)-age.recruit, lty=2, lwd=2)
  legend("topleft", fill = plot.colors, border = plot.colors, legend = c("Recruits", "SSB"))
  if(is.null(stock) & dat$n_stocks>1) title("Total SSB and Recruitment", line = 1)
  else title(paste0(mod$input$stock_names[stock], " in ", mod$input$region_names[dat$spawn_regions[stock]]), line = 1)
  par(origpar)
}  # end function
#plot.SARC.R.SSB(ssm, ssm.aux)
#revised

plot.fleet.F <- function(mod, plot.colors)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
  years = mod$years
  years_full = mod$years_full
  nyrs <- length(years_full)
	n_fleets <- mod$input$data$n_fleets
  FAA <- mod$rep$FAA
  F <- matrix(NA, length(years_full), n_fleets)
  for(y in 1:length(years_full)) for(f in 1:n_fleets) F[y,f] <- FAA[f,y,which(FAA[f,y,] == max(FAA[f,y,]))[1]]
  std <- NULL
  # if(class(mod$sdrep)[1] == "sdreport"){
  #   std = summary(mod$sdrep)
  # } else {
  #   std = mod$sdrep
  # }
  # faa.ind <- which(rownames(std) == "log_FAA_tot")
  # log.faa <- matrix(std[faa.ind,1], nyrs, mod$env$data$n_ages)

  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=n_fleets) # mypalette(n_fleets)
  plot(years_full, F[,1], xlab="Year", ylab="Full F", ylim=c(0,max(F)),	type='n', lty=1, lwd=2)
	if(n_fleets == 1){ 
    lines(years_full, F, lty=1, lwd=2, col = "black")
  } else { 
    for(i in 1:n_fleets){
  		if(i==1){
        if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
      }
      lines(years_full, F[,i],lty=1, lwd=2, col=plot.colors[i])
    }
  }
  grid(col = gray(0.7))
	leg.names <- paste0(mod$input$fleet_names, " in ", mod$input$region_names[mod$input$data$fleet_regions])
	legend('topright', legend= leg.names, col=plot.colors,lwd=2, lty=1, bty='n')
	par(origpar)
}   # end function
#plot.fleet.F(ssm,ssm.aux)

#revised
#------------------------------------
plot.cv <- function(mod)
{
  origpar <- par(no.readonly = TRUE)
	par(mfrow=c(1,1), mar=c(4,4,2,2))
  years = mod$years
  years_full = mod$years_full
  nyrs <- length(years_full)
  dat <- mod$env$data
	if(class(mod$sdrep)[1] == "sdreport"){
    std = summary(mod$sdrep)
  } else {
    std = mod$sdrep
  }
  rep_est <- as.list(mod$sdrep, "Est", report=T)
  rep_std <- as.list(mod$sdrep, "Std", report=T)
  #ssb.ind <- matrix(which(rownames(std) == "log_SSB"), nyrs, dat$n_stocks)
	log.ssb <- rep_est$log_SSB #matrix(std[ssb.ind,1], ncol = dat$n_stocks)
	ssb.cv <- rep_std$log_SSB #matrix(std[ssb.ind,2], ncol = dat$n_stocks)
  #NAA.ind = array(which(rownames(std) == "log_NAA_rep"),dim = dim(mod$rep$NAA))
  log.R <- R.cv <- matrix(NA, nyrs, dat$n_stocks)
  for(s in 1:dat$n_stocks) {
    log.R[,s] <- rep_est$log_NAA_rep[s,dat$spawn_regions[s],,1] #[NAA.ind,1], dim = dim(mod$rep$NAA))[s, dat$spawn_regions[s],,1]
    R.cv[,s] <- rep_std$log_NAA_rep[s,dat$spawn_regions[s],,1] #array(std[NAA.ind,2], dim = dim(mod$rep$NAA))[s, dat$spawn_regions[s],,1]
  }
  FAA <- exp(rep_est$log_FAA)
  FAA.cv <- rep_std$log_FAA
  #F.ind = matrix(which(rownames(std) == "log_F"), nyrs, ncol = dat$n_fleets)
  F.cv <- matrix(NA, nyrs, dat$n_fleets)
  for(i in 1:dat$n_fleets) {
    ages <- apply(FAA[i,,],1, function(x) which(x==max(x))[1])
    F.cv[,i] <- FAA.cv[i,,][cbind(1:nyrs, ages)]
  }

  any.na <- any(is.na(c(R.cv, ssb.cv, F.cv)))

  if(!any.na){
    plot.colors = viridisLite::viridis(n=dat$n_stocks + dat$n_fleets) 
  	plot(years_full, R.cv[,1], type='n', xlab="Year", ylab="CV", ylim=c(0, 1.1*max(R.cv, ssb.cv, F.cv)))
    for(s in 1:dat$n_stocks) {
      lines(years_full, R.cv[,s], lty = 1, lwd=2, col=plot.colors[s])
  	  lines(years_full, ssb.cv[,s], , lty = 2, lwd=2, col=plot.colors[s])
    }
    for(f in 1:dat$n_fleets){
  	  lines(years_full, F.cv[,f], lty = 3, lwd=2, col=plot.colors[dat$n_stocks + f])
      if(length(years_full) > length(years)) abline(v=tail(years,1), lwd=1)
    }
  	labs <- paste0(rep(mod$input$stock_names, each = dat$n_stocks), " ", rep(c("Recruits", "SSB"),dat$n_stocks))
    labs <- c(labs, paste0(mod$input$fleet_names, " F"))
    cols <- c(rep(1:dat$n_stocks, each = 2), dat$n_stocks + 1:dat$n_fleets)
    legend('topleft', legend=labs, col=plot.colors[cols], lty=c(rep(1:2, dat$n_stocks), rep(3,dat$n_fleets)), lwd=2, border = "transparent")
  }
  par(origpar)
}  # end function
#------------------------------------
#revised

#------------------------------------
plot.M <- function(mod, ages, ages.lab, alpha = 0.05, plot.colors, stock = 1, region = 1)
{
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- length(ages)
	meanMAA <- apply(mod$rep$MAA[stock,region,,],2,mean)
  years = mod$years
  years_full = mod$years_full
  n_years = length(years_full)
	if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=n_years) #mypalette(n_years)
	n.M.by.age <- lapply(1:n_ages, function(x) table(mod$rep$MAA[stock,region,,x]))
	plot(ages,meanMAA,lwd=2,xlab="Age",ylab="Natural Mortality Rate", ylim=c(0,1.1*(max(mod$rep$MAA[stock,region,,]))), type= 'n', xlim = c(min(ages)-0.5,max(ages)+1),
		axes=FALSE)
	grid(col = gray(0.7), lwd = 2)
	axis(1, at = ages, labels = ages.lab, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	sapply(1:n_ages, function(x)
	{
		y <- sort(unique(mod$rep$MAA[stock,region,,x]))
		# ind1 <- sapply(y, function(z) which(mod$rep$MAA[,x] == y)[1])
		segments(x-0.2, y, x+0.2, y, lwd = 2,col = plot.colors[1:length(y)])
		text(x+0.2, y, paste('n =', table(mod$rep$MAA[stock,region,,x])), pos = 4)
	})
  title(paste0(mod$input$stock_names[stock], " in ", mod$input$region_names[region]), line = 1)
}

# revised
#------------------------------------
#--------Data Plots------------------
plot.catch.by.fleet <- function(mod, units = "mt", plot.colors)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow = c(1,1))
  dat = mod$env$data
  years = mod$years
  nyrs = length(years)
	catch.obs <- dat$agg_catch
	n_fleets <- dat$n_fleets
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=n_fleets) #mypalette(n_fleets)
	# barplot(t(catch.obs), xlab="Year", ylab= paste0("Catch (", units, ")"), ylim=c(0,1.1*max(apply(catch.obs,1,sum))), col=plot.colors,space=0)
	# axis(side=1, at = seq(2,nyrs,2)-0.5, labels = years[seq(2,nyrs,2)], cex=0.75)
	# box(lwd = 2)
	# if (n_fleets > 1)
  # {
    #legend('top', legend=paste0("Fleet ",1:n_fleets), horiz=TRUE, pch=15, col=plot.colors)

    # do proportions only if n_fleets > 1
		catch.prop <- catch.obs/apply(catch.obs,1,sum)
		barplot(t(catch.prop), xlab="Year", ylab="Proportion of Catch", ylim=c(0,1.1), col=plot.colors, space=0)
    axis(side=1, las=2, at = seq(2,nyrs,2)-0.5, labels = years[seq(2,nyrs,2)], cex=0.75, las=2)
    box(lwd = 2)
		legend('top', legend= paste0(mod$input$fleet_names, " in ", mod$input$region_names[mod$input$data$fleet_regions]), 
      horiz=TRUE, pch=15, col=plot.colors)
	# }
	par(origpar)
}
#revised

# Bubble plots of catch age comps (set is.catch.flag to False to plot Discard age comps)
plot.catch.age.comp.bubbles <- function(mod, ages, ages.lab, bubble.col = "#8c8c8caa", i=1)
{
  dat = mod$env$data
  years = mod$years
  nyrs = length(years)
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
  n_ages = length(ages)
  n_fleets = dat$n_fleets
  acomp.obs <- dat$catch_paa[i,,]
  catch.yrs <- which(dat$use_catch_paa[,i] == 1)
  
  my.title <- paste0("Age Comps for Catch for ", mod$input$fleet_names[i], " in ", mod$input$region_names[mod$input$data$fleet_regions[i]])
  origpar <- par(no.readonly = TRUE)
  par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
  scale.catch.obs <- 5
  z3 <- as.matrix(acomp.obs) * scale.catch.obs

  plot(ages, rev(ages),  xlim = range(ages), ylim = c(years[nyrs],(years[1]-2)), xlab = "Age", ylab = "", type = "n", axes=FALSE)
  axis(1, at= ages, lab = ages.lab)
  axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
  box()
  abline(h=years, col="lightgray")
  segments(x0=ages, y0=rep(years[1],n_ages), x1=ages, y1=rep(years[nyrs],n_ages), col = "lightgray", lty = 1)
  if (length(catch.yrs)>0) for (j in 1:nyrs) points(ages, rep(years[j], n_ages), cex=z3[j,], col="black", bg = bubble.col, pch = 21)

  bubble.legend1 <- c(0.05,0.2,0.4)
  bubble.legend2 <- bubble.legend1 * scale.catch.obs
  legend("topright", xpd=TRUE, legend=bubble.legend1, pch=rep(21, 3), pt.cex=bubble.legend2, horiz=T , col='black', pt.bg = bubble.col)
  title (my.title, outer=TRUE, line=-1)
  par(origpar)
}

#revised

#------------------------------------
plot.index.input <- function(mod, plot.colors)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(2,1), mar = c(1,1,1,1), oma = c(4,4,0,15))
  years = mod$years
  nyrs = length(years)
  dat = mod$env$data
	indvals <- dat$agg_indices
	indvals[which(dat$use_indices!=1)] <- NA
	n_indices = dat$n_indices
	# rescale to mean 1 and stdev 1
	rescaled <- indvals
	my.mean <- apply(indvals,2,mean, na.rm=TRUE)
	my.std <- apply(indvals,2,sd, na.rm=TRUE)
	for (i in 1:n_indices) rescaled[,i] <- (indvals[,i] - my.mean[i]) / my.std[i]

	my.range <- range(rescaled, na.rm=TRUE)
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=n_indices) #mypalette(n_indices)
	plot(years,rescaled[,1],xlab="",ylab="", axes = FALSE, xlim = range(years), ylim=my.range,col=plot.colors[1],type='n')
	grid(col = gray(0.7))
	axis(1, labels = FALSE)
	axis(2)
	box()
	mtext(side = 2, "Rescaled Indices", outer = FALSE, line = 3)
	for (i in 1:n_indices) lines(years,rescaled[,i],col=plot.colors[i])
  leg <- paste0(mod$input$index_names, " ", mod$input$region_names[mod$input$data$index_regions])
  legend("right", legend = leg, col = plot.colors, lty = 1, horiz = FALSE, xpd = NA, inset = c(-0.5,0), bty = "n")
  #legend("top", legend = mod$input$index_names, col = plot.colors, lty = 1, horiz = TRUE, xpd = NA, inset = c(0,-0.1), bty = "n")

	# now repeat on log scale
	log.indvals <- log(indvals)
	log.rescaled <- log.indvals
	my.log.mean <- apply(log.indvals,2,mean, na.rm=T)
	my.log.std <- apply(log.indvals,2,sd, na.rm=T)
	for (i in 1:n_indices) log.rescaled[,i] <- (log.indvals[,i] - my.log.mean[i]) / my.log.std[i]

	my.log.range <- range(log.rescaled, na.rm=T)
	plot(years,log.rescaled[,1], xlab="",ylab="",ylim=my.log.range,col=plot.colors[1],type='n')
	grid(col = gray(0.7))
	for (i in 1:n_indices) lines(years,log.rescaled[,i],col=plot.colors[i])
	mtext(side = 2, "Rescaled Log(Indices)", outer = FALSE, line = 3)
	mtext(side = 1, "Year", outer = FALSE, line = 3)
	par(origpar)
}
#plot.index.input(ssm)
#revised

#------------------------------------
# Bubble plots of index age comps
plot.index.age.comp.bubbles <- function(mod, ages, ages.lab, bubble.col = "#8c8c8caa", i=1)
{
  years = mod$years
  nyrs = length(years)
  dat = mod$env$data
  n_indices = dat$n_indices
  index.yrs <- which(dat$use_index_paa[,i] == 1)
  acomp.obs <- dat$index_paa[i,,]
  my.title <- paste0("Age Comps for ", mod$input$index_names[i], " in ", mod$input$region_names[mod$input$data$index_regions[i]])
  origpar <- par(no.readonly = TRUE)
  par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
  n_ages <- length(ages)

  scale.index.obs <- 5
  z3 <- as.matrix(acomp.obs) * scale.index.obs

  plot(ages, rev(ages),  xlim = range(ages), ylim = c(years[nyrs],(years[1]-2)),
  xlab = "Age", ylab = "", type = "n", axes=FALSE)
  axis(1, at= ages, lab=ages.lab)
  axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
  box()
  abline(h=years, col="lightgray")
  segments(x0=ages, y0=rep(years[1],n_ages), x1=ages, y1=rep(years[nyrs],n_ages), col = "lightgray", lty = 1)
  if(length(index.yrs)>0) for (j in 1:nyrs) points(ages, rep(years[j], n_ages), cex=z3[j,], col="black", bg = bubble.col, pch = 21)

  bubble.legend1 <- c(0.05,0.2,0.4)
  bubble.legend2 <- bubble.legend1 * scale.index.obs
  legend("topright", xpd=TRUE, legend=bubble.legend1, pch=rep(21, 3), pt.cex=bubble.legend2, horiz=TRUE, col='black', pt.bg = bubble.col)
  title (my.title, outer=T, line=-1)
  par(origpar)
}

#revised
#------------------------------------
plot.waa <- function(mod,type="ssb",plot.colors,ind=1)
{
  origpar <- par(no.readonly = TRUE)
  if(type %in% c("ssb", "fleets")) years = mod$years_full
  else years = mod$years
  nyrs = length(years)
  dat = mod$env$data
  if(missing(plot.colors)) plot.colors = viridisLite::viridis(n=dat$n_ages) #mypalette(dat$n_ages)
  point = switch(type,
    ssb = dat$waa_pointer_ssb[ind],
    #jan1 = dat$waa_pointer_jan1,
    fleets = dat$waa_pointer_fleets[ind],
    indices = dat$waa_pointer_indices[ind]
    #totcatch = dat$waa_pointer_totcatch
  )
  waa = switch(type,
    ssb = mod$rep$waa_ssb[ind,,],
    #jan1 = dat$waa_pointer_jan1,
    fleets = mod$rep$waa_catch[ind,,],
    indices = dat$waa[dat$waa_pointer_indices[ind],,]
    #totcatch = dat$waa_pointer_totcatch
  )
  labs = switch(type,
    ssb = paste0(mod$input$stock_names[ind], "SSB"),
    #jan1 = "January 1 Biomass",
    fleets = mod$input$fleet_names[ind],# paste0("Fleet ", ind),
    indices = mod$input$index_names[ind], #paste0("Index ", ind),
    # fleets = paste0("Fleet ", 1:dat$n_fleets),
    # indices = paste0("Index ", 1:dat$n_indices),
    #totcatch = "Total Catch"
  )
  #waa = dat$waa[point,,]
  n = ifelse(length(dim(waa)) == 2, 1, dim(waa)[1])
	for(i in 1:n)
	{
		if(n>1) WAA.plot <- waa[i,,]
    else WAA.plot = waa
		plot(years,years,xlab="Year",ylab="Weight",ylim=c(0,max(WAA.plot)),type='n')
		for (a in 1:dat$n_ages)
		{
			lines(years,WAA.plot[,a],col=plot.colors[a],lwd=2)
			lines(years,rep(mean(WAA.plot[,a]),length(years)),lty=2,col=plot.colors[a])
		}
		title(main = paste0("Annual Weight-at-Age for ", labs[i]), line = 1)
	}  # end k-loop
	par(origpar)
}  # end function

#revised

#------------------------------------
plot.maturity <- function(mod, ages.lab, plot.colors, stock = 1)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  years = mod$years
  n_years = length(years)
  ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
  mature <- dat$mature[stock,,] #matrix
	meanmaturity <- apply(mature,2,mean)
	if(missing(plot.colors)) plot.colors <- viridisLite::viridis(n=n_years) #mypalette(n_years)

	plot(ages,meanmaturity,type='l',lwd=2,xlab="Age",ylab="Maturity",ylim=c(0,max(mature)), axes = FALSE)
  axis(1, at = ages, labels = ages.lab, lwd = 2)
  axis(2, lwd = 2)
  box(lwd = 2)
	if (length(unique(mature)) > length(ages))
	{
		for (i in 1:n_years) points(jitter(ages, factor=0.4), mature[i,],col=plot.colors[i])
    midi <- floor(n_years/2)
		legend('topleft', horiz=FALSE, legend=c(years[1],years[midi],years[n_years]), pch=c(1,1,1), col=c(plot.colors[1], plot.colors[midi], plot.colors[n_years]))
	}
	title(main=paste0("Maturity for ", mod$input$stock_names[stock]), outer=FALSE, line = 1)
	par(origpar)
}
#revised

get_P.fn = function(time, age, year, stock, season, fleet_seasons, fleet_regions, can_move, FAA, MAA, mu, L, mig_type)
{
  #print("in get_P")
  #print(n_regions)
  #print(age)
  #print(year)
  #print(season)
  #print(stock)
  #print(cum_n_mu)
  #print(dim(n_mu))
  #stop()
  #print(mu_row)
  #print(mu_col)
  #print(mu_pointer)
  #print(mig_type)
  #print(time)
  #print(F)
  #print(M)
  n_fleets <- dim(FAA)[1]
  n_regions <- dim(MAA)[2]
  dim <- n_regions + n_fleets + 1
  P <- matrix(0,dim,dim)
  F <- rep(0,n_fleets)
  F[which(fleet_seasons[,season]==1)] <- FAA[which(fleet_seasons[,season]==1),year, age]
  M <- Z <- MAA[stock,,year,age] #n_regions
  Z <- Z + L
  #for(f in 1:n_fleets) Z[fleet_regions[f]] = Z[fleet_regions[f]] + F[f]
  for(f in 1:n_fleets) Z[fleet_regions[f]] = Z[fleet_regions[f]] + F[f]
  if(n_regions == 1){
    if(time < 1e-15){
      P[1,1] <- 1
    } else{
      P[1,1] <- exp(-Z[1]*time)
    }
    for(f in 1:n_fleets) P[1,1+f] <- F[f] * (1 - P[1,1])/ Z[1]
    P[1,2+n_fleets] <- M[1] * (1 - P[1,1]) / Z[1]
  } else { #more than one region: movement is possible
    if(sum(can_move[stock,season,])>0){ #migration is happening
      if(mig_type == 0) { #migration is instantaneous after survival
        #probs of survival and mortality and then moving between regions
        if(time < 1e-15) { #prob of survival is 1 when interval is zero
          P[1:n_regions,1:n_regions] <-  mu
        } else { 
          for(f in 1:n_fleets) P[fleet_regions[f], n_regions + f] <- F[f] * (1 - exp(-Z[fleet_regions[f]] * time))/Z[fleet_regions[f]] #caught
          for(r in 1:n_regions) {
            #Zr <- sum(F[which(fleet_regions==r)]) + M[r]
            P[r,1:n_regions] <- exp(-Z[r] * time)  * mu[r,]
            P[r,dim] <- (M[r] + L[r]) * (1 - exp(-Z[r] * time))/Z[r] #other dead
          }
        }
        P[cbind((n_regions+1):dim,(n_regions+1):dim)] <- 1 #if dead, must stay dead
      }
      if(mig_type == 1) { #migration occurs continuously during interval
        if(time < 1e-15) {
          P <- diag(dim)
        } else {
          A <- matrix(0,dim,dim)
          A[,dim] <- M + L
          #for(i in 1:n_mu[stock,year,season,age]) A[mu_row[cum_n_mu + i],mu_col[cum_n_mu + i]] = exp(mu[mu_pointer[cum_n_mu + i]]); #log of transition intensities
          for(i in 1:n_regions) for(j in 1:n_regions) if(i!=J) A[i,j] <- mu[i,j]
          for(f in 1:n_fleets) A[fleet_regions[f],n_regions + f] <- F[f]
          A[cbind(fleet_regions, n_regions + 1:n_fleets)] = F #caught
          diag(A) = -apply(A,1,sum) #hazard
          eigs.A = eigen(A)
          P = eigs.A$vec %*% diag(exp(eigs.A$val * time)) %*% solve(eigs.A$vec)
        }
      }
    } else { #no migration
      if(time < 1e-15) P = diag(dim)
      else {
        for(f in 1:n_fleets) P[fleet_regions[f], n_regions + f] <- F[f] * (1 - exp(-Z[fleet_regions[f]] * time))/Z[fleet_regions[f]] #caught
        #probs of survival and mortality and then moving between regions
        for(r in 1:n_regions) {
          #Zr <- sum(F[which(fleet_regions==r)]) + M[r]
          P[r,r] <- exp(-Z[r] * time)
          P[r,dim] <- (M[r] + L[r]) * (1 - exp(-Z[r] * time))/Z[r] #other dead
        }
        P[cbind((n_regions+1):dim,(n_regions+1):dim)] <- 1 #if dead, must stay dead
      }
    }
  }
  return(P)
}


#-------SSB/R -----------------------------
# get_SPR <- function(mod){
#   dat <- mod$env$data

#   dat$fracyr_seasons, dat$spawn_seasons, dat$fleet_seasons, dat$fleet_regions, dat$can_move, dat$mig_type
#   FAA, MAA, mu, L, )
# }
get_SPR = function(F, M, sel, mat, waassb, fracyrssb, at.age = FALSE)
{
  n_ages = length(sel)
  SPR = numeric()
  n = 1
  F = F * sel #n_fleets x n_ages
  if(!is.null(dim(sel))) if(dim(sel)[1]>1) F = apply(F,2,sum) # total F
  Z = F + M
  for(a in 1:(n_ages-1))
  {
    SPR[a] = n[a] * mat[a] * waassb[a] * exp(-fracyrssb * Z[a])
    n[a+1] = n[a] * exp(-Z[a])
  }
  n[n_ages] = n[n_ages]/(1-exp(-Z[n_ages]))
  SPR[n_ages] = n[n_ages] * mat[n_ages] * waassb[n_ages] * exp(-fracyrssb * Z[n_ages])
  if(at.age) return(SPR)
  else return(sum(SPR))
}

#-------Y/R -----------------------------
get_YPR = function(F, M, sel, waacatch, at.age = FALSE)
{
  n_ages = length(sel)
  YPR = numeric()
  n = 1
  F = F * sel
  Z = F + M
  for(a in 1:(n_ages-1))
  {
    YPR[a] = n[a] * F[a] * waacatch[a] * (1.0 - exp(-Z[a]))/Z[a]
    n[a+1] = n[a] * exp(-Z[a])
  }
  n[n_ages] = n[n_ages]/(1 - exp(-Z[n_ages]))
  YPR[n_ages] = n[n_ages] * F[n_ages] * waacatch[n_ages] * (1.0 - exp(-Z[n_ages]))/Z[n_ages]
  if(at.age) return(YPR)
  else return(sum(YPR))
}

#-------Y/R -----------------------------
get_YPR_fleets = function(FAA,M, waacatch, at.age = FALSE)
{
  n_ages = dim(FAA)[2]
  n_fleets <- dim(FAA)[1]
  YPR = matrix(0,n_fleets, n_ages)
  #F = F * sel
  Z = apply(FAA,2,sum) + M
  for(f in 1:n_fleets) {
    n = 1
    for(a in 1:(n_ages-1)){
      YPR[f,a] <- n[a] * FAA[f,a] * waacatch[f,a] * (1.0 - exp(-Z[a]))/Z[a]
      n[a+1] = n[a] * exp(-Z[a])
    }
    n[n_ages] = n[n_ages]/(1 - exp(-Z[n_ages]))
    YPR[f,n_ages] = n[n_ages] * FAA[f,n_ages] * waacatch[f,n_ages] * (1.0 - exp(-Z[n_ages]))/Z[n_ages]
  }
  if(at.age) return(YPR) #by fleet and age
  else return(apply(YPR,1,sum)) #by fleet
}

get_SPR_BRPS_fn <- function(mod, spr_yrs, percent){
  dat = mod$env$data
  if(missing(percent)) percent <- dat$percentSPR
  if(missing(spr_yrs)) spr_yrs <- dat$avg_years_ind+1 #c++
  R_yrs <- dat$XSPR_R_avg_yrs+1
  R_type <- dat$XSPR_R_opt
  if(R_type %in% c(1,3)) Rxspr <- mean(mod$rep$NAA[Ryrs,1])
  else Rxspr <- mean(mod$rep$pred_NAA[R_yrs,1])
  n_ages<- dat$n_ages
  years <- mod$years
  n_years <- length(years)
  n_fleets <- dat$n_fleets

  mat.age <- apply(dat$mature[1,,][spr_yrs,,drop = FALSE],2,mean)
  ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb[1],,][spr_yrs,,drop = FALSE],2,mean)
  catch.waa <- apply(dat$waa[dat$waa_pointer_fleets,spr_yrs,,drop = FALSE],c(1,3),mean)
  M.age <- apply(mod$rep$MAA[1,1,,][spr_yrs,,drop = FALSE],2,mean)
  FAA <- apply(mod$rep$FAA[,spr_yrs,,drop=FALSE], c(1,3), mean) #mean FAA by fleet
  FAAtot = apply(FAA,2,sum) #sum across fleet averages 
  seltot <- FAAtot/max(FAAtot)
  selAA <- FAA/max(FAAtot) 
  spawn.time <- mean(dat$fracyr_SSB[spr_yrs])
  spr0 = get_SPR(F=0, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
  F.start <- 0.11  # starting guess for optimization routine to find F_SPR%
  spr.f <- function(F.start) {
    spr = get_SPR(F=F.start, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
    abs(100*spr/spr0 - percent)
  }
  opt <- suppressWarnings(nlminb(start=F.start, objective=spr.f, lower=0, upper=10))
  Fxspr <- opt$par
  spr_Fxspr <- get_SPR(Fxspr, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
  FAA_xspr <- Fxspr * selAA
  ypr_Fxspr <- get_YPR_fleets(FAA = FAA_xspr, M=M.age, waacatch= catch.waa)
  Y_Fxspr <- Rxspr * ypr_Fxspr
  SSB_Fxspr <- Rxspr*spr_Fxspr
  return(list(Fxspr = Fxspr, FAA_xspr = FAA_xspr, SSB_Fxspr = SSB_Fxspr, Y_Fxspr = Y_Fxspr, spr_Fxspr = spr_Fxspr, Rxspr = Rxspr))
}

#y <- get_SPR_BRPS_fn(x)


plot.SPR.table <- function(mod, nyrs.ave = 5, plot=TRUE)
{
  origpar <- par(no.readonly = TRUE)
  spr.targ.values <- seq(0.2, 0.8, 0.05)
	n.spr <- length(spr.targ.values)
  dat = mod$env$data
	n_ages<- dat$n_ages
	years <- mod$years
	n_years <- length(years)
  avg.ind <- tail(1:n_years, nyrs.ave)
  # avg.ind = (n_years-nyrs.ave+1):n_years
	#fec.age <- apply(dat$waa[dat$waa_pointer_ssb,,][avg.ind,],2,mean)
  mat.age <- apply(dat$mature[1,,][avg.ind,,drop = FALSE],2,mean)
  ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb[1],,][avg.ind,,drop = FALSE],2,mean)
  catch.waa <- apply(dat$waa[dat$waa_pointer_fleets,avg.ind,,drop = FALSE], c(1,3),mean)
	#catch.waa <- apply(dat$waa[dat$waa_pointer_totcatch,,][avg.ind,],2,mean)
  M.age <- apply(mod$rep$MAA[1,1,,][avg.ind,,drop = FALSE],2,mean)
  FAA <- apply(mod$rep$FAA[,avg.ind,,drop=FALSE], c(1,3), mean) #mean FAA by fleet
  FAAtot = apply(FAA,2,sum) #sum across fleet averages 
  seltot <- FAAtot/max(FAAtot)
  selAA <- FAA/max(FAAtot) 
  
  # sel = apply(apply(mod$rep$FAA[avg.ind,,,drop=FALSE], 2:3, mean),2,sum) #sum across fleet averages 
  #sel = apply(mod$rep$FAA_tot[avg.ind,],2,mean) #average FAA, then do selectivity
	# sel <- sel/max(sel)
	spawn.time <- mean(dat$fracyr_SSB[avg.ind,1])
  spr0 = get_SPR(F=0, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
	F.start <- 0.11  # starting guess for optimization routine to find F_SPR%

	f.spr.vals <- rep(NA, n.spr)
	ypr.spr.vals <- rep(NA, n.spr)
	conv.vals <- rep(NA, n.spr)

	for (i in 1:n.spr)
	{
		t.spr <- spr.targ.values[i]

		spr.f <- function(F.start)
		{
      spr = get_SPR(F=F.start, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
			abs(spr/spr0 - t.spr)
		}
		yyy <- suppressWarnings(nlminb(start=F.start, objective=spr.f, lower=0, upper=3))
		f.spr.vals[i] <- yyy$par
    ypr.spr.vals[i] = sum(get_YPR_fleets(FAA = f.spr.vals[i]*selAA, M=M.age, waacatch= catch.waa))
	}  #end i-loop over SPR values

	spr.target.table<- as.data.frame(cbind(spr.targ.values, f.spr.vals, ypr.spr.vals))
	colnames(spr.target.table) <- c("%SPR", "F(%SPR)", "YPR")
	par(mfrow=c(1,1), mar=c(4,4,2,4))

	if(plot){ # plot, not table
  	plot(spr.targ.values, ypr.spr.vals, type='n', xlab="% SPR Target", ylab="", lwd=2, col="blue3",
      ylim=c(0,1.2*max(ypr.spr.vals)), axes = FALSE)
    box(lwd = 2)
    axis(1, lwd = 2)
  	abline(v=seq(0.2,0.8, by=0.05), col="grey85")
  	lines(spr.targ.values, ypr.spr.vals, lwd=2, col="blue3" )
  	points(spr.targ.values, ypr.spr.vals, pch=19, col="blue3" )

    scale.f.spr <- max(f.spr.vals)/max(ypr.spr.vals)
  	axis(side=2, #at=seq(0,1,by=0.1)/scale.f.spr, lab=seq(0,1,by=0.1),
      las=2, col='blue3', col.axis="blue3", lwd = 2)
    mtext(side = 2, "Yield per Recruit", line = 3, col = "blue3")


  	lines(spr.targ.values, f.spr.vals/scale.f.spr, lwd = 2, col="red")
  	points(spr.targ.values, f.spr.vals/scale.f.spr, pch = 19, col="red")
  	axis(side=4, at=seq(0,1,by=0.1)/scale.f.spr, lab=seq(0,1,by=0.1), las=2, col='red', col.axis="red", lwd = 2)
  	mtext(side=4, "F (%SPR)", line=3, col="red")

  	title (paste("SPR Target Reference Points (Years Avg = ", nyrs.ave,")", sep=""), outer=T, line=-1)
	}

	if(!plot){ # table, not plot
  	par(mfrow=c(1,1), mar=c(2,2,2,2))
  	plot(seq(1,15), seq(1,15), type='n', axes=F, bty='n',xlab="",ylab="")

  	text(x=2,y=14, labels="% SPR", font=2, pos=4)
  	text(x=5, y=14, labels="F(%SPR)" , font=2, pos=4)
  	text(x=9, y=14, labels="YPR", font=2, pos=4 )
  	for (i in 1:n.spr)
  	{
  		text(x=2, y=seq(n.spr,1, by=-1), labels=round(spr.targ.values,2), cex=1.0, pos=4, font=1)
  		text(x=5, y=seq(n.spr,1, by=-1), labels=round(f.spr.vals,4), cex=1.0, pos=4, font=1)
  		text(x=9, y=seq(n.spr,1, by=-1), labels=round(ypr.spr.vals,4), cex=1.0, pos=4, font=1)
  	}
  	title (paste("SPR Target Reference Points (Years Avg = ", nyrs.ave,")", sep=""), outer=TRUE, line=-1)
	}
  par(origpar)
	#write.csv( spr.target.table, file=paste(od,"SPR_Target_Table.csv", sep=""),  row.names=F)
} # end function

#------------------------------------
plot.annual.SPR.targets <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  spr.targ.values <- seq(0.2, 0.5, by=0.1)
	n.spr <- length(spr.targ.values)
  dat = mod$env$data
  n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)
  years_full = mod$years_full
  n_years_full = length(years_full)

  fec.age <- mod$rep$waa_ssb[1,,]
	mat.age <- mod$rep$mature_all[1,,]
	wgt.age <- mod$rep$waa_catch
	M.age <- mod$rep$MAA[1,1,,]
	# sel.age <- mod$rep$FAA_tot/apply(mod$rep$FAA_tot,1,max)
  FAA <- mod$rep$FAA[,,,drop=FALSE] #mean FAA by fleet (nf, ny, na)
  FAAtot = apply(FAA,c(2,3),sum) #annual total FAA (ny x na)
  Ftot <- apply(FAAtot,1,max)
  seltot <- FAAtot/max(FAAtot) #by column
  selAA <- FAA
  for(f in 1:dat$n_fleets) selAA[f,,] <- selAA[f,,]/Ftot

	spawn.time <- dat$fracyr_SSB #vector length n_years

	spr0 <- numeric()
	F.start <- 0.11  # starting guess for optimization routine to find F_SPR%

	f.spr <- matrix(NA, n_years_full, n.spr)
	ypr.spr <- matrix(NA, n_years_full, n.spr)
	conv <- matrix(NA, n_years_full, n.spr)

	for (j in 1:n_years_full)
	{
    spr0[j] = get_SPR(F=0, M=M.age[j,], sel=selAA[,j,], mat=mat.age[j,], waassb=fec.age[j,], fracyrssb = spawn.time[j])
    #spr0.vals[j] <- s.per.recr(n_ages=n_ages, fec.age=fec.age[j,], mat.age=mat.age[j,], M.age= M.age[j,], F.mult=0, sel.age=sel.age[j,],
    #  spawn.time=spawn.time)
	  for (i in 1:n.spr)
	  {
			t.spr <- spr.targ.values[i]

			spr.f <- function(F.start)
			{
        spr = get_SPR(F=F.start, M=M.age[j,], sel=selAA[,j,], mat=mat.age[j,], waassb=fec.age[j,], fracyrssb = spawn.time[j])
        abs(spr/spr0[j] - t.spr)
				#abs(s.per.recr(n_ages=n_ages, fec.age=fec.age[j,], mat.age=mat.age[j,], M.age= M.age[j,], F.mult=F.start, sel.age=sel.age[j,],
				#	spawn.time=spawn.time)/spr0.vals[j] - t.spr)
			}
			yyy <- suppressWarnings(nlminb(start=F.start, objective=spr.f, lower=0, upper=3))
			f.spr[j,i] <- yyy$par
			ypr.spr[j,i] = sum(get_YPR_fleets(FAA = f.spr[j,i]*rbind(selAA[,j,]), M=M.age[j,], waacatch= rbind(wgt.age[,j,])))
      #ypr(n_ages, wgt.age=wgt.age[j,], M.age=M.age[j,],  F.mult=f.spr.vals[j,i], sel.age=sel.age[j,])
		}  # end j-loop over n_years
	}  #end i-loop over SPR values

  plot.colors = mypalette(n.spr)
  if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_annual_time.pdf")), family = fontfam, height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("FSPR_annual_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
	par(mfrow=c(2,1), mar=c(4,5,2,1))
	plot(years_full, f.spr[,1], type='n', xlab="Years", ylab=expression(paste("Full ",italic(F)["%SPR"])), lwd=2, ylim=c(0,1.2*max(f.spr)))
	for (i in 1:n.spr) lines(years_full, f.spr[,i], lwd=2, col=plot.colors[i])
  if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
	legend('top', legend= sapply(spr.targ.values, function(x) as.expression(bquote(italic(F)[paste(.(x*100),"%")]))), col=plot.colors, horiz=TRUE, lwd=2, cex=0.9)
	title (main=expression(paste("Annual ", italic(F)["%SPR"], " Reference Points")), line=1)
	#if (save.plots) savePlot(paste(od, "Annual_FSPR.", plotf, sep=''), type=plotf)
	plot(years_full, ypr.spr[,1], type='n', xlab="Years", ylab=expression(paste("YPR(",italic(F)["%SPR"],")")), lwd=2, ylim=c(0,1.2*max(ypr.spr)))
	for (i in 1:n.spr) lines(years_full, ypr.spr[,i], lwd=2, col=plot.colors[i])
	if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2, lwd=1)
  #legend('top', legend=c("YPR20%", "YPR30%", "YPR40%", "YPR50%"), col=plot.colors, horiz=TRUE, lwd=2, cex=0.9)
	title (main=expression(paste("Annual YPR(",italic(F)["%SPR"], ") Reference Points")), line=1)
	if(do.tex | do.png) dev.off() else par(origpar)

	if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_freq_annual_F.pdf")), family = fontfam, height = 10, width = 10)
	if(do.png) png(filename = file.path(od, paste0("FSPR_freq_annual_F.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
	par(mfrow=c(2,2))
  f.hist = lapply(1:length(spr.targ.values), function(x) hist(f.spr[1:n_years,x], plot = FALSE))
  par(mfrow=c(2,2), mar=c(4,1,2,2), oma=c(2,3,2,0))
  for(i in 1:length(spr.targ.values))
  {
    plot(f.hist[[i]]$mids, f.hist[[i]]$counts, xlab=bquote("Full " ~ italic(F)[paste(.(spr.targ.values[i]*100),"%")]), ylab="", type='h', lwd=2, ylim=c(0, max(f.hist[[i]]$counts)), col=plot.colors[i])
    lines(f.hist[[i]]$mids, f.hist[[i]]$counts, lwd=2, col=plot.colors[i])
  }
  mtext(side=2, outer = TRUE, "Frequency", line = 1)
	title (main=expression(paste("Frequencies of Annual ", italic(F)["%SPR"], " Reference Points")), outer=TRUE, line=0, cex.main = 2)
	if(do.tex | do.png) dev.off() else par(origpar)

	if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_freq_annual_YPR.pdf")), family = fontfam, height = 10, width = 10)
	if(do.png) png(filename = file.path(od, paste0("FSPR_freq_annual_YPR.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
	ypr.hist = lapply(1:length(spr.targ.values), function(x) hist(ypr.spr[1:n_years,x], plot = FALSE))
  par(mfrow=c(2,2), mar=c(4,1,2,2), oma=c(2,3,2,0))
  for(i in 1:length(spr.targ.values))
  {
    plot(ypr.hist[[i]]$mids, ypr.hist[[i]]$counts, xlab=bquote(paste("YPR(", italic(F)[paste(.(spr.targ.values[i]*100),"%")],")")), ylab="", type='h', lwd=2, ylim=c(0, max(ypr.hist[[i]]$counts)), col=plot.colors[i])
    lines(ypr.hist[[i]]$mids, ypr.hist[[i]]$counts, lwd=2, col=plot.colors[i])
  }
  mtext(side=2, outer = TRUE, "Frequency", line = 1)
	title (main=expression(paste("Frequencies of Annual YPR(",italic(F)["%SPR"], ") Reference Points")), outer = TRUE, line=0, cex.main = 2)
	if(do.tex | do.png) dev.off() else par(origpar)

	# par(origpar)
} # end function

plot.SR.pred.line <- function(mod, ssb.units = "mt", SR.par.year, recruits.units = "thousands", scale.ssb = 1,
	scale.recruits = 1, age.recruit = 1, plot.colors, stock = 1) {
  dat <- mod$env$data
  if(dat$recruit_model[stock] %in% (3:4)) {#B-H stock recruit function with alpha/beta, ecov effects and per-recruit inputs from last year
    if(class(mod$sdrep)[1] == "sdreport"){
      std = summary(mod$sdrep)
    } else {
      std = mod$sdrep
    }
    ssb.ind = which(rownames(std) == "log_SSB")
    years = mod$years
    nyrs = length(years)
    log.ssb <- matrix(std[ssb.ind,1], ncol = dat$n_stocks)[,stock]
    R <- mod$rep$NAA[stock, dat$spawn_regions[stock],,1]
    SR <- matrix(NA, (nyrs-age.recruit), 3)
    SR[,1] <- years[1:(nyrs-age.recruit)]
    SR[,2] <- exp(log.ssb[1:(nyrs-age.recruit)])/scale.ssb
    SR[,3] <- R[age.recruit +1:(nyrs-age.recruit)]/scale.recruits
    a.ind = which(rownames(std) == "log_SR_a")
    b.ind = which(rownames(std) == "log_SR_b")
    log_a <- matrix(std[a.ind,1], ncol = dat$n_stocks)[,stock]
    log_b <- matrix(std[b.ind,1], ncol = dat$n_stocks)[,stock]
    if(missing(SR.par.year)) SR.par.year = nyrs
    a.b.ind <- cbind(a.ind,b.ind)[SR.par.year,]
    #a.b.ind = matrix(which(rownames(std) == "mean_rec_pars"),dat$n_stocks,2)[stock,]
    l.ab = std[a.b.ind,1]
    ab.na <- which(is.na(std[a.b.ind,2]))
    a.b.cov = mod$sdrep$cov[a.b.ind,a.b.ind]
    if(length(ab.na)) a.b.cov[ab.na,] <- a.b.cov[,ab.na] <- 0
    seq.ssb <- seq(0, max(SR[,2]), length.out=300)
    if(dat$recruit_model[stock] == 3){#B-H
      lR.fn = function(la, lb, S) la  + log(S) - log(1 + exp(lb)*S)
    } else{ #Ricker
      lR.fn = function(la, lb, S) la  + log(S) - exp(lb)*S
    }
    dlR.dp = Deriv::Deriv(lR.fn, x = c("la","lb"))

    sd.pred.lR = sapply(seq.ssb, function(x) {
      d = dlR.dp(l.ab[1], l.ab[2], S= x)
      return(c(sqrt(t(d) %*% a.b.cov %*% d)))
    })
    pred.lR <- lR.fn(l.ab[1], l.ab[2], seq.ssb)
    #exp(log_a) *seq.ssb/(1 + exp(log_b)*seq.ssb)
    ci.pred.lR = pred.lR + qnorm(0.975)*cbind(-sd.pred.lR, sd.pred.lR) - log(scale.recruits)

    if(missing(plot.colors)) plot.colors = viridisLite::viridis(n = 1)
    tcol <- adjustcolor(plot.colors, alpha.f = 0.4)
    # tcol <- col2rgb(plot.colors[1])
    # tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')

    plot(SR[,2], SR[,3], type='p', col='black', pch=19,
      xlab=as.expression(substitute(paste("SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))),
      ylab= as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
        list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), ylim=c(0, max(SR[,3])), xlim=c(0,1.1*max(SR[,2])))
    lines(seq.ssb, exp(pred.lR), col=plot.colors[1], lwd=2)
    polygon(c(seq.ssb,rev(seq.ssb)), exp(c(ci.pred.lR[,1],rev(ci.pred.lR[,2]))), col = tcol, border = "transparent")
    mtext(paste0("Stock-recruit parameters, per-recruit inputs, and any environmental effects from year: ", mod$years[SR.par.year]), side = 3, outer = F)
  }
}
#revised 

kobe.plot <- function(mod, status.years=NULL, static = FALSE, msy = FALSE, single.plot = TRUE, alpha = 0.05, max.x=NULL, max.y=NULL){
  if(is.null(status.years)) {
    status.years <- length(mod$years)
    if(length(mod$years_full)> status.years) status.years <- c(status.years, length(mod$years_full))
  }
  n_stocks <- mod$env$data$n_stocks
  percentSPR = mod$env$data$percentSPR
  std <- summary(mod$sdrep, "report")
  inds <- list()
  inds$ssb <- which(rownames(std) == "log_SSB_all")
  inds$full.f <- which(rownames(std) == "log_F_tot")
  inds$F.t <- which(rownames(std) == "log_FXSPR") 
  inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR"), ncol = n_stocks+1)
  if(msy & any(rownames(std) == "log_FMSY")) {
    inds$F.t <- which(rownames(std) == "log_FMSY") 
    inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_MSY"), ncol = n_stocks+1)
  }
  if(static){
    inds$F.t <- rep(which(rownames(std) == "log_FXSPR_static"), length(mod$years_full)) #only 1 value
    inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR_static"), nrow = length(mod$years_full), ncol = n_stocks+1, byrow=TRUE)
    if(msy & any(rownames(std) == "log_FMSY_static")) {
      inds$F.t <- rep(which(rownames(std) == "log_FMSY_static"), length(mod$years_full)) #only 1 value
      inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_MSY_static"), nrow = length(mod$years_full), ncol = n_stocks+1, byrow=TRUE)
    }
  }
  rel.f.vals <- std[inds$full.f,1][status.years] - std[inds$F.t,1][status.years]
  rel.ssb.vals <- std[inds$ssb,1][status.years] - std[inds$SSB.t[,n_stocks+1],1][status.years]
  cov <- mod$sdrep$cov
  log.rel.ssb.rel.F.cov <- lapply(status.years, function(x) {
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    ind <- c(inds$ssb[x],inds$SSB.t[x,n_stocks+1],inds$full.f[x],inds$F.t[x])
    tcov <- cov[ind,ind]
    return(t(K) %*% tcov %*% K)
  })
  if(mod$env$data$n_years_proj>0) { #check whether projecting at F40/FMSY because the ratio to status those years will be 1 and variance 0, but numerical accuracy might be an issue.
    if(!msy) proj_F40 <- which(mod$env$data$proj_F_opt==3)
    if(msy) proj_F40 <- which(mod$env$data$proj_F_opt==6)
    if(length(proj_F40)){
      proj_F40 <- mod$env$data$n_years_model + proj_F40
      proj_F40 <- which(status.years %in% proj_F40)
    }
    if(length(proj_F40)){
      log.rel.ssb.rel.F.cov[proj_F40] <- lapply(log.rel.ssb.rel.F.cov[proj_F40], function(x) {
        x[cbind(c(1,2,2),c(2,2,1))] <- 0
        return(x)
      })      
    }
  }
  do.kobe <- which(sapply(log.rel.ssb.rel.F.cov, function(x) !all(!is.finite(x)))) # only if some non-infinite values for at least some status years 
  if(length(do.kobe)<length(status.years)){
    no.kobe <- which(!status.years %in% do.kobe)
    print(paste0("status confidence region not available for years: ", mod$years_full[status.years[no.kobe]]))
    if(length(do.kobe)==0) return()
  } 
  log.rel.ssb.rel.F.cr <- lapply(1:length(status.years), function(x){
    if(is.na(rel.f.vals[x]) | any(diag(log.rel.ssb.rel.F.cov[[x]])<0)) return(matrix(NA,100,2))
    else return(exp(ellipse::ellipse(log.rel.ssb.rel.F.cov[[x]], centre = c(rel.ssb.vals[x],rel.f.vals[x]), level = 1-alpha)))
  })

  x <- rep(NA,length=length(status.years))
  p.status <- cbind.data.frame(p.ssb.lo.f.lo = x, p.ssb.lo.f.hi = x, p.ssb.hi.f.lo = x, p.ssb.hi.f.hi = x)
  #p.ssb.lo.f.lo <- p.ssb.lo.f.hi <- p.ssb.hi.f.lo <- p.ssb.hi.f.hi <- rep(NA,length(status.years))
  for(i in 1:length(status.years)){
    check.bad.sd <- diag(log.rel.ssb.rel.F.cov[[i]])==0 | diag(log.rel.ssb.rel.F.cov[[i]]) < 0
    if(!any(is.na(check.bad.sd))) if(!any(check.bad.sd)){
      p.status$p.ssb.lo.f.lo[i] <- mnormt::sadmvn(lower = c(-Inf,-Inf), upper = c(-log(2), 0), mean = c(rel.ssb.vals[i],rel.f.vals[i]), 
        varcov = log.rel.ssb.rel.F.cov[[i]])
      p.status$p.ssb.lo.f.hi[i] <- mnormt::sadmvn(lower = c(-Inf,0), upper = c(-log(2), Inf), mean = c(rel.ssb.vals[i],rel.f.vals[i]), 
        varcov = log.rel.ssb.rel.F.cov[[i]])
      p.status$p.ssb.hi.f.lo[i] <- mnormt::sadmvn(lower = c(-log(2),-Inf), upper = c(Inf, 0), mean = c(rel.ssb.vals[i],rel.f.vals[i]), 
        varcov = log.rel.ssb.rel.F.cov[[i]])
      p.status$p.ssb.hi.f.hi[i] <- mnormt::sadmvn(lower = c(-log(2),0), upper = c(Inf, Inf), mean = c(rel.ssb.vals[i],rel.f.vals[i]), 
        varcov = log.rel.ssb.rel.F.cov[[i]])
    }
  }
  vals <- exp(cbind(rel.ssb.vals, rel.f.vals))
  if(is.null(max.x)) max.x <- max(sapply(log.rel.ssb.rel.F.cr, function(x) max(x[,1],na.rm = TRUE)),1.5)
  if(is.null(max.y)) max.y <- max(sapply(log.rel.ssb.rel.F.cr, function(x) max(x[,2],na.rm = TRUE)),1.5)
  if(length(do.kobe)){
    if(single.plot){
      pcols <- 1
      prows <- 1
    } else{
      pcols <- floor(sqrt(length(do.kobe)))
      prows <- ceiling(length(do.kobe)/pcols)
    }
    par(mfrow = c(prows,pcols), mar = c(1,1,0,0), oma = c(5,5,3,1))
    n.plots <- 1
    if(!single.plot) n.plots <- length(do.kobe)
    for(k in 1:n.plots){
      if(n.plots ==1)  {
        plt.which <- do.kobe
      } else plt.which <- do.kobe[k]
      st.yrs.plt <- status.years[plt.which]

      plot(vals[plt.which,1],vals[plt.which,2], ylim = c(0,max.y), xlim = c(0,max.x), xlab = "", ylab = "", type = 'n', axes = F)
      box(lwd = 2)
      if(k > (prows-1)*pcols) axis(1, lwd = 2)
      else axis(1, lwd = 2, labels = FALSE)
      if(k %in% seq(1,n.plots,pcols)) axis(2,lwd = 2)
      else axis(2, lwd = 2, labels = FALSE)
      lims = par("usr")
      tcol <- col2rgb('red')
      tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
      polygon(c(lims[1],0.5,0.5,lims[1]),c(1,1,lims[4],lims[4]), border = tcol, col = tcol)
      tcol <- col2rgb('green')
      tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
      polygon(c(0.5,lims[2],lims[2],0.5),c(lims[3],lims[3],1,1), border = tcol, col = tcol)
      tcol <- col2rgb('yellow')
      tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
      polygon(c(lims[1],0.5,0.5,lims[1]),c(lims[3],lims[3],1,1), border = tcol, col = tcol)
      polygon(c(0.5,lims[2],lims[2],0.5),c(1,1,lims[4],lims[4]), border = tcol, col = tcol)
      legend("topleft", legend = paste0("Year: ", mod$years_full[st.yrs.plt], ", Prob = ", round(p.status$p.ssb.lo.f.hi[plt.which],2)), bty = "n", text.font=1)
      legend("topright", legend = paste0("Year: ", mod$years_full[st.yrs.plt], ", Prob = ", round(p.status$p.ssb.hi.f.hi[plt.which],2)), bty = "n", text.font=1)
      legend("bottomleft", legend = paste0("Year: ", mod$years_full[st.yrs.plt], ", Prob = ", round(p.status$p.ssb.lo.f.lo[plt.which],2)), bty = "n", text.font=1)
      legend("bottomright", legend = paste0("Year: ", mod$years_full[st.yrs.plt], ", Prob = ", round(p.status$p.ssb.hi.f.lo[plt.which],2)), bty = "n", text.font=1)
      text(vals[plt.which,1],vals[plt.which,2], substr(mod$years_full[st.yrs.plt],3,4), font=1)
      for(k in plt.which) {
        polygon(log.rel.ssb.rel.F.cr[[k]][,1], log.rel.ssb.rel.F.cr[[k]][,2], lwd=1)#, border = gray(0.7))
      }
    }
    if(msy & any(rownames(std) == "log_FMSY")){
      mtext(side = 1, outer = TRUE, bquote(SSB*"/"*SSB["MSY"]), line = 2)
      mtext(side = 2, outer = TRUE, bquote(paste(italic(F),"/", italic(F)["MSY"])), line = 2)
    }
    else {
      mtext(side = 1, outer = TRUE, bquote(SSB*"/"*SSB[.(percentSPR)*"%"]), line = 2)
      mtext(side = 2, outer = TRUE, bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])), line = 2)
    }
    title.text <- paste0("Averaged inputs for per recruit calculations")
    if(!static) title.text <- paste0("Annual inputs for per recruit calculations")
    mtext(side = 3, outer = TRUE, title.text, line = 1)
  }
  return(p.status)
}

plot.FXSPR.annual <- function(mod, alpha = 0.05, status.years, max.x=NULL, max.y=NULL, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  percentSPR = mod$env$data$percentSPR
	n_ages = mod$env$data$n_ages
  years_full = mod$years_full
	n_years_full = length(years_full)
  years = mod$years
  n_years = length(years)
  if(missing(status.years)){
    status.years = n_years
    status.lwd = 1
    if(mod$env$data$n_years_proj>0){
      status.years <- c(status.years, n_years_full)
      status.lwd <- c(2,1)
    }
  } else {
    status.lwd <- rep(1, length(status.years))
  }
  all_stocks <- mod$input$data$n_stocks +1
  if(all_stocks == 2) all_stocks <- 1
  all_catch <- mod$input$data$n_fleets + mod$input$data$n_regions + 1
  if(all_catch == 2) all_catch <- 1
  std <- summary(mod$sdrep, "report")
  inds <- list()
  inds$Y.t <- matrix(which(rownames(std) == "log_Y_FXSPR"), nrow = n_years_full)[,1:all_catch, drop= F]
  inds$F.t <- which(rownames(std) == "log_FXSPR")
  inds$SSB.t <- matrix(which(rownames(std) == "log_SSB_FXSPR"), nrow = n_years_full)[,1:all_stocks, drop= F]
  # print(dim(inds$SSB.t))
  # print(all_stocks)
  inds$ssb <- which(rownames(std) == "log_SSB_all")
  inds$full.f <- which(rownames(std) == "log_F_tot")
  na.sd <- sapply(inds, function(x) any(is.na(std[x,2])))
  # log.faa <- array(std[inds$faa,1], dim = dim(mod$rep$FAA_tot))
  # age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
  # inds$full.f <- (age.full.f-1)*n_years_full + 1:n_years_full  + min(inds$faa) - 1 #cbind(1:n_years, age.full.f)
  cov <- mod$sdrep$cov
  log.rel.ssb.rel.F.cov <- lapply(1:n_years_full, function(x)
  {
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    ind <- c(inds$ssb[x],inds$SSB.t[x,all_stocks],inds$full.f[x],inds$F.t[x])
    tcov <- cov[ind,ind]
    return(t(K) %*% tcov %*% K)
  })
  if(mod$env$data$n_years_proj>0) { #check whether projecting at F40 because the ratio to status those years will be 1 and variance 0, but numerical accuracy might be an issue.
    proj_F40 <- which(mod$env$data$proj_F_opt==3)
    if(length(proj_F40)){
      proj_F40 <- mod$env$data$n_years_model + proj_F40
      log.rel.ssb.rel.F.cov[proj_F40] <- lapply(log.rel.ssb.rel.F.cov[proj_F40], function(x) {
        x[cbind(c(1,2,2),c(2,2,1))] <- 0
        return(x)
      })      
    }
  }

  ylabs <- c(
    bquote(paste('Yield(',italic(F)[paste(.(percentSPR), "%")], ')')),
    bquote(italic(F)[paste(.(percentSPR), "%")]),
    bquote(paste('SSB(', italic(F)[paste(.(percentSPR), "%")],')')))
  np <- c(all_catch,1,all_stocks)
  # FSPR absolute --------------------------------------------------
  if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_absolute.pdf")), family = fontfam, height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("FSPR_absolute.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
  par(mfrow = c(3,1), mar = c(1,5,1,1), oma = c(4,1,3,1))
  for(i in 1:3)
  {
    t.ylab <- ylabs[i]
    vals <- matrix(std[inds[[i]],1], ncol = np[i])#[1:n_years_full,,drop = F]
    cv <- matrix(std[inds[[i]],2], ncol = np[i])#[1:n_years_full,,drop = F]
    ci <-  lapply(1:np[i], function(x) vals[,x] + cbind(qnorm(1-alpha/2)*cv[,x], -qnorm(1-alpha/2)*cv[,x]))
    max_y1 <- exp(max(sapply(ci, max, na.rm = T), na.rm = T))
    plot.colors <- "black"
    if(np[i]>1) plot.colors <- viridisLite::viridis(n = np[i])
    tcol <- adjustcolor(plot.colors, alpha.f = 0.4)
    if(all(!is.nan(unique(vals)))){
  	  if(!na.sd[i]) plot(years_full, exp(vals[,1]), xlab = '', ylab = t.ylab, ylim = c(0,max_y1), type = 'n', xaxt = 'n', cex.lab = 2)
      if(na.sd[i]) plot(years_full, exp(vals[,1]), xlab = '', ylab = t.ylab, ylim = c(0,max(exp(vals),na.rm= TRUE)), type = 'n', cex.lab = 2)
  	  grid(col = gray(0.7))
      for(p in 1:np[i]){
        lines(years_full, exp(vals[,p]), col = plot.colors[p], lwd = 2)  
        polyy = exp(c(ci[[p]][,1],rev(ci[[p]][,2])))
        polyx = c(years_full,rev(years_full))
        polyx = polyx[!is.na(polyy)]
        polyy = polyy[!is.na(polyy)]
        polygon(polyx, polyy, col = tcol[p], border = tcol[p], lwd = 1)
      }
      if(mod$env$data$n_years_proj>0) abline(v=tail(years,1), lty=2, lwd=1)
      if(i == 1) legend("topright", legend = c(mod$input$fleet_names, mod$input$region_names, "Total"), lty = 1, col = plot.colors, bty = "n")
      if(i == 3) {
        if(np[i] > 1) legend("topright", legend = c(mod$input$stock_names, "Total"), lty = 1, col = plot.colors, bty = "n")
        axis(1, cex.lab= 2)
        mtext(side = 1, outer = TRUE, "Year", cex = 2, line = 2)
      } else{
        axis(1, labels = FALSE)
      }
    } else { # all nan, print error message
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.55, paste("Error in plot, all values are NaN"), cex = 1.6, col = "black")
      text(x = 0.5, y = 0.45, ylabs[[i]], cex = 1.6, col = "black")
    }
	}
  title.text <- paste0("Annual inputs used in per recruit calculations")
  mtext(side = 3, outer = TRUE, title.text, line = 1)
  if(do.tex | do.png) dev.off() else par(origpar)
  plot.colors <- "black"
  tcol <- adjustcolor(plot.colors, alpha.f = 0.4)
  # FSPR relative --------------------------------------------------
  if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_relative.pdf")), family = fontfam, height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("FSPR_relative.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
  par(mfrow=c(2,1), oma = c(4,1,3,0), mar = c(1,5,1,1))
  rel.vals <- list(ssb = std[inds$ssb,1][1:n_years_full] - std[inds$SSB.t[,all_stocks],1][1:n_years_full])
  rel.vals$f = std[inds$full.f,1][1:n_years_full] - std[inds$F.t,1][1:n_years_full]
  #rel.ssb.vals <- std[inds$ssb,1][1:n_years_full] - std[inds$SSB.t[,all_stocks],1][1:n_years_full]
  cv <- list(ssb = sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[1,1]))))
  cv$f <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[2,2])))
  ci <-  list(ssb = rel.vals$ssb + cbind(-qnorm(1-alpha/2)*cv$ssb, qnorm(1-alpha/2)*cv$ssb))
  ci$f <- rel.vals$f + cbind(-qnorm(1-alpha/2)*cv$f, qnorm(1-alpha/2)*cv$f)
  ylabs <- c(bquote(paste("SSB/", SSB[paste(.(percentSPR),"%")])), bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])))
  any_na_sd <- c(ssb = na.sd["ssb"], f = na.sd["full.f"])
  for(i in 1:2){
    if(all(!is.nan(unique(rel.vals[[i]])))){
      ymax <- max(exp(ci[[i]]),1, na.rm = TRUE)
      if(na.sd[i]) ymax <- max(exp(rel.vals[[i]]),1, na.rm = TRUE)
      plot(years_full, exp(rel.vals[[i]]), xlab = '', ylab = ylabs[i],
        ylim = c(0,ymax), xaxt = 'n', type = 'l', cex.lab = 2)
      if(i==1) axis(1, labels = FALSE)
      else {
        axis(1)
        mtext(side = 1, outer = TRUE, "Year", line = 2, cex = 2)
      }
      grid(col = gray(0.7))
      polyy = exp(c(ci[[i]][,1],rev(ci[[i]][,2])))
      polyx = c(years_full,rev(years_full))
      polyx = polyx[!is.na(polyy)]
      polyy = polyy[!is.na(polyy)]
      polygon(polyx, polyy, col = tcol, border = tcol, lwd = 1)
      abline(h=1, lty = 2)
      abline(h=0.5, lty = 2, col = 'red')
      if(mod$env$data$n_years_proj>0) abline(v=tail(mod$years,1), lty=2, lwd=1)
    } else {
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.58, paste("Error in plot, all values are NaN"), cex = 1.6, col = "black")
      text(x = 0.5, y = 0.42, ylabs[i], cex = 1.6, col = "black")    
    }
  }

  #title.text <- paste0("Averaged inputs for per recruit calculations")
  #if(!static) 
  title.text <- paste0("Annual inputs used in per recruit calculations")
  mtext(side = 3, outer = TRUE, title.text, line = 1)
  if(do.tex | do.png) dev.off() else par(origpar)

  if(do.tex) cairo_pdf(file.path(od, paste0("Kobe_status.pdf")), family = fontfam, height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("Kobe_status.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
  x <- kobe.plot(mod, status.years=status.years, static = TRUE, single.plot = TRUE, alpha = 0.05, max.x=max.x, max.y= max.y)
  if(do.tex | do.png) dev.off() else par(origpar)
  return(x)
}  # end function
#revised

plot.MSY.annual <- function(mod, alpha = 0.05, status.years, max.x=NULL, max.y=NULL, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)
  years_full = mod$years_full
  n_years_full = length(years_full)
  if(missing(status.years)){
    status.years = n_years
    status.lwd = 1
    if(mod$env$data$n_years_proj>0){
      status.years <- c(status.years, n_years_full)
      status.lwd <- c(2,1)
    }
  } else {
    status.lwd <- rep(1, length(status.years))
  }  
  std = summary(mod$sdrep, "report")
  cov <- mod$sdrep$cov
	if(dat$recruit_model %in% (3:4)) #Beverton-Holt assumed in model fit
	{ # test to make sure steepness was estimated
    tcol <- col2rgb('black')
    tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
		inds <- list(MSY = array(which(rownames(std) == "log_MSY"), dim = dim(mod$rep$log_MSY))[dat$n_fleets+1,dat$n_stocks+1,1:n_years_full])
		inds$FMSY <- which(rownames(std) == "log_FMSY")
		inds$SSBMSY <- matrix(which(rownames(std) == "log_SSB_MSY"),nrow = n_years_full, ncol = dat$n_stocks+1)[1:n_years_full,dat$n_stocks+1] #total SSBMSY
		inds$RMSY <- matrix(which(rownames(std) == "log_R_MSY"),nrow = n_years_full, ncol = dat$n_stocks+1)[1:n_years_full,dat$n_stocks+1] #total RMSY
		inds$ssb <- which(rownames(std) == "log_SSB_all")
    inds$full.f <- which(rownames(std) == "log_F_tot")
  	# inds$faa <- which(rownames(std) == "log_FAA_tot")
	  # log.faa <- matrix(std[inds$faa,1], n_years_full, n_ages)
	  # age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
	  # inds$full.f <- (age.full.f-1)*n_years_full + 1:n_years_full  + min(inds$faa) - 1 #cbind(1:n_years, age.full.f)
	  ylabs <- c(expression(MSY),expression(italic(F)[MSY]), expression(SSB[MSY]), expression(italic(R)[MSY]))
	  log.rel.ssb.rel.F.cov <- lapply(1:n_years_full, function(x)
	  {
	    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
	    ind <- c(inds$ssb[x],inds$SSBMSY[x],inds$full.f[x],inds$FMSY[x])
	    tcov <- cov[ind,ind]
	    return(t(K) %*% tcov %*% K)
	  })
    if(mod$env$data$n_years_proj>0) { #check whether projecting at Fmsy because the ratio to status those years will be 1 and variance 0, but numerical accuracy might be an issue.
      proj_Fmsy <- which(mod$env$data$proj_F_opt==6)
      if(length(proj_Fmsy)){
        proj_Fmsy <- mod$env$data$n_years_model + proj_Fmsy
        log.rel.ssb.rel.F.cov[proj_Fmsy] <- lapply(log.rel.ssb.rel.F.cov[proj_Fmsy], function(x) {
          x[cbind(c(1,2,2),c(2,2,1))] <- 0
          return(x)
        })      
      }
    }

    # 4-panel MSY plot
    if(do.tex) cairo_pdf(file.path(od, paste0("MSY_4panel_F_SSB_R.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("MSY_4panel_F_SSB_R.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mfrow=c(2,2), oma = c(4,1,3,0), mar = c(1,5,1,1))
    for(i in 1:4)
    {
      t.ind <- inds[[i]]
      t.ylab <- ylabs[i]
      vals <- std[t.ind,1]
      cv <- std[t.ind,2]
      ci <-  vals + cbind(qnorm(1-alpha/2)*cv, -qnorm(1-alpha/2)*cv)
      na.ci <- all(is.na(ci))
      # remove NaN and Inf
      # rm.rows <- which(vals < 0)
      vals[!is.finite(vals)] = NA
      # vals[rm.rows] = NA
      # rm.rows <- which(ci[,2] < 0)
      ci[!is.finite(exp(ci))] = NA
      # ci[rm.rows,] = NA
      labels <- i > 2
      if(na.ci) ylim <- c(0,max(exp(vals),na.rm=TRUE))
      else ylim <- c(0,max(exp(ci),na.rm=TRUE))
      plot(years_full, exp(vals), xlab = '', ylab = t.ylab, ylim = ylim, xaxt = "n", type = 'l')
      axis(1, labels = labels)

		  # if(!na.ci) plot(years_full, exp(vals), xlab = '', ylab = t.ylab, ylim = c(0,max(exp(ci),na.rm=TRUE)), type = 'l')
      # if(na.ci) plot(years_full, exp(vals), xlab = '', ylab = t.ylab, ylim = c(0,max(exp(vals),na.rm=TRUE)), type = 'l')
		  grid(col = gray(0.7))
      not.na.ci <- !is.na(ci[,1])
		  polygon(c(years_full[not.na.ci],rev(years_full[not.na.ci])), exp(c(ci[,1][not.na.ci],rev(ci[,2][not.na.ci]))), col = tcol, border = tcol, lwd = 1)
      # polygon(c(years_full,rev(years_full)), exp(c(ci[,1],rev(ci[,2]))), col = tcol, border = tcol, lwd = 1)
      if(mod$env$data$n_years_proj>0) abline(v=tail(years,1), lty=2, lwd=1)
		}
    title.text <- paste0("Annual inputs used in per recruit calculations")
    mtext(side = 3, outer = TRUE, title.text, line = 1)
    mtext(side = 1, outer = TRUE, "Year", line = 2, cex = 2)
    if(do.tex | do.png) dev.off() else par(origpar)
    # 2-panel SSB_MSY and F_MSY
    if(do.tex) cairo_pdf(file.path(od, paste0("MSY_relative_status.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("MSY_relative_status.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    par(mfrow=c(2,1), oma = c(4,1,3,0), mar = c(1,5,1,1))
    rel.ssb.vals <- std[inds$ssb,1] - std[inds$SSBMSY,1]
    cv <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[1,1])))
    ci <-  rel.ssb.vals + cbind(-qnorm(1-alpha/2)*cv, qnorm(1-alpha/2)*cv)
    na.ci <- all(is.na(ci))
    # remove NaN and Inf
    # rm.rows <- which(rel.ssb.vals < 0)
    rel.ssb.vals[!is.finite(rel.ssb.vals)] = NA
    # rel.ssb.vals[rm.rows] = NA
    # rm.rows <- which(ci[,2] < 0 | ci[,1] < 0)
    ci[!is.finite(exp(ci))] = NA
    # ci[rm.rows,] = NA    
    if(na.ci) ylim <- c(0,max(exp(rel.ssb.vals),na.rm=TRUE))
    else ylim <- c(0,max(exp(ci),na.rm=TRUE))
    plot(years_full, exp(rel.ssb.vals), xlab = '', xaxt = "n", ylab = expression(paste("SSB/", SSB[MSY])), ylim = ylim, type = 'l')
		# if(!na.ci) plot(years_full, exp(rel.ssb.vals), xlab = '', xaxt = "n", ylab = expression(paste("SSB/", SSB[MSY])), ylim = c(0,max(exp(ci),na.rm=TRUE)), type = 'l')
    # if(na.ci) plot(years_full, exp(rel.ssb.vals), xlab = '', xaxt = "n", ylab = expression(paste("SSB/", SSB[MSY])), ylim = c(0,max(exp(rel.ssb.vals),na.rm=TRUE)), type = 'l')
    axis(1, labels = FALSE)
	  grid(col = gray(0.7))
	  polygon(c(years_full,rev(years_full)), exp(c(ci[,1],rev(ci[,2]))), col = tcol, border = "transparent", lwd = 1)
	  abline(h=1, lty = 2)
	  abline(h=0.5, lty = 2, col = 'red')
    if(mod$env$data$n_years_proj>0) abline(v=tail(years,1), lty=2, lwd=1)

    rel.f.vals <- std[inds$full.f,1][1:n_years_full] - std[inds$FMSY,1][1:n_years_full]
    cv <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[2,2])))
    ci <-  rel.f.vals + cbind(-qnorm(1-alpha/2)*cv, qnorm(1-alpha/2)*cv)
    na.ci <- all(is.na(ci))
      # remove NaN and Inf
      # rm.rows <- which(rel.f.vals < 0)
      rel.f.vals[!is.finite(rel.f.vals)] = NA
      # rel.f.vals[rm.rows] = NA
      # rm.rows <- which(ci[,2] < 0 | ci[,1] < 0)
      ci[!is.finite(exp(ci))] = NA
      # ci[rm.rows,] = NA    
    if(na.ci) ylim <- c(0,max(exp(rel.f.vals),na.rm=TRUE))
    else ylim <- c(0,max(exp(ci),na.rm=TRUE))
    plot(years_full, exp(rel.f.vals), xlab = '', xaxt = "n",ylab = expression(paste(italic(F),"/", italic(F)[MSY])), ylim = ylim, type = 'l')
    axis(1, labels = TRUE)
    grid(col = gray(0.7))
    polygon(c(years_full,rev(years_full)), exp(c(ci[,1],rev(ci[,2]))), col = tcol, border = tcol, lwd = 1)
	  abline(h=1, lty = 2, col = 'red')
    if(mod$env$data$n_years_proj>0) abline(v=tail(years,1), lty=2, lwd=1)
    title.text <- paste0("Annual inputs used in per recruit calculations")
    mtext(side = 3, outer = TRUE, title.text, line = 1)
    mtext(side = 1, outer = TRUE, "Year", line = 2, cex = 2)
    if(do.tex | do.png) dev.off() else par(origpar)
    if(do.tex) cairo_pdf(file.path(od, paste0("Kobe_msy_status.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Kobe_msy_status.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
    x <- kobe.plot(mod, status.years=status.years, static = TRUE, msy = TRUE, single.plot = TRUE, alpha = 0.05, max.x=max.x, max.y= max.y)
    if(do.tex | do.png) dev.off() else par(origpar)
	}
}  # end function
#revised

#------------------------------------
plot.yield.curves <- function(mod, nyrs.ave = 5, plot=TRUE, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
	n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)
  avg.ind <- tail(1:n_years, nyrs.ave)
  # avg.ind = (n_years-nyrs.ave+1):n_years
  mat.age <- apply(dat$mature[1,,][avg.ind,,drop = FALSE],2,mean)
  ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb[1],,][avg.ind,,drop = FALSE],2,mean)
	catch.waa <- apply(dat$waa[dat$waa_pointer_fleets,avg.ind,,drop = FALSE], c(1,3),mean)
	M.age <- apply(mod$rep$MAA[1,1,,][avg.ind,,drop = FALSE],2,mean)
  FAA <- apply(mod$rep$FAA[,avg.ind,,drop=FALSE], c(1,3), mean) #mean FAA by fleet
  FAAtot = apply(FAA,2,sum) #sum across fleet averages 
  seltot <- FAAtot/max(FAAtot)
  selAA <- FAA/max(FAAtot) 


	spawn.time <- mean(dat$fracyr_SSB[avg.ind,1])
  spr0 = get_SPR(F=0, M=M.age, sel=selAA, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)

	F.range <- seq(0,2.0, by=0.01)
	nF <- length(F.range)

  spr = suppressWarnings(sapply(F.range, function(x) get_SPR(F=x, M=M.age, sel=selAA, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)))
  ypr = suppressWarnings(sapply(F.range, function(x) sum(get_YPR_fleets(F=x*selAA, M=M.age, waacatch=catch.waa))))
  pr = spr/spr0
  if(plot & all(!is.na(pr))){
    if(do.tex) cairo_pdf(file.path(od, paste0("YPR_F_curve_plot.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("YPR_F_curve_plot.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
  	par(mfrow=c(1,1), mar=c(4,4,2,4))
  	plot(F.range, ypr, type='n', xlab="Full F", ylab="", lwd=2, col="blue3", ylim=c(0,max(ypr)), axes = FALSE)
    box()
    axis(1)
    axis(2, col='blue3', col.axis="blue3")
  	abline(v=seq(0.1,2.0, by=0.1), col="grey85")
  	lines(F.range, ypr, lwd=2, col="blue3" )
  	points(F.range, ypr, pch=19, col="blue3")
    mtext(side =2, "Yield per Recruit", col = "blue3", line = 2)
  	scale.pr <- max(pr)/max(ypr)

  	lines(F.range, pr/scale.pr, col="red", lwd=2)
  	points(F.range, pr/scale.pr, pch=19, col="red")
  	axis(side=4, at=seq(0,1,by=0.1)/scale.pr, lab=seq(0,1,by=0.1), las=2, col='red', col.axis="red")
  	mtext(side=4, "% SPR", line=3, col="red")
  	title (paste("YPR-SPR Reference Points (Years Avg = ", nyrs.ave,")", sep=""), outer=T, line=-1)
  	if(do.tex | do.png) dev.off() else par(origpar)
  }

  if(!plot){
    if(do.tex) cairo_pdf(file.path(od, paste0("YPR_F_curve_table.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("YPR_F_curve_table.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)

  	par(mfrow=c(1,1), mar=c(2,2,2,2))
  	plot(seq(1,42), seq(-3,38), type='n', axes=F, bty='n',xlab="",ylab="")
  	text(x=0,y=36, labels="F", font=2, pos=4)
  	text(x=4, y=36, labels="YPR" , font=2, pos=4)
  	text(x=9, y=36, labels="SPR", font=2, pos=4)
  	text(x=15,y=36, labels="F", font=2, pos=4)
  	text(x=19, y=36, labels="YPR" , font=2, pos=4)
  	text(x=24, y=36, labels="SPR", font=2, pos=4)
  	text(x=30,y=36, labels="F", font=2, pos=4)
  	text(x=34, y=36, labels="YPR" , font=2, pos=4)
  	text(x=39, y=36, labels="SPR", font=2, pos=4)

    text(x=0, y=35:1, labels=F.range[1:35], cex=0.82, pos=4, font=1)
    text(x=4, y=35:1, labels=round(ypr[1:35],4), cex=0.82, pos=4, font=1)
    text(x=9, y=35:1, labels=round(pr[1:35],4), cex=0.82, pos=4, font=1)

    text(x=15, y=35:1, labels=F.range[36:70], cex=0.82, pos=4, font=1)
    text(x=19, y=35:1, labels=round(ypr[36:70],4), cex=0.82, pos=4, font=1)
    text(x=24, y=35:1, labels=round(pr[36:70],4), cex=0.82, pos=4, font=1)

    text(x=30, y=35:1, labels=F.range[71:105], cex=0.82, pos=4, font=1)
    text(x=34, y=35:1, labels=round(ypr[71:105],4), cex=0.82, pos=4, font=1)
    text(x=39, y=35:1, labels=round(pr[71:105],4), cex=0.82, pos=4, font=1)

  	title (paste("YPR-SPR Reference Points (Years Avg = ", nyrs.ave,")", sep=""), outer=T, line=-1)
  	if(do.tex | do.png) dev.off() else par(origpar)
  }
  ypr.table = cbind.data.frame(F.range, ypr, pr)
  colnames(ypr.table) <- c("Full.F", "YPR", "SPR")
	#write.csv(ypr.table, file=paste(od,"YPR_Table.csv", sep=""), row.names=F)
  # par(origpar)
} # end function

#------------------------------------
#------------------------------------

plot.retro <- function(mod,y.lab,y.range1,y.range2, alpha = 0.05, what = "SSB", age=NULL, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od) {
  
  origpar <- par(no.readonly = TRUE)
  years = mod$years
	n_years <- length(years) # don't use projections
  sfns <- chartr(" ", "_", mod$input$stock_names)
  rfns <- chartr(" ", "_", mod$input$region_names)

  npeels = length(mod$peels)
  if(npeels) {
    repwhat <- what
    if(what == "NAA_age") {
      age.ind <- (1:mod$input$data$n_ages)[age] 
      #age.ind <- match(age, mod$ages.lab)
      what.print <- paste0(what, "_", age)
      repwhat <- "NAA"
    } else {
      what.print <- what
    }
    n_ages <- dim(mod$rep$NAA)[4]

    # standard retro plot
    plot.colors <- viridisLite::viridis(n=npeels+1)
    tcol <- adjustcolor(plot.colors, alpha.f=0.4)
    # tcol = col2rgb(plot.colors)
    # tcol = rgb(tcol[1,],tcol[2,],tcol[3,], maxColorValue = 255, alpha = 200)
    # if(what %in% c("NAA","NAA_age")){
    #   res = list(mod$rep$NAA[,,1:n_years,, drop = FALSE])
    #   res[2:(npeels+1)] = lapply(mod$peels, function(x) x$rep$NAA[,,1:n_years,, drop = FALSE])
    # } else {
      # res = list(head(mod$rep[[what]],n_years))
      res = list(mod$rep[[repwhat]])
      res[2:(npeels+1)] = lapply(mod$peels, function(x) x$rep[[repwhat]])
    # }

    if(what %in% c("NAA","NAA_age")) for(s in 1:mod$input$data$n_stocks) {
      regions <- 1:mod$input$data$n_regions
      #if(sum(mod$input$data$NAA_where[s,,,])== 0) regions <- mod$input$data$spawn_regions[s]

      for(r in regions) {
        if(what == "NAA") if(sum(mod$input$data$NAA_where[s,r,]) > 0){
          what.print.s <- paste0(sfns[s], "_", rfns[r], "_", what.print)
          if(do.tex) cairo_pdf(file.path(od, paste0(what.print.s,"_retro.pdf")), family = fontfam, height = 10, width = 10)
          if(do.png) png(filename = file.path(od, paste0(what.print.s,"_retro.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
          #n_ages <- dim(res[[1]])[4]
          par(mfcol = c(ceiling(n_ages/2),2), mar = c(4,4,1,1), oma = c(0,0,4,0))
          for(i in 1:n_ages) if(mod$input$data$NAA_where[s,r,i]>0) {
            y.range1 <- range(sapply(res, function(x) range(x[s,r,,i])))
            plt.type = 'l'
            if(mod$input$data$NAA_where[s,r,i]==0) {
              plt.type <- 'n'
              #y.range1 <- c(-1,1)
            }
            plot(years,res[[1]][s,r,1:n_years,i],lwd=1,col=plot.colors[1],type=plt.type,xlab="Year",ylab=paste0("Numbers at age ", mod$ages.lab[i]),ylim=y.range1)
            grid(col = gray(0.7), lty = 2)
            for (j in 1:npeels) {
              lines(years[1:(n_years-j)],res[[j+1]][s,r,1:(n_years-j),i], col = tcol[j+1])
              points(years[n_years-j],res[[j+1]][s,r,n_years-j,i],pch=16,col=plot.colors[j+1])
            }
          }
          title(paste0(mod$input$stock_names[s], " in ", mod$input$region_names[r]), line = 1, outer = TRUE)
          if(do.tex | do.png) dev.off() else par(origpar)
        }
        if(what == "NAA_age") if(mod$input$data$NAA_where[s,r,age.ind] > 0){# only specified age
          what.print.s <- paste0(sfns[s], "_", rfns[r], "_", what.print)
          if(do.tex) cairo_pdf(file.path(od, paste0(what.print.s,"_retro.pdf")), family = fontfam, height = 10, width = 10)
          if(do.png) png(filename = file.path(od, paste0(what.print.s,"_retro.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
          par(mfrow = c(1,1))
          par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,4,0))
          y.range1 <- range(sapply(res, function(x) range(x[s,r,,age.ind])))
          plot(years,res[[1]][s,r,1:n_years,age.ind],lwd=1,col=plot.colors[1],type='l',xlab="Year",ylab=paste0("Numbers at age ", mod$ages.lab[age.ind]),ylim=y.range1)
          grid(col = gray(0.7), lty = 2)
          for (j in 1:npeels) {
            lines(years[1:(n_years-j)],res[[j+1]][s,r,1:(n_years-j),age.ind], col = tcol[j+1])
            points(years[n_years-j],res[[j+1]][s,r,n_years-j,age.ind],pch=16,col=plot.colors[j+1])
          }
          title(paste0(mod$input$stock_names[s], " in ", mod$input$region_names[r]), line = 1, outer = TRUE)
          if(do.tex | do.png) dev.off() else par(origpar)
        }
      }
    }
    if(what %in% c("SSB","Fbar")) {
      if(what == "SSB") names.plot <- mod$input$stock_names
      if(what == "Fbar") names.plot <- mod$input$region_names
      n <- length(names.plot)

      for(p in 1:n){
        what.print.p <- paste0(names.plot[p], "_", what.print)
        if(do.tex) cairo_pdf(file.path(od, paste0(what.print.p,"_retro.pdf")), family = fontfam, height = 10, width = 10)
        if(do.png) png(filename = file.path(od, paste0(what.print.p,"_retro.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
        if(missing(y.range1)) y.range1 <- range(sapply(res, function(x) range(x[,p])))
        par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,4,0))
        plot(years,res[[1]][1:n_years,p],lwd=1,col=plot.colors[1],type='l',xlab="Year",ylab=what,ylim=y.range1)
        grid(col = gray(0.7), lty = 2)
        for (i in 1:npeels)
        {
          lines(years[1:(n_years-i)],res[[i+1]][1:(n_years-i),p], col = tcol[i+1])
          points(years[n_years-i],res[[i+1]][n_years-i,p],pch=16,col=plot.colors[i+1])
        }
        title(names.plot[p], line = 1, outer = TRUE)
        if(do.tex | do.png) dev.off() else par(origpar)
      }
    }

    # relative retro plot
    if(missing(y.lab)) y.lab = bquote(paste("Mohn's ", rho, "(",.(what),")"))
    #if(what %in% c("SSB","Fbar")) rel.res = lapply(1:length(res), function(x) res[[x]][1:(n_years - x + 1),]/res[[1]][1:(n_years - x + 1),] - 1)
    rho.vals = mohns_rho(mod)

    if(what %in% c("NAA","NAA_age")) for(s in 1:mod$input$data$n_stocks){
      regions <- 1:mod$input$data$n_regions
      #if(sum(mod$input$data$can_move[s,,,])== 0) regions <- mod$input$data$spawn_regions[s]
      for(r in regions) {
        rel.res = lapply(1:length(res), function(x) res[[x]][s,r,1:(n_years - x + 1),]/res[[1]][s,r,1:(n_years - x + 1),] - 1)
        if(what == "NAA") if(sum(mod$input$data$NAA_where[s,r,]) > 0) {
          what.print.s <- paste0(sfns[s], "_", rfns[r], "_", what.print)
          if(do.tex) cairo_pdf(file.path(od, paste0(what.print.s,"_retro_relative.pdf")), family = fontfam, height = 10, width = 10)
          if(do.png) png(filename = file.path(od, paste0(what.print.s,"_retro_relative.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
          par(mfcol = c(ceiling(n_ages/2),2), mar = c(4,4,1,1), oma = c(0,0,4,0))
          for(i in 1:n_ages) {
            y.range2 <- c(-1,max(sapply(rel.res, function(x) range(x[,i]))))
            plt.type = 'l'
            if(mod$input$data$NAA_where[s,r,i]==0) {
              plt.type <- 'n'
              y.range2 <- c(-1,1)
            }
            # print(r)
            # print(s)
            # print(i)
            # print(y.range2)
            # print(file.path(od, paste0(what.print.s,"_retro.pdf")))
            # print(do.tex)
            # print(do.png)
            plot(years,rel.res[[1]][,i],lwd=1,col=plot.colors[1],type=plt.type,xlab="Year",ylab=bquote(paste("Mohn's ", rho, "(Numbers at age ", .(mod$ages.lab[i]), ")")),ylim = y.range2)
            grid(col = gray(0.7), lty = 2)
            for (j in 1:npeels) {
              lines(years[1:(n_years-j)],rel.res[[j+1]][1:(n_years-j),i], col = tcol[j+1])
              points(years[n_years-j],rel.res[[j+1]][n_years-j,i],pch=16,col=plot.colors[j+1])
            }
            rho.plot <- round(rho.vals$naa[s,r,i],3)
          }
          legend("bottomleft", legend = bquote(rho == .(rho.plot)), bty = "n")
          title(paste0(mod$input$stock_names[s], " in ", mod$input$region_names[r]), line = 1, outer = TRUE)
          if(do.tex | do.png) dev.off() else par(origpar)
        }
        if(what == "NAA_age") if(mod$input$data$NAA_where[s,r,age.ind] > 0) {
          what.print.s <- paste0(sfns[s], "_", rfns[r], "_", what.print)
          if(do.tex) cairo_pdf(file.path(od, paste0(what.print.s,"_retro_relative.pdf")), family = fontfam, height = 10, width = 10)
          if(do.png) png(filename = file.path(od, paste0(what.print.s,"_retro_relative.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
          par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,4,0))
          y.range2 <- c(-1,max(sapply(rel.res, function(x) range(x[,age.ind]))))
          plot(years,rel.res[[1]][,age.ind],lwd=1,col=plot.colors[1],type='l',xlab="Year",ylab=bquote(paste("Mohn's ", rho, "(Numbers at age ", .(mod$ages.lab[age.ind]), ")")),ylim = y.range2)
          grid(col = gray(0.7), lty = 2)
          for (j in 1:npeels) {
            lines(years[1:(n_years-j)],rel.res[[j+1]][1:(n_years-j),age.ind], col = tcol[j+1])
            points(years[n_years-j],rel.res[[j+1]][n_years-j,age.ind],pch=16,col=plot.colors[j+1])
          }
          rho.plot <- round(rho.vals$naa[s,r,age.ind],3)
          legend("bottomleft", legend = bquote(rho == .(rho.plot)), bty = "n")
          title(paste0(mod$input$stock_names[s], " in ", mod$input$region_names[r]), line = 1, outer = TRUE)
          if(do.tex | do.png) dev.off() else par(origpar)
        }
      }
    }

    if(what %in% c("SSB","Fbar")) {
      rel.res = lapply(1:length(res), function(x) res[[x]][1:(n_years - x + 1),,drop=FALSE]/res[[1]][1:(n_years - x + 1),,drop=FALSE] - 1)
      if(what == "SSB") names.plot <- mod$input$stock_names
      if(what == "Fbar") names.plot <- mod$input$region_names
      n <- length(names.plot)

      for(p in 1:n) {
        what.print.p <- paste0(names.plot[p], "_", what.print)
        if(do.tex) cairo_pdf(file.path(od, paste0(what.print.p,"_retro_relative.pdf")), family = fontfam, height = 10, width = 10)
        if(do.png) png(filename = file.path(od, paste0(what.print.p,"_retro_relative.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
        if(missing(y.range2)) y.range2 <- c(-1,max(sapply(rel.res, function(x) range(x[,p]))))
        par(mfrow = c(1,1), mar = c(4,4,1,1), oma = c(0,0,4,0))
        plot(years,rel.res[[1]][1:n_years,p],lwd=1,col="black",type='l',xlab="Year",ylab=y.lab,ylim=y.range2)
        grid(col = gray(0.7), lty = 2)
        for (i in 1:npeels) {
          lines(years[1:(n_years-i)],rel.res[[i+1]][1:(n_years-i),p], col = tcol[i+1])
          points(years[n_years-i],rel.res[[i+1]][n_years-i,p],pch=16,col=plot.colors[i+1])
        }
        rho.plot <- round(rho.vals[[what]][p],3)
        legend("bottomleft", legend = bquote(rho == .(rho.plot)), bty = "n")
        title(names.plot[p], line = 1, outer = TRUE)
        if(do.tex | do.png) dev.off() else par(origpar)
      }
    }
    # par(origpar)
    return(rho.vals)
  }
}

#revised

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-----Catch curve and cohort correspondence plots-------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# function replace zeros and take log
rep0log <- function(mat)
{
  mat[mat==0] = NA
	return(log(mat))
}
#-------------------------------------------------------------------------------

# function make cohorts
makecohorts <- function(mat)
{
	NAmat <- matrix(NA, NCOL(mat), NCOL(mat))
	matcoh1 <- rbind(NAmat,mat,NAmat)
	nr <- NROW(matcoh1)
	nc <- NCOL(mat)
	matcoh <- matrix(NA, nrow=nr, ncol=nc)
	for (i in 1:nc) matcoh[1:(nr-i),i] <- matcoh1[i:(nr-1),i]
	return(matcoh)
}
#-------------------------------------------------------------------------------

# function to plot all ages by all ages of data by cohorts
# assumes matrix has ages as columns and that cohorts are in rows
plotcoh <- function(matcoh,mytitle="",mylabels=NA, save.plots = FALSE, mod)
{
	origpar <- par(no.readonly = TRUE)
  nc <- NCOL(matcoh)
	my.cor <- cor(matcoh,use="pairwise.complete.obs")
	my.cor.round <- round(my.cor,2)
	par(mfcol=c(nc,nc))
	par(oma=c(0,0,3,1),mar=c(1,1,0,0))
	for (i in 1:nc)
	{
		for (j in nc:1)
		{
			if (i == j)
			{
				plot(1:10,1:10,type='n',axes=FALSE)
				text(5,5,paste0("age-",mod$ages.lab[i]),cex=1.4)
			}
			if (i < j)
			{
				if (!is.na(my.cor[i,j]))
				{
					plot(matcoh[,i],matcoh[,j],axes=FALSE) # make sure have some data to plot
					xx <- matcoh[,i]
					yy <- matcoh[,j]
					my.fit <- lm(yy~xx)
					if (!is.na(my.fit$coefficients[2])) abline(my.fit,col="red")
					xrng <- data.frame(xx = seq(min(xx,na.rm=T),max(xx,na.rm=T),length.out=100))
					zz <- suppressWarnings(predict(my.fit,xrng,interval="confidence"))
					lines(xrng[,1],zz[,2],col="blue")
					lines(xrng[,1],zz[,3],col="blue")
				}
				if (is.na(my.cor[i,j]))
				{  # if not data, just make empty box
					plot(1:10,1:10,type='n',axes=FALSE)
				}
        box()
			}
			if (i > j)
			{
				plot(1:10,1:10,type='n',axes=FALSE)
				txt <- format(my.cor.round[i,j], nsmall=2)
				text(5,5,txt)
				box()
			}
		}
	}
	title(mytitle, outer=TRUE)
	# par(origpar)
	return(my.cor)
}
#-------------------------------------------------------------------------------

plot_catch_at_age_consistency <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
	# create plots of ages vs each other (correctly lagged) on log scale for catch by fleet
	cat.corr <- list()
  dat = mod$env$data
  rep = mod$rep
  n_years = length(mod$years)
	for (i in 1:dat$n_fleets)
	{
		title1 = paste("Catch for ",mod$input$fleet_names[i], " in ", mod$input$region_names[mod$input$data$fleet_regions[i]], sep="")

		# get catch at age
    catchob = dat$catch_paa[i,,] * dat$agg_catch[,i]/apply(dat$catch_paa[i,,] * dat$waa[dat$waa_pointer_fleets[i],1:n_years,],1,sum)
    catchpr = rep$pred_catch_paa[i,1:n_years,] * exp(rep$pred_log_catch[1:n_years,i])/apply(rep$pred_catch_paa[i,1:n_years,] * dat$waa[dat$waa_pointer_fleets[i],1:n_years,],1,sum)
		# replace zeros with NA and take logs
		cob <- rep0log(catchob)
		cpr <- rep0log(catchpr)

		# make cohorts
		cob.coh <- makecohorts(cob)
		cpr.coh <- makecohorts(cpr)

		# make the plots
    fn <- paste0(chartr(" ", "_", mod$input$fleet_names[i]), "_", chartr(" ", "_", mod$input$region_names[mod$input$data$fleet_regions[i]]))
		if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_obs_",fn, ".pdf")), family = fontfam, height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_obs_",fn, ".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
		cob.cor <- plotcoh(cob.coh,mytitle=paste(title1," Observed", sep=""),mod=mod)
		if(do.tex | do.png) dev.off()

		if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_pred_",fn, ".pdf")), family = fontfam, height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_pred_",fn, ".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
		cpr.cor <- plotcoh(cpr.coh,mytitle=paste(title1," Predicted", sep=""),mod=mod)
		if(do.tex | do.png) dev.off()

		cat.corr[[i]] <- list(cob.cor,cpr.cor)
	}
	return(cat.corr)
}
#revised

#-------------------------------------------------------------------------------

convert_survey_to_at_age <- function(mod)
{
	# takes West Coast style surveys and converts them to catch at age matrices
  dat = mod$env$data
  rep = mod$rep
	index.mats <- list(ob = list(), pr = list())
	weight.mat.counter <- 0
	for (i in 1:dat$n_indices)
	{
		if (sum(dat$use_index_paa[,i])>0)
		{  # used age composition for the index
			# get the aggregate index observed and predicted time series
			agg.ob <- dat$agg_indices[which(dat$use_index_paa[,i]==1),i]
			agg.pr <- exp(rep$pred_log_indices[which(dat$use_index_paa[,i]==1),i]) # bias corrected

			# get proportions for correct years and ages only
			props.ob <- dat$index_paa[i,which(dat$use_index_paa[,i]==1),]
			props.pr <- rep$pred_index_paa[i,which(dat$use_index_paa[,i]==1),]

      # figure out units for aggregate and proportions
			agg.units <- dat$units_indices[i]
			prp.units <- dat$units_index_paa[i]

      waa <- dat$waa[dat$waa_pointer_indices[i],which(dat$use_index_paa[,i]==1),]
			# get weight (matrix if necessary)
			if (agg.units==1 | prp.units==1)
			{  # either in weight
			}

			# create index.obs and pred based on which of the four possible combinations of units is used for this index
			if (agg.units==1 & prp.units==1)
			{  # both in weight
				index.ob <- agg.ob * props.ob / waa
				index.pr <- agg.pr * props.pr / waa
			}
			if (agg.units==1 & prp.units==2)
			{  # agg in weight, props in numbers
        index.ob = props.ob * agg.ob/apply(props.ob * waa,1,sum)
        index.pr = props.pr * agg.pr/apply(props.pr * waa,1,sum)
			}
			if (agg.units==2 & prp.units==1)
			{  # agg in numbers, props in weight
				# need to search for correct agg total in weight to result in observed agg total in number
				# for now just use simple approximation that agg.wt = sum(waa*prop) *ctot and then solve using both in weight approach
				agg.wt.ob <- apply((waa * props.ob),1,sum) * agg.ob
				agg.wt.pr <- apply((waa * props.pr),1,sum) * agg.pr
				index.ob <- agg.wt.ob * props.ob / waa
				index.pr <- agg.wt.pr * props.pr / waa
			}
			if (agg.units==2 & prp.units==2)
			{  # both in numbers
				index.ob <- agg.ob * props.ob
				index.pr <- agg.pr * props.pr
			}

			# put matrices into full year matrix (ages only for selected ages though)
			# need to do this to account for missing years of data interspersed in time series
			index.ob.full <- index.pr.full <- matrix(NA, nrow=length(mod$years), ncol=dat$n_ages)
			index.ob.full[dat$use_index_paa[,i]==1,] <- index.ob
			index.pr.full[dat$use_index_paa[,i]==1,] <- index.pr

			# save the results for this index
			index.mats$ob[[i]] <- index.ob.full
			index.mats$pr[[i]] <- index.pr.full
		}
		else
		{  # cannot use this index
			index.mats$ob[[i]] <- NA
			index.mats$pr[[i]] <- NA
		}
	}
	return(index.mats)
}
#revised 
#-------------------------------------------------------------------------------

plot_index_at_age_consistency <- function(mod, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
  dat = mod$env$data
  n_years = length(mod$years)
	# now loop through indices, check to make sure actually estimating proportions at age
	# also need to use only the age range selected in the proportions
	index.corr <- list()

	# convert the west coast style indices to catch at age matrices
	index.mats <- convert_survey_to_at_age(mod)

	# loop through all the indices
	for (ind in 1:dat$n_indices)
	{
		if (sum(dat$use_index_paa[,ind]) >0)
		{  # used age composition for the index
		  title1 <- paste(mod$input$index_names[ind], " in ", mod$input$region_names[mod$input$data$index_regions[ind]], sep="")

			# replace zeros with NA and take logs
			iob <- rep0log(index.mats$ob[[ind]])
			ipr <- rep0log(index.mats$pr[[ind]])

			# make cohorts
			iob.coh <- makecohorts(iob)
			ipr.coh <- makecohorts(ipr)


			# create labels for plot (if necessary)
			mylabels <- paste0("age-",mod$ages.lab, sep="")

			# make the plots
      fn <- paste0(chartr(" ", "_", mod$input$index_names[ind]), "_", chartr(" ", "_", mod$input$region_names[mod$input$data$index_regions[ind]]))

			if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_obs_",fn, ".pdf")), family = fontfam, height = 10, width = 10)
			if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_obs_",fn, ".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
			iob.cor <- plotcoh(iob.coh,mytitle=paste(title1," Observed", sep=""),mylabels,mod=mod)
			if(do.tex | do.png) dev.off()

			if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_pred_",fn, ".pdf")), family = fontfam, height = 10, width = 10)
			if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_pred_",fn, ".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
			ipr.cor <- plotcoh(ipr.coh,mytitle=paste(title1," Predicted", sep=""),mylabels,mod=mod)
			if(do.tex | do.png) dev.off()

			index.corr[[ind]] <- list(iob.cor,ipr.cor)
		}
	}
	return(index.corr)
}
#revised
#plot_index_at_age_consistency(ssm)

#-------------------------------------------------------------------------------

find_peak_age <- function(cohmat)
{
	# determines peak age within each cohort and replaces younger ages with NA
	ages <- seq(1,NCOL(cohmat))
	temp.sum <- apply(cohmat,1,sum, na.rm = TRUE)
	for (i in 1:NROW(cohmat))
	{
		if (temp.sum[i] > 0)
		{  # necessary to avoid cohorts that are all NA
			age.peak <- max(ages[cohmat[i,]==max(cohmat[i,],na.rm=TRUE)], na.rm=TRUE)  # find the peak age
			if (age.peak > 1) cohmat[i,1:(age.peak-1)] <- NA  # replace ages < peak.age with NA
		}
	}
	return(cohmat)
}
#-------------------------------------------------------------------------------

calc_Z_cohort <- function(cohmat)
{
	# calculate Z along cohort using linear regression and return point estimate and 80% confidence interval
	ages <- seq(1,NCOL(cohmat))
	z <- matrix(NA, nrow=NROW(cohmat), ncol=3)
	for (i in 1:NROW(cohmat))
	{
		if (length(cohmat[i,which(!is.na(cohmat[i,]))]) >= 2)
		{  # make sure there are at least 2 data points for point estimate
			z.fit <- lm(cohmat[i,]~ages)  # linear regression of cohort abundance vs age
			z[i,1] <- -1 * z.fit$coefficients[2]  # note change in sign for Z
			if (length(cohmat[i,!is.na(cohmat[i,])]) >= 3)
			{  # need at least 3 data points for CI
				z[i,2:3] <- -1 * rev(confint(z.fit, "ages", level=0.80))  # note change in sign and order for Z CI
			}
		}
	}
	return(z)
}
#-------------------------------------------------------------------------------

plot_catch_curves_for_catch <- function(mod, first.age=-999, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
	# create catch curve plots for catch by fleet
  origpar <- par(no.readonly = TRUE)
	lastyr <- max(mod$years)
  dat = mod$env$data
  rep = mod$rep
  n_years = length(mod$years)
	cohort <- (min(mod$years) - dat$n_ages-1):(lastyr+dat$n_ages-1)
	ages <- 1:dat$n_ages
  my.col = rep(viridisLite::viridis(n=5),50)
	#my.col <- rep(c("blue","red","green","orange","gray50"),50)
	for (i in 1:dat$n_fleets)
	{
		title1 = paste("Catch for " , mod$input$fleet_names[i], " in ", mod$input$region_names[mod$input$data$fleet_regions[i]], sep="")

		# get catch at age
    catchob = dat$catch_paa[i,,] * dat$agg_catch[,i]/apply(dat$catch_paa[i,,] * dat$waa[dat$waa_pointer_fleets[i],1:n_years,],1,sum)
    catchpr = rep$pred_catch_paa[i, 1:n_years,] * exp(rep$pred_log_catch[1:n_years,i])/apply(rep$pred_catch_paa[i,1:n_years,] * dat$waa[dat$waa_pointer_fleets[i],1:n_years,],1,sum)

		# replace zeros with NA and take logs
		cob <- rep0log(catchob)
		cpr <- rep0log(catchpr)

		# make cohorts
		cob.coh <- makecohorts(cob)
		cpr.coh <- makecohorts(cpr)

		# drop plus group
		cob.coh[,dat$n_ages] <- NA
		cpr.coh[,dat$n_ages] <- NA

		first.age.label <- 1
		if (first.age==1) title1 <- paste(title1," First Age = 1", sep="")

		# determine which ages to use for each cohort (default)
		if (first.age == -999)
    {
      cob.coh <- find_peak_age(cob.coh)
      cpr.coh <- find_peak_age(cpr.coh)
      first.age.label <- "find_peak"
      title1 <- paste0(title1," (Peak Age)")
		}

		# or drop youngest ages based on user control
		if (first.age > 1)
    {
      cob.coh[,1:(first.age-1)] <- NA
      cpr.coh[,1:(first.age-1)] <- NA
      title1 <- paste0(title1," First Age = ",first.age)
      first.age.label <- first.age
		}

		# compute Z by cohort
		z.ob <- calc_Z_cohort(cob.coh)
		z.pr <- calc_Z_cohort(cpr.coh)

		# make the plots
    fn <- paste0(chartr(" ", "_", mod$input$fleet_names[i]), "_", chartr(" ", "_", mod$input$region_names[mod$input$data$fleet_regions[i]]))

		if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_", fn,"_obs.pdf")), family = fontfam, height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_curves_", fn,"_obs.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
		par(oma=c(1,1,1,1),mar=c(4,4,1,0.5),mfrow=c(2,1))
		plot(cohort,cohort,type='n',ylim=range(c(cob.coh,cpr.coh),na.rm=TRUE),xlab="",ylab="Log(Catch)", main=paste0(title1," Observed"))
		grid(col = gray(0.7))
		for (j in 1:NROW(cob.coh))
		{
			lines(cohort[j]:(cohort[j]+dat$n_ages-1),cob.coh[j,],type='p',lty=1,pch=1:dat$n_ages,col="gray50")
			lines(cohort[j]:(cohort[j]+dat$n_ages-1),cob.coh[j,],type='l',lty=1,col=my.col[j])
		}
		Hmisc::errbar(cohort,z.ob[,1],z.ob[,3],z.ob[,2],xlab="Year Class",ylab="Z",ylim=range(c(z.ob,z.pr),na.rm=T))
		grid(col = gray(0.7))
		if(do.tex | do.png) dev.off()

		if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_", fn,"_pred.pdf")), family = fontfam, height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_curves_", fn,"_pred.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
		par(oma=c(1,1,1,1),mar=c(4,4,1,0.5),mfrow=c(2,1))
		plot(cohort,cohort,type='n',ylim=range(c(cob.coh,cpr.coh),na.rm=TRUE),xlab="",ylab="Log(Catch)", main=paste0(title1," Predicted"))
		grid(col = gray(0.7))
		for (j in 1:length(cob.coh[,1]))
		{
			lines(cohort[j]:(cohort[j]+dat$n_ages-1),cpr.coh[j,],type='p',lty=1,pch=1:dat$n_ages,col="gray50")
			lines(cohort[j]:(cohort[j]+dat$n_ages-1),cpr.coh[j,],type='l',lty=1,col=my.col[j])
		}
		Hmisc::errbar(cohort,z.pr[,1],z.pr[,3],z.pr[,2],xlab="Year Class",ylab="Z",ylim=range(c(z.ob,z.pr),na.rm=T))
		grid(col = gray(0.7))
		if(do.tex | do.png) dev.off()

    # write out .csv files for Z, one file for each fleet
		colnames(z.ob) <-c("Z.obs","low.80%", "high.80%")
		#write.csv(z.ob, file=paste(od,"Z.Ob.Fleet.",i,".csv", sep=""), row.names=cohort)

		colnames(z.pr) <-c("Z.pred","low.80%", "high.80%")
		#write.csv(z.pr, file=paste(od,"Z.Pr.Fleet.",i,".csv", sep=""), row.names=cohort)
	}
  if(!(do.tex | do.png)) par(origpar)
}
#plot_catch_curves_for_catch(ssm)

#-------------------------------------------------------------------------------

plot_catch_curves_for_index <- function(mod, first.age=-999, do.tex = FALSE, do.png = FALSE, fontfam="", res = 72, od)
{
	# create catch curve plots for each west coast style index
  origpar <- par(no.readonly = TRUE)
  lastyr <- max(mod$years)
  dat = mod$env$data
  rep = mod$rep
  ages = 1:dat$n_ages
  n_ages = dat$n_ages
  cohort <- (min(mod$years)-n_ages-min(ages)):(lastyr+n_ages-min(ages))
	my.col <- rep(c("blue","red","green","orange","gray50"),50)

	# convert the west coast style indices to catch at age matrices
	index.mats <- convert_survey_to_at_age(mod)

	# loop through all the indices
	for (ind in 1:dat$n_indices)
	{
		if (sum(dat$use_index_paa[,ind]) > 0)
		{  # used age composition for the index
		  title1 <- paste0(mod$input$index_names[ind], " in ", mod$input$region_names[mod$input$data$index_regions[ind]])
			# replace zeros with NA and take logs
			iob <- rep0log(index.mats$ob[[ind]])
			ipr <- rep0log(index.mats$pr[[ind]])

			# make cohorts
			iob.coh <- makecohorts(iob)
			ipr.coh <- makecohorts(ipr)

			# drop plus group
			iob.coh[,NCOL(iob.coh)] <- NA
			ipr.coh[,NCOL(iob.coh)] <- NA

			first.age.label <- 1
			if (first.age==1) title1 <- paste0(title1," First Age = 1")

			# determine which ages to use for each cohort (default)
			if (first.age == -999)
			{
				iob.coh <- find_peak_age(iob.coh)
				ipr.coh <- find_peak_age(ipr.coh)
				first.age.label <- "find_peak"
				title1 <- paste0(title1," (Peak Age)")
			}

			# or drop youngest ages based on user control
			if (first.age > min(ages))
			{
				iob.coh[,1:(first.age-min(ages))] <- NA
				ipr.coh[,1:(first.age-min(ages))] <- NA
				title1 <- paste0(title1," First Age = ",first.age)
				first.age.label <- first.age
			}

			# compute Z by cohort
			z.ob <- calc_Z_cohort(iob.coh)
			z.pr <- calc_Z_cohort(ipr.coh)

			# make the plots
		  fn <- paste0(chartr(" ", "_", mod$input$index_names[ind]), "_", chartr(" ", "_", mod$input$region_names[mod$input$data$index_regions[ind]]))
			if(!(all(is.na(iob.coh)) & all(is.na(ipr.coh))))
			{
			  if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_", fn ,"_obs.pdf")), family = fontfam, height = 10, width = 10)
			  if(do.png) png(filename = file.path(od, paste0("catch_curves_", fn ,"_obs.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
			  par(oma=c(1,1,1,1),mar=c(4,4,1,0.5),mfrow=c(2,1))
			  plot(cohort,cohort,type='n',ylim=range(c(iob.coh,ipr.coh),na.rm=T),xlab="",ylab="Log(Index)",
          main=paste0(title1," Observed"))
				grid(col = gray(0.7))
				for (i in 1:length(iob.coh[,1]))
				{
					lines(cohort[i]:(cohort[i]+n_ages-1),iob.coh[i,],type='p',lty=1,pch=1:n_ages,col="gray50")
					lines(cohort[i]:(cohort[i]+n_ages-1),iob.coh[i,],type='l',lty=1,col=my.col[i])
				}
				Hmisc::errbar(cohort,z.ob[,1],z.ob[,3],z.ob[,2],xlab="Year Class",ylab="Z",ylim=range(c(z.ob,z.pr),na.rm=TRUE))
				grid(col = gray(0.7))
				if(do.tex | do.png) dev.off()

				if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_", fn ,"_pred.pdf")), family = fontfam, height = 10, width = 10)
				if(do.png) png(filename = file.path(od, paste0("catch_curves_", fn ,"_pred.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
				par(oma=c(1,1,1,1),mar=c(4,4,1,0.5),mfrow=c(2,1))
				plot(cohort,cohort,type='n',ylim=range(c(iob.coh,ipr.coh),na.rm=T),xlab="",ylab="Log(Index)", main=paste0(title1," Predicted"))
				grid(col = gray(0.7))
				for (i in 1:NROW(iob.coh))
				{
					lines(cohort[i]:(cohort[i]+n_ages-1),ipr.coh[i,],type='p',lty=1,pch=1:n_ages,col="gray50")
					lines(cohort[i]:(cohort[i]+n_ages-1),ipr.coh[i,],type='l',lty=1,col=my.col[i])
				}
				Hmisc::errbar(cohort,z.pr[,1],z.pr[,3],z.pr[,2],xlab="Year Class",ylab="Z",ylim=range(c(z.ob,z.pr),na.rm=TRUE))
				grid(col = gray(0.7))
				if(do.tex | do.png) dev.off()
			}

			# write out .csv files for Z, one file for each fleet
			colnames(z.ob) <-c("Z.obs","low.80%", "high.80%")
			#write.csv(z.ob, file=paste(od,"Z.Ob.Index.",ind,".csv", sep=""), row.names=cohort)

			colnames(z.pr) <-c("Z.pred","low.80%", "high.80%")
			#write.csv(z.pr, file=paste(od,"Z.Pr.Index.",ind,".csv", sep=""), row.names=cohort)
		}
	}   # end loop over n_indices
	if(!(do.tex | do.png)) par(origpar)
}
#revised
#plot_catch_curves_for_index(ssm)
#-------------------------------------------------------------------------------
plot.ecov.diagnostic <- function(mod, use.i, plot.pad = FALSE, do.tex = FALSE, do.png = FALSE, fontfam="", od){
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  ecov.pred = mod$rep$Ecov_x
  ecov.obs = dat$Ecov_obs[1:dat$n_years_Ecov,,drop=F]
  
  years <- seq(from=mod$input$years_Ecov[1], by=1, length.out=NROW(ecov.obs))
  years_full <- seq(from=mod$input$years_Ecov[1], by=1, length.out=NROW(ecov.pred))#dat$n_years_Ecov+dat$n_years_proj_Ecov)

  # ecov.obs.sig = mod$rep$Ecov_obs_sigma # Ecov_obs_sigma now a derived quantity in sdrep
  if(class(mod$sdrep)[1] == "sdreport"){
    sdrep = summary(mod$sdrep)
  } else {
    sdrep = mod$sdrep
  }
  ecov.obs.sig = mod$rep$Ecov_obs_sigma # Ecov_obs_sigma is filled with fixed, or estimated values (fe or re) for each covariate depending on the respective options
  # if("Ecov_obs_logsigma" %in% names(mod$env$par)){
  #   ecov.obs.sig = matrix(exp(sdrep[rownames(sdrep) %in% "Ecov_obs_logsigma",1]), ncol=dat$n_Ecov) # all the same bc obs_sig_var --> 0
  #   if(dim(ecov.obs.sig)[1] == 1) ecov.obs.sig = matrix(rep(ecov.obs.sig, dim(ecov.obs)[1]), ncol=dat$n_Ecov, byrow=T)
  # } else {
  #   ecov.obs.sig = exp(mod$input$par$Ecov_obs_logsigma)
  # }
  ecov.use = dat$Ecov_use_obs[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig = ecov.obs.sig[1:dat$n_years_Ecov,,drop=F]
  ecov.obs.sig[ecov.use == 0] <- NA
  ecov.pred.se = matrix(sdrep[rownames(sdrep) %in% "Ecov_x",2], ncol=dat$n_Ecov)

  # default: don't plot the padded entries that weren't used in ecov likelihood
  if(!plot.pad) ecov.obs[ecov.use == 0] <- NA

  ecov.res = (ecov.obs - ecov.pred[1:dat$n_years_Ecov,]) / ecov.obs.sig # standard residual (obs - pred)

  if(!missing(use.i)) ecovs <- use.i
  else ecovs <- 1:dat$n_Ecov
  plot.colors = viridisLite::viridis(n = dat$n_Ecov) #mypalette(dat$n_Ecov)
  efns <- chartr(" ", "_", mod$input$Ecov_names)
  for (i in ecovs)
  {
    if(do.tex) cairo_pdf(file.path(od, paste0(efns[i],"_diagnostic.pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0(efns[i],'_diagnostic.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)

    m <- rbind(c(1,1), c(2,3))
    layout(m)
    par(mar=c(4,4,2,0), oma=c(0,0,0.5,0.5))

    ecov.pred.low <- ecov.pred[,i] - 1.96 * ecov.pred.se[,i]
    ecov.pred.high <- ecov.pred[,i] + 1.96 * ecov.pred.se[,i]
    ecov.low <- ecov.obs[,i] - 1.96 * ecov.obs.sig[,i]
    ecov.high <- ecov.obs[,i] + 1.96 * ecov.obs.sig[,i]
    y.min <- ifelse(min(ecov.low,na.rm=T) < 0, 1.1*min(ecov.low,na.rm=T), 0.9*min(ecov.low,na.rm=T))
    y.max <- ifelse(max(ecov.high,na.rm=T) < 0, 0.9*max(ecov.high,na.rm=T), 1.1*max(ecov.high,na.rm=T))
    if(max(ecov.pred[,i],na.rm=T) > y.max) y.max <- max(ecov.pred[,i],na.rm=T)
    if(min(ecov.pred[,i],na.rm=T) < y.min) y.min <- min(ecov.pred[,i],na.rm=T)
    plot(years_full, ecov.pred[,i], type='n', xlab="Year", ylab=mod$input$Ecov_names[i],
         ylim=c(y.min, y.max))
    polygon(c(years_full,rev(years_full)), c(ecov.pred.low, rev(ecov.pred.high)), col=adjustcolor(plot.colors[i], alpha.f=0.4), border = "transparent")
    arrows(years, ecov.low, years, ecov.high, length=0)
    points(years, ecov.obs[,i], pch=19)
    lines(years_full, ecov.pred[,i], col=plot.colors[i], lwd=3)
    title (paste0("Ecov ",i, ": ",mod$input$Ecov_names[i]), outer=T, line=-1)
    if(length(years_full) > length(years)) abline(v=tail(years,1), lty=2)

    plot(years, ecov.res[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Std. Residual")
    abline(h=0)

    hist(ecov.res[,i], breaks=10, plot=T, xlab="Std. Residual", ylab="Probability Density", freq=F, main=NULL)
    if(do.tex | do.png) dev.off() else par(origpar)
  }
}
#revised

#-------------------------------------------------------------------------------
# 2D tile plot by age and year (e.g. selAA, MAA)
plot.tile.age.year <- function(mod, type="selAA", do.tex = FALSE, do.png = FALSE, fontfam="", od){
  dat = mod$env$data
  rep = mod$rep
  years = mod$years
  n_years = length(years)
  n_ages = dat$n_ages
  ages <- 1:n_ages
  ages.lab = 1:n_ages
  if(!is.null(mod$ages.lab)) ages.lab = mod$ages.lab

  # selAA for all blocks using facet_wrap
  if(type=="selAA"){ 
    n_selblocks <- length(rep$selAA)
    sel_mod <- c("age-specific","logistic","double-logistic","decreasing-logistic")[dat$selblock_models]
    sel_re <- c("no","IID","AR1","AR1_y","2D AR1")[dat$selblock_models_re]
    df.selAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
    colnames(df.selAA) <- c(paste0("Age_",1:n_ages),"Year","Block")
    block.names <- paste0("Block ",1:n_selblocks,": ", sel_mod,"\n(",sel_re," random effects)")
    block.fleets.indices <- lapply(1:n_selblocks, function(x){
      y <- dat$selblock_pointer_fleets
      z <- matrix(as.integer(y == x), NROW(y), NCOL(y))
      fleet_ind <- apply(z,2,any)
      out <- mod$input$fleet_names[which(fleet_ind)]
      y <- dat$selblock_pointer_indices
      z <- matrix(as.integer(y == x), NROW(y), NCOL(y))
      index_ind <- apply(z,2,any)
      out <- c(out, mod$input$index_names[which(index_ind)])
    })
    include.selblock <- sapply(block.fleets.indices, length) > 0
    for(i in 1:n_selblocks) if(include.selblock[i]){
      block.names[i] <- paste0(block.names[i], "\n", paste(block.fleets.indices[[i]], collapse = ", "))
    }
    for(i in 1:n_selblocks) if(include.selblock[i]){
      tmp = as.data.frame(rep$selAA[[i]])
      tmp$Year <- years
      colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
      tmp$Block = block.names[i]
      df.selAA <- rbind(df.selAA, tmp)
    }
    df.plot <- df.selAA %>% tidyr::pivot_longer(-c(Year,Block),
              names_to = "Age", 
              names_prefix = "Age_",
              names_ptypes = list(Age = character()),
              values_to = "Selectivity")
    df.plot$Age <- as.factor(as.integer(df.plot$Age))
    levels(df.plot$Age) = ages.lab
    df.plot$Block <- factor(as.character(df.plot$Block), levels=block.names[include.selblock])
    fn <- "SelAA_tile"
    if(do.tex) cairo_pdf(file.path(od, paste0(fn, ".pdf")), family = fontfam, height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0(fn, ".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
      print(ggplot2::ggplot(df.plot, ggplot2::aes(x=Year, y=Age, fill=Selectivity)) + 
        ggplot2::geom_tile() +
        ggplot2::scale_x_continuous(expand=c(0,0)) +
        ggplot2::scale_y_discrete(expand=c(0,0)) + #, breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +        
#        ggplot2::scale_y_continuous(expand=c(0,0), breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +        
        ggplot2::theme_bw() + 
        ggplot2::facet_wrap(~Block, dir="v") +
        viridis::scale_fill_viridis())
    if(do.tex | do.png) dev.off()
  }

  # MAA
  if(type=="MAA"){ 
    if(mod$env$data$n_years_proj>0){
      years_full = mod$years_full
    } else {
      years_full = years
    }
    sfns <- chartr(" ", "_", mod$input$stock_names)
    rfns <- chartr(" ", "_", mod$input$region_names)
    df.plot <- data.frame(Stock = character(), Region = character(), Year = integer(), Age = integer(), Deviation = numeric())
    for(s in 1:mod$env$data$n_stocks) for(r in 1:mod$env$data$n_regions){
      df.MAA <- as.data.frame(rep$MAA[s,r,,])
      colnames(df.MAA) <- paste0("Age_",1:n_ages)
      df.MAA <- cbind.data.frame(Stock = mod$input$stock_names[s], Region = mod$input$region_names[r], Year = years_full, df.MAA)
      temp <- df.MAA %>% tidyr::pivot_longer(tidyr::starts_with("Age"),
                names_to = "Age", 
                names_prefix = "Age_",
                names_ptypes = list(Age = character()),
                values_to = "M") %>% as.data.frame()
      df.plot <- rbind(df.plot, temp)
    }
    df.plot$Age <- as.factor(as.integer(df.plot$Age))
    levels(df.plot$Age) = ages.lab
    fn <- paste0("MAA_tile")

    if(do.tex) cairo_pdf(file.path(od, paste0(fn, ".pdf")), family = fontfam, height = 5, width = 10)
    if(do.png) png(filename = file.path(od, paste0(fn, ".png")), width = 10*144, height = 5*144, res = 144, pointsize = 12, family = fontfam)
      plt <- ggplot2::ggplot(df.plot, ggplot2::aes(x=Year, y=Age, fill=M)) + 
      ggplot2::geom_tile() +
      ggplot2::scale_x_continuous(expand=c(0,0)) +
      ggplot2::scale_y_discrete(expand=c(0,0)) + 
      ggplot2::theme_bw() + 
      ggplot2::facet_grid(Stock~Region, labeller = ggplot2::label_both) +
      ggplot2::ggtitle("MAA") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      viridis::scale_fill_viridis()
      if(mod$env$data$n_years_proj>0) plt  <- plt + ggplot2::geom_vline(xintercept = tail(years,1), linetype = 2)
      print(plt)
    if(do.tex | do.png) dev.off()
  }

  # NAA_devs
  if(type=="NAA_devs"){ 
    if(mod$env$data$n_years_proj>0){
      years_full = mod$years_full
    } else {
      years_full = years
    }
    df.plot <- data.frame(Stock = character(), Region = character(), Year = integer(), Age = integer(), Deviation = numeric())
    for(s in 1:mod$env$data$n_stocks) for(r in 1:mod$env$data$n_regions){
      df.NAA <- as.data.frame(rep$NAA_devs[s,r,,])
      colnames(df.NAA) <- paste0("Age_",1:n_ages)
      df.NAA <- cbind.data.frame(Stock = mod$input$stock_names[s], Region = mod$input$region_names[r], Year = years_full, df.NAA)
      temp <- df.NAA %>% tidyr::pivot_longer(tidyr::starts_with("Age"),
                names_to = "Age", 
                names_prefix = "Age_",
                names_ptypes = list(Age = character()),
                values_to = "Deviation") %>% as.data.frame()
      df.plot <- rbind(df.plot, temp)
    }
    df.plot$Age <- as.factor(as.integer(df.plot$Age))
    levels(df.plot$Age) = ages.lab
    fn <- paste0("NAA_dev_tile")

    if(do.tex) cairo_pdf(file.path(od, paste0(fn, ".pdf")), family = fontfam, height = 5, width = 10)
    if(do.png) png(filename = file.path(od, paste0(fn, ".png")), width = 10*144, height = 5*144, res = 144, pointsize = 12, family = fontfam)
      plt <- ggplot2::ggplot(df.plot, ggplot2::aes(x=Year, y=Age, fill=Deviation)) + 
        ggplot2::geom_tile() +
        ggplot2::scale_x_continuous(expand=c(0,0)) +
        ggplot2::scale_y_discrete(expand=c(0,0)) + 
        ggplot2::theme_bw() + 
        ggplot2::facet_grid(Stock~Region, labeller = ggplot2::label_both) +
        ggplot2::ggtitle("NAA deviations") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        viridis::scale_fill_viridis()
      if(mod$env$data$n_years_proj>0) plt  <- plt + ggplot2::geom_vline(xintercept = tail(years,1), linetype = 2)
      print(plt)
    if(do.tex | do.png) dev.off()
  }
}  

#revised

#pdf of a univariate logit-normal with any min and max
dlogitnorm = function(p,mu,sd,min,max) {
  logitp = log((p-min)/(max-p))
  (exp(-(logitp - mu)^2/(2 * sd^2))/(sd * sqrt(2*pi))) * (max-p)/((p-min) * (max-p))
}

plot_q_prior_post = function(mod, do.tex = F, do.png = F, fontfam="", od){
  origpar <- par(no.readonly = TRUE)
  ind = which(mod$input$data$use_q_prior == 1)
  if(length(ind) & "sdrep" %in% names(mod)){
    logit_q = cbind(as.list(mod$sdrep, "Est")$q_prior_re, as.list(mod$sdrep, "Std")$q_prior_re)
    ht = 10
    wd = 10*length(ind)
    priorq = approx_postq = list()
    if(do.tex) cairo_pdf(file.path(od, "q_prior_post.pdf"), family = fontfam, height = ht, width = wd)
    if(do.png) png(filename = file.path(od, "q_prior_post.png"), width = wd*144, height = ht*144, res = 144, pointsize = 12, family = fontfam)
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(1,length(ind)))
    pal = viridisLite::viridis(n=2)
    for(i in ind) {
      qmax = mod$input$data$q_upper[i]
      x = seq(0.001,qmax,0.001)
      y = log(x-0) - log((qmax-x))
      approx_postq[[i]] = dnorm(y, logit_q[i,1], logit_q[i,2])
      priorq[[i]] = dlogitnorm(x, mod$parList$logit_q[i], 0.3, 0, 1000) 
      maxx = max(x[which(approx_postq[[i]] > 1e-5)], x[which(priorq[[i]] > 1e-5)], na.rm = T)
      plot(x,priorq[[i]], type = 'l', xlab = "q", ylab = "pdf", col = pal[1], lwd = 2, ylim = c(0,max(approx_postq[[i]],priorq[[i]], na.rm =T)), xlim = c(0,maxx))
      lines(x,approx_postq[[i]], col = pal[2], lwd = 2)
      ci = qmax/(1+ exp(-(logit_q[i,1] + c(-1,1)*qnorm(0.975) * logit_q[i,2])))
      abline(v= ci, lty = 2)
      legend("topright", legend = c("prior", "approx. posterior", "95% CI"), col = c(pal, "black"), lty = c(1,1,2), lwd = 2)
      mtext(paste0("Index ", ind), side = 3, line = 1, outer = F, cex = 1.5)
    }
    if(do.tex | do.png) dev.off() else par(origpar)
  }  
}

plot_q = function(mod, do.tex = F, do.png = F, fontfam = '', od){
  origpar <- par(no.readonly = TRUE)
  yrs <- 1:length(mod$years)
  q = t(mod$input$data$q_lower + (mod$input$data$q_upper - mod$input$data$q_lower)/(1+exp(-t(mod$rep$logit_q_mat[yrs,,drop = FALSE]))))
  if(do.tex) cairo_pdf(file.path(od, "q_time_series.pdf"), family = fontfam, height = 10, width = 10)
  if(do.png) png(filename = file.path(od, "q_time_series.png"), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
  par(mar=c(4,4,3,1), oma=c(1,1,1,15))
  pal = viridisLite::viridis(n=mod$input$data$n_indices)
  ymax = max(q, na.rm = TRUE)
  if("sdrep" %in% names(mod)){
    if("q_re" %in% mod$input$random) se = as.list(mod$sdrep, "Std. Error", report=TRUE)$logit_q_mat[yrs,,drop = FALSE]
    else se = t(matrix(as.list(mod$sdrep, "Std. Error")$logit_q, nrow = NCOL(mod$rep$logit_q_mat), 
      ncol = NROW(mod$rep$logit_q_mat[yrs,,drop = FALSE])))
    logit_q_lo = mod$rep$logit_q_mat[yrs,,drop = FALSE] - qnorm(0.975)*se
    logit_q_hi = mod$rep$logit_q_mat[yrs,,drop = FALSE] + qnorm(0.975)*se
    q_lo = t(mod$input$data$q_lower + (mod$input$data$q_upper - mod$input$data$q_lower)/(1+exp(-t(logit_q_lo))))
    q_hi = t(mod$input$data$q_lower + (mod$input$data$q_upper - mod$input$data$q_lower)/(1+exp(-t(logit_q_hi))))
    ymax = max(c(q,q_hi), na.rm = TRUE)
  }
  plot(mod$years, q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,ymax), ylab = "q", xlab = "Year")
  for( i in 1:mod$input$data$n_indices){
    lines(mod$years, q[,i], lwd = 2, col = pal[i])
    if("sdrep" %in% names(mod)){
      polygon(c(mod$years,rev(mod$years)), c(q_lo[,i],rev(q_hi[,i])), col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
    }
  }
  leg <- paste0(mod$input$index_names, " ", mod$input$region_names[mod$input$data$index_regions])
  legend("right", legend = leg, col = pal, lty = 1, xpd = NA, inset = c(-0.5,0), bty = "n", lwd = 2)
  if(do.tex | do.png) dev.off() else par(origpar)
}

#NOT DONE YET
plot_mu = function(mod, do.tex = F, do.png = F, fontfam = '', od){
  #only call if n_regions=2
  origpar <- par(no.readonly = TRUE)
  dat <- mod$input$data
  ymax = max(mod$rep$mu, na.rm = TRUE)
  if(data$n_regions == 2) if(sum(dat$can_move)>0){
    if("sdrep" %in% names(mod)){
      se <- as.list(mod$sdrep, "Std. Error", report=TRUE)$trans_mu_sdrep
    } else se <- NULL

    for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) {
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
      if(data$mu_model[r,rr] %in% c(2,4,6,8,10,12,14,16)) na <- data$n_ages
      else na <- 1
      modify <- matrix(modify, nt, ns)
      if(do.tex) cairo_pdf(file.path(od, "mu_", r,"_", k, "_time_series.pdf"), family = fontfam, height = 10, width = 10)
      if(do.png) png(filename = file.path(od, "mu_", r,"_", k, "_time_series.png"), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = fontfam)
      par(mar=c(4,4,3,1), oma=c(1,1,1,15))
      pal = viridisLite::viridis(n=nt+ns)
      plot(mod$years, q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,ymax), ylab = "q", xlab = "Year")
      for(s in 1:ns) for(a in 1:na) for(t in 1:nt) if(dat$can_move[s,t,r,rr]) {
        lines(mod$years, mod$rep$mu[s,a,t,r,rr], lwd = 2, col = pal[t + (s-1)*nt])
        if(!is.null(se)) {
          mu_lo <- mu_hi <- NULL #NEED TO SEE IF YOU CAN transform logit CIs for n_regions>2
          polygon(c(mod$years,rev(mod$years)), c(q_lo[,i],rev(q_hi[,i])), col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
        }
      }
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }
  # q = t(mod$input$data$q_lower + (mod$input$data$q_upper - mod$input$data$q_lower)/(1+exp(-t(mod$rep$logit_q_mat))))
  # ymax = max(q, na.rm = TRUE)
  # if("sdrep" %in% names(mod)){
  #   if("q_re" %in% mod$input$random) se = as.list(mod$sdrep, "Std. Error", report=TRUE)$logit_q_mat
  #   else se = t(matrix(as.list(mod$sdrep, "Std. Error")$logit_q, nrow = NCOL(mod$rep$logit_q_mat), 
  #     ncol = NROW(mod$rep$logit_q_mat)))
  #   logit_q_lo = mod$rep$logit_q_mat - qnorm(0.975)*se
  #   logit_q_hi = mod$rep$logit_q_mat + qnorm(0.975)*se
  #   q_lo = t(mod$input$data$q_lower + (mod$input$data$q_upper - mod$input$data$q_lower)/(1+exp(-t(logit_q_lo))))
  #   q_hi = t(mod$input$data$q_lower + (mod$input$data$q_upper - mod$input$data$q_lower)/(1+exp(-t(logit_q_hi))))
  #   ymax = max(q_hi, na.rm = TRUE)
  # }
  # plot(mod$years, q[,1], type = 'n', lwd = 2, col = pal[1], ylim = c(0,ymax), ylab = "q", xlab = "Year")
  # for( i in 1:mod$input$data$n_indices){
  #   lines(mod$years, q[,i], lwd = 2, col = pal[i])
  #   if("sdrep" %in% names(mod)){
  #     polygon(c(mod$years,rev(mod$years)), c(q_lo[,i],rev(q_hi[,i])), col=adjustcolor(pal[i], alpha.f=0.4), border = "transparent")
  #   }
  # }
  # leg <- paste0(mod$input$index_names, " ", mod$input$region_names[mod$input$data$index_regions])
  # legend("right", legend = leg, col = pal, lty = 1, xpd = NA, inset = c(-0.5,0), bty = "n", lwd = 2)
}

sci_note <- function(x, cols=1:4){
  for(i in cols){
    temp <- formatC(x[,i], format = "e")
    temp <- strsplit(temp, "e")
    exps <- suppressWarnings(as.integer(sapply(temp, function(x) x[2])))
    coefs <- suppressWarnings(as.numeric(sapply(temp, function(x) x[1])))
    coefs <- formatC(coefs, format = "f", digits = 3)
    do.notation <- which(exps < -3)
    NA.ind <- which(is.na(x[,i]))
    x[,i] <- formatC(round(x[,i],3), format = "f", digits = 3)
    x[do.notation,i] <- paste0(coefs[do.notation], "\\times 10^{", exps[do.notation], "}") 
    x[,i] <- paste0("$", x[,i], "$")
    x[NA.ind,i] <- NA
  }
  return(x)
}

#revised

