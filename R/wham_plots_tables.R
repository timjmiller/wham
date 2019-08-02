plot.osa.residuals <- function(mod, do.tex=FALSE, do.png=FALSE, res=72, od){
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  if("logcatch" %in% mod$osa$type){
    dat <- subset(mod$osa, type=="logcatch")
    n.fleets <- length(table(dat$fleet))
    plot.colors = mypalette(n.fleets)
    for(f in 1:n.fleets){
      tmp <- subset(dat, fleet==names(table(dat$fleet))[f])
      if(do.tex) cairo_pdf(file.path(od, paste0("OSAresid_catch_4panel_fleet",f,".pdf")), family = "Times", height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSAresid_catch_4panel_fleet",f,'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = "Times")
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(mod$osa$residual, na.rm=TRUE)))
      ylims <- c(-ylim.max, ylim.max)

      # 1. trend vs. year
      plot(years, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 2. trend vs. fitted val
      plot(log(mod$rep$pred_catch[,f]), tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Log(Predicted Catch)", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 3. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(tmp$residual, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)

      # 4. QQ plot modified from car:::qqPlot.default
      ord.x <- tmp$residual[order(tmp$residual)]
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

      title (paste0("OSA residual diagnostics: Fleet ",f), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }

  if("logindex" %in% mod$osa$type){
    dat <- subset(mod$osa, type=="logindex")
    n.fleets <- length(table(dat$fleet))
    plot.colors = mypalette(n.fleets)
    for(f in 1:n.fleets){
      tmp <- subset(dat, fleet==names(table(dat$fleet))[f])
      if(do.tex) cairo_pdf(file.path(od, paste0("OSAresid_catch_4panel_index",f,".pdf")), family = "Times", height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSAresid_catch_4panel_index",f,'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = "Times")
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(mod$osa$residual, na.rm=TRUE)))
      ylims <- c(-ylim.max, ylim.max)

      # 1. trend vs. year
      plot(years, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 2. trend vs. fitted val
      plot(log(mod$rep$pred_indices[,f]), tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Log(Predicted Index)", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 3. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(tmp$residual, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)

      # 4. QQ plot modified from car:::qqPlot.default
      ord.x <- tmp$residual[order(tmp$residual)]
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

      title (paste0("OSA residual diagnostics: Index ",f), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }

  if(!all(mod$env$data$Ecov_model == 0)){
    dat <- subset(mod$osa, type=="ecov")
    n.fleets <- length(table(dat$fleet))
    plot.colors = mypalette(n.fleets)
    for(f in 1:n.fleets){
      tmp <- subset(dat, fleet==names(table(dat$fleet))[f])
      tmp$year <- seq(mod$env$data$year1_Ecov, by=1, length.out=mod$env$data$n_years_Ecov)
      tmp$pred <- mod$rep$Ecov_x[,f]
      tmp <- subset(tmp, !is.nan(dat$residual))
      if(do.tex) cairo_pdf(file.path(od, paste0("OSAresid_ecov_4panel_",f,".pdf")), family = "Times", height = 10, width = 10)
      if(do.png) png(filename = file.path(od, paste0("OSAresid_ecov_4panel_",f,'.png')), width = 10*res, height = 10*res, res = res, pointsize = 12, family = "Times")
      par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))

      # set plot lims using max residual for any component (easier to compare if all the same)
      ylim.max <- max(abs(range(mod$osa$residual, na.rm=TRUE)))
      ylims <- c(-ylim.max, ylim.max)

      # 1. trend vs. year
      plot(tmp$year, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab="Year", ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 2. trend vs. fitted val
      plot(tmp$pred, tmp$residual, type='p', col=plot.colors[f], pch=19, xlab=paste0("Predicted ", mod$env$data$Ecov_label[f]), ylab="OSA Residuals",
           ylim=ylims)
      abline(h=0, col=plot.colors[f], lwd=2)

      # 3. histogram
      xfit<-seq(-ylim.max, ylim.max, length=100)
      yfit<-dnorm(xfit)
      hist(tmp$residual, ylim=c(0,1.05*max(yfit)), xlim=ylims, plot=T, xlab="OSA Residuals", ylab="Probability Density", col=plot.colors[f], freq=F, main=NULL, breaks="scott")
      lines(xfit, yfit)

      # 4. QQ plot modified from car:::qqPlot.default
      ord.x <- tmp$residual[order(tmp$residual)]
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

      title (paste0("OSA residual diagnostics: Ecov ",f," (",mod$env$data$Ecov_label[f],")"), outer=T, line=-1)
      if(do.tex | do.png) dev.off() else par(origpar)
    }
  }
}

mypalette = function(n){
  palette.fn <- colorRampPalette(c("dodgerblue","green","red"), space = "Lab")
  palette.fn(n)
}

fit.summary.text.plot.fn <- function(mod){
  acm = c("Multinomial", "Dirichlet-multinomial", "Dirichlet", "ZI-logistic normal(1)","logistic normal(1)","ZI-logistic normal(2)","logistic normal(2)")
  selmods = c("Age-specific", "Logistic(+)", "Double-Logistic", "Logistic(-)")
  recs <- c("Random walk","Random about mean","Bev-Holt","Ricker")
  env.mod <- c("RW", "AR1")
  env.where <- c('Recruitment','Growth','Mortality')
  env.how <- c("controlling", "limiting", "lethal", "masking", "directive")
  fleet_selblocks = lapply(1:mod$env$data$n_fleets, function(x) unique(mod$env$data$selblock_pointer_fleets[,x]))
  index_selblocks = lapply(1:mod$env$data$n_indices, function(x) unique(mod$env$data$selblock_pointer_indices[,x]))
  plot(1:10,1:10,type='n',axes=F,xlab="",ylab="")
  nl = 10
  text(5,nl <- nl-0.5,mod$model_name)
  text(5,nl <- nl-0.5,paste0("Model years: ", min(mod$years), "-", max(mod$years)))
  text(5,nl <- nl-0.5,paste0("Number of fleets: ", mod$env$data$n_fleets))
  text(5,nl <- nl-0.5, paste0("Fleet Age Comp Models: ", paste(acm[mod$env$data$age_comp_model_fleets], collapse = ", ")))
  text(5,nl <- nl-0.5,paste0("Number of indices: ", mod$env$data$n_indices))
  text(5,nl <- nl-0.5, paste0("Index Age Comp Models: ", paste(acm[mod$env$data$age_comp_model_indices], collapse = ", ")))
  text(5,nl <- nl-0.5,paste0("Recruitment model: ", recs[mod$env$data$recruit_model]))
  if(!all(mod$env$data$Ecov_model == 0)){
    for(ec in 1:mod$env$data$n_Ecov) text(5,nl <- nl-0.5, paste0("Environmental effect ", ec,": ", mod$env$data$Ecov_label[ec]," (",env.mod[mod$env$data$Ecov_model[ec]],") on ",env.where[mod$env$data$Ecov_where[ec]], " (", env.how[mod$env$data$Ecov_how[ec]],")"))
  } else {
    text(5,nl <- nl-0.5, "Environmental effects: none")
  }
  text(5,nl <- nl-0.5,paste0("Number of Selectivity blocks: ", mod$env$data$n_selblocks))
  text(5,nl <- nl-0.5, paste0("Selectivity Block Types: ", paste(selmods[mod$env$data$selblock_models], collapse = ", ")))
  for(i in 1:length(fleet_selblocks)) text(5,nl <- nl-0.5, paste0("Fleet ", i, " Selectivity Blocks: ", paste(fleet_selblocks[i], collapse = ", ")))
  for(i in 1:length(index_selblocks)) text(5,nl <- nl-0.5, paste0("Index ", i, " Selectivity Blocks: ", paste(index_selblocks[i], collapse = ", ")))
  text(5,nl <- nl-0.5,paste0("WHAM run on ",format(mod$date, usetz = TRUE), " in directory ", mod$dir))

  if(mod$is_sdrep)
  {
    text(5,nl <- nl-0.5,paste0("sdreport() performed",
      ifelse(mod$na_sdrep, ", but with NAs for some variance estimates.", " with all variance estimates provided.")))
  }
  else
  {
    text(5,nl <- nl-0.5,"Warning: run did not provide pos-def Hessian or sdreport() not performed", col="red")
  }
  mgind = which(abs(mod$final_gradient) == max(abs(mod$final_gradient)))
  text(5,nl <- nl-0.5,paste0("Maximum absolute gradient: ", names(mod$opt$par)[mgind]," ", format(mod$final_gradient[mgind],digits =4)))
  text(5,nl <- nl-0.5,paste0("Number of fixed effects = ",length(mod$opt$par),", Number of random effects = ", length(mod$env$random)))

  return()
}

plot.ll.table.fn <- function(mod,plot.colors){
  par(mfrow=c(1,1) )

  npar <- length(mod$opt$par)
  lls = mod$rep[c(grep("nll",names(mod$rep)), grep("lprior_b", names(mod$rep)))]
  ll.names = names(lls)
  #n.like <- length(lls)
  n_fleets = mod$env$data$n_fleets
  n_indices = mod$env$data$n_indices
  obs.lls = lls[names(lls) %in% c("nll_agg_catch", "nll_catch_acomp", "nll_agg_indices", "nll_index_acomp")]
  obs.lls = lapply(obs.lls, function(x) apply(x,2,sum))
  names(obs.lls$nll_agg_catch) = paste0("Fleet ", 1:n_fleets, " Catch")
  names(obs.lls$nll_catch_acomp) = paste0("Fleet ", 1:n_fleets, " Age Comp")
  names(obs.lls$nll_agg_indices) = paste0("Index ", 1:n_indices, " Catch")
  names(obs.lls$nll_index_acomp) = paste0("Index ", 1:n_indices, " Age Comp")
  names(obs.lls) = NULL
  obs.lls = unlist(obs.lls)
  n.obs.ll = length(obs.lls)
  obs.dists = character(n.obs.ll)
  names(obs.dists) = names(obs.lls)
  obs.dists[grep("Catch", names(obs.lls))] = "log(x) ~ Gaussian"
  acm = c("Multinomial", "Dirichlet-multinomial", "Dirichlet", "ZI-logistic normal(1)","logistic normal(1)","ZI-logistic normal(2)","logistic normal(2)")
  obs.dists[paste0("Fleet ", 1:n_fleets, " Age Comp")] = paste0("x ~ ", acm[mod$env$data$age_comp_model_fleets])
  obs.dists[paste0("Index ", 1:n_indices, " Age Comp")] = paste0("x ~ ", acm[mod$env$data$age_comp_model_fleets])

  proc.lls = lls[names(lls) %in% c("nll_M", "nll_NAA", "nll_recruit", "lprior_b")]
  names(proc.lls) = c("M", "NAA", "recruit", "W_b_M")[match(names(proc.lls),c("nll_M", "nll_NAA", "nll_recruit", "lprior_b"))]
  proc.lls = unlist(lapply(proc.lls, sum))
  n.proc.ll = length(proc.lls)
  proc.dists = rep("log(x) ~ Gaussian", n.proc.ll)
  likes = -c(obs.lls, proc.lls)
  n.likes = length(likes)
  my.range <- 1.2*range(likes)#c(min(likes), 1.2*max(likes))
  par(mar=c(5,10,1,1), oma=c(1,0,0,0))
  if(missing(plot.colors)) plot.colors = mypalette(n.likes)
  barplot(horiz=TRUE, likes, beside=FALSE, col=plot.colors, xlab="Conditional log-likelihood components",  axisnames=FALSE,  axes=FALSE,  space=0,
    xlim=my.range)
  axis(side=1, at=pretty(seq(my.range[1],my.range[2]), n=10), labels=pretty(seq(my.range[1],my.range[2]), n=10), cex=.75 )
  axis(side=2, at=seq(0.5,(n.likes-0.5)), labels= names(likes), las=2)
  axis(side=2, at=seq(0.5,(n.likes-0.5))-0.2, labels = c(obs.dists,proc.dists), las = 2, tick = FALSE)
  text(x= likes, y=seq(0.5,(n.likes-0.5)), labels=round(likes,0), cex=0.8, pos=ifelse(likes>0, 4, 2))
  box()
  #title(paste0("Components of Obj. Function (", round(as.numeric(asap$like[1]),0), "), npar=", npar), cex=0.9 )
  title(sub=paste0("Model: ", mod$model_name, "     ", mod$date))
}

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
  sdrep = summary(model$sdrep)
  temp = sdrep[rownames(sdrep) %in% "log_catch_resid",]
  catch_stdresid <- matrix(temp[,1]/temp[,2], model$env$data$n_years_model, model$env$data$n_fleets)
  temp = sdrep[rownames(sdrep) %in% "log_index_resid",]
  index_stdresid <- matrix(temp[,1]/temp[,2], model$env$data$n_years_model, model$env$data$n_indices)
  temp = model$env$data$catch_paa - aperm(model$rep$pred_catch_paa,c(2,1,3))
  temp = sdrep[rownames(sdrep) %in% "log_index_resid",]
  index_stdresid <- matrix(temp[,1]/temp[,2], model$env$data$n_years_model, model$env$data$n_indices)
  #temp = model$env$data$catch_paa - aperm(model$rep$pred_catch_paa,c(2,1,3))

  out$RMSE$catch <- sqrt(apply(catch_stdresid^2,2,mean, na.rm = TRUE))
  out$RMSE_n$catch = apply(catch_stdresid^2,2,function(x) sum(!is.na(x)))
  out$RMSE$catch_total = sqrt(mean(catch_stdresid^2, na.rm = TRUE))
  out$RMSE_n$catch_total <- sum(!is.na(catch_stdresid^2))
  out$RMSE$index <- sqrt(apply(index_stdresid^2,2,mean, na.rm = TRUE))
  out$RMSE_n$index = apply(index_stdresid^2,2,function(x) sum(!is.na(x)))
  out$RMSE$index_total = sqrt(mean(index_stdresid^2, na.rm = TRUE))
  out$RMSE_n$index_total <- sum(!is.na(index_stdresid^2))
  return(out)
}
RMSE.table.fn <- function(model)
{
  origpar <- par(no.readonly = TRUE)
  RMSEs = get.RMSEs.fn(model)
	par(mfrow=c(1,1), oma=rep(2,4), mar=c(0,0,1,0))
	max.txt<-16
	rmses <- which(unlist(RMSEs$RMSE)>0)
	n.rmses<-length(rmses)
	plot(seq(1,15), seq(1,15), type='n', axes=F, xlab="", ylab="", xlim=c(1,max.txt+2+ 6+2+ 8+2), ylim=c(n.rmses+4, 1))
	text(rep(1, n.rmses), seq(3,n.rmses+2), labels=names(unlist(RMSEs$RMSE)[rmses]), pos=4)
	text(rep(max.txt+2, n.rmses), seq(3,n.rmses+2), labels=unlist(RMSEs$RMSE_n)[rmses], pos=4)
	text(rep(max.txt+2+ 6+2, n.rmses), seq(3,n.rmses+2), labels=signif(as.numeric(unlist(RMSEs$RMSE)[rmses]),3), pos=4)
	text(c(1, max.txt+2, max.txt+2+ 6+2), rep( 2, 3), labels=c("Component","# resids","RMSE"), font=2, pos=4)
	title(main="Root Mean Square Error computed from Standardized Residuals", outer=T, cex=0.85)
  par(origpar)
	#if (save.plots) savePlot(paste(od, "RMSE_Comp_Table.",plotf, sep=""), type=plotf)
}
#RMSE.table.fn(mod)
#plot.RMSE.table(x)

get.wham.results.fn = function(mod, out.dir, do.tex = FALSE, do.png = FALSE)
{
  years = mod$years #not used in cpp
  ages = mod$ages #not used in cpp, assumes last age is a plus group
  ny = length(years)
  na = length(ages)
  ni = mod$env$data$n_indices
  nf = mod$env$data$n_fleets

  res = list()
  res$ll = -mod$opt$obj
  res$np = length(mod$opt$par)
  res$aic = 2*(mod$opt$obj + res$np)
  #rho = mohns_rho(mod) #if no peels, then this will be NULL
  tcol <- col2rgb('black')
  black.poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  tcol <- col2rgb('red')
  red.poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  tcol <- col2rgb('blue')
  blue.poly <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')

  if(do.tex | do.png)
  {
    use_outer = TRUE
    x_line = y_line = 2.5
  }
  else
  {
    use_outer = FALSE
    x_line = -1
    y_line = 3
  }
  origpar <- par(no.readonly = TRUE)
  par(mar = c(5,5,1,1), oma = c(1,1,1,1), mfrow = c(3,1))
  if(do.tex | do.png)
  {
    fn = paste0(out.dir, '/SSB')
    if(do.tex) cairo_pdf(paste0(fn,'.pdf'), family = "Times", height = 10, width = 10)
    else png(filename = paste0(fn, '.png'), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar = c(0,0,0,0), oma = c(4,4,1,1), mfrow = c(1,1))
  }
  temp = summary(mod$sdrep)
  temp = temp[rownames(temp) == "log_SSB",]
  temp = exp(cbind(temp, temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2])))/1000
  max.y <- max(temp[,4])
  na.lim <- is.na(max.y)
  if(!na.lim){
    plot(years,temp[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    grid(col = gray(0.7))
    lines(years,temp[,1], lwd = 2)
    polygon(c(years,rev(years)), c(temp[,3],rev(temp[,4])), col = black.poly, border = "transparent")
    box(lwd = 2)
    mtext(side = 1, "Year", cex = 2, outer = TRUE, line = x_line)
    mtext(side = 2, "SSB (kmt)", cex = 2, outer = use_outer, line = y_line)
  } else {
    max.y <- max(temp[,1])
    plot(years,temp[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    grid(col = gray(0.7))
    lines(years,temp[,1], lwd = 2)
    # polygon(c(years,rev(years)), c(temp[,3],rev(temp[,4])), col = black.poly, border = "transparent")
    box(lwd = 2)
    mtext(side = 1, "Year", cex = 2, outer = TRUE, line = x_line)
    mtext(side = 2, "SSB (kmt)", cex = 2, outer = use_outer, line = y_line)    
  }
  if(do.tex | do.png) dev.off() else par(origpar)

  if(do.tex | do.png)
  {
    fn = paste0(out.dir, '/recruits')
    if(do.tex) cairo_pdf(paste0(fn,'.pdf'), family = "Times", height = 10, width = 10)
    else png(filename = paste0(fn, '.png'), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar = c(0,0,0,0), oma = c(4,4,1,1), mfrow = c(1,1))
  }
  temp = summary(mod$sdrep)
  ind = rownames(temp) == "log_NAA_rep"
  templo = exp(array(temp[ind,1] - qnorm(0.975)*temp[ind,2], dim = c(ny, na)))
  temphi = exp(array(temp[ind,1] + qnorm(0.975)*temp[ind,2], dim = c(ny, na)))
  temp = exp(array(temp[ind,1], dim = c(ny, 8)))
  max.y <- max(temphi[,1])
  na.lim <- is.na(max.y)
  if(!na.lim){
    plot(years,temp[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    grid(col = gray(0.7))
    lines(years,temp[,1], lwd = 2)
    polygon(c(years,rev(years)), c(templo[,1],rev(temphi[,1])), col = black.poly, border = "transparent")
    box(lwd = 2)
    if(use_outer) mtext(side = 1, "Year", cex = 2, outer = TRUE, line = x_line)
    mtext(side = 2, "Recruits (1000s)", cex = 2, outer = use_outer, line = y_line)
  } else {
    max.y <- max(temp[,1])
    plot(years,temp[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    grid(col = gray(0.7))
    lines(years,temp[,1], lwd = 2)
    # polygon(c(years,rev(years)), c(templo[,1],rev(temphi[,1])), col = black.poly, border = "transparent")
    box(lwd = 2)
    if(use_outer) mtext(side = 1, "Year", cex = 2, outer = TRUE, line = x_line)
    mtext(side = 2, "Recruits (1000s)", cex = 2, outer = use_outer, line = y_line)    
  }
  if(do.tex | do.png) dev.off() else par(origpar)

  if(do.tex | do.png)
  {
    fn = paste0(out.dir, '/Fbar')
    if(do.tex) cairo_pdf(paste0(fn,'.pdf'), family = "Times", height = 10, width = 10)
    else png(filename = paste0(fn, '.png'), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar = c(0,0,0,0), oma = c(4,4,1,1), mfrow = c(1,1))
  }
  temp = summary(mod$sdrep)
  temp = temp[rownames(temp) == "log_Fbar",]
  temp = exp(cbind(temp, temp[,1] + qnorm(0.975)*cbind(-temp[,2],temp[,2])))
  max.y <- max(temp[,4])
  na.lim <- is.na(max.y)
  if(!na.lim){
    plot(years,temp[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    grid(col = gray(0.7))
    lines(years,temp[,1], lwd = 2)
    polygon(c(years,rev(years)), c(temp[,3],rev(temp[,4])), col = black.poly, border = "transparent")
    box(lwd = 2)
    if(use_outer) mtext(side = 1, "Year", cex = 2, outer = TRUE,line = x_line)
    ar = c(min(mod$env$data$Fbar_ages), max(mod$env$data$Fbar_ages))
    mtext(side = 2, paste0("Average F (",ar[1],"-",ar[2],")"), cex = 2, outer = use_outer, line = y_line)
  } else {
    max.y <- max(temp[,1])
    plot(years,temp[,1], type = 'n', ylim = c(0,max.y), xlab = "", ylab = '', axes = FALSE)
    axis(1, lwd = 2, cex.axis = 1.5)
    axis(2, lwd = 2, cex.axis = 1.5)
    grid(col = gray(0.7))
    lines(years,temp[,1], lwd = 2)
    # polygon(c(years,rev(years)), c(temp[,3],rev(temp[,4])), col = black.poly, border = "transparent")
    box(lwd = 2)
    if(use_outer) mtext(side = 1, "Year", cex = 2, outer = TRUE,line = x_line)
    ar = c(min(mod$env$data$Fbar_ages), max(mod$env$data$Fbar_ages))
    mtext(side = 2, paste0("Average F (",ar[1],"-",ar[2],")"), cex = 2, outer = use_outer, line = y_line)
  }
  if(do.tex | do.png) dev.off() else par(origpar)
  # par(origpar)
}

plot.all.stdresids.fn = function(mod, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
  # load Ecov residuals
  xe <- NULL
  if(!all(mod$env$data$Ecov_model == 0)){
    ny = mod$env$data$n_years_Ecov
    ni = mod$env$data$n_Ecov
    years = mod$years
    temp = summary(mod$sdrep)
    ind = rownames(temp) == "Ecov_resid"
    templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
    temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
    temp = matrix(temp[ind,1], ny, ni)
    xe = data.frame(Label = integer(),
      Year = numeric(),
      Stdres = numeric(),
      lo = numeric(),
      hi = numeric())
    for(i in 1:ni)
    {
      ind = which(mod$env$data$Ecov_use_obs[,i] == 1)
      td = data.frame(Label = rep(i,length(ind)),
        Year = years[ind],
        Stdres = temp[ind,i],
        lo = templo[ind,i],
        hi = temphi[ind,i])
      xe <- rbind(xe, td)
    }
    xe$row = xe$Label
    xe$Label = factor(xe$Label)
    levels(xe$Label) = mod$env$data$Ecov_label
    xe$type = "Ecov"
  }

  # load Index residuals
  ny = mod$env$data$n_years_model
  ni = mod$env$data$n_indices
  temp = summary(mod$sdrep)
  ind = rownames(temp) == "log_index_resid"
  templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
  temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
  temp = matrix(temp[ind,1], ny, ni)
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
      Stdres = temp[ind,i],
      lo = templo[ind,i],
      hi = temphi[ind,i])
    xi <- rbind(xi, td)
  }
  xi$row = xi$Label
  xi$Label = factor(xi$Label)
  levels(xi$Label) = paste0("Index ",1:length(table(xi$Label)))
  xi$type = "Index"
  # if(!is.null(index.names)) levels(x$Index) = index.names

  # load catch data (fleet)
  ny = mod$env$data$n_years_model
  ni = mod$env$data$n_fleets
  if(missing(years)) years = mod$years
  temp = summary(mod$sdrep)
  ind = rownames(temp) == "log_catch_resid"
  templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
  temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
  temp = matrix(temp[ind,1], ny, ni)
  xc = data.frame(Label = integer(),
    Year = numeric(),
    Stdres = numeric(),
    lo = numeric(),
    hi = numeric())
  for(i in 1:ni)
  {
    td = data.frame(Label = rep(i,ny),
      Year = years,
      Stdres = temp[,i],
      lo = templo[,i],
      hi = temphi[,i])
    xc <- rbind(xc, td)
  }
  xc$row = xc$Label
  xc$Label = factor(xc$Label)
  levels(xc$Label) = paste0("Fleet ",1:length(table(xc$Label)))
  xc$type = "Catch"
  # if(!is.null(fleet.names)) levels(xc$Fleet) = fleet.names

  x <- rbind(xe, xi, xc)
  x$row = factor(x$row)

  ggp = ggplot2::ggplot(x, ggplot2::aes(x=Year, y = Stdres, color=type)) +
    ggplot2::geom_line(size=1.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=type), alpha=0.3, linetype = 0) +
    ggplot2::ylab("Standardized residual") +
    # expand_limits(y=0) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    # ggplot2::scale_color_manual(values=plot.colors) +
    # ggplot2::scale_fill_manual(values=plot.colors) +
    # ggplot2::facet_wrap(~Ecov, ncol=1)
    ggplot2::facet_grid(type ~ row)
  if(do.tex) cairo_pdf(file.path(od, paste0("Residuals_time.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("Residuals_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
  print(ggp)
  if(do.tex | do.png) dev.off()
  # return(ggp)
}

plot.ecov.stdresids.fn = function(mod, years, do.tex = FALSE, do.png = FALSE, res = 72, plot.colors, od)
{
  ny = mod$env$data$n_years_Ecov
  ni = mod$env$data$n_Ecov
  if(missing(years)) years = mod$years
  if(missing(plot.colors)) plot.colors = mypalette(ni)
  temp = summary(mod$sdrep)
  ind = rownames(temp) == "Ecov_resid"
  templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
  temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
  temp = matrix(temp[ind,1], ny, ni)
  x = data.frame(Ecov = integer(),
    Year = numeric(),
    Stdres = numeric(),
    lo = numeric(),
    hi = numeric())
  for(i in 1:ni)
  {
    ind = which(mod$env$data$Ecov_use_obs[,i] == 1)
    td = data.frame(Ecov = rep(i,length(ind)),
      Year = years[ind],
      Stdres = temp[ind,i],
      lo = templo[ind,i],
      hi = temphi[ind,i])
    x <- rbind(x, td)
  }
  x$Ecov = factor(x$Ecov)
  levels(x$Ecov) = mod$env$data$Ecov_label
  names(plot.colors) = levels(x$Ecov)
  ggp = ggplot2::ggplot(x, ggplot2::aes(x=Year, y = Stdres, color=Ecov)) +
    ggplot2::geom_line(size=1.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=Ecov), alpha=0.3, linetype = 0) +
    ggplot2::ylab("Standardized residual") +
    # expand_limits(y=0) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=plot.colors) +
    ggplot2::scale_fill_manual(values=plot.colors) +
    ggplot2::facet_wrap(~Ecov, ncol=1)
  if(do.tex) cairo_pdf(file.path(od, paste0("Residuals_ecov_time.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("Residuals_ecov_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
  print(ggp)
  if(do.tex | do.png) dev.off()
  # return(ggp)
}

plot.index.stdresids.fn = function(mod, years, index.names = NULL, do.tex = FALSE, do.png = FALSE, res = 72, plot.colors, od)
{
  ny = mod$env$data$n_years_model
  ni = mod$env$data$n_indices
  if(missing(years)) years = mod$years
  if(missing(plot.colors)) plot.colors = mypalette(ni)
  temp = summary(mod$sdrep)
  ind = rownames(temp) == "log_index_resid"
  templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
  temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
  temp = matrix(temp[ind,1], ny, ni)
  x = data.frame(Index = integer(),
    Year = numeric(),
    Stdres = numeric(),
    lo = numeric(),
    hi = numeric())
  for(i in 1:ni)
  {
    ind = which(mod$env$data$use_indices[,i] == 1)
    td = data.frame(Index = rep(i,length(ind)),
      Year = years[ind],
      Stdres = temp[ind,i],
      lo = templo[ind,i],
      hi = temphi[ind,i])
    x <- rbind(x, td)
  }
  x$Index = factor(x$Index)
  names(plot.colors)= levels(x$Index)
  if(!is.null(index.names)) levels(x$Index) = index.names
  ggp = ggplot2::ggplot(x, ggplot2::aes(x=Year, y = Stdres, color=Index)) +
    ggplot2::geom_line(size=1.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=Index), alpha=0.3, linetype = 0) +
    ggplot2::ylab("Standardized residual") +
    ggplot2::expand_limits(y=0) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values=plot.colors) +
    ggplot2::scale_fill_manual(values=plot.colors) +
    ggplot2::facet_wrap(~Index, ncol=1)
  if(do.tex) cairo_pdf(file.path(od, paste0("Residuals_log_index_time.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("Residuals_log_index_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
  print(ggp)
  if(do.tex | do.png) dev.off()
  # return(ggp)
}

plot.fleet.stdresids.fn = function(mod, years, fleet.names = NULL, do.tex = FALSE, do.png = FALSE, res = 72, plot.colors, od)
{
  ny = mod$env$data$n_years_model
  ni = mod$env$data$n_fleets
  if(missing(years)) years = mod$years
  if(missing(plot.colors)) plot.colors = mypalette(ni)
  temp = summary(mod$sdrep)
  ind = rownames(temp) == "log_catch_resid"
  templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
  temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
  temp = matrix(temp[ind,1], ny, ni)
  x = data.frame(Fleet = integer(),
    Year = numeric(),
    Stdres = numeric(),
    lo = numeric(),
    hi = numeric())
  for(i in 1:ni)
  {
    td = data.frame(Fleet = rep(i,ny),
      Year = years,
      Stdres = temp[,i],
      lo = templo[,i],
      hi = temphi[,i])
    x <- rbind(x, td)
  }
  x$Fleet = factor(x$Fleet)
  if(!is.null(fleet.names)) levels(x$Fleet) = fleet.names
  names(plot.colors)= levels(x$Fleet)
  ggp = ggplot2::ggplot(x, ggplot2::aes(x=Year, y = Stdres, color=Fleet)) +
    ggplot2::geom_line(size=1.1) +
    ggplot2::scale_color_manual(values=plot.colors) +
    ggplot2::scale_fill_manual(values=plot.colors) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lo, ymax=hi, fill=Fleet), alpha=0.3, linetype = 0) +
    ggplot2::ylab("Standardized residual") +
    ggplot2::expand_limits(y=0) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~Fleet, ncol=1)
  if(do.tex) cairo_pdf(file.path(od, paste0("Residuals_log_catch_time.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("Residuals_log_catch_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
  print(ggp)
  if(do.tex | do.png) dev.off()
  # return(ggp)
}

plot.catch.4.panel <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  dat = mod$env$data
  pred_catch = mod$rep$pred_catch
  catch = dat$agg_catch
  sigma = dat$agg_catch_sigma
  log_stdres = (log(catch)-log(pred_catch))/sigma
  if(!missing(use.i)) fleets <- use.i
  else fleets <- 1:dat$n_fleets
  if(missing(plot.colors)) plot.colors = mypalette(dat$n_fleets)
	for (i in fleets)
	{
		if(do.tex) cairo_pdf(file.path(od, paste0("Catch_4panel_fleet",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_4panel_fleet",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
		plot(years, catch[,i], type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Total Catch",
			ylim=c(0, 1.1*max(catch[,i])))
		lines(years, pred_catch[,i], col=plot.colors[i], lwd=2)
		log.ob.min <- log(catch[,i])-1.96*sigma[,i]
		log.ob.max <- log(catch[,i])+1.96*sigma[,i]
		plot(years, log(catch[,i]), type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Ln(Total Catch)",
			ylim=c(min(log.ob.min,log(pred_catch[,i])), 1.1*max(log.ob.max,log(pred_catch[,i]))))
		lines(years, log(pred_catch[,i]), col=plot.colors[i], lwd=2)
		arrows(years, log.ob.min, years, log.ob.max, length=0)
		title (paste0("Fleet ",i, " Catch"), outer=T, line=-1)
		plot(years, log_stdres[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Log-scale Std. Residual")
		abline(h=0)
		hist(log_stdres[,i], plot=T, xlab="Std. Residual", ylab="Probability Density", freq=F, main=NULL)
		if(do.tex | do.png) dev.off() else par(origpar)
	}
}

plot.index.4.panel <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  dat = mod$env$data
  pred_index = aperm(mod$rep$pred_IAA, c(2,1,3))
  pred_index = sapply(1:dat$n_indices, function(x)
  {
      if(dat$units_indices[x] == 2) apply(pred_index[x,,],1,sum)
      else apply(pred_index[x,,] * waa[dat$waa_pointer_indices[x],,],1,sum)
  })
  index = dat$agg_indices
  sigma = dat$agg_index_sigma
  log_stdres = (log(index)-log(pred_index))/sigma
  if(!missing(use.i)) indices <- use.i
  else indices <- 1:dat$n_indices
  if(missing(plot.colors)) plot.colors = mypalette(dat$n_indices)
	for (i in indices)
	{
		if(do.tex) cairo_pdf(file.path(od, paste0("Index_4panel_",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Index_4panel_",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
		plot(years, index[,i], type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Index",
			ylim=c(0, 1.1*max(index[,i])))
		lines(years, pred_index[,i], col=plot.colors[i], lwd=2)
		log.ob.min <- log(index[,i])-1.96*sigma[,i]
		log.ob.max <- log(index[,i])+1.96*sigma[,i]
		plot(years, log(index[,i]), type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Ln(Index)",
			ylim=c(min(log.ob.min,log(pred_index[,i])), 1.1*max(log.ob.max,log(pred_index[,i]))))
		lines(years, log(pred_index[,i]), col=plot.colors[i], lwd=2)
		arrows(years, log.ob.min, years, log.ob.max, length=0)
		title (paste0("Index ",i), outer=T, line=-1)
		plot(years, log_stdres[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Log-scale Std. Residual")
		abline(h=0)
		hist(log_stdres[,i], plot=T, xlab="Std. Residual", ylab="Probability Density", freq=F, main=NULL)
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

plot.NAA.4.panel <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
  years <- mod$years
  dat = mod$env$data
  pred_NAA = mod$rep$pred_NAA
  NAA = mod$rep$NAA
  sigma = exp(mod$parList$log_NAA_sigma[mod$env$data$NAA_sigma_pointers])
  sigma = matrix(sigma, length(years), dat$n_ages, byrow = TRUE)
  log_stdres = (log(NAA)-log(pred_NAA))/sigma
  if(!missing(use.i)) ages <- use.i
  else ages <- 1:dat$n_ages
  if(missing(plot.colors)) plot.colors = mypalette(dat$n_ages)
	for (i in ages)
	{
		if(do.tex) cairo_pdf(file.path(od, paste0("NAA_4panel_",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("NAA_4panel_",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(4,4,3,2), oma=c(1,1,1,1), mfrow=c(2,2))
		plot(years, NAA[,i], type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Abundance (1000s)",
			ylim=c(0, 1.1*max(NAA[,i])))
		lines(years, pred_NAA[,i], col=plot.colors[i], lwd=2)
		log.ob.min <- log(NAA[,i])-1.96*sigma[,i]
		log.ob.max <- log(NAA[,i])+1.96*sigma[,i]
		plot(years, log(NAA[,i]), type='p', col=plot.colors[i], pch=1, xlab="Year", ylab="Ln(Abundance)",
			ylim=c(min(log.ob.min,log(pred_NAA[,i])), 1.1*max(log.ob.max,log(pred_NAA[,i]))))
		lines(years, log(pred_NAA[,i]), col=plot.colors[i], lwd=2)
		arrows(years, log.ob.min, years, log.ob.max, length=0)
		title (paste0("Conditional Expected and Posterior Estimates of Age ",i, " Abundance "), outer=T, line=-1)
		plot(years, log_stdres[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Log-scale (Conditional) Std. Residual")
		abline(h=0)
		hist(log_stdres[,i], plot=T, xlab="(Conditional) Std. Residual", ylab="Probability Density", freq=F, main=NULL)
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

plot.NAA.res <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  ages = mod$ages
	n.ages <- length(ages)
  if(missing(plot.colors)) plot.colors = mypalette(n.ages)
  years = mod$years
	n.yrs <- length(years)
  pred_NAA = mod$rep$pred_NAA
  NAA = mod$rep$NAA
  sigma = exp(mod$parList$log_NAA_sigma[mod$env$data$NAA_sigma_pointers])
  sigma = matrix(sigma, length(years), n.ages, byrow = TRUE)
  log_stdres = (log(NAA)-log(pred_NAA))/sigma

  ymin <- min(apply(log_stdres,1, function(x) sum(x[which(x<0)])))
  ymax <- max(apply(log_stdres,1, function(x) sum(x[which(x>0)])))
  dat <- as.data.frame.table(log_stdres)
  x.at <- seq(1,n.yrs,5)
  x.lab <- years[x.at]
  if(do.tex) cairo_pdf(file.path(od, paste0("NAA_res_barplot_stacked.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("NAA_res_barplot_stacked.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
	par(mfrow=c(1,1), mar=c(5,5,1,1), oma = c(0,0,0,0))
	lattice::barchart(Freq ~ Var1, data = dat, groups = Var2, stack = TRUE, col = plot.colors, xlab = "Year", ylab = "Std. Abundance Residuals", box.ratio = 10, reference = TRUE,
	  scales = list(x=list(at = x.at, labels = x.lab), alternating = FALSE),
	  key = list(text = list(lab = as.character(ages)), rectangles = list(col = plot.colors), columns = n.ages, title = "Age"))
  if(do.tex | do.png) dev.off() else par(origpar)
  # par(origpar)
}

plot.catch.age.comp <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years = mod$years
  ages = 1:mod$env$data$n_ages
  ages.lab = mod$ages.lab
  if(!missing(use.i)) fleets <- use.i
  else fleets <- 1:mod$env$data$n_fleets
  if(missing(plot.colors)) plot.colors = mypalette(mod$env$data$n_fleets)

	for (i in fleets)
	{
    acomp.obs = mod$env$data$catch_paa[i,,]
    acomp.pred = aperm(mod$rep$pred_catch_paa, c(2,1,3))[i,,]
    if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_fleet",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_fleet",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
    my.title <- paste0("Fleet ", i)
    for (j in 1:length(years))
    {
      plot(1:length(ages), acomp.obs[j,], type='p', col=plot.colors[i], pch=1, xlab="", ylab="",
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
      lines(1:length(ages), acomp.pred[j,], col=plot.colors[i],  lwd=2)
      title(paste("Year = ", years[j], sep=""), outer = FALSE)

      # if 5x3 multipanel is full, save png and open new one
      if((j %% 15 == 0) & (do.tex | do.png) & (j < length(years))){
        dev.off()
        if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_fleet",i,"_",letters[j/15],".pdf")), family = "Times", height = 10, width = 10)
        if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_fleet",i,"_",letters[j/15],".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
        par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
      }
    }  #end loop on n_years
    if(length(years) %% 15 != 0) frame()
    if(do.tex | do.png) dev.off() else par(origpar)
	}  #end loop on n_fleets
  # par(origpar)
}

plot.index.age.comp <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, use.i, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  years = mod$years
  ages = 1:mod$env$data$n_ages
  ages.lab = mod$ages.lab
  if(!missing(use.i)) indices <- use.i
  else indices <- 1:mod$env$data$n_indices
  if(missing(plot.colors)) plot.colors = mypalette(mod$env$data$n_indices)

	for (i in indices)
	{
    acomp.obs = mod$env$data$index_paa[i,,]
    acomp.pred = aperm(mod$rep$pred_IAA, c(2,1,3))[i,,]
    if(mod$env$data$units_index_paa[i] == 1) acomp.pred = acomp.pred * waa[mod$env$data$waa_pointer_indices[i],,]
    acomp.pred = acomp.pred/apply(acomp.pred,1,sum)
    if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_index",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_index",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
    my.title <- paste0("Index ", i)
    for (j in 1:length(years))
    {
      plot(1:length(ages), acomp.obs[j,], type='p', col=plot.colors[i], pch=1, xlab="", ylab="",
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
      lines(1:length(ages), acomp.pred[j,], col=plot.colors[i],  lwd=2)
      title(paste("Year = ", years[j], sep=""), outer = FALSE)

      # if 5x3 multipanel is full, save png and open new one
      if((j %% 15 == 0) & (do.tex | do.png) & (j < length(years))){
        dev.off()
        if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_index",i,"_",letters[j/15],".pdf")), family = "Times", height = 10, width = 10)
        if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_index",i,"_",letters[j/15],".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
        par(mar=c(1,1,2,1), oma=c(4,4,2,1), mfcol=c(5,3))
      }
    }  #end loop on n_years
    if(length(years) %% 15 != 0) frame()
		if(do.tex | do.png) dev.off() else par(origpar)
	}  #end loop on n_indices
  # par(origpar)
}

multinomial.pearson.fn = function(mod, ind = 1)
{
  dat = mod$env$data
  rep = mod$rep
  x = dat$index_paa[ind,,] - rep$pred_index_paa[,ind,]
  temp = dat$index_Neff[,ind]
  temp = rep$pred_index_paa[,ind,]*(1-rep$pred_index_paa[,ind,])/temp
  x = x/sqrt(temp)
  x[which(dat$use_index_paa[,ind] == 0),] = NA
  return(x)
}
#mean(x < 0, na.rm = TRUE)

plot.catch.age.comp.resids <- function(mod, ages, ages.lab, scale.catch.bubble2 = 2, pos.resid.col = "#ffffffaa", neg.resid.col = "#8c8c8caa",
                                       do.tex = FALSE, do.png = FALSE, res = 72, use.i, od)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- dat$n_ages
  years = mod$years
	nyrs <- length(years)
	if(!missing(use.i)) fleets <- use.i
	else fleets <- 1:mod$env$data$n_fleets

	for (i in fleets)
	{
    acomp.obs = dat$catch_paa[i,,]
    acomp.pred = aperm(mod$rep$pred_catch_paa, c(2,1,3))[i,,]
    if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_resids_fleet",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_resids_fleet",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
		my.title <- "Age Comp Residuals for Catch by Fleet "
		#my.save <- "catch_resid_bubble_plots_"
    resids <- acomp.obs - acomp.pred  # NOTE obs-pred
    tylab <- "Residuals (Observed-Predicted)"
    z1 <- resids
    range.resids<-range(abs((as.vector(z1))), na.rm=T)
    scale.resid.bubble.catch <- 25

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
    title (paste0(my.title,i), outer=T, line=-1)
    if(do.tex | do.png) dev.off() else par(origpar)
	}   #end loop n_fleets
  # par(origpar)
}

plot.index.age.comp.resids <- function(mod, ages, ages.lab, scale.catch.bubble2 = 2, pos.resid.col = "#ffffffaa", neg.resid.col = "#8c8c8caa",
                                       do.tex = FALSE, do.png = FALSE, res = 72, use.i, od)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- dat$n_ages
  years = mod$years
	nyrs <- length(years)
	if(!missing(use.i)) indices <- use.i
	else indices <- 1:mod$env$data$n_indices

	for (i in indices)
	{
    acomp.obs = mod$env$data$index_paa[i,,]
    acomp.pred = aperm(mod$rep$pred_IAA, c(2,1,3))[i,,]
    if(mod$env$data$units_index_paa[i] == 1) acomp.pred = acomp.pred * waa[mod$env$data$waa_pointer_indices[i],,]
    acomp.pred = acomp.pred/apply(acomp.pred,1,sum)
    if(do.tex) cairo_pdf(file.path(od, paste0("Catch_age_comp_resids_index",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Catch_age_comp_resids_index",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
    my.title <- "Age Comp Residuals for Index "
    resids <- acomp.obs - acomp.pred  # NOTE obs-pred
    tylab <- "Residuals (Observed-Predicted)"
    z1 <- resids
    range.resids<-range(abs((as.vector(z1))), na.rm=T)
    scale.resid.bubble.catch <- 25

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
    title (paste0(my.title,i), outer=T, line=-1)
    if(do.tex | do.png) dev.off() else par(origpar)
	}   #end loop n_fleets
  # par(origpar)
}

plot.fleet.sel.blocks <- function(mod, ages, ages.lab, plot.colors, do.tex = FALSE, do.png = FALSE, res = 72, use.i, od)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
	cc<-0
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = ages
  sb_p = dat$selblock_pointer_fleets #selblock pointer by year and fleet
  if(missing(plot.colors)) plot.colors = mypalette(length(unique(sb_p)))
	years <- mod$years
	if(!missing(use.i)) fleets <- use.i
	else fleets <- 1:mod$env$data$n_fleets

	for (i in fleets)
	{
	  if(do.tex) cairo_pdf(file.path(od, paste0("Selectivity_fleet",i,".pdf")), family = "Times", height = 10, width = 10)
	  if(do.png) png(filename = file.path(od, paste0("Selectivity_fleet",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
	  blocks = unique(sb_p[,i])
		n.blocks <- length(blocks)
    sel = rbind(mod$rep$selblocks[blocks,])
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
		title(paste0("Fleet ",i))
		legend("topright", col=my.col, legend=paste0(minyr, " - ", maxyr), lwd=2, bg = "white")
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

plot.index.sel.blocks <- function(mod, ages, ages.lab, plot.colors, do.tex = FALSE, do.png = FALSE, res = 72, use.i, od)
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
	  if(do.tex) cairo_pdf(file.path(od, paste0("Selectivity_index",i,".pdf")), family = "Times", height = 10, width = 10)
	  if(do.png) png(filename = file.path(od, paste0("Selectivity_index",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
	  blocks = unique(sb_p[,i])
		n.blocks <- length(blocks)
    sel = rbind(mod$rep$selblocks[blocks,])
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
		title(paste0("Index ",i))
		legend("topright", col=my.col, legend=paste0(minyr, " - ", maxyr), lwd=2, bg = "white")
		if(do.tex | do.png) dev.off() else par(origpar)
	}
  # par(origpar)
}

plot.SSB.F.trend<-function(mod, alpha = 0.05)
{
  origpar <- par(no.readonly = TRUE)
  years <- mod$years
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  std = summary(mod$sdrep)
	par(mfrow=c(2,1), mar=c(1,1,1,1), oma = c(4,4,0,0))

	ssb.ind <- which(rownames(std) == "log_SSB")
	log.ssb <- std[ssb.ind,1]
  ssb = exp(log.ssb)/1000
	ssb.cv <- std[ssb.ind,2]
  log.ssb.ci <- log.ssb + cbind(qnorm(1-alpha/2)*ssb.cv, -qnorm(1-alpha/2)*ssb.cv)
  ssb.ci = exp(log.ssb.ci)/1000
  no.ssb.ci <- all(is.na(ssb.ci))
  if(!no.ssb.ci){ # have CI
  	plot(years, ssb, type='l', lwd=2, xlab="", ylab="", ylim=c(0,max(ssb.ci)), axes = FALSE)
  	axis(1, labels = FALSE)
  	axis(2)
  	box()
  	mtext(side = 2, "SSB (kmt)", outer = FALSE, line = 3)
  	grid(col = gray(0.7))
  	polygon(c(years,rev(years)), c(ssb.ci[,1],rev(ssb.ci[,2])), col = tcol, border = tcol, lwd = 1)
	} else { # no CI but plot SSB trend
    plot(years, ssb, type='l', lwd=2, xlab="", ylab="", ylim=c(0,max(ssb)), axes = FALSE)
    axis(1, labels = FALSE)
    axis(2)
    box()
    mtext(side = 2, "SSB (kmt)", outer = FALSE, line = 3)
    grid(col = gray(0.7))
    # polygon(c(years,rev(years)), c(ssb.ci[,1],rev(ssb.ci[,2])), col = tcol, border = tcol, lwd = 1)
  }
  # F trend
  n_ages = mod$env$data$n_ages
	faa.ind <- which(rownames(std) == "log_FAA_tot")
	log.faa <- matrix(std[faa.ind,1], length(years), n_ages)
	faa.cv <- matrix(std[faa.ind,2], length(years), n_ages)
	age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
  full.f.ind = cbind(1:length(years), age.full.f)
	#full.f.ind <- c(age.full.f[1], n_ages + cumsum(age.full.f[-1]))
  log.full.f <- log.faa[full.f.ind]
  full.f.cv <- faa.cv[full.f.ind]
  log.f.ci <- log.full.f + cbind(qnorm(1-alpha/2)*full.f.cv, -qnorm(1-alpha/2)*full.f.cv)
  full.f = exp(log.full.f)
  no.f.ci <- all(is.na(log.f.ci))
  if(!no.f.ci){ # have CI
  	plot(years, full.f, type='l', lwd=2, col='black', xlab="", ylab="", ylim=c(0,max(exp(log.f.ci))), axes = FALSE)
  	axis(1)
  	axis(2)
  	box()
  	mtext(side = 1, "Year", outer = FALSE, line = 3)
  	mtext(side = 2, "Fully-selected F", outer = FALSE, line = 3)
  	grid(col = gray(0.7))
  	polygon(c(years,rev(years)), exp(c(log.f.ci[,1],rev(log.f.ci[,2]))), col = tcol, border = tcol, lwd = 1)
  } else { # CI all NA
    plot(years, full.f, type='l', lwd=2, col='black', xlab="", ylab="", ylim=c(0,max(full.f)), axes = FALSE)
    axis(1)
    axis(2)
    box()
    mtext(side = 1, "Year", outer = FALSE, line = 3)
    mtext(side = 2, "Fully-selected F", outer = FALSE, line = 3)
    grid(col = gray(0.7))
    # polygon(c(years,rev(years)), exp(c(log.f.ci[,1],rev(log.f.ci[,2]))), col = tcol, border = tcol, lwd = 1)
  }
  par(origpar)
}  #end function

plot.SSB.AA <- function(mod, ages, ages.lab, plot.colors, prop=FALSE)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n.ages <- length(ages)
  if(missing(plot.colors)) plot.colors = mypalette(n.ages)
	years<-  mod$years
	n.yrs <- length(years)
  ssbfrac = dat$fracyr_SSB
	ssb.aa <- (mod$rep$NAA * exp(-ssbfrac * (mod$rep$FAA_tot + mod$rep$MAA)) * dat$waa[dat$waa_pointer_ssb,,] * dat$mature)/1000
	ssb.max <- max(apply(ssb.aa,1,sum))

	par(mfrow=c(1,1), mar=c(5,5,1,1), oma = c(0,0,0,0))
	if(!prop){ # plot SSB at age
  	barplot(t(ssb.aa), beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year', ylab =paste('SSB at age (', "kmt", ')', sep = ''),
  		ylim = c(0,1.15*ssb.max), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8)
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years[seq(5,n.yrs,5)])
    box()
	}
	if(prop){ # plot *proportion* SSB at age
  	barplot(t(ssb.aa/apply(ssb.aa,1,sum)), beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year', ylab ='Proportion SSB at age',
  		ylim = c(0,1.1), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8)
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years[seq(5,n.yrs,5)])
    box()
	}
  par(origpar)
}  #end funciton
#plot.SSB.AA(ssm)

#------------------------------------
plot.NAA <- function(mod, ages, ages.lab, plot.colors, scale = 1000, units = expression(10^6), prop=FALSE)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
  dat = mod$env$data
	## stacked barplot of NAA
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n.ages <- length(ages)
  if(missing(plot.colors)) plot.colors = mypalette(n.ages)
  years = mod$years
	n.yrs <-length(years)
	NAA <- mod$rep$NAA
	N.max=max(apply(NAA,1,sum))/scale

	par(mfrow=c(1,1), mar=c(5,5,1,1), oma = c(0,0,0,0))
	if(!prop){ # plot numbers at age
  	barplot(t(NAA)/scale, beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year',
  		ylab =as.expression(substitute(paste("January 1 numbers at age (", units, ")", sep = ''), list(units = units[[1]]))),
  		ylim = c(0,1.15*N.max), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8)
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years[seq(5,n.yrs,5)])
    box()
	}
	if(prop){ # plot *proportion* of numbers at age
  	barplot(t(NAA/apply(NAA,1,sum)), beside=F, cex.names=0.75, width=1, space=rep(0,n.yrs), xlab = 'Year',
  		ylab ='January 1 proportions at age', ylim = c(0,1.1), xlim=c(0.5,n.yrs+1-0.5), col=plot.colors)
  	legend('top', horiz=TRUE, legend=ages.lab, pch=15, col=plot.colors, cex=0.8 )
    axis(1, at = seq(5,n.yrs,5)-0.5, labels = years[seq(5,n.yrs,5)])
    box()
	}
  par(origpar)
} # end function
#plot.NAA(ssm)
#------------------------------------
#__annual estimates of Recruitment deviations (2 panel plot)____
plot.recruitment.devs <- function(mod, age.recruit = 1, units = expression(10^3), save.plots = FALSE)
{
  origpar <- par(no.readonly = TRUE)
	par(mfrow=c(2,1), mar=c(1,5,1,1), oma = c(4,1,1,0))
  dat = mod$env$data
	years <- mod$years
  nyrs = length(years)
	R <- mod$rep$NAA[,1]
	R.pred <- mod$rep$pred_NAA[,1]
	R.resids <- (log(R)-log(R.pred))/exp(mod$parList$log_NAA_sigma)[dat$NAA_sigma_pointers[1]]

	plot(years, R.pred, type='l', col='#114466', lty=1, lwd=2, xlab="",
	ylab= as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
		list(age.recruit = age.recruit, units = units[[1]]))), xlim = range(years), ylim=c(0, 1.1*max(R, R.pred)), axes = FALSE)
	lines(years, R, col="grey35", lwd=2, pch = 19, type = 'b')
	axis(1, labels = FALSE)
	axis(2)
	box()

	plot(years, R.resids, type='h', col='black', xlab="", ylab="Standardized (Log) Residuals", lwd=2,
	  xlim = range(years), ylim=1.1*range(R.resids))
	abline(h=0, lwd=1)
	mtext(side = 1, outer = TRUE, "Year", line = 2)
	par(origpar)
} # end function
#plot.recruitment.devs(ssm)

#scatter plot of SSB, R with 2-digit year as symbol (lag by 1 year)
plot.recr.ssb.yr <- function(mod, ssb.units = "kmt", recruits.units = expression(10^6), alpha = 0.05, 
  scale.ssb = 1000, scale.recruits = 1000, age.recruit = 1, plot.colors, loglog=FALSE)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1), mar = c(4,5,1,1), oma = c(1,1,1,1))

  std <- summary(mod$sdrep, "report")
  cov <- mod$sdrep$cov
  dat = mod$env$data
  years = mod$years
  nyrs <- length(years)
  nages <- dat$n_ages
  #nstates <- model$dimensions$n_states
  #nprojyrs <- model$dimensions$n_proj_years
	ssb.ind <- which(rownames(std) == "log_SSB")
	log.ssb <- std[ssb.ind,1]
	ssb.cv <- std[ssb.ind,2]
  log.ssb.ci <- log.ssb + cbind(qnorm(1-alpha/2)*ssb.cv, -qnorm(1-alpha/2)*ssb.cv)
  R.ind = which(rownames(std) == "log_NAA_rep")[1:nyrs]
  log.R = std[R.ind,1]
  R.cv = std[R.ind,2]
	ssb.R.cov <- cov[c(ssb.ind,R.ind),c(ssb.ind,R.ind)]

  if(missing(plot.colors)) plot.colors = mypalette(nyrs-age.recruit)
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
  par(origpar)
}  #end function

#------------------------------------
plot.SARC.R.SSB <- function(mod, scale.ssb=1, scale.recruits=1, age.recruit = 1, ssb.units = 'mt', recruits.units = expression(10^3))
{
  origpar <- par(no.readonly = TRUE)
  par(mar = c(5,5,1,5), oma = c(0,0,0,1), family='serif')
  years = mod$years
  nyrs <- length(years)
  std = summary(mod$sdrep)
	ssb.ind <- which(rownames(std) == "log_SSB")
	log.ssb <- std[ssb.ind,1]
	R <- mod$rep$NAA[,1]
	ssb.plot <- exp(log.ssb[1:(nyrs-age.recruit)])/scale.ssb
	recr.plot <- R[age.recruit + 1:(nyrs-age.recruit)]/scale.recruits
	yr.text <- substr(years,3,4)

	max.r <- max(recr.plot)
	max.ssb <- max(ssb.plot)
	scale.r <- max(ssb.plot)/max(recr.plot)
	ylimr <- c(0,1.1*max(recr.plot))
	barplot(recr.plot/scale.recruits, axisnames=FALSE, width=1, space=rep(0,nyrs-age.recruit), offset=rep(-0.5,nyrs-age.recruit), axes=FALSE, xpd=FALSE,
		xlab = '', ylab ='', ylim = ylimr, xlim=c(0.5,nyrs-age.recruit - 0.5), col="lightcyan2")
	xr <-pretty(c(0,recr.plot/scale.recruits))
	axis(2, at = xr, lab = xr )
	axis(side=1, las=2, at=seq(0.5,nyrs-age.recruit-0.5, by=2),
	labels=as.character(seq(years[1],years[nyrs-age.recruit], by=2)), cex=0.75, las=2)

	y.ssb <- (ssb.plot)*max.r/max.ssb
	lines(seq(0.5,nyrs-age.recruit-0.5, by=1), y.ssb, lwd=2, col = 'navyblue')
	x <- pretty(c(0,ssb.plot))
	axis(4, at = c(0,x*max.r/max.ssb), lab = c(0,x), col='navyblue', col.axis="navyblue")
	box()
	mtext(side = 1, 'Year', line = 3)
	mtext(side = 4, as.expression(substitute(paste("SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))), line = 3, col='navyblue')
	mtext(side = 2, as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
		list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), line = 3)
	par(origpar)
}  # end function
#plot.SARC.R.SSB(ssm, ssm.aux)

plot.fleet.F <- function(mod, plot.colors)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(1,1))
  years = mod$years
  nyrs <- length(years)
	n_fleets <- mod$env$data$n_fleets
  if(missing(plot.colors)) plot.colors = mypalette(n_fleets)
	for (i in 1:n_fleets)
	{
		if (i==1)
		{
		  plot(years, mod$rep$F[,i], xlab="Year", ylab="Full F", ylim=c(0,max(mod$rep$F)),	type='l', lty=1, lwd=2, col = plot.colors[i])
		  grid(col = gray(0.7))
		}
		if (i>1) lines(years, mod$rep$F[,i],lty=i, lwd=2, col=plot.colors[i])
	}

	leg.names <- paste0("Fleet ",i)
	legend('topleft', legend=leg.names, col=plot.colors,lwd=rep(2, n_fleets), lty=seq(1, n_fleets), horiz=TRUE, bty='n')
	par(origpar)
}   # end function
#plot.fleet.F(ssm,ssm.aux)

#------------------------------------
plot.cv <- function(mod)
{
  origpar <- par(no.readonly = TRUE)
	par(mfrow=c(1,1), mar=c(4,4,2,2))
  years = mod$years
  nyrs <- length(years)
	std <- summary(mod$sdrep)
	ssb.ind <- which(rownames(std) == "log_SSB")
	log.ssb <- std[ssb.ind,1]
	ssb.cv <- std[ssb.ind,2]
	R.ind <- which(rownames(std) == "log_NAA_rep")[1:nyrs]
	log.R <- std[R.ind,1]
  R.cv <- std[R.ind,2]
	faa.ind <- which(rownames(std) == "log_FAA_tot")
	log.faa <- matrix(std[faa.ind,1], nyrs, mod$env$data$n_ages)
	faa.cv <- matrix(std[faa.ind,2], nyrs, mod$env$data$n_ages)

	full.f.ind <- cbind(1:nyrs, apply(log.faa,1, function(x) max(which(x == max(x)))))
  log.full.f <- log.faa[full.f.ind]
  full.f.cv <- faa.cv[full.f.ind]

  any.na <- any(is.na(c(R.cv, faa.cv, full.f.cv)))

  if(!any.na){
  	plot(years, R.cv, type='l', lwd=2, col='black', xlab="Year", ylab="CV", ylim=c(0, 1.1*max(R.cv, ssb.cv, full.f.cv)))
  	lines(years, ssb.cv, lwd=2, col="blue")
  	lines(years, full.f.cv, lwd=2, col="green3")
  	legend('top', legend=c("Recruits", "SSB", "Full F"), col=c("black", "blue", "green3"), lty=rep(1,3), lwd=rep(2,3), horiz=T)
  }
  par(origpar)
}  # end function
#------------------------------------

#------------------------------------
plot.M <- function(mod, ages, ages.lab, alpha = 0.05, plot.colors)
{
  dat = mod$env$data
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- length(ages)
	meanMAA <- apply(mod$rep$MAA,2,mean)
  years = mod$years
  n_years = length(years)
	if(missing(plot.colors)) plot.colors = mypalette(n_years)
	n.M.by.age <- lapply(1:n_ages, function(x) table(mod$rep$MAA[,x]))
	plot(ages,meanMAA,lwd=2,xlab="Age",ylab="Natural Mortality Rate", ylim=c(0,1.1*(max(mod$rep$MAA))), type= 'n', xlim = c(min(ages)-0.5,max(ages)+1),
		axes=FALSE)
	grid(col = gray(0.7), lwd = 2)
	axis(1, at = ages, labels = ages.lab, lwd = 2)
	axis(2, lwd = 2)
	box(lwd = 2)
	sapply(1:n_ages, function(x)
	{
		y <- sort(unique(mod$rep$MAA[,x]))
		ind1 <- sapply(y, function(z) which(mod$rep$MAA[,x] == y)[1])
		segments(x-0.2, y, x+0.2, y, lwd = 2,col = plot.colors[1:length(y)])
		text(x+0.2, y, paste('n =', table(mod$rep$MAA[,x])), pos = 4)
	})
}

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
  if(missing(plot.colors)) plot.colors = mypalette(n_fleets)
	barplot(t(catch.obs), xlab="Year", ylab= paste0("Catch (", units, ")"), ylim=c(0,1.1*max(apply(catch.obs,1,sum))), col=plot.colors,space=0)
	axis(side=1, at = seq(2,nyrs,2)-0.5, labels = years[seq(2,nyrs,2)], cex=0.75)
	box(lwd = 2)
	if (n_fleets > 1)
  {
    legend('top', legend=paste0("Fleet ",1:n_fleets), horiz=TRUE, pch=15, col=plot.colors)

    # do proportions only if n_fleets > 1
		catch.prop <- catch.obs/apply(catch.obs,1,sum)
		barplot(t(catch.prop), xlab="Year", ylab="Proportion of Catch", ylim=c(0,1.1), col=plot.colors, space=0)
    axis(side=1, las=2, at = seq(2,nyrs,2)-0.5, labels = years[seq(2,nyrs,2)], cex=0.75, las=2)
    box(lwd = 2)
		legend('top', legend=paste0("Fleet ",1:n_fleets), horiz=TRUE, pch=15, col=plot.colors)
	}
	par(origpar)
}

# Bubble plots of catch age comps (set is.catch.flag to False to plot Discard age comps)
plot.catch.age.comp.bubbles <- function(mod, ages, ages.lab, bubble.col = "#8c8c8caa", i=1)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  years = mod$years
  nyrs = length(years)
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
  n_ages = length(ages)
	par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
  n_fleets = dat$n_fleets
	# for (i in 1:n_fleets)
	# {
		acomp.obs <- dat$catch_paa[i,,]
		catch.yrs <- which(dat$use_catch_paa[,i] == 1)
		my.title <- "Age Comps for Catch for Fleet "
		if (length(catch.yrs)>0)
		{
			scale.catch.obs <- 5
			z3 <- as.matrix(acomp.obs) * scale.catch.obs

			plot(ages, rev(ages),  xlim = range(ages), ylim = c(years[nyrs],(years[1]-2)), xlab = "Age", ylab = "", type = "n", axes=FALSE)
			axis(1, at= ages, lab = ages.lab)
			axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
			box()
			abline(h=years, col="lightgray")
			segments(x0=ages, y0=rep(years[1],n_ages), x1=ages, y1=rep(years[nyrs],n_ages), col = "lightgray", lty = 1)
			for (j in 1:nyrs) points(ages, rep(years[j], n_ages), cex=z3[j,], col="black", bg = bubble.col, pch = 21)

			bubble.legend1 <- c(0.05,0.2,0.4)
			bubble.legend2 <- bubble.legend1 * scale.catch.obs
			legend("topright", xpd=TRUE, legend=bubble.legend1, pch=rep(21, 3), pt.cex=bubble.legend2, horiz=T , col='black', pt.bg = bubble.col)
			title (paste0(my.title,i), outer=TRUE, line=-1)
		} # end catch.yrs test
	# }   #end loop n_fleets
  par(origpar)
}

#------------------------------------
plot.index.input <- function(mod, plot.colors)
{
  origpar <- par(no.readonly = TRUE)
  par(mfrow=c(2,1), mar = c(1,1,1,1), oma = c(4,4,2,0))
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
  if(missing(plot.colors)) plot.colors = mypalette(n_indices)
	plot(years,rescaled[,1],xlab="",ylab="", axes = FALSE, xlim = range(years), ylim=my.range,col=plot.colors[1],type='n')
	grid(col = gray(0.7))
	axis(1, labels = FALSE)
	axis(2)
	box()
	mtext(side = 2, "Rescaled Indices", outer = FALSE, line = 3)
	for (i in 1:n_indices) lines(years,rescaled[,i],col=plot.colors[i])
  legend("top", legend = paste0("Index " , 1:n_indices), col = plot.colors, lty = 1, horiz = TRUE, xpd = NA, inset = c(0,-0.1), bty = "n")

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

#------------------------------------
# Bubble plots of index age comps
plot.index.age.comp.bubbles <- function(mod, ages, ages.lab, bubble.col = "#8c8c8caa", i=1)
{
  origpar <- par(no.readonly = TRUE)
  par(mar=c(4,4,2,2), oma=c(1,1,1,1), mfrow=c(1,1))
  years = mod$years
  nyrs = length(years)
  dat = mod$env$data
  n_indices = dat$n_indices
  if(missing(ages)) ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	n_ages <- length(ages)

	# for (i in 1:n_indices)
	# {
		acomp.obs <- dat$index_paa[i,,]
		index.yrs <- which(dat$use_index_paa[,i] == 1)
		my.title <- "Age Comps for Index "
		if (length(index.yrs)>0)
		{
			scale.index.obs <- 5
			z3 <- as.matrix(acomp.obs) * scale.index.obs

			plot(ages, rev(ages),  xlim = range(ages), ylim = c(years[nyrs],(years[1]-2)),
			xlab = "Age", ylab = "", type = "n", axes=FALSE)
			axis(1, at= ages, lab=ages.lab)
			axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
			box()
			abline(h=years, col="lightgray")
			segments(x0=ages, y0=rep(years[1],n_ages), x1=ages, y1=rep(years[nyrs],n_ages), col = "lightgray", lty = 1)
			for (j in 1:nyrs) points(ages, rep(years[j], n_ages), cex=z3[j,], col="black", bg = bubble.col, pch = 21)

			bubble.legend1 <- c(0.05,0.2,0.4)
			bubble.legend2 <- bubble.legend1 * scale.index.obs
			legend("topright", xpd=TRUE, legend=bubble.legend1, pch=rep(21, 3), pt.cex=bubble.legend2, horiz=TRUE, col='black', pt.bg = bubble.col)
			title (paste0(my.title,i), outer=T, line=-1)
		} # end index.yrs test
	# }   #end loop n_fleets
	par(origpar)
}

#------------------------------------
plot.waa <- function(mod,type="ssb",plot.colors,ind=1)
{
  origpar <- par(no.readonly = TRUE)
  years = mod$years
  nyrs = length(years)
  dat = mod$env$data
  if(missing(plot.colors)) plot.colors = mypalette(dat$n_ages)
  point = switch(type,
    ssb = dat$waa_pointer_ssb,
    jan1 = dat$waa_pointer_jan1,
    fleets = dat$waa_pointer_fleets[ind],
    indices = dat$waa_pointer_indices[ind],
    totcatch = dat$waa_pointer_totcatch
  )
  labs = switch(type,
    ssb = "SSB",
    jan1 = "January 1 Biomass",
    fleets = paste0("Fleet ", ind),
    indices = paste0("Index ", ind),
    # fleets = paste0("Fleet ", 1:dat$n_fleets),
    # indices = paste0("Index ", 1:dat$n_indices),
    totcatch = "Total Catch"
  )
  waa = dat$waa[point,,]
  n = ifelse(length(dim(waa)) == 2, 1, dim(waa)[1])
	for(i in 1:n)
	{
		if(n>1) WAA.plot <- dat$waa[i,,]
    else WAA.plot = waa
		plot(years,years,xlab="Year",ylab="Weight",ylim=c(0,max(WAA.plot)),type='n')
		for (a in 1:dat$n_ages)
		{
			lines(years,WAA.plot[,a],col=plot.colors[a],lwd=2)
			lines(years,rep(mean(WAA.plot[,a]),length(years)),lty=2,col=plot.colors[a])
		}
		title(main = paste0("Annual Weight-at-Age for ", labs[i]))
	}  # end k-loop
	par(origpar)
}  # end function

#------------------------------------
plot.maturity <- function(mod, ages.lab, plot.colors)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  years = mod$years
  n_years = length(years)
  ages = 1:dat$n_ages
  if(missing(ages.lab)) ages.lab = mod$ages.lab
	meanmaturity <- apply(dat$mature,2,mean)
	if(missing(plot.colors)) plot.colors <- mypalette(n_years)

	plot(ages,meanmaturity,type='l',lwd=2,xlab="Age",ylab="Maturity",ylim=c(0,max(dat$mature)), axes = FALSE)
  axis(1, at = ages, labels = ages.lab, lwd = 2)
  axis(2, lwd = 2)
  box(lwd = 2)
	if (length(unique(dat$mature)) > length(ages))
	{
		for (i in 1:n_years) points(jitter(ages, factor=0.25), dat$mature[i,],col=yr.col[i])
		legend('topleft', horiz=FALSE, legend=c("First Year", "Last Year"), pch=c(1,1), col=c(plot.colors[1], plot.colors[n_years]))
	}
	title(main="Maturity", outer=FALSE)
	par(origpar)
}

#-------SSB/R -----------------------------
get_SPR = function(F, M, sel, mat, waassb, fracyrssb, at.age = FALSE)
{
  n_ages = length(sel)
  SPR = numeric()
  n = 1
  F = F * sel
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

plot.SPR.table <- function(mod, nyrs.ave = 5, plot=TRUE)
{
  origpar <- par(no.readonly = TRUE)
  spr.targ.values <- seq(0.2, 0.8, 0.05)
	n.spr <- length(spr.targ.values)
  dat = mod$env$data
	n_ages<- dat$n_ages
	years <- mod$years
	n_years <- length(years)
  avg.ind = (n_years-nyrs.ave+1):n_years
	#fec.age <- apply(dat$waa[dat$waa_pointer_ssb,,][avg.ind,],2,mean)
	mat.age <- apply(dat$mature[avg.ind,],2,mean)
	ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb,,][avg.ind,],2,mean)
	catch.waa <- apply(dat$waa[dat$waa_pointer_totcatch,,][avg.ind,],2,mean)
	M.age <- apply(mod$rep$MAA[avg.ind,],2,mean)
  sel = apply(mod$rep$FAA_tot[avg.ind,],2,mean) #average FAA, then do selectivity
	sel <- sel/max(sel)
	spawn.time <- mean(dat$fracyr_SSB[avg.ind])
  spr0 = get_SPR(F=0, M=M.age, sel=sel, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
	F.start <- 0.11  # starting guess for optimization routine to find F_SPR%

	f.spr.vals <- rep(NA, n.spr)
	ypr.spr.vals <- rep(NA, n.spr)
	conv.vals <- rep(NA, n.spr)

	for (i in 1:n.spr)
	{
		t.spr <- spr.targ.values[i]

		spr.f <- function(F.start)
		{
      spr = get_SPR(F=F.start, M=M.age, sel=sel, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
			abs(spr/spr0 - t.spr)
		}
		yyy <- nlminb(start=F.start, objective=spr.f, lower=0, upper=3)
		f.spr.vals[i] <- yyy$par
    ypr.spr.vals[i] = get_YPR(F = f.spr.vals[i], M=M.age, sel = sel, waacatch= catch.waa)
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
plot.annual.SPR.targets <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, plot.colors, od)
{
  origpar <- par(no.readonly = TRUE)
  spr.targ.values <- seq(0.2, 0.5, by=0.1)
	n.spr <- length(spr.targ.values)
  dat = mod$env$data
  n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)

  fec.age <- dat$waa[dat$waa_pointer_ssb,,]
	mat.age <- dat$mature
	wgt.age <- dat$waa[dat$waa_pointer_totcatch,,]
	M.age <- mod$rep$MAA
	sel.age <- mod$rep$FAA_tot/apply(mod$rep$FAA_tot,1,max)

	spawn.time <- dat$fracyr_SSB #vector length n_years

	spr0 <- numeric()
	F.start <- 0.11  # starting guess for optimization routine to find F_SPR%

	f.spr <- matrix(NA, n_years, n.spr)
	ypr.spr <- matrix(NA, n_years, n.spr)
	conv <- matrix(NA, n_years, n.spr)

	for (j in 1:n_years)
	{
    spr0[j] = get_SPR(F=0, M=M.age[j,], sel=sel.age[j,], mat=mat.age[j,], waassb=fec.age[j,], fracyrssb = spawn.time[j])
    #spr0.vals[j] <- s.per.recr(n_ages=n_ages, fec.age=fec.age[j,], mat.age=mat.age[j,], M.age= M.age[j,], F.mult=0, sel.age=sel.age[j,],
    #  spawn.time=spawn.time)
	  for (i in 1:n.spr)
	  {
			t.spr <- spr.targ.values[i]

			spr.f <- function(F.start)
			{
        spr = get_SPR(F=F.start, M=M.age[j,], sel=sel.age[j,], mat=mat.age[j,], waassb=fec.age[j,], fracyrssb = spawn.time[j])
        abs(spr/spr0[j] - t.spr)
				#abs(s.per.recr(n_ages=n_ages, fec.age=fec.age[j,], mat.age=mat.age[j,], M.age= M.age[j,], F.mult=F.start, sel.age=sel.age[j,],
				#	spawn.time=spawn.time)/spr0.vals[j] - t.spr)
			}
			yyy <- nlminb(start=F.start, objective=spr.f, lower=0, upper=3)
			f.spr[j,i] <- yyy$par
			ypr.spr[j,i] = get_YPR(F = f.spr[j,i], M=M.age[j,], sel = sel.age[j,], waacatch= wgt.age[j,])
      #ypr(n_ages, wgt.age=wgt.age[j,], M.age=M.age[j,],  F.mult=f.spr.vals[j,i], sel.age=sel.age[j,])
		}  # end j-loop over n_years
	}  #end i-loop over SPR values

  plot.colors = mypalette(n.spr)
  if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_annual_time.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("FSPR_annual_time.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
	par(mfrow=c(2,1), mar=c(4,5,2,1))
	plot(years, f.spr[,1], type='n', xlab="Years", ylab=expression(paste("Full ",italic(F)["%SPR"])), lwd=2, ylim=c(0,1.2*max(f.spr)))
	for (i in 1:n.spr) lines(years, f.spr[,i], lwd=2, col=plot.colors[i])
	legend('top', legend= sapply(spr.targ.values, function(x) as.expression(bquote(italic(F)[paste(.(x*100),"%")]))), col=plot.colors, horiz=TRUE, lwd=2, cex=0.9)
	title (main=expression(paste("Annual ", italic(F)["%SPR"], " Reference Points")), line=1)
	#if (save.plots) savePlot(paste(od, "Annual_FSPR.", plotf, sep=''), type=plotf)
	plot(years, ypr.spr[,1], type='n', xlab="Years", ylab=expression(paste("YPR(",italic(F)["%SPR"],")")), lwd=2, ylim=c(0,1.2*max(ypr.spr)))
	for (i in 1:n.spr) lines(years, ypr.spr[,i], lwd=2, col=plot.colors[i])
	#legend('top', legend=c("YPR20%", "YPR30%", "YPR40%", "YPR50%"), col=plot.colors, horiz=TRUE, lwd=2, cex=0.9)
	title (main=expression(paste("Annual YPR(",italic(F)["%SPR"], ") Reference Points")), line=1)
	if(do.tex | do.png) dev.off() else par(origpar)

	if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_freq_annual_F.pdf")), family = "Times", height = 10, width = 10)
	if(do.png) png(filename = file.path(od, paste0("FSPR_freq_annual_F.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
	par(mfrow=c(2,2))
  f.hist = lapply(1:length(spr.targ.values), function(x) hist(f.spr[,x], plot = FALSE))
  par(mfrow=c(2,2), mar=c(4,1,2,2), oma=c(2,3,2,0))
  for(i in 1:length(spr.targ.values))
  {
    plot(f.hist[[i]]$mids, f.hist[[i]]$counts, xlab=bquote("Full " ~ italic(F)[paste(.(spr.targ.values[i]*100),"%")]), ylab="", type='h', lwd=2, ylim=c(0, max(f.hist[[i]]$counts)), col=plot.colors[i])
    lines(f.hist[[i]]$mids, f.hist[[i]]$counts, lwd=2, col=plot.colors[i])
  }
  mtext(side=2, outer = TRUE, "Frequency", line = 1)
	title (main=expression(paste("Frequencies of Annual ", italic(F)["%SPR"], " Reference Points")), outer=TRUE, line=0, cex.main = 2)
	if(do.tex | do.png) dev.off() else par(origpar)

	if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_freq_annual_YPR.pdf")), family = "Times", height = 10, width = 10)
	if(do.png) png(filename = file.path(od, paste0("FSPR_freq_annual_YPR.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
	ypr.hist = lapply(1:length(spr.targ.values), function(x) hist(ypr.spr[,x], plot = FALSE))
  par(mfrow=c(2,2), mar=c(4,1,2,2), oma=c(2,3,2,0))
  for(i in 1:length(spr.targ.values))
  {
    plot(ypr.hist[[i]]$mids, ypr.hist[[i]]$counts, xlab=bquote(paste("YPR(", italic(F)[paste(.(spr.targ.values[i]*100),"%")],")")), ylab="", type='h', lwd=2, ylim=c(0, max(f.hist[[i]]$counts)), col=plot.colors[i])
    lines(ypr.hist[[i]]$mids, ypr.hist[[i]]$counts, lwd=2, col=plot.colors[i])
  }
  mtext(side=2, outer = TRUE, "Frequency", line = 1)
	title (main=expression(paste("Frequencies of Annual YPR(",italic(F)["%SPR"], ") Reference Points")), outer = TRUE, line=0, cex.main = 2)
	if(do.tex | do.png) dev.off() else par(origpar)

	# par(origpar)
} # end function

plot.SR.pred.line <- function(mod, ssb.units = "mt", SR.par.year, recruits.units = "thousands", scale.ssb = 1,
	scale.recruits = 1, age.recruit = 1, plot.colors)
{
  if(mod$env$data$recruit_model == 3 & mod$env$data$use_steepness != 1) #B-H stock recruit function with alpha/beta
  {
    std = summary(mod$sdrep)
    ssb.ind = which(rownames(std) == "log_SSB")
    years = mod$years
    nyrs = length(years)
    log.ssb <- std[ssb.ind,1]
    R <- mod$rep$NAA[,1]
    SR <- matrix(NA, (nyrs-age.recruit), 3)
    SR[,1] <- years[1:(nyrs-age.recruit)]/scale.ssb
    SR[,2] <- exp(log.ssb[1:(nyrs-age.recruit)])
    SR[,3] <- R[age.recruit +1:(nyrs-age.recruit)]/scale.recruits
    log_a <- mod$parList$mean_rec_pars[1]
    log_b <- mod$parList$mean_rec_pars[2]
    if(missing(SR.par.year)) SR.par.year = nyrs
    a.b.ind = which(rownames(std) == "mean_rec_pars")
    l.ab = std[a.b.ind,1]
    a.b.cov = mod$sdrep$cov.fixed[a.b.ind,a.b.ind]
    lR.fn = function(la, lb, S) la  + log(S) - log(1 + exp(lb)*S)
    dlR.dp = Deriv::Deriv(lR.fn, x = c("la","lb"))
    seq.ssb <- seq(0, max(SR[,2]), length.out=300)
    sd.pred.lR = sapply(seq.ssb, function(x)
    {
      d = dlR.dp(l.ab[1], l.ab[2], S= x)
      return(c(sqrt(t(d) %*% a.b.cov %*% d)))
    })
    pred.lR <- lR.fn(l.ab[1], l.ab[2], seq.ssb)
    #exp(log_a) *seq.ssb/(1 + exp(log_b)*seq.ssb)
    ci.pred.lR = pred.lR + qnorm(0.975)*cbind(-sd.pred.lR, sd.pred.lR) - log(scale.recruits)

    if(missing(plot.colors)) plot.colors = mypalette(nyrs-age.recruit)
    tcol <- col2rgb(plot.colors[1])
    tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')

    plot(SR[,2], SR[,3], type='p', col='black', pch=19,
      xlab=as.expression(substitute(paste("SSB (", ssb.units, ")", sep = ""), list(ssb.units = ssb.units[[1]]))),
      ylab= as.expression(substitute(paste("Age-", age.recruit, " Recruits (", units, ")", sep = ''),
        list(age.recruit = age.recruit[[1]], units = recruits.units[[1]]))), ylim=c(0, max(SR[,3])), xlim=c(0,1.1*max(SR[,2])))
    lines(seq.ssb, exp(pred.lR), col=plot.colors[1], lwd=2)
    polygon(c(seq.ssb,rev(seq.ssb)), exp(c(ci.pred.lR[,1],rev(ci.pred.lR[,2]))), col = tcol, border = "transparent")
  }
}

plot.FXSPR.annual <- function(mod, alpha = 0.05, status.years, max.x, max.y, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  percentSPR = mod$env$data$percentSPR
	n_ages = mod$env$data$n_ages
  years = mod$years
	n_years = length(years)
  if(missing(status.years)) status.years = n_years
  tcol <- col2rgb('black')
  tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
  std <- summary(mod$sdrep, "report")
	inds <- list(Y.t = which(rownames(std) == "log_Y_FXSPR"))
	inds$F.t <- which(rownames(std) == "log_FXSPR")
	inds$SSB.t <- which(rownames(std) == "log_SSB_FXSPR")
	inds$ssb <- which(rownames(std) == "log_SSB")
	inds$faa <- which(rownames(std) == "log_FAA_tot")
  log.faa <- matrix(std[inds$faa,1], n_years, n_ages)
  age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
  inds$full.f <- (age.full.f-1)*n_years + 1:n_years  + min(inds$faa) - 1 #cbind(1:n_years, age.full.f)
  na.sd <- sapply(inds, function(x) any(is.na(std[x,2])))
  ylabs <- c(
    bquote(paste('Yield(',italic(F)[paste(.(percentSPR), "%")], ')')),
    bquote(italic(F)[paste(.(percentSPR), "%")]),
    bquote(paste('SSB(', italic(F)[paste(.(percentSPR), "%")],')')))
	cov <- mod$sdrep$cov
  log.rel.ssb.rel.F.cov <- lapply(1:n_years, function(x)
  {
    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
    ind <- c(inds$ssb[x],inds$SSB.t[x],inds$full.f[x],inds$F.t[x])
    tcov <- cov[ind,ind]
    return(t(K) %*% tcov %*% K)
  })

  # FSPR absolute --------------------------------------------------
  if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_absolute.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("FSPR_absolute.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
  par(mfrow = c(3,1), mar = c(2,5,1,1), oma = c(4,2,1,1))
  for(i in 1:3)
  {
    t.ind <- inds[[i]]
    t.ylab <- ylabs[i]
    vals <- std[t.ind,1][1:n_years]
    cv <- std[t.ind,2][1:n_years]
    ci <-  vals + cbind(qnorm(1-alpha/2)*cv, -qnorm(1-alpha/2)*cv)
	  if(!na.sd[i]) plot(years, exp(vals), xlab = '', ylab = t.ylab, ylim = c(0,max(exp(ci),na.rm= TRUE)), type = 'l', cex.lab = 2)
    if(na.sd[i]) plot(years, exp(vals), xlab = '', ylab = t.ylab, ylim = c(0,max(exp(vals),na.rm= TRUE)), type = 'l', cex.lab = 2)
	  grid(col = gray(0.7))
    polyy = exp(c(ci[,1],rev(ci[,2])))
    polyx = c(years,rev(years))
    polyx = polyx[!is.na(polyy)]
    polyy = polyy[!is.na(polyy)]
    polygon(polyx, polyy, col = tcol, border = tcol, lwd = 1)
	}
  mtext(side = 1, outer = TRUE, "Year", cex = 2, line = 2)
  if(do.tex | do.png) dev.off() else par(origpar)

  # FSPR relative --------------------------------------------------
  if(do.tex) cairo_pdf(file.path(od, paste0("FSPR_relative.pdf")), family = "Times", height = 10, width = 10)
  if(do.png) png(filename = file.path(od, paste0("FSPR_relative.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
  par(mfrow=c(2,1))
  rel.ssb.vals <- std[inds$ssb,1][1:n_years] - std[inds$SSB.t,1][1:n_years]
  cv <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[1,1])))
  ci <-  rel.ssb.vals + cbind(-qnorm(1-alpha/2)*cv, qnorm(1-alpha/2)*cv)
  plot(years, exp(rel.ssb.vals), xlab = '', ylab = bquote(paste("SSB/", SSB[paste(.(percentSPR),"%")])), ylim = c(0,5), type = 'l')
  grid(col = gray(0.7))
  polyy = exp(c(ci[,1],rev(ci[,2])))
  polyx = c(years,rev(years))
  polyx = polyx[!is.na(polyy)]
  polyy = polyy[!is.na(polyy)]
  polygon(polyx, polyy, col = tcol, border = tcol, lwd = 1)
  abline(h=1, lty = 2)
  abline(h=0.5, lty = 2, col = 'red')

  rel.f.vals <- std[inds$full.f,1][1:n_years] - std[inds$F.t,1][1:n_years]
  cv <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[2,2])))
  ci <-  rel.f.vals + cbind(-qnorm(1-alpha/2)*cv, qnorm(1-alpha/2)*cv)
  if(!na.sd["full.f"]) plot(years, exp(rel.f.vals), xlab = '', ylab = bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])),
    ylim = c(0,max(exp(ci),1, na.rm = TRUE)), type = 'l')
  if(na.sd["full.f"]) plot(years, exp(rel.f.vals), xlab = '', ylab = bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])),
    ylim = c(0,max(exp(rel.f.vals),1, na.rm = TRUE)), type = 'l')
  grid(col = gray(0.7))
  polyy = exp(c(ci[,1],rev(ci[,2])))
  polyx = c(years,rev(years))
  polyx = polyx[!is.na(polyy)]
  polyy = polyy[!is.na(polyy)]
  polygon(polyx, polyy, col = tcol, border = tcol, lwd = 1)
  abline(h=1, lty = 2, col = 'red')
  mtext(side =1, "Year", outer = TRUE, line = 2, cex = 1.5)
  if(do.tex | do.png) dev.off() else par(origpar)

  # Kobe plot - only if sdreport was successful ---------------------------
  if(!mod$na_sdrep){
    log.rel.ssb.rel.F.ci.regs <- lapply(status.years, function(x){
      if(is.na(rel.f.vals[x])) return(maxtrix(NA,100,2))
      else return(exp(ellipse::ellipse(log.rel.ssb.rel.F.cov[[x]], centre = c(rel.ssb.vals[x],rel.f.vals[x]), level = 1-alpha)))
      })
    p.ssb.lo.f.lo <- sapply(status.years, function(x)
      mnormt::sadmvn(lower = c(-Inf,-Inf), upper = c(-log(2), 0), mean = c(rel.ssb.vals[x],rel.f.vals[x]), varcov = log.rel.ssb.rel.F.cov[[x]]))
    p.ssb.lo.f.hi <- sapply(status.years, function(x)
      mnormt::sadmvn(lower = c(-Inf,0), upper = c(-log(2), Inf), mean = c(rel.ssb.vals[x],rel.f.vals[x]), varcov = log.rel.ssb.rel.F.cov[[x]]))
    p.ssb.hi.f.lo <- sapply(status.years, function(x)
      mnormt::sadmvn(lower = c(-log(2),-Inf), upper = c(Inf, 0), mean = c(rel.ssb.vals[x],rel.f.vals[x]), varcov = log.rel.ssb.rel.F.cov[[x]]))
    p.ssb.hi.f.hi <- sapply(status.years, function(x)
      mnormt::sadmvn(lower = c(-log(2),0), upper = c(Inf, Inf), mean = c(rel.ssb.vals[x],rel.f.vals[x]), varcov = log.rel.ssb.rel.F.cov[[x]]))

    vals <- exp(cbind(rel.ssb.vals, rel.f.vals))
    if(missing(max.x)) max.x <- max(sapply(log.rel.ssb.rel.F.ci.regs, function(x) max(x[,1],na.rm = TRUE)),2)
    if(missing(max.y)) max.y <- max(sapply(log.rel.ssb.rel.F.ci.regs, function(x) max(x[,2],na.rm = TRUE)),2)

    if(do.tex) cairo_pdf(file.path(od, paste0("Kobe_status.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Kobe_status.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mfrow = c(1,1))
    plot(vals[status.years,1],vals[status.years,2], ylim = c(0,max.y), xlim = c(0,max.x), xlab = bquote(paste("SSB/", SSB[paste(.(percentSPR),"%")])),
      ylab = bquote(paste(italic(F),"/", italic(F)[paste(.(percentSPR),"%")])),type = 'n')
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
    legend("topleft", legend = paste0("Prob = ", round(p.ssb.lo.f.hi,2)), bty = "n")
    legend("topright", legend = paste0("Prob = ", round(p.ssb.hi.f.hi,2)), bty = "n")
    legend("bottomleft", legend = paste0("Prob = ", round(p.ssb.lo.f.lo,2)), bty = "n")
    legend("bottomright", legend = paste0("Prob = ", round(p.ssb.hi.f.lo,2)), bty = "n")
    text(vals[status.years,1],vals[status.years,2], substr(years[status.years],3,4))
    for(i in 1:length(status.years)) polygon(log.rel.ssb.rel.F.ci.regs[[i]][,1],log.rel.ssb.rel.F.ci.regs[[i]][,2])#, border = gray(0.7))
    if(do.tex | do.png) dev.off() else par(origpar)
    return(list(p.ssb.lo.f.lo = p.ssb.lo.f.lo, p.ssb.hi.f.lo = p.ssb.hi.f.lo, p.ssb.hi.f.hi = p.ssb.hi.f.hi, p.ssb.lo.f.hi = p.ssb.lo.f.hi))
  } else { return(NULL) }
}  # end function

plot.MSY.annual <- function(mod, alpha = 0.05, status.years, max.x, max.y, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)
  if(missing(status.years)) status.years = n_years
  std = summary(mod$sdrep, "report")
  cov <- mod$sdrep$cov
	if(dat$recruit_model == 3) #Beverton-Holt assumed in model fit
	{ # test to make sure steepness was estimated
    tcol <- col2rgb('black')
    tcol <- paste(rgb(tcol[1,],tcol[2,], tcol[3,], maxColorValue = 255), "55", sep = '')
		inds <- list(MSY = which(rownames(std) == "log_MSY"))
		inds$FMSY <- which(rownames(std) == "log_FMSY")
		inds$SSBMSY <- which(rownames(std) == "log_SSB_MSY")
		inds$RMSY <- which(rownames(std) == "log_R_MSY")
		inds$ssb <- which(rownames(std) == "log_SSB")
  	inds$faa <- which(rownames(std) == "log_FAA_tot")
	  log.faa <- matrix(std[inds$faa,1], n_years, n_ages)
	  age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
	  inds$full.f <- (age.full.f-1)*n_years + 1:n_years  + min(inds$faa) - 1 #cbind(1:n_years, age.full.f)
	  ylabs <- c(expression(MSY),expression(italic(F)[MSY]), expression(SSB[MSY]), expression(italic(R)[MSY]))
	  log.rel.ssb.rel.F.cov <- lapply(1:n_years, function(x)
	  {
	    K <- cbind(c(1,-1,0,0),c(0,0,1,-1))
	    ind <- c(inds$ssb[x],inds$SSBMSY[x],inds$full.f[x],inds$FMSY[x])
	    tcov <- cov[ind,ind]
	    return(t(K) %*% tcov %*% K)
	  })

    # 4-panel MSY plot
    if(do.tex) cairo_pdf(file.path(od, paste0("MSY_4panel_F_SSB_R.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("MSY_4panel_F_SSB_R.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mfrow=c(2,2))
    for(i in 1:4)
    {
      t.ind <- inds[[i]]
      t.ylab <- ylabs[i]
      vals <- std[t.ind,1][1:n_years]
  	  cv <- std[t.ind,2][1:n_years]
      ci <-  vals + cbind(qnorm(1-alpha/2)*cv, -qnorm(1-alpha/2)*cv)
      na.ci <- any(is.na(ci))
		  if(!na.ci) plot(years, exp(vals), xlab = 'Year', ylab = t.ylab, ylim = c(0,max(exp(ci))), type = 'l')
      if(na.ci) plot(years, exp(vals), xlab = 'Year', ylab = t.ylab, ylim = c(0,max(exp(vals))), type = 'l')
		  grid(col = gray(0.7))
		  polygon(c(years,rev(years)), exp(c(ci[,1],rev(ci[,2]))), col = tcol, border = tcol, lwd = 1)
		}
    if(do.tex | do.png) dev.off() else par(origpar)

    # 2-panel SSB_MSY and F_MSY
    if(do.tex) cairo_pdf(file.path(od, paste0("MSY_2panel_SSB_F.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("MSY_2panel_SSB_F.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    par(mfrow=c(2,1))
    rel.ssb.vals <- std[inds$ssb,1][1:n_years] - std[inds$SSBMSY,1][1:n_years]
    cv <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[1,1])))
    ci <-  rel.ssb.vals + cbind(-qnorm(1-alpha/2)*cv, qnorm(1-alpha/2)*cv)
		plot(years, exp(rel.ssb.vals), xlab = 'Year', ylab = expression(paste("SSB/", SSB[MSY])), ylim = c(0,5), type = 'l')
	  grid(col = gray(0.7))
	  polygon(c(years,rev(years)), exp(c(ci[,1],rev(ci[,2]))), col = tcol, border = "transparent", lwd = 1)
	  abline(h=1, lty = 2)
	  abline(h=0.5, lty = 2, col = 'red')

    rel.f.vals <- std[inds$full.f,1][1:n_years] - std[inds$FMSY,1][1:n_years]
    cv <- sapply(log.rel.ssb.rel.F.cov, function(x) return(sqrt(x[2,2])))
    ci <-  rel.f.vals + cbind(-qnorm(1-alpha/2)*cv, qnorm(1-alpha/2)*cv)
    na.ci <- any(is.na(ci))
		if(!na.ci) plot(years, exp(rel.f.vals), xlab = 'Year', ylab = expression(paste(italic(F),"/", italic(F)[MSY])),
      ylim = c(0,max(exp(ci),1)), type = 'l')
    if(na.ci) plot(years, exp(rel.f.vals), xlab = 'Year', ylab = expression(paste(italic(F),"/", italic(F)[MSY])),
      ylim = c(0,max(exp(rel.f.vals),1)), type = 'l')
	  grid(col = gray(0.7))
	  polygon(c(years,rev(years)), exp(c(ci[,1],rev(ci[,2]))), col = tcol, border = tcol, lwd = 1)
	  abline(h=1, lty = 2, col = 'red')
    if(do.tex | do.png) dev.off() else par(origpar)
	}
}  # end function

#------------------------------------
plot.yield.curves <- function(mod, nyrs.ave = 5, plot=TRUE, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
	n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)
  avg.ind = (n_years-nyrs.ave+1):n_years
	mat.age <- apply(dat$mature[avg.ind,],2,mean)
	ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb,,][avg.ind,],2,mean)
	catch.waa <- apply(dat$waa[dat$waa_pointer_totcatch,,][avg.ind,],2,mean)
	M.age <- apply(mod$rep$MAA[avg.ind,],2,mean)
  sel = apply(mod$rep$FAA_tot[avg.ind,],2,mean) #average FAA, then do selectivity
	sel <- sel/max(sel)
	spawn.time <- mean(dat$fracyr_SSB[avg.ind])
  spr0 = get_SPR(F=0, M=M.age, sel=sel, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)

	F.range <- seq(0,2.0, by=0.01)
	nF <- length(F.range)

  spr = sapply(F.range, function(x) get_SPR(F=x, M=M.age, sel=sel, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time))
  ypr = sapply(F.range, function(x) get_YPR(F=x, M=M.age, sel=sel, waacatch=catch.waa))
  pr = spr/spr0

  if(plot){
    if(do.tex) cairo_pdf(file.path(od, paste0("YPR_F_curve_plot.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("YPR_F_curve_plot.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
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
    if(do.tex) cairo_pdf(file.path(od, paste0("YPR_F_curve_table.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("YPR_F_curve_table.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")

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
plot.exp.spawn <- function(mod, nyrs.ave = 5)
{
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  n_ages = dat$n_ages
  years = mod$years
  n_years = length(years)
  avg.ind = (n_years-nyrs.ave+1):n_years
	mat.age <- apply(dat$mature[avg.ind,],2,mean)
	ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb,,][avg.ind,],2,mean)
	catch.waa <- apply(dat$waa[dat$waa_pointer_totcatch,,][avg.ind,],2,mean)
	M.age <- apply(mod$rep$MAA[avg.ind,],2,mean)
  sel = apply(mod$rep$FAA_tot[avg.ind,],2,mean) #average FAA, then do selectivity
	sel <- sel/max(sel)
	spawn.time <- mean(dat$fracyr_SSB[avg.ind])
  spr0 = get_SPR(F=0, M=M.age, sel=sel, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)

	F.range <- seq(0,2.0, by=0.01)
	nF <- length(F.range)
	# plot maturity vs selectivity based on "nyrs.ave" average

	par(mfrow=c(1,1), mar=c(4,4,2,4))
	plot( seq(1,n_ages), mat.age, type='b', pch = 16, col='black', lwd=2, ylim=c(0,1.1), xlab="Age", ylab="Selectivity or Maturity at age")
	lines( seq(1,n_ages), sel, col='orange2', lwd=1)
	points( seq(1,n_ages), sel, col='orange2', pch=16)
	legend('topleft', legend=c("Maturity", "Selectivity"), col=c("black", "orange2"), lwd=c(2,1), pch=c(NA, 16))
  spr = sapply(F.range, function(x) get_SPR(F=x, M=M.age, sel=sel, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time))
  expected.spawners = sapply(F.range, function(x) get_SPR(F=x, M=M.age, sel=sel, mat=mat.age, waassb=rep(1,n_ages), fracyrssb = spawn.time))
  ypr = sapply(F.range, function(x) get_YPR(F=x, M=M.age, sel=sel, waacatch=catch.waa))
  pr = spr/spr0

	par(mfrow=c(1,1), mar=c(4,4,2,4))

	plot(F.range, expected.spawners, type='n', xlab="Full F", ylab="", lwd=2, col="skyblue3", ylim=c(0,max(expected.spawners)), axes = FALSE)
  box()
  axis(1)
  axis(2, col='skyblue3', col.axis="skyblue3")
  mtext(side = 2, col = "skyblue3", "Expected Spawnings", line = 2)
	abline(v=seq(0.1,2.0, by=0.1), col="grey85")
	abline(h=c(1,2,3), col="skyblue1")
	lines(F.range, expected.spawners, lwd=2, col="skyblue3")
	points(F.range, expected.spawners, pch=19, col="skyblue3")
	scale.pr <- max(pr)/max(expected.spawners)

	lines(F.range, pr/scale.pr, col="red", lwd=2)
	points(F.range, pr/scale.pr, col="red", pch=19)
	axis(side=4, at=seq(0,1,by=0.1)/scale.pr, lab=seq(0,1,by=0.1), las=2, col='red', col.axis="red")
	mtext(side=4, "% SPR", line=3, col="red")
	title (paste("Expected Spawnings and SPR Reference Points (Years Avg = ", nyrs.ave,")",	sep=""), cex=0.9, outer=T, line=-1)

	par(mfrow=c(1,1), mar=c(2,2,2,2))
	plot(seq(1,42), seq(-3,38), type='n', axes=F, bty='n',xlab="",ylab="")
	text(x=0,y=36, labels="F", font=2, pos=4)
	text(x=4, y=36, labels="E[Sp]" , font=2, pos=4,cex=0.9)
	text(x=9, y=36, labels="SPR", font=2, pos=4)
	text(x=15,y=36, labels="F", font=2, pos=4)
	text(x=19, y=36, labels="E[Sp]" , font=2, pos=4,cex=0.9)
	text(x=24, y=36, labels="SPR", font=2, pos=4)
	text(x=30,y=36, labels="F", font=2, pos=4)
	text(x=34, y=36, labels="E[Sp]" , font=2, pos=4,cex=0.9)
	text(x=39, y=36, labels="SPR", font=2, pos=4)

  text(x=0, y=seq(35,1, by=-1), labels=F.range[1:35], cex=0.82, pos=4, font=1)
  text(x=4, y=seq(35,1, by=-1), labels=round(expected.spawners[1:35],4), cex=0.82, pos=4, font=1)
  text(x=9, y=seq(35,1, by=-1), labels=round(pr[1:35],4), cex=0.82, pos=4, font=1)

  text(x=15, y=seq(35,1, by=-1), labels=F.range[36:70], cex=0.82, pos=4, font=1)
  text(x=19, y=seq(35,1, by=-1), labels=round(expected.spawners[36:70],4), cex=0.82, pos=4, font=1)
  text(x=24, y=seq(35,1, by=-1), labels=round(pr[36:70],4), cex=0.82, pos=4, font=1)

  text(x=30, y=seq(35,1, by=-1), labels=F.range[71:105], cex=0.82, pos=4, font=1)
  text(x=34, y=seq(35,1, by=-1), labels=round(expected.spawners[71:105],4), cex=0.82, pos=4, font=1)
  text(x=39, y=seq(35,1, by=-1), labels=round(pr[71:105],4), cex=0.82, pos=4, font=1)

	title (paste("Expected Spawnings & SPR Reference Points (Years Avg = ", nyrs.ave,")", sep=""), cex=0.9, outer=T, line=-1)


	exp.spawn.table<- as.data.frame(cbind(F.range, expected.spawners, pr))
	colnames(exp.spawn.table) <- c("Full.F", "Exp.Spawn", "SPR")
	#write.csv(exp.spawn.table, file=paste(od,"Exp.Spawn.Table.csv", sep=""), row.names=F)
  par(origpar)
} # end function

plot.retro <- function(mod,y.lab,y.range1,y.range2, alpha = 0.05, what = "SSB", do.tex = FALSE, do.png = FALSE, res = 72, od)
{
  origpar <- par(no.readonly = TRUE)
  years = mod$years
	n_years <- length(years)
  npeels = length(mod$peels)
  if(npeels)
  {
    # standard retro plot
    if(do.tex) cairo_pdf(file.path(od, paste0(what,"_retro.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0(what,"_retro.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    plot.colors = mypalette(npeels+1)
    tcol = col2rgb(plot.colors)
    tcol = rgb(tcol[1,],tcol[2,],tcol[3,], maxColorValue = 255, alpha = 200)
    res = list(mod$rep[[what]])
    res[2:(npeels+1)] = lapply(mod$peels, function(x) x$rep[[what]])
    if(what == "NAA")
    {
      par(mfcol = c(ceiling(NCOL(res[[1]])/2),2))
      for(i in 1:NCOL(res[[1]]))
      {
        y.range1 <- range(sapply(res, function(x) range(x[,i])))
        plot(years,res[[1]][,i],lwd=1,col=plot.colors[1],type='l',xlab="Year",ylab=paste0("Numbers at age ", mod$ages.lab[i]),ylim=y.range1)
        grid(col = gray(0.7), lty = 2)
        for (j in 1:npeels)
        {
          lines(years[1:(n_years-j)],res[[j+1]][,i], col = tcol[j+1])
          points(years[n_years-j],res[[j+1]][n_years-j,i],pch=16,col=plot.colors[j+1])
        }
      }
    }
    else
    {
      if(missing(y.range1)) y.range1 <- range(sapply(res, function(x) range(x)))
      par(mfrow = c(1,1))
      plot(years,res[[1]],lwd=1,col=plot.colors[1],type='l',xlab="Year",ylab=what,ylim=y.range1)
      grid(col = gray(0.7), lty = 2)
      for (i in 1:npeels)
      {
        lines(years[1:(n_years-i)],res[[i+1]], col = tcol[i+1])
        points(years[n_years-i],res[[i+1]][n_years-i],pch=16,col=plot.colors[i+1])
      }
    }
    if(do.tex | do.png) dev.off() else par(origpar)

    # relative retro plot
    if(do.tex) cairo_pdf(file.path(od, paste0(what,"_retro_relative.pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0(what,"_retro_relative.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
    if(missing(y.lab)) y.lab = bquote(paste("Mohn's ", rho, "(",.(what),")"))
    if(what == "NAA") rel.res = lapply(1:length(res), function(x) res[[x]]/res[[1]][1:(n_years - x + 1),] - 1)
    else rel.res = lapply(1:length(res), function(x) res[[x]]/res[[1]][1:(n_years - x + 1)] - 1)
    rho.vals = mohns_rho(mod)

    if(what == "NAA")
    {
      par(mfcol = c(ceiling(NCOL(rel.res[[1]])/2),2))
      for(i in 1:NCOL(res[[1]]))
      {
        y.range2 <- c(-1,max(sapply(rel.res, function(x) range(x[,i]))))
        plot(years,rel.res[[1]][,i],lwd=1,col=plot.colors[1],type='l',xlab="Year",ylab=bquote(paste("Mohn's ", rho, "(Numbers at age ", .(mod$ages.lab[i]), ")")),ylim = y.range2)
        grid(col = gray(0.7), lty = 2)
        for (j in 1:npeels)
        {
          lines(years[1:(n_years-j)],rel.res[[j+1]][,i], col = tcol[j+1])
          points(years[n_years-j],rel.res[[j+1]][n_years-j,i],pch=16,col=plot.colors[j+1])
        }
        rho.nm = ifelse(i==1, "R",paste0("N",mod$ages.lab[i]))
        rho.plot <- round(rho.vals[rho.nm],3)
        legend("bottomleft", legend = bquote(rho == .(rho.plot)), bty = "n")
      }
    }
    else
    {
      if(missing(y.range2)) y.range2 <- c(-1,max(sapply(rel.res, function(x) range(x))))
      par(mfrow = c(1,1))
      plot(years,rel.res[[1]],lwd=1,col="black",type='l',xlab="Year",ylab=y.lab,ylim=y.range2)
      grid(col = gray(0.7), lty = 2)
      for (i in 1:npeels)
      {
        lines(years[1:(n_years-i)],rel.res[[i+1]], col = tcol[i+1])
        points(years[n_years-i],rel.res[[i+1]][n_years-i],pch=16,col=plot.colors[i+1])
      }
      rho.plot <- round(rho.vals[what],3)
      legend("bottomleft", legend = bquote(rho == .(rho.plot)), bty = "n")
    }
    if(do.tex | do.png) dev.off() else par(origpar)
    # par(origpar)
    return(rho.vals)
  }
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-----Catch curve and cohort correspondence plots-------------------------------
#-------------------------------------------------------------------------------
# function to compute catch at age (matrix)
# from proportions at age (matrix), weight at age (matrix) and total weight (vector)
wtprop2caa <- function(fleet)#totwt,waa,props)
{

  return(paa * totwt/apply(paa*waa,1,sum))
	caa <- mod$rep$pred_caa[,fleet,] * mod$env$data
	for (i in 1:length(totwt))
	{
		if (sum(props[i,]) == 0) caa[i,] <- 0
		if (sum(props[i,]) > 0) caa[i,] <- props[i,] * (totwt[i] / sum(props[i,] * waa[i,]))
	}
	return(caa)
}

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
					zz <- predict(my.fit,xrng,interval="confidence")
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

plot_catch_at_age_consistency <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
	# create plots of ages vs each other (correctly lagged) on log scale for catch by fleet
	cat.corr <- list()
  dat = mod$env$data
  rep = mod$rep
	for (i in 1:dat$n_fleets)
	{
		if (dat$n_fleets == 1) title1 = "Catch"
		if (dat$n_fleets >= 2) title1 = paste("Catch for Fleet ",i, sep="")

		# get catch at age
    catchob = dat$catch_paa[i,,] * dat$agg_catch[,i]/apply(dat$catch_paa[i,,] * dat$waa[dat$waa_pointer_fleets[i],,],1,sum)
    catchpr = rep$pred_catch_paa[,i,] * rep$pred_catch[,i]/apply(rep$pred_catch_paa[,i,] * dat$waa[dat$waa_pointer_fleets[i],,],1,sum)

		# replace zeros with NA and take logs
		cob <- rep0log(catchob)
		cpr <- rep0log(catchpr)

		# make cohorts
		cob.coh <- makecohorts(cob)
		cpr.coh <- makecohorts(cpr)

		# make the plots
		if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_obs_fleet",i,".pdf")), family = "Times", height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_obs_fleet",i,".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
		cob.cor <- plotcoh(cob.coh,mytitle=paste(title1," Observed", sep=""),mod=mod)
		if(do.tex | do.png) dev.off()

		if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_pred_fleet",i,".pdf")), family = "Times", height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_pred_fleet",i,".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
		cpr.cor <- plotcoh(cpr.coh,mytitle=paste(title1," Predicted", sep=""),mod=mod)
		if(do.tex | do.png) dev.off()

		cat.corr[[i]] <- list(cob.cor,cpr.cor)
	}
	return(cat.corr)
}

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
			agg.ob <- dat$agg_indices[dat$use_index_paa[,i]==1,i]
			agg.pr <- rep$pred_indices[dat$use_index_paa[,i]==1,i]

			# get proportions for correct years and ages only
			props.ob <- dat$index_paa[i,dat$use_index_paa[,i]==1,]
			props.pr <- rep$pred_index_paa[dat$use_index_paa[,i]==1,i,]

      # figure out units for aggregate and proportions
			agg.units <- dat$units_indices[i]
			prp.units <- dat$units_index_paa[i]

      waa <- dat$waa[dat$waa_pointer_indices[i],dat$use_index_paa[,i]==1,]
			# get weight (matrix if necessary)
			if (agg.units==1 | prp.units==1)
			{  # either in weight
			}

			# create index.obs and pred based on which of the four possible combinations of units is used for this index
			if (agg.units==1 && prp.units==1)
			{  # both in weight
				index.ob <- agg.ob * props.ob / waa
				index.pr <- agg.pr * props.pr / waa
			}
			if (agg.units==1 && prp.units==2)
			{  # agg in weight, props in numbers
        index.ob = props.obs * agg.obs/apply(props.obs * waa,1,sum)
        index.pr = props.pr * agg.pr/apply(props.pr * waa,1,sum)
				#index.ob <- wtprop2caa(agg.ob,waa,props.ob)  # use catch function
				#index.pr <- wtprop2caa(agg.pr,waa,props.pr)
			}
			if (agg.units==2 && prp.units==1)
			{  # agg in numbers, props in weight
				# need to search for correct agg total in weight to result in observed agg total in number
				# for now just use simple approximation that agg.wt = sum(waa*prop) *ctot and then solve using both in weight approach
				agg.wt.ob <- apply((waa * props.ob),1,sum) * agg.ob
				agg.wt.pr <- apply((waa * props.pr),1,sum) * agg.pr
				index.ob <- agg.wt.ob * props.ob / waa
				index.pr <- agg.wt.pr * props.pr / waa
			}
			if (agg.units==2 && prp.units==2)
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
#-------------------------------------------------------------------------------

plot_index_at_age_consistency <- function(mod, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
  dat = mod$env$data
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
			title1 <- paste("Index ",ind, sep="")

			# replace zeros with NA and take logs
			iob <- rep0log(index.mats$ob[[ind]])
			ipr <- rep0log(index.mats$pr[[ind]])

			# make cohorts
			iob.coh <- makecohorts(iob)
			ipr.coh <- makecohorts(ipr)


			# create labels for plot (if necessary)
			mylabels <- paste0("age-",mod$ages.lab, sep="")

			# make the plots
			if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_obs_index",ind,".pdf")), family = "Times", height = 10, width = 10)
			if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_obs_index",ind,".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
			iob.cor <- plotcoh(iob.coh,mytitle=paste(title1," Observed", sep=""),mylabels,mod=mod)
			if(do.tex | do.png) dev.off()

			if(do.tex) cairo_pdf(file.path(od, paste0("catch_at_age_consistency_pred_index",ind,".pdf")), family = "Times", height = 10, width = 10)
			if(do.png) png(filename = file.path(od, paste0("catch_at_age_consistency_pred_index",ind,".png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
			ipr.cor <- plotcoh(ipr.coh,mytitle=paste(title1," Predicted", sep=""),mylabels,mod=mod)
			if(do.tex | do.png) dev.off()

			index.corr[[ind]] <- list(iob.cor,ipr.cor)
		}
	}
	return(index.corr)
}
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

plot_catch_curves_for_catch <- function(mod, first.age=-999, do.tex = FALSE, do.png = FALSE, res = 72, od)
{
	# create catch curve plots for catch by fleet
  origpar <- par(no.readonly = TRUE)
	lastyr <- max(mod$years)
  dat = mod$env$data
  rep = mod$rep
	cohort <- (min(mod$years) - dat$n_ages-1):(lastyr+dat$n_ages-1)
	ages <- 1:dat$n_ages
  my.col = rep(mypalette(5),50)
	#my.col <- rep(c("blue","red","green","orange","gray50"),50)
	for (i in 1:dat$n_fleets)
	{
		if (dat$n_fleets == 1) title1 = "Catch"
		if (dat$n_fleets >= 2) title1 = paste0("Catch for Fleet ",i)

		# get catch at age
    catchob = dat$catch_paa[i,,] * dat$agg_catch[,i]/apply(dat$catch_paa[i,,] * dat$waa[dat$waa_pointer_fleets[i],,],1,sum)
    catchpr = rep$pred_catch_paa[,i,] * rep$pred_catch[,i]/apply(rep$pred_catch_paa[,i,] * dat$waa[dat$waa_pointer_fleets[i],,],1,sum)

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
		if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_fleet",i,"_obs.pdf")), family = "Times", height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_curves_fleet",i,"_obs.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
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

		if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_fleet",i,"_pred.pdf")), family = "Times", height = 10, width = 10)
		if(do.png) png(filename = file.path(od, paste0("catch_curves_fleet",i,"_pred.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
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

plot_catch_curves_for_index <- function(mod, first.age=-999, do.tex = FALSE, do.png = FALSE, res = 72, od)
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
			title1 <- paste0("Index ",ind)
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
			if(!(all(is.na(iob.coh)) & all(is.na(ipr.coh))))
			{
			  if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_index",ind,"_obs.pdf")), family = "Times", height = 10, width = 10)
			  if(do.png) png(filename = file.path(od, paste0("catch_curves_index",ind,"_obs.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
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

				if(do.tex) cairo_pdf(file.path(od, paste0("catch_curves_index",ind,"_pred.pdf")), family = "Times", height = 10, width = 10)
				if(do.png) png(filename = file.path(od, paste0("catch_curves_index",ind,"_pred.png")), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")
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
#plot_catch_curves_for_index(ssm)

#-------------------------------------------------------------------------------
plot.ecov <- function(mod, use.i, plot.pad = FALSE, do.tex = FALSE, do.png = FALSE, od){
  origpar <- par(no.readonly = TRUE)
  dat = mod$env$data
  years <- seq(from=dat$year1_Ecov, by=1, length.out=dat$n_years_Ecov)

  ecov.pred = mod$rep$Ecov_x
  ecov.obs = dat$Ecov_obs
  ecov.obs.sig = dat$Ecov_obs_sigma

  # default: don't plot the padded entries that weren't used in ecov likelihood
  if(!plot.pad){
    ecov.use = dat$Ecov_use_obs
    ecov.obs[ecov.use == 0] <- NA
    ecov.obs.sig[ecov.use == 0] <- NA
  }
  ecov.res = (ecov.obs - ecov.pred) / ecov.obs.sig # standard residual (obs - pred)

  if(!missing(use.i)) ecovs <- use.i
  else ecovs <- 1:dat$n_Ecov
  plot.colors = mypalette(dat$n_Ecov)
  for (i in ecovs)
  {
    if(do.tex) cairo_pdf(file.path(od, paste0("Ecov_",i,".pdf")), family = "Times", height = 10, width = 10)
    if(do.png) png(filename = file.path(od, paste0("Ecov_",i,'.png')), width = 10*144, height = 10*144, res = 144, pointsize = 12, family = "Times")

    m <- rbind(c(1,1), c(2,3))
    layout(m)
    par(mar=c(4,4,2,0), oma=c(0,0,0.5,0.5))

    ecov.low <- ecov.obs[,i] - 1.96 * ecov.obs.sig[,i]
    ecov.high <- ecov.obs[,i] + 1.96 * ecov.obs.sig[,i]
    y.min <- ifelse(min(ecov.low,na.rm=T) < 0, 1.1*min(ecov.low,na.rm=T), 0.9*min(ecov.low,na.rm=T))
    y.max <- ifelse(max(ecov.high,na.rm=T) < 0, 0.9*max(ecov.high,na.rm=T), 1.1*max(ecov.high,na.rm=T))
    if(max(ecov.pred[,i],na.rm=T) > y.max) y.max <- max(ecov.pred[,i],na.rm=T)
    if(min(ecov.pred[,i],na.rm=T) < y.min) y.min <- min(ecov.pred[,i],na.rm=T)
    plot(years, ecov.obs[,i], type='p', col=plot.colors[i], pch=1, xlab="Year", ylab=dat$Ecov_label[i],
         ylim=c(y.min, y.max))
    lines(years, ecov.pred[,i], col=plot.colors[i], lwd=3)
    arrows(years, ecov.low, years, ecov.high, length=0)
    title (paste0("Ecov ",i, ": ",dat$Ecov_label[i]), outer=T, line=-1)
    plot(years, ecov.res[,i], type='h', lwd=2, col=plot.colors[i], xlab="Year", ylab="Std. Residual")
    abline(h=0)
    hist(ecov.res[,i], breaks=10, plot=T, xlab="Std. Residual", ylab="Probability Density", freq=F, main=NULL)
    if(do.tex | do.png) dev.off() else par(origpar)
  }
}

#-------------------------------------------------------------------------------
