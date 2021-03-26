read_asap3_fit <- function(wd,asap.name, pSPR=40)  {
  #' @param wd directory where ASAP run is located (ex: "C:/MY/file/directories/model/" )
  #' @param asap.name Base name of original dat file (without the .dat extension)
  #' @param pSPR is user-specified SPR reference point, 40 by default (expressed as 100*SSBPR(Fspr)/SSBPR(F=0) ); used in ASAPplots function PlotAnnualSPRtargets

  #' @return years: numeric vector, model years only
  #' @return years_full: numeric vector, model + proj years (for ASAP this will be the same as $years)
  #' @return ages: numeric vector, e.g. 1:9
  #' @return selAA: list of length(n_selblocks), first the fleet blocks then indices; e.g. if 4 fleet blocks and 3 indices, selAA[[5]] is for index 1. Each element is a matrix, years (rows) x ages (cols)
  #' @return selblock_pointer_fleets: years x fleets, indices of selAA used by each fleet in each year
  #' @return selblock_pointer_indices: years x indices, indices of selAA used by each index in each year
  #' @return MAA: matrix, years x ages
  #' @return log_SSB: matrix nyears x 2; first column is log_SSB, second column is SE (from .std file in ADMB)
  #' @return log_F: matrix nyears x 2; first column is log_Fmult, second column is SE (from .std file in ADMB)
  #' @return log_NAA: matrix, years x ages
  #' @return NAA_CV: matrix, years x ages
  #' @return percentSPR: scalar, 40 by default
  #' @return log_Y_FXSPR: matrix nyears x 2; first column is log_yield at FXSPR, second column is SE log_yield at FXSPR
  #' @return log_FXSPR: matrix nyears x 2; first column is log_F at FXSPR, second column is SE log_F at FXSPR
  #' @return log_SSB_FXSPR: matrix nyears x 2; first column is log_SSB_FXSPR (annual numerator of SPR(FXSPR)*Recruits), second column is SE log_SSB_FXSPR
  #' @return log_rel_ssb_F_cov: list of length n_years, each element is a 2x2 covariance matrix, SSB/SSB_FXSPR first, F/F_FXSPR second

  
  # !!! need to check '\\' vs '/' defaults in ASAPplots  (made some checks below)
  # !!! check selAA and pointers on multifleet case  (selAA works; log_F is more complicated... needs more work; fine for 1 fleet)
  # !!! should check that user selected SRscalar is R0 not SSB0 (did not start this yet)
  # !!! need to load library(ASAPplots)  (brian stock will take care of this)
  
  # call some ASAPplots functions to get needed pieces
  #make sure wd has trailing '/' or '\\'
  wd.len <- nchar(wd)
  test=0
  if(substr(wd,wd.len, wd.len)=='/' ) test=1
  if(substr(wd,(wd.len-1), wd.len)=='\\' ) test=2
  if(test==0) wd <- ifelse(substr(wd, 3,3)=='/', paste0(wd, '/'), paste0(wd, '\\') )
  asap.rdat <- ifelse(  file.exists(paste0(wd, asap.name, '.RDAT')), paste0(wd, asap.name, '.RDAT'), paste0(wd, asap.name, '.rdat')  )
  asap <- dget(asap.rdat)
  
  gn <- GrabNames(wd,asap.name,asap)
  fleet.names <- gn$fleet.names
  index.names <- gn$index.names
  a1 <- GrabAuxFiles(wd,asap.name,asap,fleet.names,index.names)
  years <- seq( asap$parms$styr, asap$parms$endyr  )
  nyears <- length(years)
  years_full <- years
  ages <- seq(1, asap$parms$nages)
  nages <- length(ages)
  n.fleet <- asap$parms$nfleets
  n.fleet.sel <- max(asap$fleet.sel.blocks)
  n.index.sel <- length(asap$index.sel[,1])
  total.sel <- n.fleet.sel+n.index.sel
  selAA <- list()
   for (s in 1:total.sel) {
     if (s<=n.fleet.sel) {
       if(n.fleet==1) tmp.mat <- asap$fleet.sel.mats[[1]]
       if(n.fleet>1) tmp.mat <- matrix(unlist(asap$fleet.sel.mats), nrow=nyears*n.fleet, ncol=nages, byrow=T)
       selAA[[s]] <- tmp.mat[ which(asap$fleet.sel.blocks==s), ]
     } # end case s<=n.fleet.sel
     if (s>n.fleet.sel) {
       selAA[[s]] <- matrix(rep(asap$index.sel[(s-n.fleet.sel),], nyears), nrow=nyears, ncol=nages, byrow=T)
       rownames(selAA[[s]]) <- years
     }
     
   } #end s loop over selectivities
  
  
  selblock_pointer_fleets <- t(asap$fleet.sel.blocks) 
  selblock_pointer_indices <- matrix(rep(seq(1, asap$parms$nindices)+ n.fleet.sel, nyears ), nrow=nyears, ncol=asap$parms$nindices, byrow=T)
  rownames(selblock_pointer_indices) <- years
  
  MAA <- asap$M.age
  
  SSB.rows <- which(a1$asap.std$name=="SSB")
  log_SSB <- matrix(NA, nrow=nyears, ncol=2)
  log_SSB[,1] <- log(a1$asap.std$value[SSB.rows])
  SSB.CV <- a1$asap.std[SSB.rows, 4]/a1$asap.st[SSB.rows, 3]
  log_SSB[,2] <- sqrt( log(SSB.CV*SSB.CV +1 )   )
  
  logFmult1  <- unlist(a1$asap.std[ which(a1$asap.std$name=="log_Fmult_year1")  ,3] )
  logFmult_devs <- a1$asap.std[which(a1$asap.std$name=="log_Fmult_devs")  ,3]
  log_F <- matrix(NA, nrow=nyears*n.fleet, ncol=2)
  log_F[1,1] <- logFmult1
   
  log.fmult.rows <- which(substr(a1$asap.cor.names, 1, 9)=="log_Fmult")  #this is nyears*n.fleet
  log.fmult.cor <- a1$asap.cor.mat[log.fmult.rows, log.fmult.rows]  #this has dimension nyears*n.fleet x nyears*n.fleet
  diag(log.fmult.cor) <- rep(1, nyears*n.fleet)  
  log.fmult.std <- a1$asap.std[log.fmult.rows, 4]
  log.fmult.cov <- log.fmult.cor*(log.fmult.std %o% log.fmult.std) #creates var-cov matrix
  var.logFmult <- rep(NA, nyears*n.fleet)
  sigma.logFmult <- rep(NA, nyears*n.fleet)
  var.logFmult[1] <- log.fmult.cov[1,1] #need to modify for multifleet
  

  for (y in 2:nyears ) {
    
    log_F[y,1] <- log_F[y-1] + logFmult_devs[y-1]
    var.logFmult[y] <- var.logFmult[(y-1)] + a1$asap.std[log.fmult.rows[y],4]*a1$asap.std[log.fmult.rows[y],4]+2*(log.fmult.cov[y,(y-1)])
  }
  
  sigma.logFmult <- sqrt(var.logFmult)
  log_F[,2] <- sigma.logFmult
  
  log_NAA <- log(asap$N.age)
  NAA_CV <- matrix(NA, nrow=nyears, ncol=nages) # not sure how to derive this, but I can get it for age 1 recruits
  
  recr.rows <- which(a1$asap.std$name=="recruits")
  log_recruits <- log(a1$asap.std$value[recr.rows])
  recruits <- a1$asap.std$value[recr.rows]
  recruits.std <- a1$asap.std$stdev[recr.rows]
  log_recruits.std <- sqrt(log((recruits.std/recruits)*(recruits.std/recruits) +1)  )
  log_NAA <- log(asap$N.age)
  NAA_CV <- matrix(NA, nrow=nyears, ncol=nages) # not sure how to derive this, but I can get it for age 1 recruits
  NAA_CV[,1] <- log_recruits.std

  percentSPR <- pSPR   
 
  asap.spr <- PlotAnnualSPRtargets(asap, pspr=(pSPR/100), save.plots=F,od=paste0(wd,"plots"),plotf='png' )
 
  log_Y_FXSPR <- matrix(NA, nyears, 2)
  log_FXSPR <- matrix(NA, nyears, 2)
  log_SSB_FXSPR <- matrix(NA, nyears, 2)

  log_Y_FXSPR[,1] <- log(asap.spr$ypr.spr.vals)
  log_FXSPR[,1] <- log(asap.spr$f.spr.vals)
  

  log_SSB_FXSPR[,1] <- log(recruits*(pSPR/100)*asap$SR.annual.parms$s.per.r.vec)  
  log_SSB_FXSPR[,2] <- sqrt( (log_recruits.std*log_recruits.std)*(ssb_pr*ssb_pr) )
  
  log_rel_ssb_F_cov <- rep(list(matrix(NA, 2,2)), nyears ) 

  
ret.list <-list(years=years, years_full=years_full, ages=ages, selAA=selAA, selblock_pointer_fleets=selblock_pointer_fleets,
                selblock_pointer_indices=selblock_pointer_indices, MAA=MAA, log_SSB=log_SSB, log_F=log_F,
                log_NAA=log_NAA, NAA_CV=NAA_CV, percentSPR=percentSPR, log_Y_FXSPR=log_Y_FXSPR, log_FXSPR=log_FXSPR, 
                log_SSB_FXSPR=log_SSB_FXSPR, log_rel_ssb_F_cov=log_rel_ssb_F_cov )  
}  #end function read_asap3_fit
