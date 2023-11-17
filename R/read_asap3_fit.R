#' Read ASAP3 fit
#'
#' Gets output from a fit ASAP3 model for plotting with WHAM models.
#'
#' @param wd character, directory where ASAP3 output files are located (ex: 'C:/MY/file/directories/model/'). 5 files are needed: \code{.rdat}, \code{.dat}, \code{.std}, \code{.cor}, and \code{.par}.
#' @param asap.name character, base name of original .dat file (i.e. without the .dat extension)
#' @param pSPR scalar, user-specified percent SPR to use for reference points, expressed as 100*SSBPR(Fspr)/SSBPR(F=0). Default = 40.
#'
#' @return a named list with the following elements:
#'   \describe{
#'     \item{\code{$years}}{numeric vector, model years only, e.g. \code{1972:2020}}
#'     \item{\code{$years_full}}{numeric vector, model + proj years, e.g. \code{1972:2022}. For ASAP this will be the same as \code{$years}.}
#'     \item{\code{$selAA}}{list of length(n_selblocks), first the fleet blocks then indices, i.e. if 4 fleet blocks and 3 indices, \code{selAA[[5]]} is for index 1. Each element is a matrix, years (rows) x ages (cols), selectivity at age}
#'     \item{\code{$selblock_pointer_fleets}}{matrix, n_years x n_fleets, indices of selAA used by each fleet in each year}
#'     \item{\code{$selblock_pointer_indices}}{matrix, n_years x n_indices, indices of selAA used by each index in each year}
#'     \item{\code{$MAA}}{matrix, n_years x n_ages, natural mortality}
#'     \item{\code{$log_SSB}}{matrix, n_years x 2, log-scale spawning stock biomass. 1st col = MLE, 2nd col = SE (from .std file in ADMB).}
#'     \item{\code{$log_F}}{matrix, n_years x 2, log-scale fully-selected F. 1st col = MLE, 2nd col = SE (from .std file in ADMB).}
#'     \item{\code{$log_NAA}}{matrix, n_years x n_ages, numbers at age}
#'     \item{\code{$NAA_CV}}{matrix, n_years x n_ages, CV of numbers at age}
#'     \item{\code{$percentSPR}}{scalar, X\% SPR used to calculate reference points, default = 40}
#'     \item{\code{$log_Y_FXSPR}}{matrix, n_years x 2, log-scale yield at FXSPR. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_FXSPR}}{matrix, n_years x 2, log-scale FXSPR. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_SSB_FXSPR}}{matrix, n_years x 2, log-scale SSB at FXSPR, i.e. annual numerator of SPR(FXSPR) * Recruits. 1st col = MLE, 2nd col = SE.}
#'     \item{\code{$log_rel_ssb_F_cov}}{list, length n_years, each element is a 2x2 covariance matrix with SSB / SSB_FXSPR first and F / F_FXSPR second}
#'  }
#'
#' @export
#'
#' @seealso \code{\link{compare_wham_models}}, \code{\link{read_wham_fit}}
#'
#' @examples
#' \dontrun{
#' base <- read_asap3_fit(wd=file.path(getwd(),'asap_results'), asap.name='BASE_5C.DAT', pSPR=40)
#' m1 <- fit_wham(input1)
#' m2 <- fit_wham(input2)
#' mods <- list(base=base, m1=m1, m2=m2)
#' res <- compare_wham_models(mods)
#' }
#'
read_asap3_fit <- function(wd, asap.name, pSPR=40)  {
  # !!! check selAA and pointers on multifleet case  (selAA works; log_F is more complicated... needs more work; fine for 1 fleet)
  # !!! should check that user selected SRscalar is R0 not SSB0 (did not start this yet)
  asap.rdat <- list.files(wd, pattern=paste0(asap.name,".rdat"), ignore.case=TRUE, full.names=TRUE)
  if(length(asap.rdat)==0){
    stop("ASAP .rdat (or .RDAT) file not found in 'wd'")
  } else {
    asap <- dget(asap.rdat)
  }

  gn <- asap_grab_names(wd, asap.name, asap)
  fleet.names <- gn$fleet.names
  index.names <- gn$index.names
  a1 <- asap_grab_aux(wd, asap.name, asap, fleet.names, index.names)
  years <- seq( asap$parms$styr, asap$parms$endyr  )
  nyears <- length(years)
  years_full <- years
  ages <- seq(1, asap$parms$nages)
  nages <- length(ages)
  n.fleet <- asap$parms$nfleets
  n.fleet.sel <- max(asap$fleet.sel.blocks)
  n.index.sel <- length(asap$index.sel[,1])
  total.sel <- n.fleet.sel+n.index.sel
  selblock_pointer_fleets <- t(asap$fleet.sel.blocks)
  selblock_pointer_indices <- matrix(rep(seq(1, asap$parms$nindices)+ n.fleet.sel, nyears ), nrow=nyears, ncol=asap$parms$nindices, byrow=T)
  rownames(selblock_pointer_indices) <- years
  selAA <- list()
  for (s in 1:total.sel) {
    if (s<=n.fleet.sel) {
      selAA[[s]] <- matrix(NA, nrow=nyears*n.fleet, ncol=nages)
      if(n.fleet==1) selAA[[s]][which(asap$fleet.sel.blocks==s), ] <- asap$fleet.sel.mats[[1]][which(asap$fleet.sel.blocks==s), ]
      if(n.fleet>1) selAA[[s]] <- matrix(unlist(asap$fleet.sel.mats), nrow=nyears*n.fleet, ncol=nages, byrow=T)
    } # end case s<=n.fleet.sel
    if (s>n.fleet.sel) {
     selAA[[s]] <- matrix(rep(asap$index.sel[(s-n.fleet.sel),], nyears), nrow=nyears, ncol=nages, byrow=T)
     rownames(selAA[[s]]) <- years
    }
  } #end s loop over selectivities

  MAA <- asap$M.age
  SSB.rows <- which(a1$asap.std$name=="SSB")
  log_SSB <- matrix(NA, nrow=nyears, ncol=2)
  log_SSB[,1] <- log(a1$asap.std$value[SSB.rows])
  SSB.CV <- a1$asap.std[SSB.rows, 4]/a1$asap.st[SSB.rows, 3]
  log_SSB[,2] <- sqrt( log(SSB.CV*SSB.CV +1 )   )

  log.fmult1.rows <- which(a1$asap.std$name=="log_Fmult_year1")
  log.fmultdev.rows <- which(a1$asap.std$name=="log_Fmult_devs")
  log.fmult.rows <- integer()
  #put the indexing for F1, Fdevs in order by fleet
  log.fmult.rows[1+nyears*(0:(n.fleet-1))] <- log.fmult1.rows
  for(f in 1:n.fleet) log.fmult.rows[(f-1)*nyears + 2:nyears] <- log.fmultdev.rows[1:(nyears-1) + (f-1)*(nyears-1)]
  
  logFmult  <- unlist(a1$asap.std[log.fmult.rows,3] )
  log.fmult.cor <- a1$asap.cor.mat[log.fmult.rows, log.fmult.rows]  #this has dimension nyears*n.fleet x nyears*n.fleet
  log.fmult.cor[is.na(log.fmult.cor)] <- 0
  log.fmult.cor <- log.fmult.cor + t(log.fmult.cor)
  diag(log.fmult.cor) <- 1
  log.fmult.std <- a1$asap.std[log.fmult.rows, 4]
  log.fmult.cov <- log.fmult.cor*(log.fmult.std %o% log.fmult.std) #creates var-cov matrix
  var.logFmult <- numeric()
  for(f in 1:n.fleet) {
    #estimates and var are just cumulative sums
    logFmult[(f-1)*nyears + 1:nyears] <- cumsum(logFmult[(f-1)*nyears + 1:nyears])
    for(y in 1:nyears){
      var.logFmult[(f-1)*nyears + y] <- sum(log.fmult.cov[(f-1)*nyears + 1:y,(f-1)*nyears + 1:y]) 
    }
  }
  log_F <- cbind(logFmult,sqrt(var.logFmult))

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

  asap.spr <- asap_get_SPR(asap, pspr=(pSPR/100))

  log_Y_FXSPR <- matrix(NA, nyears, 2)
  log_FXSPR <- matrix(NA, nyears, 2)
  log_SSB_FXSPR <- matrix(NA, nyears, 2)

  log_SR_scalar.row <- which(a1$asap.std$name=="log_SR_scaler")
  log_R0.std <- a1$asap.std$stdev[log_SR_scalar.row]
  R0 <- asap$SR.annual.parms$R0.vec[nyears]

  # log_Y_FXSPR[,1] <- log(R0*asap.spr$ypr.spr.vals)
  log_Y_FXSPR[,1] <- log(mean(recruits)*asap.spr$ypr.spr.vals) # match wham default to use average estimated recruitment (all years)
  log_FXSPR[,1] <- log(asap.spr$f.spr.vals)

  # log_SSB_FXSPR[,1] <- log(R0*(pSPR/100)*asap$SR.annual.parms$s.per.r.vec)
  # log_SSB_FXSPR[,2] <- log_R0.std
  log_SSB_FXSPR[,1] <- log(mean(recruits)*(pSPR/100)*asap$SR.annual.parms$s.per.r.vec) # match wham default to use average estimated recruitment (all years)
  # log_SSB_FXSPR[,2] <- sd(log(recruits*(pSPR/100)*asap$SR.annual.parms$s.per.r.vec))/sqrt(length(recruits))

  #log_rel_ssb_F_cov <- rep(list(matrix(-99, 2,2)), nyears )
  log_rel_ssb_F_cov <- lapply(1:nyears, function(x) {
    list(c(log_SSB[x,1]-log_SSB_FXSPR[x,1], log_F[x,1]-log_FXSPR[x,1]), matrix(-99, 2,2))
  })


  return(list(years=years, years_full=years_full, selAA=selAA, selblock_pointer_fleets=selblock_pointer_fleets,
                selblock_pointer_indices=selblock_pointer_indices, MAA=MAA, log_SSB=log_SSB, log_F=log_F,
                log_NAA=log_NAA, NAA_CV=NAA_CV, percentSPR=percentSPR, log_Y_FXSPR=log_Y_FXSPR, log_FXSPR=log_FXSPR,
                log_SSB_FXSPR=log_SSB_FXSPR, log_rel_ssb_F_cov=log_rel_ssb_F_cov))
}  #end function read_asap3_fit

# helper functions from ASAPplots: https://github.com/cmlegault/ASAPplots
ypr <- function(nages, wgt.age, M.age, F.mult, sel.age ) {
  yield=0.0
  cum.survive=1.0
  z=0.0
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    yield=yield + wgt.age[i]*F.mult*sel.age[i]*(1-exp(-z) )*cum.survive/z
    cum.survive=cum.survive*exp(-z)
  }
  z= M.age[nages] + F.mult*sel.age[nages]
  yield=yield + wgt.age[nages]*F.mult*sel.age[nages]*cum.survive/z
  return(yield)
}

s.per.recr<-function(nages,fec.age,mat.age,M.age, F.mult, sel.age, spawn.time ) {
  spr=0.0
  cum.survive=1.0
  z=0.0
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    z.ts=(M.age[i]+F.mult*sel.age[i])*spawn.time
    spr=spr+cum.survive*fec.age[i]*mat.age[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
  }
  z= M.age[nages] + F.mult*sel.age[nages]
  z.ts=(M.age[nages]+F.mult*sel.age[nages])*spawn.time
  spr=spr + fec.age[nages]*mat.age[nages]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  return(spr)
}

asap_get_SPR <- function(asap, pspr){
  spr.targ.values <- pspr
  n.spr <- length(spr.targ.values)
  nages<- asap$parms$nages
  nyears <- asap$parms$nyears
  years <- seq(asap$parms$styr,asap$parms$endyr)
  fec.age <- asap$WAA.mats$WAA.ssb
  mat.age <- asap$maturity
  wgt.age <- asap$WAA.mats$WAA.catch.all
  M.age <- asap$M.age
  sel.age <- asap$F.age/apply(asap$F.age,1,max)
  spawn.time <- asap$options$frac.yr.spawn
  spr0.vals<- asap$SR.annual.parms$s.per.r.vec
  F.start <-0.11  # starting guess for optimization routine to find F_SPR%

  f.spr.vals <- matrix(NA, nyears, n.spr)
  ypr.spr.vals <- matrix(NA, nyears, n.spr)
  conv.vals <- matrix(NA, nyears, n.spr)

  for (i in 1:n.spr) {
    for (j in 1:nyears) {
      t.spr <- spr.targ.values[i]
      spr.f <- function(F.start) {
        abs(s.per.recr(nages=nages, fec.age=fec.age[j,], mat.age=mat.age[j,], M.age= M.age[j,], F.mult=F.start, sel.age=sel.age[j,], spawn.time=spawn.time)/spr0.vals[j] - t.spr )
      }
      yyy <- nlminb(start=F.start, objective=spr.f, lower=0, upper=3)
      f.spr.vals[j,i] <- yyy$par
      ypr.spr.vals[j,i] <- ypr(nages, wgt.age=wgt.age[j,], M.age=M.age[j,],  F.mult=f.spr.vals[j,i], sel.age=sel.age[j,] )
      conv.vals[j,i] <- ifelse(yyy$convergence==0, TRUE, FALSE)
    }  # end j-loop over nyears
  }  #end i-loop over SPR values

  return(list(f.spr.vals=f.spr.vals, ypr.spr.vals=ypr.spr.vals, conv.vals=conv.vals))
} # end get_asap_SPR

asap_grab_names <- function(wd, asap.name, asap){
  my.names <- list()
  my.names$fleet.names <- paste0("FLEET-",1:asap$parms$nfleets)
  my.names$index.names <- paste0("INDEX-",1:asap$parms$nindices)
  my.file.name <- list.files(wd, pattern=paste0(asap.name,".dat"), ignore.case=TRUE, full.names=TRUE)
  if (file.exists(my.file.name)){
    datfile <- readLines(con = my.file.name)
    nlines <- length(datfile)
    nfinis <- nlines-asap$parms$nfleets-asap$parms$navailindices-3
    if (datfile[nfinis] == "###### FINIS ######"){
      my.names$fleet.names <- substr(datfile[(nfinis+2):(nfinis+2+asap$parms$nfleets-1)],3,100)
      avail.index.names <- substr(datfile[(nfinis+3+asap$parms$nfleets):(nlines-1)],3,100)
      my.names$index.names <- avail.index.names[asap$initial.guesses$index.use.flag==1]
    }  # end if-test for nfinis
  } else {
    stop(paste0(asap.name,".dat or ",asap.name,".DAT file not found in 'wd'"))
  }
  return(my.names)
}

asap_grab_aux <- function(wd, asap.name, asap, fleet.names, index.names) {
  aux.list <-  list("npar"=-999, "asap.cor.names"=NA, "asap.cor.mat"=NA,
                    "asap.std"=NA, "max.grad"=-999, "F.rep"=NA, "Tot.B"=NA , "SSB"=NA,
                    "Expl.B"=NA, "recr"=NA, "asap.name"=asap.name )
  std.fname <- list.files(wd, pattern=paste0(asap.name,".std"), ignore.case=TRUE, full.names=TRUE)
  cor.fname <- list.files(wd, pattern=paste0(asap.name,".cor"), ignore.case=TRUE, full.names=TRUE)
  par.fname <- list.files(wd, pattern=paste0(asap.name,".par"), ignore.case=TRUE, full.names=TRUE)
  if(all(file.exists(std.fname), file.exists(cor.fname), file.exists(par.fname))){
    # Read in std file from admb
    asap.std <- read.table(std.fname, header = F, skip=1, stringsAsFactors = FALSE)
    names(asap.std) <- c("index", "name", "value", "stdev" )
    years <- seq(asap$parms$styr, asap$parms$endyr)

    # Read in cor file from admb
    ncol.cor <- dim(asap.std) [1]
    asap.cor <- read.table(cor.fname, header = F, skip=2, col.names=c("index", "name", "value", "stdev", seq(1,ncol.cor)), fill=T)
    asap.cor.mat <- as.matrix(asap.cor[,5:length(asap.cor[1,])])
    asap.cor.names <- as.character(as.vector(asap.cor[,2]) )
    levels.cor.names <- unique(asap.cor.names)
    F.rep <- which(asap.cor.names=="Freport")
    asap.cor.names[F.rep] <- paste(asap.cor.names[F.rep],years, sep=".")
    Tot.B <- which(asap.cor.names=="TotJan1B")
    asap.cor.names[Tot.B] <- paste(asap.cor.names[Tot.B],years, sep=".")
    SSB <- which(asap.cor.names=="SSB")
    asap.cor.names[SSB] <- paste(asap.cor.names[SSB],years, sep=".")
    Expl.B <- which(asap.cor.names=="ExploitableB")
    asap.cor.names[Expl.B] <- paste(asap.cor.names[Expl.B],years, sep=".")
    recr <- which(asap.cor.names=="recruits")
    asap.cor.names[recr] <- paste(asap.cor.names[recr],years, sep=".")

    N.yr1 <- which(asap.cor.names=="log_N_year1_devs")
    asap.cor.names[N.yr1] <- paste(asap.cor.names[N.yr1],"Age",seq(2,asap$parms$nages), sep=".")

    recr.devs <- which(asap.cor.names=="log_recruit_devs")
    asap.cor.names[recr.devs] <- paste(asap.cor.names[recr.devs],years, sep=".")
    Fmult.yr1 <- which(asap.cor.names=="log_Fmult_year1")
    asap.cor.names[Fmult.yr1] <- paste(asap.cor.names[Fmult.yr1],fleet.names, sep=".")
    Fmult.devs <- which(asap.cor.names=="log_Fmult_devs")
    asap.cor.names[Fmult.devs] <- paste(asap.cor.names[Fmult.devs],paste(fleet.names, rep(years[2:length(years)], each=asap$parms$nfleets),  sep="."), sep=".")
    q.yr1 <- which(asap.cor.names=="log_q_year1")
    asap.cor.names[q.yr1] <- paste(asap.cor.names[q.yr1],index.names, sep=".")

    q.devs <- which(asap.cor.names=="log_q_devs")
    #asap.cor.names[q.devs] <- paste(asap.cor.names[q.devs],index.names, sep=".")
    if (length(q.devs)>0) asap.cor.names[q.devs] <- paste(asap.cor.names[q.devs],years, sep=".")
    diag(asap.cor.mat) <- rep(NA, ncol.cor)

    # Read in par file from admb
    asap.grad <- readLines(par.fname, n=1)
    par.split<- unlist(strsplit(asap.grad, " "))
    max.grad <- par.split[length(par.split)]
    npar <- as.numeric(par.split[ 6])
    aux.list <-  list("npar"=npar, "asap.cor.names"=asap.cor.names, "asap.cor.mat"=asap.cor.mat,
                      "asap.std"=asap.std, "max.grad"=max.grad, "F.rep"=F.rep[1], "Tot.B"=Tot.B[1],
                      "SSB"=SSB[1], "Expl.B"=Expl.B[1] , "recr"=recr[1], "asap.name"=asap.name )
  } else {
    stop(paste0(asap.name,".std or ",asap.name,".STD file not found in 'wd'"))
  }
  return(aux.list)
}


