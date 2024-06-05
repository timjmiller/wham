#' Compare multiple WHAM (or ASAP) models
#'
#' After fitting multiple WHAM (or ASAP) models, \code{compare_wham_models} produces plots
#' and a table of AIC and Mohn's rho to aid model comparison.
#'
#' \code{plot.opts$which} specifies which plots to make:
#'   \describe{
#'     \item{1}{3-panel of SSB (spawning stock biomass), F (fully-selected fishing mortality), and Recruitment}
#'     \item{2}{CV (coefficient of variation) for SSB, F, and Recruitment}
#'     \item{3}{Fleet selectivity (by block, averaged across years)}
#'     \item{4}{Index selectivity (by block, averaged across years)}
#'     \item{5}{Selectivity tile (fleets + indices, useful for time-varying random effects)}
#'     \item{6}{M time series (natural mortality, can specify which age with plot.opts$M.age)}
#'     \item{7}{M tile (useful for time-varying random effects)}
#'     \item{8}{3-panel of F X\% SPR, SSB at F_X\%SPR, and yield at F_X\%SPR}
#'     \item{9}{2-panel of relative status (SSB / SSB at F_X\%SPR and F / F_X\%SPR)}
#'     \item{10}{Kobe status (relative SSB vs. relative F)}
#'   }
#' If \code{plot.opts$return.ggplot = TRUE}, a list \code{g} is returned holding the above ggplot2 objects for later modification.
#' \code{g[[i]]} holds the plot corresponding to \code{i} above, e.g. \code{g[[2]]} is the CV plot.
#'
#' @param mods (named) list of fit WHAM/ASAP models. To read in ASAP model output, use \code{\link{read_asap3_fit}}. If no names are given, 'm1', 'm2', ...
#' will be used.
#' @param do.table T/F, produce table of AIC and/or Mohn's rho? Default = TRUE.
#' @param do.plot T/F, produce plots? Default = TRUE.
#' @param fdir character, path to directory to save table and/or plots. Default = \code{getwd()}.
#' @param compare.opts list of options to generate comparison results:
#'   \describe{
#'     \item{\code{$stock}}{integer, which stock to include in results. Default = 1.}
#'     \item{\code{$region}}{integer, which region to include in results. Default = 1.}
#'   }
#' @param table.opts list of options for AIC/rho table:
#'   \describe{
#'     \item{\code{$fname}}{character, filename to save CSV results table (.csv will be appended). Default = \code{'model_comparison'}.}
#'     \item{\code{$sort}}{T/F, sort by AIC? Default = TRUE.}
#'     \item{\code{$calc.rho}}{T/F, calculate Mohn's rho? Retrospective analysis must have been run for all modes. Default = TRUE.}
#'     \item{\code{$calc.aic}}{T/F, calculate AIC? Default = TRUE.}
#'     \item{\code{$print}}{T/F, print table to console? Default = TRUE.}
#'     \item{\code{$save.csv}}{T/F, save table as a CSV file? Default = FALSE.}
#'   }
#' @param plot.opts list of options for plots:
#'   \describe{
#'     \item{\code{$out.type}}{character, either \code{'pdf'} or \code{'png'} (default = \code{'png'} because I am not sure \code{system('pdftk')} will work across platforms.)}
#'     \item{\code{$ci}}{vector of T/F, length = 1 (applied to all models) or number of models}
#'     \item{\code{$years}}{vector, which years to plot? Default = all (model and projection years).}
#'     \item{\code{$which}}{vector, which plots to make? Default = all. See details.}
#'     \item{\code{$relative.to}}{character, name of "base" model to plot differences relative to.}
#'     \item{\code{$alpha}}{scalar, (1-alpha)\% confidence intervals will be plotted. Default = 0.05 for 95\% CI.}
#'     \item{\code{$ages.lab}}{vector, overwrite model age labels.}
#'     \item{\code{$kobe.yr}}{integer, which year to use in Kobe plot (relative status). Default = terminal model year.}
#'     \item{\code{$M.age}}{integer, which age to use in M time-series plot. Default = \code{max(data$which_F_age)} (max age of F to use for full total F).}
#'     \item{\code{$return.ggplot}}{T/F, return a list of ggplot2 objects for later modification? Default = TRUE.}
#'     \item{\code{$kobe.prob}}{T/F, print probabilities for each model in each quadrant of Kobe plot? Default = TRUE.}
#'     \item{\code{$refpt}}{"XSPR" or "MSY", which reference point to use. Default = "XSPR".}
#'   }
#'
#' @return a list with the following components:
#'   \describe{
#'     \item{\code{daic}}{Vector of delta-AIC by model (if \code{do.table=T} and \code{table.opts$calc.aic=T})}
#'     \item{\code{aic}}{Vector of AIC by model (if \code{do.table=T} and \code{table.opts$calc.aic=T})}
#'     \item{\code{rho}}{Matrix of Mohn's rho by model (if \code{do.table=T} and \code{table.opts$calc.rho=T})}
#'     \item{\code{best}}{Name of best model (lowest AIC) (if \code{do.table=T} and \code{table.opts$calc.aic=T})}
#'     \item{\code{tab}}{Results table of AIC and Mohn's rho (if \code{do.table=T})}
#'     \item{\code{g}}{List of ggplot2 objects for later modification (if \code{do.plot=T} and \code{plot.opts$return.ggplot=T})}
#'   }
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{read_asap3_fit}, \code{\link{read_wham_fit}}}
#'
#' @examples
#' \dontrun{
#' base <- read_asap3_fit()
#' m1 <- fit_wham(input1)
#' m2 <- fit_wham(input2)
#' mods <- list(base=base, m1=m1, m2=m2)
#' res <- compare_wham_models(mods)
#' }
#'
#' @export
compare_wham_models <- function(mods, do.table=TRUE, do.plot=TRUE, fdir=getwd(), compare.opts=NULL, table.opts=NULL, plot.opts=NULL, fname=NULL, sort=NULL, calc.rho=NULL, calc.aic=NULL, do.print=NULL){
  if(is.null(names(mods))) names(mods) <- paste0("m",1:length(mods))
  if(!any(c(missing("fname"),missing("sort"),missing("calc.rho"),missing("calc.aic"),missing("do.print")))){
    warning("Argument deprecated, use table.opts=list(fname, sort, calc.rho, calc.aic, print, save.csv) instead.
      Using default: table.opts=list(fname='model_comparison', sort=TRUE, calc.rho=TRUE, calc.aic=TRUE, print=TRUE, save.csv=TRUE)")
  }
  out <- list() # to return
  wham.mods.ind <- which(sapply(mods, function(x) "wham_version" %in% names(x))) # otherwise asap model output read in from read_asap3_fit
  asap.mods.ind <- c(1:length(mods))[-wham.mods.ind] # otherwise asap model output read in from read_asap3_fit
  wham.mods <- mods[wham.mods.ind] # get wham models only
  asap.mods <- mods[asap.mods.ind] # get asap models only
  all.wham <- ifelse(length(wham.mods.ind)==length(mods), TRUE, FALSE)
  no.wham <- ifelse(length(wham.mods.ind)==0, TRUE, FALSE)

  if(is.null(compare.opts)) compare.opts <- list(stock = 1, region = 1)
  if(is.null(compare.opts$stock)) compare.opts$stock = 1
  if(is.null(compare.opts$region)) compare.opts$region = 1
  stock <- compare.opts$stock
  region <- compare.opts$region

  if(do.table){
    if(is.null(table.opts)) table.opts=list(fname = "model_comparison", sort = TRUE, calc.rho = TRUE, calc.aic = TRUE, print=TRUE, save.csv=TRUE)
    if(is.null(table.opts$fname)) table.opts$fname = "model_comparison"
    if(is.null(table.opts$sort)) table.opts$sort = TRUE
    if(is.null(table.opts$calc.rho)) table.opts$calc.rho = TRUE
    if(is.null(table.opts$calc.aic)) table.opts$calc.aic = TRUE
    if(is.null(table.opts$print)) table.opts$print = TRUE
    if(is.null(table.opts$save.csv)) table.opts$save.csv = FALSE
    aic.tab <- aic <- daic <- NULL
    if(no.wham) warning("No WHAM models. Consider setting do.table = FALSE.\n")
    if(!all.wham){
      cat("Cannot calculate AIC or Mohn's rho for ASAP3 models.\n
        Returning AIC/rho table for WHAM models only.\n")
    }
    if(table.opts$calc.aic & !no.wham) if(length(wham.mods)>1) {
      if(sum(mapply(function(x) x$env$data$Ecov_model==0, wham.mods)) %in% c(0,length(wham.mods))){
        ecov.obs <- lapply(wham.mods, function(x) x$env$data$Ecov_use_obs)
        all.identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))
        if(!all.identical(ecov.obs)){ #covariate obs in each model, but not the same covariates.
          warning("Inclusion of environmental covariate observations differs among models, cannot compare models' AIC. \n
            Setting 'calc.aic = FALSE' to make other comparisons. \n\n")
          table.opts$calc.aic = FALSE
        }
      } else { 
        warning("Inclusion of environmental covariate observations differs among models, cannot compare models' AIC. \n
          Setting 'calc.aic = FALSE' to make other comparisons. \n\n")
          table.opts$calc.aic = FALSE
      }
    }
    tab <- NULL
    if(table.opts$calc.aic){
      aic <- sapply(wham.mods, function(x){
        k = length(x$opt$par)
        2*(x$opt$obj + k) # AIC
        # 2*(x$opt$obj + k + k*(k+1)/(n-k-1)) # AICc
      })
      aic <- round(aic, 1)
      daic <- round(aic - min(aic), 1)
      aic.tab <- cbind(daic, aic)
      colnames(aic.tab) <- c("dAIC","AIC")
      tab <- cbind(aic.tab)
    }
    rho <- NULL
    if(table.opts$calc.rho){
      if(any(mapply(function(x) is.null(x$peels), wham.mods))){
        warning("Not all models have peels --> Cannot compare Mohn's rho. \n
             Setting 'calc.rho = FALSE' to make other comparisons. \n")
        table.opts$calc.rho = FALSE
      } else{
        rho <- t(sapply(wham.mods, function(x){
          rho <- mohns_rho(x)
          rho <- c(rho$naa[stock,region,1], rho$SSB[stock], rho$Fbar[region])
        }))#[ ,c("R","SSB","Fbar")]
        rho <- round(rho, 4)
        colnames(rho) <- paste0("rho_",c("R","SSB","Fbar"))
        # apply(rho, 1, function(y) mean(abs(y)))
        if(!is.null(tab)) tab <- cbind(tab, rho)
        else tab <- cbind(rho)
      }
    }

    best <- NULL
    if(table.opts$calc.aic) best <- names(wham.mods)[which(aic == min(aic))]
    if(table.opts$sort){
      ord <- order(aic)
      daic <- daic[ord]
      aic <- aic[ord]
      rho <- rho[ord,]
      if(!is.null(tab)) tab <- tab[ord,]
    }
    if(!is.null(tab)){ 
      if(table.opts$save.csv) write.csv(tab, file = paste0(file.path(fdir, table.opts$fname),".csv"))
      if(table.opts$print) print(tab) # print to console
    }
    out[c("daic","aic","rho","best","tab")] <- list(daic,aic,rho,best,tab)
  }
  if(do.plot){
    if(is.null(plot.opts)) plot.opts=list(out.type='png', ci=TRUE, years=NULL, which=1:10, relative.to=NULL, alpha=0.05, ages.lab=mods[[1]]$ages.lab, kobe.yr=NULL, M.age=NULL, return.ggplot=TRUE, kobe.prob=TRUE)
    if(is.null(plot.opts$out.type)) plot.opts$out.type <- 'png'
    if(!plot.opts$out.type %in% c("pdf","png")) stop("plot.opts$out.type must be 'pdf' or 'png' (default)")
    if(is.null(plot.opts$ci)) plot.opts$ci <- TRUE
    if(!length(plot.opts$ci) %in% c(1,length(mods))) stop("plot.opts$ci must have length = 1 or n.models")
    if(length(plot.opts$ci) == 1) plot.opts$ci <- rep(plot.opts$ci, length(mods))
    all.yrs <- unique(unlist(lapply(mods, function(x) x$years_full)))
    if(is.null(plot.opts$years)) plot.opts$years <- all.yrs
    if(!all(plot.opts$years %in% all.yrs)) stop("plot.opts$years must be a subset of $years_full in model fits")
    if(!is.null(plot.opts$relative.to)) if(!plot.opts$relative.to %in% names(mods)) stop("plot.opts$relative.to must match a model name, e.g. 'm1'.")
    if(is.null(plot.opts$which)) plot.opts$which <- 1:10
    if(!all(plot.opts$which %in% 1:10)) stop("All elements of plot.opts$which must be in 1:10. See ?compare_wham_models for available plots.")
    if(is.null(plot.opts$alpha)) plot.opts$alpha <- 0.05
    if(is.null(plot.opts[["kobe.yr"]])) plot.opts$kobe.yr <- tail(mods[[1]]$years, 1)
    if(!(plot.opts$kobe.yr %in% all.yrs)) stop("plot.opts$kobe.yr must be a year in $years_full from model fits.")
    if(is.null(plot.opts[["return.ggplot"]])) plot.opts$return.ggplot <- TRUE
    if(is.null(plot.opts[["kobe.prob"]])) plot.opts$kobe.prob <- TRUE
    if(is.null(plot.opts[["refpt"]])) plot.opts$refpt <- "XSPR"
    if(!plot.opts$refpt %in% c("XSPR","MSY")) stop("plot.opts$refpt must be either 'XSPR' or 'MSY'.")

    x <- list()
    for(i in 1:length(mods)){
      if(i %in% wham.mods.ind){
        x[[i]] <- read_wham_fit(mods[[i]])
        x[[i]]$is.wham = TRUE
      } else {
        x[[i]] <- mods[[i]] # already processed by read_asap3_fit
        x[[i]]$is.wham <- FALSE
      }
    }
    names(x) <- names(mods)
    if(is.null(plot.opts$ages.lab)){
      if(!no.wham) plot.opts$ages.lab <- wham.mods[[1]]$ages.lab
      if(is.null(plot.opts$ages.lab)) plot.opts$ages.lab <- paste0(1:dim(x[[1]]$MAA)[2], c(rep("",dim(x[[1]]$MAA)[2]-1),"+"))
    }
    if(is.null(plot.opts$M.age)){
      if(!no.wham) plot.opts$M.age <- max(wham.mods[[1]]$env$data$which_F_age)
      if(is.null(plot.opts$M.age)) plot.opts$M.age <- dim(x[[1]]$MAA)[2]
    }

    gg_facet_dims <- function(p){
      if(!is.null(p)){
        nrows <- p %>%
         ggplot2::ggplot_build() %>%
         magrittr::extract2('layout') %>%
         magrittr::extract2('layout') %>%
         magrittr::extract2('ROW') %>%
         unique() %>%
         length()
        ncols <- p %>%
         ggplot2::ggplot_build() %>%
         magrittr::extract2('layout') %>%
         magrittr::extract2('layout') %>%
         magrittr::extract2('COL') %>%
         unique() %>%
         length()
        return(c(2*nrows+2, 3*ncols+2))
      } else return(NULL)
    }
    g <-  vector("list", 10)
    no.SPR <- sapply(x, function(x) !("log_FXSPR" %in% names(x)))
    no.MSY <- sapply(x, function(x) !("log_FMSY" %in% names(x)))
    if((any(no.SPR) & plot.opts$refpt == "XSPR") | (any(no.MSY) & plot.opts$refpt == "MSY")){
      if(plot.opts$refpt == "XSPR"){
        warning("Not plotting SPR-based reference points and status because at least one model does not have them reported.")
      }
      if(plot.opts$refpt == "MSY"){
        warning("Not plotting MSY-based reference points and status because at least one model does not have them reported.")
      }
      plot.opts$which <- plot.opts$which[which(!(plot.opts$which %in% 8:10))]
    }
    for(i in plot.opts$which){
      if(i==1) g[[i]] <- suppressWarnings(plot.SSB.F.R.compare(x, compare.opts, plot.opts))
      if(i==2) g[[i]] <- suppressWarnings(plot.cv.compare(x, compare.opts, plot.opts))
      if(i==3) g[[i]] <- suppressWarnings(plot.selectivity.compare(x, plot.opts, type="fleet"))
      if(i==4) g[[i]] <- suppressWarnings(plot.selectivity.compare(x, plot.opts, type="indices"))
      if(i==5) g[[i]] <- suppressWarnings(plot.tile.compare(x, compare.opts, plot.opts, type="selAA"))
      if(i==6) g[[i]] <- suppressWarnings(plot.M.compare(x, compare.opts, plot.opts))
      if(i==7) g[[i]] <- suppressWarnings(plot.tile.compare(x, compare.opts, plot.opts, type="MAA"))
      if(i==8) g[[i]] <- suppressWarnings(plot.refpt.compare(x, plot.opts))
      if(i==9) g[[i]] <- suppressWarnings(plot.rel.compare(x, plot.opts))
    }
    pdims <- lapply(g, gg_facet_dims)
    pdims[[10]] <- c(7,7)
    plabs <- c("SSB_F_R","CV","sel_fleets","sel_indices","sel_tile",paste0("M_age",plot.opts$M.age),"M_tile","ref_pts","rel_status_timeseries","rel_status_kobe")
    if(plot.opts$out.type == 'pdf'){
      pnames <- file.path(fdir,"compare_pdf",paste0("compare_",plabs,".pdf"))
      if(!dir.exists(file.path(fdir,"compare_pdf"))) dir.create(file.path(fdir,"compare_pdf"))
    }
    if(plot.opts$out.type == 'png'){
      pnames <- file.path(fdir,"compare_png",paste0("compare_",plabs,".png"))
      if(!dir.exists(file.path(fdir,"compare_png"))) dir.create(file.path(fdir,"compare_png"))
    }
    for(i in plot.opts$which){
      if(plot.opts$out.type == 'pdf') grDevices::cairo_pdf(filename=pnames[i], height = pdims[[i]][1], width = pdims[[i]][2])
      if(plot.opts$out.type == 'png') png(pnames[i], width=pdims[[i]][2], height=pdims[[i]][1], units="in", res=300)
      if(i < 10) suppressWarnings(print(g[[i]]))
      if(i == 10) g[[i]] <- suppressWarnings(plot.kobe.compare(x, plot.opts))
      dev.off()
    }
    if(plot.opts$return.ggplot) out$g <- g
  }

  return(out)
}
fancy_scientific <- function(l) {
  if(max(l, na.rm=T) < 100){
    l <- format(l, scientific = FALSE, digits=2)
  } else {
    l <- format(l, scientific = TRUE, digits=2)
    l <- gsub("0e\\+00","0",l)
    l <- gsub("^(.*)e", "'\\1'e", l)
    l <- gsub("e\\+","e",l)
    l <- gsub("e", "%*%10^", l)
    l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  }
  parse(text=l)
}
# get.ci <- function(x, alpha=0.05, getci=TRUE){
#   if(!getci) x[,2] <- 0
#   ci <- exp(x[,1] + cbind(0, -qnorm(1-alpha/2)*x[,2], qnorm(1-alpha/2)*x[,2]))
#   ci[is.nan(ci[,2]),2] = ci[is.nan(ci[,2]),1]
#   ci[is.nan(ci[,3]),3] = ci[is.nan(ci[,3]),1]
#   return(ci)
# }
get.ci = function(par,se, p=0.975, lo = 0, hi = 1, type = "I", k = 1, alpha.ci=0.05, getci=TRUE, asap = FALSE){
  if(!getci) {
    ci <- list(lo = par, hi = par)
    ci$lo[] <- 0
    ci$hi[] <- 0
  } else {
    if(!asap) ci <- list(lo = par - qnorm(1-alpha.ci/2) * se, hi = par + qnorm(1-alpha.ci/2) * se)
    else ci <- list(lo = cbind(par[,1] - qnorm(1-alpha.ci/2)* par[,2]), hi = cbind(par[,1] + qnorm(1-alpha.ci/2)* par[,2]))
  }
  if(type == "I") {
  }
  if(type == "exp") {
    ci <- lapply(ci, exp)
  }
  if(type == "expit") { #Delta-method: V(lo + (hi-lo)/(1 + exp(-x))) ~ ((hi-lo) * p * (1-p))^2 * V(x)
    ci <- lapply(ci, function(x) lo + (hi-lo)/(1+ exp(-k * x)))
  }
  return(ci)
}
plot.timeseries.compare <- function(df, x, plot.opts){
  df <- df[df$Year %in% plot.opts$years,]
  df$Model <- factor(df$Model, levels=names(x), labels=names(x))
  df$Year <- as.integer(df$Year)
  if(any(plot.opts$ci)){
    df <- df %>% dplyr::group_by(var) %>% dplyr::mutate(y_max = 1.2*max(val)) %>% as.data.frame   # trim large CIs to zoom to 120% of max MLE
    ind <- which(df$hi > df$y_max)
    df$hi[ind] = df$y_max[ind]
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=val, color=Model, group=Model))
  if(any(plot.opts$ci)){
    g <- g + ggplot2::geom_ribbon(ggplot2::aes(x=Year, ymin=lo, ymax=hi, fill=Model), color=NA, alpha=.15) +
          ggplot2::scale_fill_viridis_d()
  }
  g <- g + ggplot2::geom_line(size=.8) +
          ggplot2::facet_wrap(ggplot2::vars(var), scales="free_y", ncol=1, strip.position = "left") +
          ggplot2::ylab(NULL) +
          ggplot2::scale_x_continuous(expand=c(0.01,0.01)) + # breaks=scales::breaks_extended(5)
          ggplot2::scale_colour_viridis_d() +
          ggplot2::theme_bw() +
          ggplot2::theme(strip.background = ggplot2::element_blank(), strip.placement = "outside",
                legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
  # if not relative, force y min to 0
  if(is.null(plot.opts$relative.to)){
    g <- g + ggplot2::scale_y_continuous(expand=c(0.01,0.01), limits = c(0,NA), labels=fancy_scientific)
  } else {
    g <- g + ggplot2::scale_y_continuous(expand=c(0.01,0.01), labels=fancy_scientific)
  }
  # if projections, add vline at terminal year
  last_proj <- sapply(x, function(x) tail(x$years_full,1))
  last_mod <- sapply(x, function(x) tail(x$years,1))
  v.yr <- unique(last_mod[last_proj > last_mod])
  g <- g + ggplot2::geom_vline(xintercept = v.yr, linetype=2, size=.4)

  return(g)
}
plot.SSB.F.R.compare <- function(x, compare.opts, plot.opts){
  stock <- compare.opts$stock
  region <- compare.opts$region
  df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(df) <- c("Year","var","val","lo","hi","Model")
  if(!is.null(plot.opts$relative.to)){
    base.i <- which(names(x) == plot.opts$relative.to)
    plot.opts$ci <- rep(FALSE, length(x)) #no confidence intervals
    for(i in 1:length(x)){
      if(x[[i]]$is.wham){
        if(x[[base.i]]$is.wham) ratio <- exp(x[[i]]$log_SSB_all$est)/exp(x[[base.i]]$log_SSB_all$est)
        else ratio <- exp(x[[i]]$log_SSB_all$est)/exp(x[[base.i]]$log_SSB[,1])
      } else {
        if(x[[base.i]]$is.wham) ratio <- exp(x[[i]]$log_SSB[,1])/exp(x[[base.i]]$log_SSB_all$est)
        else ratio <- exp(x[[i]]$log_SSB[,1])/exp(x[[base.i]]$log_SSB[,1])
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("SSB relative to ",plot.opts$relative.to),
        val=ratio, lo=ratio, hi=ratio, Model=names(x)[i]))
      if(x[[i]]$is.wham){
        if(x[[base.i]]$is.wham) ratio <- exp(x[[i]]$log_F_tot$est)/exp(x[[base.i]]$log_F_tot$est)
        else ratio <- exp(x[[i]]$log_F_tot$est)/exp(x[[base.i]]$log_F[,1])
      } else {
        if(x[[base.i]]$is.wham) ratio <- exp(x[[i]]$log_SSB[,1])/exp(x[[base.i]]$log_SSB_all$est)
        else ratio <- exp(x[[i]]$log_F[,1])/exp(x[[base.i]]$log_F[,1])
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("F relative to ",plot.opts$relative.to),
        val=ratio, lo=ratio, hi=ratio, Model=names(x)[i]))
      if(x[[i]]$is.wham){
        if(x[[base.i]]$is.wham) ratio <- exp(x[[i]]$log_NAA_rep$est[stock,region,,1])/exp(x[[base.i]]$log_NAA_rep$est[stock,region,,1])
        else ratio <- exp(x[[i]]$log_NAA_rep$est[stock,region,,1])/exp(x[[base.i]]$log_NAA[,1])
      } else {
        if(x[[base.i]]$is.wham) ratio <- exp(x[[i]]$log_NAA[,1])/exp(x[[base.i]]$log_NAA_rep$est[stock,region,,1])
        else ratio <- exp(x[[i]]$log_NAA[,1])/exp(x[[base.i]]$log_NAA[,1])
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("Recruitment relative to ",plot.opts$relative.to),
        val=ratio, lo=ratio, hi=ratio, Model=names(x)[i]))
    }
  } else {
    for(i in 1:length(x)){
      if(x[[i]]$is.wham){
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="SSB", val=exp(x[[i]]$log_SSB_all$est), lo=x[[i]]$log_SSB_all$ci$lo, 
          hi=x[[i]]$log_SSB_all$ci$hi, Model=names(x)[i]))
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="F", val=exp(x[[i]]$log_F_tot$est), lo=x[[i]]$log_F_tot$ci$lo, 
          hi=x[[i]]$log_F_tot$ci$hi, Model=names(x)[i]))
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="Recruitment", val=exp(x[[i]]$log_NAA_rep$est[stock,region,,1]), 
          lo=x[[i]]$log_NAA_rep$ci$lo[stock,region,,1], hi=x[[i]]$log_NAA_rep$ci$hi[stock,region,,1], Model=names(x)[i]))
      } else{
        ci <- get.ci(x[[i]][["log_SSB"]], plot.opts$alpha, plot.opts$ci[i], asap = TRUE)
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="SSB", val=exp(x[[i]][["log_SSB"]][,1]), lo=ci$lo[,1], hi=ci$hi[,1], Model=names(x)[i]))
        ci <- get.ci(x[[i]][["log_F"]], plot.opts$alpha, plot.opts$ci[i], asap = TRUE)
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="F", val=exp(x[[i]][["log_F"]][,1]), lo=ci$lo[,1], hi=ci$hi[,1], Model=names(x)[i]))
        ci <- get.ci(cbind(x[[i]]$log_NAA[,1], x[[i]]$NAA_CV[,1]), plot.opts$alpha, plot.opts$ci[i], asap = TRUE)
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="Recruitment", val=exp(x[[i]]$log_NAA[,1]), lo=ci$lo[,1], hi=ci$hi[,1], Model=names(x)[i]))
      }
    }
  }
  g <- plot.timeseries.compare(df, x, plot.opts)
  return(g)
}
plot.cv.compare <- function(x, compare.opts, plot.opts){
  stock <- compare.opts$stock
  region <- compare.opts$region
  plot.opts$ci <- rep(FALSE, length(x))
  df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(df) <- c("Year","var","val","lo","hi","Model")
  if(!is.null(plot.opts$relative.to)){
    base.i <- which(names(x) == plot.opts$relative.to)
    for(i in 1:length(x)){
      if(x[[i]]$is.wham){
        if(x[[base.i]]$is.wham) ratio <- x[[i]]$log_SSB_all$se/x[[base.i]]$log_SSB_all$se
        else ratio <- x[[i]]$log_SSB_all$se/x[[base.i]]$log_SSB[,2]
      } else {
        if(x[[base.i]]$is.wham) ratio <- x[[i]]$log_SSB[,2]/x[[base.i]]$log_SSB_all$se
        else ratio <- x[[i]]$log_SSB[,2]/x[[base.i]]$log_SSB[,2]
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("CV(SSB) relative to ",plot.opts$relative.to),
        val=ratio, lo=NA, hi=NA, Model=names(x)[i]))
      if(x[[i]]$is.wham){
        if(x[[base.i]]$is.wham) ratio <- x[[i]]$log_F_tot$se/x[[base.i]]$log_F_tot$se
        else ratio <- x[[i]]$log_F_tot$se/x[[base.i]]$log_F[,2]
      } else {
        if(x[[base.i]]$is.wham) ratio <- x[[i]]$log_F[,2]/x[[base.i]]$log_F_tot$se
        else ratio <- x[[i]]$log_F[,2]/x[[base.i]]$log_F[,2]
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("CV(F) relative to ",plot.opts$relative.to),
        val=ratio, lo=NA, hi=NA, Model=names(x)[i]))
      if(x[[i]]$is.wham){
        if(x[[base.i]]$is.wham) ratio <- x[[i]]$log_NAA_rep$se[stock,region,,1]/x[[base.i]]$log_NAA_rep$se[stock,region,,1]
        else ratio <- x[[i]]$log_NAA_rep$se[stock,region,,1]/x[[base.i]]$NAA_CV[,1]
      } else {
        if(x[[base.i]]$is.wham) ratio <- x[[i]]$NAA_CV[,1]/x[[base.i]]$log_NAA_rep$se[stock,region,,1]
        else ratio <- x[[i]]$NAA_CV[,1]/x[[base.i]]$NAA_CV[,1]
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("CV(Rec) relative to ",plot.opts$relative.to),
        val=ratio, lo=NA, hi=NA, Model=names(x)[i]))
    }
  } else {
    for(i in 1:length(x)){
      if(x[[i]]$is.wham){
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="CV(SSB)", val=x[[i]][["log_SSB_all"]]$se, lo=NA, hi=NA, Model=names(x)[i]))
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="CV(F)", val=x[[i]][["log_F_tot"]]$se, lo=NA, hi=NA, Model=names(x)[i]))
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="CV(Rec)", val=x[[i]][["log_NAA_rep"]]$se[stock,region,,1], lo=NA, hi=NA, Model=names(x)[i]))
      } else {
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="CV(SSB)", val=x[[i]][["log_SSB"]][,2], lo=NA, hi=NA, Model=names(x)[i]))
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="CV(F)", val=x[[i]][["log_F"]][,2], lo=NA, hi=NA, Model=names(x)[i]))
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var="CV(Rec)", val=x[[i]][["NAA_CV"]][,1], lo=NA, hi=NA, Model=names(x)[i]))
      }
    }
  }
  g <- plot.timeseries.compare(df, x, plot.opts)
  return(g)
}
plot.refpt.compare <- function(x, plot.opts){
  df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(df) <- c("Year","var","val","lo","hi","Model")
  nms.df <- c("F_refpt","SSB_refpt","Y_refpt")
  if(plot.opts$refpt == "XSPR"){
    nms <- c("log_FXSPR","log_SSB_FXSPR","log_Y_FXSPR")
  } else{ # MSY ref pt
    nms <- c("log_FMSY","log_SSB_MSY","log_MSY")
  }
  if(!is.null(plot.opts$relative.to)){
    plot.opts$ci <- rep(FALSE, length(x))
    base.i <- which(names(x) == plot.opts$relative.to)
    for(i in 1:length(x)){
      for(j in 1:length(nms)){
        if(x[[base.i]]$is.wham) {
          denom <- cbind(exp(x[[base.i]][[nms[j]]]$est[,]))
          denom <- denom[,NCOL(denom)]
        } else denom <- exp(x[[base.i]][[nms[j]]][,1])
        if(x[[i]]$is.wham) {
          num <- cbind(exp(x[[i]][[nms[j]]]$est))
          num <- num[,NCOL(num)]
        } else num <- exp(x[[i]][[nms[j]]][,1])
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=nms.df[j], #changed later
          val=num/denom, lo=num/denom, hi=num/denom, Model=names(x)[i]))
      }
    }
  } else {
    for(i in 1:length(x)){
      for(j in 1:length(nms)){
        if(x[[i]]$is.wham) {
          val <- cbind(exp(x[[i]][[nms[j]]]$est))
          col <- NCOL(val)
          val <- val[,col]
          lo <- cbind(x[[i]][[nms[j]]]$ci$lo)[,col]
          hi <- cbind(x[[i]][[nms[j]]]$ci$hi)[,col]
        }
        else {
          val <- exp(x[[i]][[nms[j]]][,1])
          ci <- get.ci(x[[i]][[nms[j]]], plot.opts$alpha, plot.opts$ci[i], asap = TRUE)
          lo <- ci$lo
          hi <- ci$hi
        }
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=nms.df[j], val=val, lo=lo, hi=hi, Model=names(x)[i]))
      }
    }
  }
  if(plot.opts$refpt == "XSPR"){
    pSPR <- sapply(x, function(y) y$percentSPR)
    if(length(unique(pSPR)) != 1) stop("FXSPR does not make sense to compare because percent SPR is not equal for all models.")
    pSPR <- as.character(pSPR[1])
  }
  if(!is.null(plot.opts$relative.to)) {
    if(plot.opts$refpt == "XSPR"){
      df$var <- factor(df$var, levels=nms.df, #c("F_refpt","SSB_refpt","Y_refpt"),
                           labels=c(bquote(italic(F)[paste(.(pSPR), "%")] ~relative~to~.(plot.opts$relative.to)),
                                    bquote(paste('SSB (', italic(F)[paste(.(pSPR), "%")],')')~relative~to~.(plot.opts$relative.to)),
                                    bquote(paste('Yield (',italic(F)[paste(.(pSPR), "%")], ')')~relative~to~.(plot.opts$relative.to))))
    } else {
      df$var <- factor(df$var, levels=nms.df,#c("F_refpt","SSB_refpt","Y_refpt"),
                           labels=c(bquote(italic(F)[MSY] ~relative~to~.(plot.opts$relative.to)),
                                    bquote(SSB[MSY]~relative~to~.(plot.opts$relative.to)),
                                    bquote(MSY~relative~to~.(plot.opts$relative.to))))      
    }
  } else {
    if(plot.opts$refpt == "XSPR"){
      df$var <- factor(df$var, levels=nms.df, #c("F_refpt","SSB_refpt","Y_refpt"),
                           labels=c(bquote(italic(F)[paste(.(pSPR), "%")]),
                                    bquote(paste('SSB (', italic(F)[paste(.(pSPR), "%")],')')),
                                    bquote(paste('Yield (',italic(F)[paste(.(pSPR), "%")], ')'))))
    } else {
      df$var <- factor(df$var, levels=nms.df, #c("F_refpt","SSB_refpt","Y_refpt"),
                           labels=c(bquote(italic(F)[MSY]),
                                    bquote(SSB[MSY]),
                                    bquote(MSY)))      
    }
  }
  g <- plot.timeseries.compare(df, x, plot.opts)
  g <- g + ggplot2::facet_wrap(ggplot2::vars(var), scales="free_y", ncol=1, strip.position = "left", labeller = ggplot2::label_parsed)
  return(g)
}
plot.rel.compare <- function(x, plot.opts){
  df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(df) <- c("Year","var","val","lo","hi","Model")
  if(!is.null(plot.opts$relative.to)){
    plot.opts$ci <- rep(FALSE, length(x))
    base.i <- which(names(x) == plot.opts$relative.to)
    vars <- c("relSSB", "relF")
    if(plot.opts$refpt == "XSPR"){
      for(i in 1:length(x)){
        for(j in 1:length(vars)){
          vals <- t(sapply(x[[i]][["log_rel_ssb_F_cov"]], function(x) c(x[[1]][j], sqrt(x[[2]][j,j]))))
          vals2 <- t(sapply(x[[base.i]][["log_rel_ssb_F_cov"]], function(x) c(x[[1]][j], sqrt(x[[2]][j,j]))))
          df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=vars[j],
            val=exp(vals[,1])/exp(vals2[,1]), lo=NA, hi=NA, Model=names(x)[i]))
        }
      }
    } else { # msy ref pt
      for(i in 1:length(x)){
        for(j in 1:length(vars)){
          vals <- t(sapply(x[[i]][["log_rel_ssb_F_cov_msy"]], function(x) c(x[[1]][j], sqrt(x[[2]][j,j]))))
          vals2 <- t(sapply(x[[base.i]][["log_rel_ssb_F_cov_msy"]], function(x) c(x[[1]][j], sqrt(x[[2]][j,j]))))
          df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=vars[j],
            val=exp(vals[,1])/exp(vals2[,1]), lo=NA, hi=NA, Model=names(x)[i]))
        }
      }      
    }
    df$lo = df$val
    df$hi = df$val
  } else {
    vars <- c("relSSB", "relF")
    if(plot.opts$refpt == "XSPR"){
      for(i in 1:length(x)){
        for(j in 1:length(vars)){
          vals <- t(sapply(x[[i]][["log_rel_ssb_F_cov"]], function(x) c(x[[1]][j], sqrt(x[[2]][j,j]))))
          ci <- exp(vals[,1] +qnorm(1-plot.opts$alpha/2) * cbind(-vals[,2],vals[,2]))
          df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=vars[j],
            val=exp(vals[,1]), lo=ci[,1], hi=ci[,2], Model=names(x)[i]))
        }
      }
    } else { # msy ref pt
      for(i in 1:length(x)){
        for(j in 1:length(vars)){
          vals <- t(sapply(x[[i]][["log_rel_ssb_F_cov_msy"]], function(x) c(x[[1]][j], sqrt(x[[2]][j,j]))))
          ci <- exp(vals[,1] +qnorm(1-plot.opts$alpha/2) * cbind(-vals[,2],vals[,2]))
          df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=vars[j],
            val=exp(vals[,1]), lo=ci[,1], hi=ci[,2], Model=names(x)[i]))
        }
      }
    }
  }
  if(plot.opts$refpt == "XSPR"){
    pSPR <- sapply(x, function(y) y$percentSPR)
    if(length(unique(pSPR)) != 1) stop("Percent SPR is not equal for all models.")
    pSPR <- as.character(pSPR[1])
  }
  if(!is.null(plot.opts$relative.to)) {
    if(plot.opts$refpt == "XSPR"){
      df$var <- factor(df$var, levels=c("relF","relSSB"),
                           labels=c(bquote(italic(F) / italic(F)[paste(.(pSPR), "%")]~relative~to~.(plot.opts$relative.to)),
                                    bquote(SSB / SSB[paste(.(pSPR), "%")]~relative~to~.(plot.opts$relative.to))))
    } else { # msy ref pt
      df$var <- factor(df$var, levels=c("relF","relSSB"),
                           labels=c(bquote(italic(F) / italic(F)[MSY]~relative~to~.(plot.opts$relative.to)),
                                    bquote(SSB / SSB[MSY]~relative~to~.(plot.opts$relative.to))))
    }
  } else {
    if(plot.opts$refpt == "XSPR"){
      df$var <- factor(df$var, levels=c("relF","relSSB"),
                           labels=c(bquote(italic(F) / italic(F)[paste(.(pSPR), "%")]),
                                    bquote(SSB / SSB[paste(.(pSPR), "%")])))
    } else { # msy ref pt
      df$var <- factor(df$var, levels=c("relF","relSSB"),
                           labels=c(bquote(italic(F) / italic(F)[MSY]),
                                    bquote(SSB / SSB[MSY])))
    }
  }
  g <- plot.timeseries.compare(df, x, plot.opts)
  g <- g + ggplot2::facet_wrap(ggplot2::vars(var), scales="free_y", ncol=1, strip.position = "left", labeller = ggplot2::label_parsed)
  if(is.null(plot.opts$relative.to)){
    g <- g + ggplot2::geom_hline(yintercept = 1, linetype=2, size=.4) +
           ggplot2::geom_hline(data = df %>% dplyr::filter(var == levels(df$var)[2]) %>% dplyr::mutate(half=0.5), mapping=ggplot2::aes(yintercept = half), color="red", linetype=2, size=.4)
  }
  return(g)
}
#-------------------------------------------------------------------------------
# 2D tile plot by age and year (e.g. selAA, MAA)
plot.tile.compare <- function(x, compare.opts, plot.opts, type="selAA"){
  n_ages = length(plot.opts$ages.lab)
  stock <- compare.opts$stock
  region <- compare.opts$region  

  if(type=="selAA"){
    n_years = length(x[[1]]$years)
    df <- data.frame(matrix(NA, nrow=0, ncol=5))
    colnames(df) <- c("Year","Age","Selectivity","Fleet","Model")
    for(j in 1:length(x)){
      selblock_pointer_all <- cbind(x[[j]]$selblock_pointer_fleets, x[[j]]$selblock_pointer_indices)
      n_selblocks <- length(x[[j]]$selAA)
      n_f <- dim(x[[j]]$selblock_pointer_fleets)[2]
      n_i <- dim(x[[j]]$selblock_pointer_indices)[2]
      selblocks_f <- as.numeric(unique(x[[j]]$selblock_pointer_fleets))
      selAA.byfleet <- vector("list", n_f+n_i)
      for(i in 1:(n_f+n_i)){
        selAA.byfleet[[i]] <- matrix(NA, nrow=n_years, ncol=n_ages)
        for(y in 1:n_years){
          selAA.byfleet[[i]][y,] <- x[[j]]$selAA[[selblock_pointer_all[y,i]]][y,]
        }
      }
      names(selAA.byfleet) <- c(paste0("Fleet ",1:n_f), paste0("Index ",1:n_i))

      df.selAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
      colnames(df.selAA) <- c(paste0("Age_",1:n_ages),"Year","Fleet")
      for(i in 1:(n_f+n_i)){
        tmp = as.data.frame(selAA.byfleet[[i]])
        tmp$Year <- x[[j]]$years
        colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
        tmp$Fleet = names(selAA.byfleet)[i]
        df.selAA <- rbind(df.selAA, tmp)
      }
      df <- rbind(df, df.selAA %>% tidyr::pivot_longer(-c(Year,Fleet),
                                    names_to = "Age",
                                    names_prefix = "Age_",
                                    names_ptypes = list(Age = character()),
                                    values_to = "Selectivity") %>%
                                  dplyr::mutate(Model = names(x)[j]))
    }
    df$Age <- as.integer(df$Age)
    df$Fleet <- factor(as.character(df$Fleet), levels=names(table(df$Fleet)))
    g <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=Age, fill=Selectivity)) +
        ggplot2::facet_grid(cols=ggplot2::vars(Fleet), rows=ggplot2::vars(Model))
  }

  if(type=="MAA"){
    n_years = length(x[[1]]$years_full)
    df <- data.frame(matrix(NA, nrow=0, ncol=4))
    colnames(df) <- c("Year","Age","M","Model")
    for(j in 1:length(x)){
      if(x[[j]]$is.wham){
        MAA <- as.data.frame(x[[j]][["MAA"]][stock,region,,])
      } else {
        MAA <- as.data.frame(x[[j]][["MAA"]])
      }
      MAA$Year <- x[[j]]$years_full
      colnames(MAA) <- c(paste0("Age_",1:n_ages),"Year")
      df <- rbind(df, MAA %>% tidyr::pivot_longer(-Year,
                names_to = "Age",
                names_prefix = "Age_",
                names_ptypes = list(Age = character()),
                values_to = "M") %>%
            dplyr::mutate(Model = names(x)[j]))
    }
    df$Age <- as.integer(df$Age)
    g <- ggplot2::ggplot(df, ggplot2::aes(x=Year, y=Age, fill=M)) +
         ggplot2::facet_wrap(ggplot2::vars(Model), dir="v")
  }
  g <- g + ggplot2::geom_tile() +
        ggplot2::scale_x_continuous(expand=c(0,0)) +
        ggplot2::scale_y_continuous(expand=c(0,0), breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
        ggplot2::theme_bw() +
        viridis::scale_fill_viridis()
  return(g)
}
plot.kobe.compare <- function(x, plot.opts){
  status.year.ind <- sapply(x, function(x) which(x$years_full == plot.opts$kobe.yr))
  #if any of the models have non-finite values for cov, then don't do that status years
  do.kobe <- sapply(status.year.ind, function(y) !any(sapply(x, function(z) any(!is.finite(z$log_rel_ssb_F_cov[[y]][[2]])))))
  #do.kobe <- unlist(sapply(mapply(function(x,i) x$log_rel_ssb_F_cov[i], x, status.year.ind), function(y) !all(!is.finite(y)))) # only if some non-infinite values for at least some status years
  if(any(do.kobe)){ #length(do.kobe) == length(x)
    #fxn <- function(y,i) y[["log_F"]][i,1]-y[["log_FXSPR"]][i,1]
    fxn <- function(y,i,j) y$log_rel_ssb_F_cov[[i]][[1]][j] #i: yr, j: rel ssb or f 
    rel.f.vals <- mapply(fxn, x, status.year.ind, 2)
    rel.ssb.vals <- mapply(fxn, x, status.year.ind, 1)
    fxn <- function(y,i) y$log_rel_ssb_F_cov[[i]][[2]] #i: yr
    log.rel.ssb.rel.F.cov <- mapply(fxn, x, status.year.ind, SIMPLIFY=F)

    log.rel.ssb.rel.F.ci.regs <- lapply(1:length(x), function(x){
      if(is.na(rel.f.vals[x]) | any(diag(log.rel.ssb.rel.F.cov[[x]])<0)) return(matrix(NA,100,2))
      else return(exp(ellipse::ellipse(log.rel.ssb.rel.F.cov[[x]], centre = c(rel.ssb.vals[x],rel.f.vals[x]), level = 1-plot.opts$alpha)))
      })
    p.ssb.lo.f.lo <- p.ssb.lo.f.hi <- p.ssb.hi.f.lo <- p.ssb.hi.f.hi <- rep(NA,length(status.year.ind))
    for(i in 1:length(status.year.ind)){
      check.zero.sd <- diag(log.rel.ssb.rel.F.cov[[i]])==0 | diag(log.rel.ssb.rel.F.cov[[i]]) < 0
      if(!any(is.na(check.zero.sd))) if(!any(check.zero.sd)){
        p.ssb.lo.f.lo[i] <- mnormt::sadmvn(lower = c(-Inf,-Inf), upper = c(-log(2), 0), mean = c(rel.ssb.vals[i],rel.f.vals[i]), varcov = log.rel.ssb.rel.F.cov[[i]])
        p.ssb.lo.f.hi[i] <- mnormt::sadmvn(lower = c(-Inf,0), upper = c(-log(2), Inf), mean = c(rel.ssb.vals[i],rel.f.vals[i]), varcov = log.rel.ssb.rel.F.cov[[i]])
        p.ssb.hi.f.lo[i] <- mnormt::sadmvn(lower = c(-log(2),-Inf), upper = c(Inf, 0), mean = c(rel.ssb.vals[i],rel.f.vals[i]), varcov = log.rel.ssb.rel.F.cov[[i]])
        p.ssb.hi.f.hi[i] <- mnormt::sadmvn(lower = c(-log(2),0), upper = c(Inf, Inf), mean = c(rel.ssb.vals[i],rel.f.vals[i]), varcov = log.rel.ssb.rel.F.cov[[i]])
      }
    }

    vals <- exp(cbind(rel.ssb.vals, rel.f.vals))
    max.x <- max(sapply(log.rel.ssb.rel.F.ci.regs, function(x) max(x[,1],na.rm = TRUE)),1.25, vals[,1])
    max.y <- max(sapply(log.rel.ssb.rel.F.ci.regs, function(x) max(x[,2],na.rm = TRUE)),1.25, vals[,2])

    plot(vals[,1],vals[,2], ylim = c(0,1.05*max.y), xlim = c(0,1.05*max.x), xlab = bquote(paste("SSB / ", SSB[paste(.(x[[1]]$percentSPR),"%")])),
      ylab = bquote(paste(italic(F)," / ", italic(F)[paste(.(x[[1]]$percentSPR),"%")])),type = 'n')
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
    if(plot.opts$kobe.prob){
      legend("topleft", legend = paste0(names(x), ": Prob = ", round(p.ssb.lo.f.hi,2)), bty = "n", cex=0.7)
      legend("topright", legend = paste0(names(x), ": Prob = ", round(p.ssb.hi.f.hi,2)), bty = "n", cex=0.7)
      legend("bottomleft", legend = paste0(names(x), ": Prob = ", round(p.ssb.lo.f.lo,2)), bty = "n", cex=0.7)
      legend("bottomright", legend = paste0(names(x), ": Prob = ", round(p.ssb.hi.f.lo,2)), bty = "n", cex=0.7)
    }
    text(vals[,1],vals[,2], paste0(rownames(vals)," (",plot.opts$kobe.yr,")"), cex=0.7)
    for(i in 1:length(status.year.ind)) polygon(log.rel.ssb.rel.F.ci.regs[[i]][,1], log.rel.ssb.rel.F.ci.regs[[i]][,2], lty=i)#, border = gray(0.7))
    return(list(rel.status = vals, p.ssb.lo.f.lo = p.ssb.lo.f.lo, p.ssb.hi.f.lo = p.ssb.hi.f.lo, p.ssb.hi.f.hi = p.ssb.hi.f.hi, p.ssb.lo.f.hi = p.ssb.lo.f.hi))
  } else {
    return(NULL)
  }
}
plot.M.compare <- function(x, compare.opts, plot.opts){
  stock <- compare.opts$stock
  region <- compare.opts$region  
  plot.opts$ci <- rep(FALSE, length(x))
  df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(df) <- c("Year","var","val","lo","hi","Model")
  if(!is.null(plot.opts$relative.to)){
    base.i <- which(names(x) == plot.opts$relative.to)
    for(i in 1:length(x)){
      if(x[[base.i]]$is.wham){
        if(x[[i]]$is.wham){
          ratio <- x[[i]][["MAA"]][stock,region,,plot.opts$M.age]/x[[i]][["MAA"]][stock,region,,plot.opts$M.age]
        } else{
          ratio <- x[[i]][["MAA"]][,plot.opts$M.age]/x[[i]][["MAA"]][stock,region,,plot.opts$M.age]
        }
      } else{
        if(x[[i]]$is.wham){
          ratio <- x[[i]][["MAA"]][stock,region,,plot.opts$M.age]/x[[base.i]][["MAA"]][,plot.opts$M.age]
        } else{
          ratio <- x[[i]][["MAA"]][,plot.opts$M.age]/x[[base.i]][["MAA"]][,plot.opts$M.age]
        }        
      }
      df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("M at age ",plot.opts$M.age," relative to ",plot.opts$relative.to),
                                 val=ratio, lo=NA, hi=NA, Model=names(x)[i]))
    }
  } else {
    for(i in 1:length(x)){
      if(x[[i]]$is.wham){
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("M at age ",plot.opts$M.age), val=x[[i]][["MAA"]][stock,region,,plot.opts$M.age], lo=NA, hi=NA, Model=names(x)[i]))
      } else {
        df <- rbind(df, data.frame(Year=x[[i]]$years_full, var=paste0("M at age ",plot.opts$M.age), val=x[[i]][["MAA"]][,plot.opts$M.age], lo=NA, hi=NA, Model=names(x)[i]))
      }
    }
  }
  g <- plot.timeseries.compare(df, x, plot.opts)
  return(g)
}
plot.selectivity.compare <- function(x, plot.opts, type="fleet"){
  plot.opts$ci <- rep(FALSE, length(x))
  n_ages = length(plot.opts$ages.lab)
  allSame <- function(x) length(unique(x)) == 1
  if(type == 'fleet'){
    if(!allSame(lapply(x, function(y) unname(y$selblock_pointer_fleets)))) {
      cat("Fleet selectivity blocks not identical, cannot produce comparison plot.")
      return(NULL)
    } else {
      selblocks <- unique(as.integer(x[[1]]$selblock_pointer_fleets))
      yrs <- lapply(selblocks, function(y){
                                  tmp <- which(x[[1]]$selblock_pointer_fleets == y);
                                  tmp <- tmp %% length(x[[1]]$years);
                                  tmp[tmp == 0] = length(x[[1]]$years);
                                  return(tmp)})
    }
  }
  if(type == 'indices'){
    if(!allSame(lapply(x, function(y) unname(y$selblock_pointer_indices)))) stop("Index selectivity blocks not identical, cannot produce comparison plot")
    selblocks <- unique(as.integer(x[[1]]$selblock_pointer_indices))
    yrs <- lapply(selblocks, function(y){
                                tmp <- which(x[[1]]$selblock_pointer_indices == y);
                                tmp <- tmp %% length(x[[1]]$years);
                                tmp[tmp == 0] = length(x[[1]]$years);
                                return(tmp)})
  }
  df <- data.frame(matrix(NA, nrow=0, ncol=6))
  colnames(df) <- c("Age","Selectivity","Block","Model")
  for(i in 1:length(x)){
    for(j in selblocks){
      sel <- apply(x[[i]][["selAA"]][[j]][yrs[[which(selblocks==j)]],], 2, mean, na.rm=T)
      df <- rbind(df, data.frame(Age=1:n_ages, Selectivity=sel, Block=paste0("Block ",j), Model=names(x)[i]))
    }
  }
  df$Model <- factor(df$Model, levels=names(x), labels=names(x))
  df$Block <- as.factor(df$Block)
  g <- ggplot2::ggplot(df, ggplot2::aes(x=Age, y=Selectivity, color=Model, group=Model)) +
          ggplot2::geom_line(size=.8) +
          ggplot2::facet_wrap(ggplot2::vars(Block), ncol=1, strip.position = 'right') +
          ggplot2::ylab("Selectivity") +
          ggplot2::scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01), labels=scales::number_format(accuracy=0.1)) +
          ggplot2::scale_x_continuous(expand=c(0.01,0.01), breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
          ggplot2::scale_colour_viridis_d() +
          ggplot2::theme_bw() +
          ggplot2::theme(legend.position="top", legend.box.margin = ggplot2::margin(0,0,0,0), legend.margin = ggplot2::margin(0,0,0,0))
  return(g)
}

