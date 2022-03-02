#' Plot WHAM output
#'
#' Generates many output plots and tables for a fit WHAM model.
#'
#' \code{out.type = 'pdf'} makes one pdf file of all plots. \code{out.type = 'png'} (default)
#' creates a subdirectory `plots_png`` in \code{dir.main} and saves .png files within.
#' \code{out.type = 'html'} makes a html files for viewing plot .png files and html tables of parameter estimates in a browser.
#' \code{out.type = 'pdf' or 'png'} makes LaTeX and pdf files of tables of parameter estimates.
#' (tabs: 'input data', 'diagnostics', 'results', 'ref_points', 'retro', and 'misc').
#' 
#' \code{plot.opts} holds optional arguments to modify plots:
#'   \describe{
#'     \item{\code{$ages.lab}}{Character vector, will change age labels in plots (default is \code{1:n.ages}).}
#'     \item{\code{$font.family}}{Font family, e.g. \code{"Times"}.}
#'   }
#'
#' Plot functions are located in \code{wham_plots_tables.R}
#' Table function is located in \code{par_tables_fn.R}
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param dir.main character, directory to save plots to (default = \code{getwd()})
#' @param out.type character, either \code{'html'}, \code{'pdf'}, or \code{'png'} (default = \code{'png'})
#' @param res resolution to save .png files (dpi)
#' @param plot.opts (optional) list of plot modifications
#'
#' @return NULL
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{\link{wham_html}}, \code{wham_plots_tables}
#'
#' @examples
#' \dontrun{
#' data("input4_SNEMAYT") # load fit wham model
#' mod <- fit_wham(input4_SNEMAYT)
#' plot_wham_output(mod)
#' }
plot_wham_output <- function(mod, dir.main = getwd(), out.type = 'png', res = 72, plot.opts = NULL){
  # if sdreport succeeded but didn't save full sdreport object in mod, recalculate it here
  if(mod$is_sdrep & class(mod$sdrep)[1] != "sdreport"){
    mod$sdrep <- TMB::sdreport(mod)
  }  

  # allow overwrite of default ages.lab = 1:n.ages
  if(!is.null(plot.opts)){
    if(!is.null(plot.opts[["ages.lab"]])) mod$ages.lab = plot.opts$ages.lab
    if(!is.null(plot.opts[["font.family"]])) fontfam = plot.opts$font.family
  } else {
    fontfam = ""
  }

  if(!out.type %in% c('html', 'pdf', 'png')){
    stop("out.type must be one of 'html', 'pdf', or 'png'. See ?plot_wham_output")
  }
  if(out.type %in% c("png", "pdf")) {
    table.type = "pdf"
  } else {
    table.type = "html"
  }

  # if(dir.exists(dir.main)){
  #   warning("Output directory already exists. Potentially overwriting previously saved output...")
  # }
  if(!dir.exists(dir.main)){
    dir.create(dir.main, showWarnings = FALSE)
    # stop("Output directory does not exist. Check 'dir.main' and try again.")
  }
  dir.res.tables <- file.path(dir.main, "res_tables")
  if(dir.exists(dir.res.tables)) {
    file.remove(file.path(dir.res.tables, list.files(dir.res.tables, recursive = TRUE)))
  }
  dir.create(dir.res.tables, showWarnings = FALSE)

  cat("Generating plot files... Please wait ~30 seconds...\n")
  if(out.type == "html") cat("html output works best with Google Chrome browser\n")
  graphics.off() # close any open windows

  if(out.type == 'pdf'){
    # PDF input_data -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main,"input_data.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    plot.catch.by.fleet(mod)
    for(i in 1:mod$env$data$n_fleets) plot.catch.age.comp.bubbles(mod, i=i)
    plot.index.input(mod)
    for(i in 1:mod$env$data$n_indices) plot.index.age.comp.bubbles(mod, i=i)
    plot.waa(mod,"ssb")
    plot.waa(mod,"jan1")
    plot.waa(mod,"totcatch")
    for(i in 1:mod$env$data$n_fleets) plot.waa(mod,"fleets", ind=i)
    for(i in 1:mod$env$data$n_indices) plot.waa(mod,"indices", ind=i)
    plot.maturity(mod)
    dev.off()

    # PDF diagnostics -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main,"diagnostics.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    fit.summary.text.plot.fn(mod)
    plot.ll.table.fn(mod)
    plot.catch.4.panel(mod)
    plot.index.4.panel(mod)
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov.diagnostic(mod)
    plot.NAA.4.panel(mod)
    plot.catch.age.comp(mod)
    plot.catch.age.comp.resids(mod)
    plot.index.age.comp(mod)
    plot.index.age.comp.resids(mod)
    plot.NAA.res(mod)
    if(!is.null(mod$osa)) plot.osa.residuals(mod)
    if(mod$is_sdrep) plot.all.stdresids.fn(mod)
    # plot.fleet.stdresids.fn(mod)
    # plot.index.stdresids.fn(mod)
    # if(!all(mod$env$data$Ecov_model == 0)) plot.ecov.stdresids.fn(mod)
    #plot.recruitment.devs(mod)
    dev.off()

    # PDF results -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main,"results.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    for(i in 1:mod$env$data$n_fleets) plot.fleet.sel.blocks(mod, use.i=i)
    for(i in 1:mod$env$data$n_indices) plot.index.sel.blocks(mod, use.i=i)
    if(mod$is_sdrep) plot.SSB.F.trend(mod)
    plot.SSB.AA(mod, prop=FALSE)
    plot.SSB.AA(mod, prop=TRUE)
    plot.NAA(mod, prop=FALSE)
    plot.NAA(mod, prop=TRUE)
    if(mod$env$data$recruit_model == 3 & mod$is_sdrep){ # these only work if Bev-Holt S-R was fit
      plot.SR.pred.line(mod)
    }
    if(mod$is_sdrep){
      plot.recr.ssb.yr(mod, loglog=FALSE)
      plot.recr.ssb.yr(mod, loglog=TRUE)
      plot.SARC.R.SSB(mod)
      plot.cv(mod)
    }
    plot.fleet.F(mod)
    plot.M(mod)
    plot.tile.age.year(mod, type="selAA")
    plot.tile.age.year(mod, type="MAA")
    plot_q_prior_post(mod) #flag inside to plot if prior is being used. 
    plot_q(mod)
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov(mod)
    dev.off()

    # PDF reference points -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main, "ref_points.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    plot.SPR.table(mod, plot=TRUE)
    plot.SPR.table(mod, plot=FALSE)
    plot.annual.SPR.targets(mod)
    if(mod$is_sdrep) plot.FXSPR.annual(mod)
    if(mod$env$data$recruit_model == 3 & mod$is_sdrep){ # these only work if Bev-Holt S-R was fit
      plot.MSY.annual(mod)
    }
    plot.yield.curves(mod, plot=TRUE)
    plot.yield.curves(mod, plot=FALSE)
    dev.off()

    # PDF retrospective -----------------
    if(!is.null(mod$peels)){
      grDevices::cairo_pdf(filename=file.path(dir.main, "retro.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
      plot.retro(mod, what = "SSB")
      plot.retro(mod, what = "Fbar")
      plot.retro(mod, what = "NAA")
      plot.retro(mod, what = "NAA_age", age=1)
      dev.off()
    }

    # PDF misc -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main, "misc.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    plot_catch_at_age_consistency(mod)
    plot_index_at_age_consistency(mod)
    plot_catch_curves_for_catch(mod)
    plot_catch_curves_for_index(mod)
    dev.off()

  } # end PDF section =============================================================

  if(out.type %in% c('png','html')){
    dir.plots <- file.path(dir.main, "plots_png")
    #This removes any plots from previous models. Some plots that are left in the directories may not pertain to the current model.
    if(dir.exists(dir.plots)) {
      file.remove(file.path(dir.plots, list.files(dir.plots, recursive = TRUE)))
    }
    dir.create(dir.plots, showWarnings = FALSE)

    # PNG input_data -----------------
    dir.data <- file.path(dir.plots, "input_data")
    dir.create(dir.data, showWarnings = FALSE)
    png(file.path(dir.data,"catch_by_fleet.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.catch.by.fleet(mod)
    dev.off()
    for(i in 1:mod$env$data$n_fleets){
      png(file.path(dir.data, paste0("catch_age_comp_fleet",i,".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.catch.age.comp.bubbles(mod, i=i)
      dev.off()
    }
    png(file.path(dir.data,"index.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.index.input(mod)
    dev.off()
    for(i in 1:mod$env$data$n_indices){
      png(file.path(dir.data, paste0("index",i,"_age_comp.png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.index.age.comp.bubbles(mod, i=i)
      dev.off()
    }
    png(file.path(dir.data,"weight_at_age_SSB.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.waa(mod,"ssb")
    dev.off()
    png(file.path(dir.data,"weight_at_age_Jan1.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.waa(mod,"jan1")
    dev.off()
    png(file.path(dir.data,"weight_at_age_catch.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.waa(mod,"totcatch")
    dev.off()
    for(i in 1:mod$env$data$n_fleets){
      png(file.path(dir.data, paste0("weight_at_age_fleet",i,".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.waa(mod,"fleets", ind=i)
      dev.off()
    }
    for(i in 1:mod$env$data$n_indices){
      png(file.path(dir.data, paste0("weight_at_age_index",i,".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.waa(mod,"indices", ind=i)
      dev.off()
    }
    png(file.path(dir.data,"maturity.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.maturity(mod)
    dev.off()

    # PNG diagnostics -----------------
    dir.diag <- file.path(dir.plots, "diagnostics")
    dir.create(dir.diag, showWarnings = FALSE)
    png(file.path(dir.diag,"summary_text.png"),width=10,height=10,units="in",res=res,family=fontfam)
    fit.summary.text.plot.fn(mod)
    dev.off()
    png(file.path(dir.diag,"likelihood.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.ll.table.fn(mod)
    dev.off()
    for(i in 1:mod$env$data$n_fleets){
      plot.catch.4.panel(mod, do.png = TRUE, fontfam=fontfam, use.i=i, od=dir.diag)
      plot.catch.age.comp(mod, do.png = TRUE, fontfam=fontfam, use.i=i, od=dir.diag)
      plot.catch.age.comp.resids(mod, do.png = TRUE, fontfam=fontfam, use.i=i, od=dir.diag)
    }
    for(i in 1:mod$env$data$n_indices){
      plot.index.4.panel(mod, do.png = TRUE, fontfam=fontfam, use.i=i, od=dir.diag)
      plot.index.age.comp(mod, do.png = TRUE, fontfam=fontfam, use.i=i, od=dir.diag)
      plot.index.age.comp.resids(mod, do.png = TRUE, fontfam=fontfam, use.i=i, od=dir.diag)
    }
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov.diagnostic(mod, do.png = TRUE, fontfam=fontfam, od=dir.diag)
    plot.NAA.4.panel(mod, do.png = TRUE, fontfam=fontfam, od=dir.diag)
    plot.NAA.res(mod, do.png = TRUE, fontfam=fontfam, od=dir.diag)
    if(!is.null(mod$osa)) plot.osa.residuals(mod, do.png=TRUE, fontfam=fontfam, res=res, od=dir.diag)
    if(mod$is_sdrep) plot.all.stdresids.fn(mod, do.png=TRUE, fontfam=fontfam, res=res, od=dir.diag)
    # plot.fleet.stdresids.fn(mod, do.png=TRUE, fontfam=fontfam, res=res,family=fontfam, od=dir.diag)
    # plot.index.stdresids.fn(mod, do.png=TRUE, fontfam=fontfam, res=res,family=fontfam, od=dir.diag)
    # if(!all(mod$env$data$Ecov_model == 0)) plot.ecov.stdresids.fn(mod, do.png=TRUE, fontfam=fontfam, res=res,family=fontfam, od=dir.diag)
    #plot.recruitment.devs(mod)

    # PNG results -----------------
    dir.res <- file.path(dir.plots, "results")
    dir.create(dir.res, showWarnings = FALSE)
    for(i in 1:mod$env$data$n_fleets){
      plot.fleet.sel.blocks(mod, ages.lab = mod$ages.lab, do.png=TRUE, fontfam=fontfam, use.i=i, od=dir.res)
    }
    for(i in 1:mod$env$data$n_indices){
      plot.index.sel.blocks(mod, ages.lab = mod$ages.lab, do.png=TRUE, fontfam=fontfam, use.i=i, od=dir.res)
    }
    if(mod$is_sdrep){
      png(file.path(dir.res,"SSB_F_trend.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SSB.F.trend(mod)
      dev.off()
    }
    png(file.path(dir.res,"SSB_at_age.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.SSB.AA(mod, prop=FALSE)
    dev.off()
    png(file.path(dir.res,"SSB_at_age_proportion.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.SSB.AA(mod, prop=TRUE)
    dev.off()
    png(file.path(dir.res,"Numbers_at_age.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.NAA(mod, prop=FALSE)
    dev.off()
    png(file.path(dir.res,"Numbers_at_age_proportion.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.NAA(mod, prop=TRUE)
    dev.off()
    if(mod$env$data$recruit_model == 3 & mod$is_sdrep){ # these only work if Bev-Holt S-R was fit
      png(file.path(dir.res,"SSB_Rec_fit.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SR.pred.line(mod)
      dev.off()
    }
    if(mod$is_sdrep){
      png(file.path(dir.res,"SSB_Rec.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.recr.ssb.yr(mod, loglog=FALSE)
      dev.off()
      png(file.path(dir.res,"SSB_Rec_loglog.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.recr.ssb.yr(mod, loglog=TRUE)
      dev.off()
      png(file.path(dir.res,"SSB_Rec_time.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SARC.R.SSB(mod)
      dev.off()
      png(file.path(dir.res,"CV_SSB_Rec_F.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.cv(mod)
      dev.off()
    }
    png(file.path(dir.res,"F_byfleet.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.fleet.F(mod)
    dev.off()
    png(file.path(dir.res,"M_at_age.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.M(mod)
    dev.off()
    plot.tile.age.year(mod, type="selAA", do.png=TRUE, fontfam=fontfam, od=dir.res)
    plot.tile.age.year(mod, type="MAA", do.png=TRUE, fontfam=fontfam, od=dir.res)
    plot_q_prior_post(mod, do.png=TRUE, fontfam=fontfam, od=dir.res) #flag inside to plot if prior is being used. 
    plot_q(mod, do.png=TRUE, fontfam=fontfam, od=dir.res)
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov(mod, do.png=TRUE, fontfam=fontfam, od=dir.res, res=res)

    # PNG reference points -----------------
    dir.refpts <- file.path(dir.plots, "ref_points")
    dir.create(dir.refpts, showWarnings = FALSE)

    png(file.path(dir.refpts,"SPR_targets_ave_plot.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.SPR.table(mod, plot=TRUE)
    dev.off()
    png(file.path(dir.refpts,"SPR_targets_ave_table.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.SPR.table(mod, plot=FALSE)
    dev.off()
    plot.annual.SPR.targets(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam)
    if(mod$is_sdrep) plot.FXSPR.annual(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam)
    plot.yield.curves(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam, plot=TRUE)
    plot.yield.curves(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam, plot=FALSE)
    if(mod$env$data$recruit_model == 3 & mod$is_sdrep) plot.MSY.annual(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam)

    # PNG retrospective -----------------
    if(!is.null(mod$peels)){
      dir.retro <- file.path(dir.plots, "retro")
      dir.create(dir.retro, showWarnings = FALSE)
      plot.retro(mod, what = "SSB", od=dir.retro, do.png=TRUE, fontfam=fontfam)
      plot.retro(mod, what = "Fbar", od=dir.retro, do.png=TRUE, fontfam=fontfam)
      plot.retro(mod, what = "NAA", od=dir.retro, do.png=TRUE, fontfam=fontfam)
      plot.retro(mod, what = "NAA_age", age=1, od=dir.retro, do.png=TRUE, fontfam=fontfam)
      dev.off()
    }

    # PNG misc -----------------
    dir.misc <- file.path(dir.plots, "misc")
    dir.create(dir.misc, showWarnings = FALSE)
    plot_catch_at_age_consistency(mod, od=dir.misc, do.png=TRUE, fontfam=fontfam)
    plot_index_at_age_consistency(mod, od=dir.misc, do.png=TRUE, fontfam=fontfam)
    plot_catch_curves_for_catch(mod, od=dir.misc, do.png=TRUE, fontfam=fontfam)
    plot_catch_curves_for_index(mod, od=dir.misc, do.png=TRUE, fontfam=fontfam)
    dev.off()
  } # end PNG section =====================================================


  # uses png output, automatically opens in browser
  if(out.type == 'html'){
    wham_html(dir.main = dir.main)
  }
  if(table.type == "pdf"){
    cat(paste0("Making LaTeX/pdf tables: ", file.path(dir.res.tables, "wham_par_tables.pdf"), "\n"))
    par_tables_fn(mod, od = dir.res.tables, do.tex = T)
  }
  if(table.type == "html"){
    cat("Making HTML tables.\n")
    par_tables_fn(mod, od = dir.res.tables, do.html = T)
    cat("Opening HTML tables in your default web-browser.\n")
    browseURL(file.path(dir.res.tables, "wham_par_tables.html"))
  }
}


