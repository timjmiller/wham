#' Make WHAM plots
#'
#' Generates many output plots for a fit WHAM model.
#'
#' \code{out.type = 'pdf'} makes one pdf file of all plots (default). In the future,
#' \code{out.type = 'html'} will be implemented, making an html file for viewing in a
#' browser (tabs: 'input data', 'diagnostics', 'results', 'ref_points', 'retro',
#' and 'misc'). \code{out.type = 'png'} creates a subdirectory 'plots_png' in \code{dir} and
#' saves .png files within.
#'
#' Plot functions are located in \code{wham_plots_tables.R}
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param dir.main character, directory to save plots to (default = \code{getwd()})
#' @param out.type character, either 'html', 'pdf', or 'png' (only 'pdf' currently supported)
#'
#' @return NULL
#'
#' @export
#'
#' @seealso \code{\link{fit_wham}}, \code{wham_plots_tables}
#'
#' @examples
#' \dontrun{
#' data("SNEMA_ytl") # load SNEMA yellowtail flounder data and parameter settings
#' mod = fit_wham(input)
#' make_wham_plots(mod)
#' }
make_wham_plots <- function(mod, dir.main = getwd(), out.type = 'pdf'){
  if(!out.type %in% c('html', 'pdf', 'png')){
    stop("out.type must be one of 'html', 'pdf', or 'png'. See ?make_wham_plots")
  }
  if(!dir.exists(dir.main)){
    stop("Output directory does not exist. Check 'dir' and try again.")
  }
  # if(out.type != 'pdf'){
  #   stop("Only pdf output is currently supported. Stay tuned for 'html' or 'png' plots")
  # }

  graphics.off()     # close any open windows
  origpar = par()

  if(out.type == 'pdf'){
    # PDF input_data
    grDevices::cairo_pdf(filename=file.path(dir.main,"input_data.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
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

    # PDF diagnostics
    grDevices::cairo_pdf(filename=file.path(dir.main,"diagnostics.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    fit.summary.text.plot.fn(mod)
    plot.ll.table.fn(mod)
    plot.catch.4.panel(mod)
    plot.index.4.panel(mod)
    plot.NAA.4.panel(mod)
    plot.catch.age.comp(mod)
    plot.catch.age.comp.resids(mod)
    plot.index.age.comp(mod)
    plot.index.age.comp.resids(mod)
    plot.NAA.res(mod)
    #plot.recruitment.devs(mod)
    dev.off()

    # PDF results
    grDevices::cairo_pdf(filename=file.path(dir.main,"results.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    for(i in 1:mod$env$data$n_fleets) plot.fleet.sel.blocks(mod, use.i=i)
    for(i in 1:mod$env$data$n_indices) plot.index.sel.blocks(mod, use.i=i)
    plot.SSB.F.trend(mod)
    plot.SSB.AA(mod, prop=FALSE)
    plot.SSB.AA(mod, prop=TRUE)
    plot.NAA(mod, prop=FALSE)
    plot.NAA(mod, prop=TRUE)
    plot.recr.ssb.yr(mod, loglog=FALSE)
    plot.recr.ssb.yr(mod, loglog=TRUE)
    plot.SARC.R.SSB(mod)
    plot.fleet.F(mod)
    plot.cv(mod)
    plot.M(mod)
    dev.off()

    # PDF reference points
    grDevices::cairo_pdf(filename=file.path(dir.main, "ref_points.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.SPR.table(mod, plot=TRUE)
    plot.SPR.table(mod, plot=FALSE)
    plot.annual.SPR.targets(mod)
    plot.FXSPR.annual(mod)
    plot.SR.pred.line(mod)
    plot.MSY.annual(mod)
    plot.yield.curves(mod, plot=TRUE)
    plot.yield.curves(mod, plot=FALSE)
    dev.off()

    # PDF retrospective
    grDevices::cairo_pdf(filename=file.path(dir.main, "retro.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.retro(mod, what = "SSB")
    plot.retro(mod, what = "Fbar")
    plot.retro(mod, what = "NAA")
    dev.off()

    # PDF misc
    grDevices::cairo_pdf(filename=file.path(dir.main, "misc.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot_catch_at_age_consistency(mod)
    plot_index_at_age_consistency(mod)
    plot_catch_curves_for_catch(mod)
    plot_catch_curves_for_index(mod)
    dev.off()
  } # end PDF section

  if(out.type %in% c('png','html')){
    dir.plots <- file.path(dir.main, "plots_png")
    dir.create(dir.plots)

    # PNG input_data
    dir.data <- file.path(dir.plots, "input_data")
    dir.create(dir.data)
    png(file.path(dir.data,"catch_by_fleet.png"),width=10,height=10,units="in",res=300)
    plot.catch.by.fleet(mod)
    dev.off()
    for(i in 1:mod$env$data$n_fleets){
      png(file.path(dir.data, paste0("catch_age_comp_fleet",i,".png")),width=10,height=10,units="in",res=300)
      plot.catch.age.comp.bubbles(mod, i=i)
      dev.off()
    }
    png(file.path(dir.data,"index.png"),width=10,height=10,units="in",res=300)
    plot.index.input(mod)
    dev.off()
    for(i in 1:mod$env$data$n_indices){
      png(file.path(dir.data, paste0("index",i,"_age_comp.png")),width=10,height=10,units="in",res=300)
      plot.index.age.comp.bubbles(mod, i=i)
      dev.off()
    }
    png(file.path(dir.data,"weight_at_age_SSB.png"),width=10,height=10,units="in",res=300)
    plot.waa(mod,"ssb")
    dev.off()
    png(file.path(dir.data,"weight_at_age_Jan1.png"),width=10,height=10,units="in",res=300)
    plot.waa(mod,"jan1")
    dev.off()
    png(file.path(dir.data,"weight_at_age_catch.png"),width=10,height=10,units="in",res=300)
    plot.waa(mod,"totcatch")
    dev.off()
    for(i in 1:mod$env$data$n_fleets){
      png(file.path(dir.data, paste0("weight_at_age_fleet",i,".png")),width=10,height=10,units="in",res=300)
      plot.waa(mod,"fleets", ind=i)
      dev.off()
    }
    for(i in 1:mod$env$data$n_indices){
      png(file.path(dir.data, paste0("weight_at_age_index",i,".png")),width=10,height=10,units="in",res=300)
      plot.waa(mod,"indices", ind=i)
      dev.off()
    }
    png(file.path(dir.data,"maturity.png"),width=10,height=10,units="in",res=300)
    plot.maturity(mod)
    dev.off()

    # PNG diagnostics
    dir.diag <- file.path(dir.plots, "diagnostics")
    dir.create(dir.diag)
    png(file.path(dir.diag,"summary_text.png"),width=10,height=10,units="in",res=300)
    fit.summary.text.plot.fn(mod)
    dev.off()
    png(file.path(dir.diag,"likelihood.png"),width=10,height=10,units="in",res=300)
    plot.ll.table.fn(mod)
    dev.off()
    for(i in 1:mod$env$data$n_fleets){
      plot.catch.4.panel(mod, do.png = TRUE, use.i=i, od=dir.diag)
      plot.catch.age.comp(mod, do.png = TRUE, use.i=i, od=dir.diag)
      plot.catch.age.comp.resids(mod, do.png = TRUE, use.i=i, od=dir.diag)
    }
    for(i in 1:mod$env$data$n_indices){
      plot.index.4.panel(mod, do.png = TRUE, use.i=i, od=dir.diag)
      plot.index.age.comp(mod, do.png = TRUE, use.i=i, od=dir.diag)
      plot.index.age.comp.resids(mod, do.png = TRUE, use.i=i, od=dir.diag)
    }
    plot.NAA.4.panel(mod, do.png = TRUE, od=dir.diag)
    plot.NAA.res(mod, do.png = TRUE, od=dir.diag)
    #plot.recruitment.devs(mod)

    # PNG results
    dir.res <- file.path(dir.plots, "results")
    dir.create(dir.res)
    for(i in 1:mod$env$data$n_fleets){
      plot.fleet.sel.blocks(mod, do.png=TRUE, use.i=i, od=dir.res)
    }
    for(i in 1:mod$env$data$n_indices){
      plot.index.sel.blocks(mod, do.png=TRUE, use.i=i, od=dir.res)
    }
    png(file.path(dir.res,"SSB_F_trend.png"),width=10,height=10,units="in",res=300)
    plot.SSB.F.trend(mod)
    dev.off()
    png(file.path(dir.res,"SSB_at_age.png"),width=10,height=10,units="in",res=300)
    plot.SSB.AA(mod, prop=FALSE)
    dev.off()
    png(file.path(dir.res,"SSB_at_age_proportion.png"),width=10,height=10,units="in",res=300)
    plot.SSB.AA(mod, prop=TRUE)
    dev.off()
    png(file.path(dir.res,"Numbers_at_age.png"),width=10,height=10,units="in",res=300)
    plot.NAA(mod, prop=FALSE)
    dev.off()
    png(file.path(dir.res,"Numbers_at_age_proportion.png"),width=10,height=10,units="in",res=300)
    plot.NAA(mod, prop=TRUE)
    dev.off()
    png(file.path(dir.res,"SSB_Rec.png"),width=10,height=10,units="in",res=300)
    plot.recr.ssb.yr(mod, loglog=FALSE)
    dev.off()
    png(file.path(dir.res,"SSB_Rec_loglog.png"),width=10,height=10,units="in",res=300)
    plot.recr.ssb.yr(mod, loglog=TRUE)
    dev.off()
    png(file.path(dir.res,"SSB_Rec_time.png"),width=10,height=10,units="in",res=300)
    plot.SARC.R.SSB(mod)
    dev.off()
    png(file.path(dir.res,"F_byfleet.png"),width=10,height=10,units="in",res=300)
    plot.fleet.F(mod)
    dev.off()
    png(file.path(dir.res,"CV_SSB_Rec_F.png"),width=10,height=10,units="in",res=300)
    plot.cv(mod)
    dev.off()
    png(file.path(dir.res,"M_at_age.png"),width=10,height=10,units="in",res=300)
    plot.M(mod)
    dev.off()

    # PNG reference points
    dir.refpts <- file.path(dir.plots, "ref_points")
    dir.create(dir.refpts)

    png(file.path(dir.refpts,"SPR_targets_ave_plot.png"),width=10,height=10,units="in",res=300)
    plot.SPR.table(mod, plot=TRUE)
    dev.off()
    png(file.path(dir.refpts,"SPR_targets_ave_table.png"),width=10,height=10,units="in",res=300)
    plot.SPR.table(mod, plot=FALSE)
    dev.off()
    plot.annual.SPR.targets(mod, od=dir.refpts, do.png=TRUE)
    plot.FXSPR.annual(mod, od=dir.refpts, do.png=TRUE)
    plot.yield.curves(mod, od=dir.refpts, do.png=TRUE, plot=TRUE)
    plot.yield.curves(mod, od=dir.refpts, do.png=TRUE, plot=FALSE)
    if(mod$env$data$recruit_model == 3){ # these only work if Bev-Holt S-R was fit
      png(file.path(dir.refpts,"SSB_Rec_fit.png"),width=10,height=10,units="in",res=300)
      plot.SR.pred.line(mod)
      dev.off()
      png(file.path(dir.refpts,"MSY_annual.png"),width=10,height=10,units="in",res=300)
      plot.MSY.annual(mod)
      dev.off()
    }

    # PNG retrospective
    dir.retro <- file.path(dir.plots, "retro")
    dir.create(dir.retro)
    plot.retro(mod, what = "SSB", od=dir.retro, do.png=TRUE)
    plot.retro(mod, what = "Fbar", od=dir.retro, do.png=TRUE)
    plot.retro(mod, what = "NAA", od=dir.retro, do.png=TRUE)
    dev.off()

    # PNG misc
    dir.misc <- file.path(dir.plots, "misc")
    dir.create(dir.misc)
    plot_catch_at_age_consistency(mod, od=dir.misc, do.png=TRUE)
    plot_index_at_age_consistency(mod, od=dir.misc, do.png=TRUE)
    plot_catch_curves_for_catch(mod, od=dir.misc, do.png=TRUE)
    plot_catch_curves_for_index(mod, od=dir.misc, do.png=TRUE)
    dev.off()
  } # end PNG section
}


