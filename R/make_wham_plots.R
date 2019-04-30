#' Make WHAM plots
#'
#' Generates many output plots for a fit WHAM model.
#'
#' \code{out.type = 'html'} is the default, makes an html file viewable in any
#' browser (tabs: 'input data', 'diagnostics', 'results', 'SPR / MSY', 'retro',
#' and 'misc'). \code{out.type = 'pdf'} makes one pdf file of all plots.
#' \code{out.type = 'png'} creates a subdirectory 'plots_png' in \code{dir} and
#' saves .png files within.
#'
#' Plot functions are located in \code{wham_plots_tables.R}
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param dir character, directory to save plots to (default = \code{getwd()})
#' @param out.type character, either 'html', 'pdf', or 'png' (default = 'html')
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
make_wham_plots <- function(mod, dir = getwd(), out.type = 'html'){
  if(!out.type %in% c('html', 'pdf', 'png')){
    stop("out.type must be one of 'html', 'pdf', or 'png'. See ?make_wham_plots")
  }
  if(!dir.exists(dir)){
    stop("Output directory does not exist. Check 'dir' and try again.")
  }

  graphics.off()     # close any open windows
  origpar = par()

  if(out.type == 'pdf'){
    # PDF input_data
    grDevices::cairo_pdf(filename=file.path(dir,"input_data.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.catch.by.fleet(mod)
    plot.catch.age.comp.bubbles(mod)
    plot.index.input(mod)
    plot.index.age.comp.bubbles(mod)
    plot.waa(mod,"ssb")
    plot.waa(mod,"jan1")
    plot.waa(mod,"totcatch")
    plot.waa(mod,"fleets")
    plot.waa(mod,"indices")
    plot.maturity(mod)
    dev.off()

    # PDF diagnostics
    grDevices::cairo_pdf(filename=file.path(dir,"diagnostics.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
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
    grDevices::cairo_pdf(filename=file.path(dir,"results.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.fleet.sel.blocks(mod)
    plot.index.sel.blocks(mod)
    plot.SSB.F.trend(mod)
    plot.SSB.AA(mod)
    plot.NAA(mod)
    plot.recr.ssb.yr(mod)
    plot.SARC.R.SSB(mod)
    plot.fleet.F(mod)
    plot.cv(mod)
    plot.M(mod)
    dev.off()

    # PDF SPR & MSY
    grDevices::cairo_pdf(filename=file.path(dir, "spr_msy.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.SPR.table(mod)
    plot.annual.SPR.targets(mod)
    plot.FXSPR.annual(mod)
    plot.SR.pred.line(mod)
    plot.MSY.annual(mod)
    plot.yield.curves(mod)
    dev.off()

    # PDF retrospective
    grDevices::cairo_pdf(filename=file.path(dir, "retro.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.retro(mod, what = "SSB")
    plot.retro(mod, what = "Fbar")
    plot.retro(mod, what = "NAA")
    dev.off()

    # PDF misc
    grDevices::cairo_pdf(filename=file.path(dir, "misc.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot_catch_at_age_consistency(mod)
    plot_index_at_age_consistency(mod)
    plot_catch_curves_for_catch(mod)
    plot_catch_curves_for_index(mod)
    dev.off()
  } # end PDF section

  if(out.type == 'png'){
    dir.plots <- file.path(dir, "plots_png")
    dir.create(dir.plots)

    # PNG input_data
    dir.data <- file.path(dir.plots, "input_data")
    dir.create(dir.data)
    grDevices::cairo_PNG(filename=file.path(dir,"input_data.PNG"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.catch.by.fleet(mod)
    plot.catch.age.comp.bubbles(mod)
    plot.index.input(mod)
    plot.index.age.comp.bubbles(mod)
    plot.waa(mod,"ssb")
    plot.waa(mod,"jan1")
    plot.waa(mod,"totcatch")
    plot.waa(mod,"fleets")
    plot.waa(mod,"indices")
    plot.maturity(mod)
    dev.off()

    # PNG diagnostics
    grDevices::cairo_PNG(filename=file.path(dir,"diagnostics.PNG"), family = "Times", height = 10, width = 10, onefile = TRUE)
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

    # PNG results
    grDevices::cairo_PNG(filename=file.path(dir,"results.PNG"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.fleet.sel.blocks(mod)
    plot.index.sel.blocks(mod)
    plot.SSB.F.trend(mod)
    plot.SSB.AA(mod)
    plot.NAA(mod)
    plot.recr.ssb.yr(mod)
    plot.SARC.R.SSB(mod)
    plot.fleet.F(mod)
    plot.cv(mod)
    plot.M(mod)
    dev.off()

    # PNG SPR & MSY
    grDevices::cairo_PNG(filename=file.path(dir, "spr_msy.PNG"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.SPR.table(mod)
    plot.annual.SPR.targets(mod)
    plot.FXSPR.annual(mod)
    plot.SR.pred.line(mod)
    plot.MSY.annual(mod)
    plot.yield.curves(mod)
    dev.off()

    # PNG retrospective
    grDevices::cairo_PNG(filename=file.path(dir, "retro.PNG"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot.retro(mod, what = "SSB")
    plot.retro(mod, what = "Fbar")
    plot.retro(mod, what = "NAA")
    dev.off()

    # PNG misc
    grDevices::cairo_PNG(filename=file.path(dir, "misc.PNG"), family = "Times", height = 10, width = 10, onefile = TRUE)
    plot_catch_at_age_consistency(mod)
    plot_index_at_age_consistency(mod)
    plot_catch_curves_for_catch(mod)
    plot_catch_curves_for_index(mod)
    dev.off()
  } # end PNG section
}


