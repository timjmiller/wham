#' Plot WHAM output
#'
#' Generates many output plots and tables for a fit WHAM model.
#'
#' \code{out.type = 'html'} (default) makes a html file for viewing plot .png files and html tables of parameter estimates in a browser.
#' \code{out.type = 'pdf'} makes one pdf file of all plots and tables. 
#' \code{out.type = 'png'} creates a subdirectory `plots_png`` in \code{dir.main} and saves .png files within.
#' \code{out.type = 'pdf' or 'png'} makes LaTeX and pdf files of tables of parameter estimates.
#' (tabs: 'input data', 'diagnostics', 'results', 'ref_points', 'retro', and 'misc').
#' 
#' \code{plot.opts} holds optional arguments to modify plots:
#'   \describe{
#'     \item{\code{$ages.lab}}{Character vector, will change age labels in plots (default is \code{1:n.ages}).}
#'     \item{\code{$font.family}}{Font family, e.g. \code{"Times"}.}
#'     \item{\code{$browse}}{T/F whether to open the html file in a browser. Default = T.}
#'   }
#'
#' Plot functions are located in \code{wham_plots_tables.R}
#' Table function is located in \code{par_tables_fn.R}
#'
#' @param mod output from \code{\link{fit_wham}}
#' @param dir.main character, directory to save plots to (default = \code{getwd()})
#' @param out.type character, either \code{'html'}, \code{'pdf'}, or \code{'png'} (default = \code{'html'})
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
plot_wham_output <- function(mod, dir.main = getwd(), out.type = 'html', res = 72, plot.opts = NULL){
  # if sdreport succeeded but didn't save full sdreport object in mod, recalculate it here
  if(mod$is_sdrep & class(mod$sdrep)[1] != "sdreport"){
    mod$sdrep <- TMB::sdreport(mod)
  }  
  fslash <- function(fp) chartr('\\','/',fp)
  dir.main = fslash(dir.main)

# allow overwrite of default ages.lab = 1:n.ages
  fontfam = ""
  browse = TRUE
  if(!is.null(plot.opts)){
    if(!is.null(plot.opts[["ages.lab"]])) mod$ages.lab = plot.opts$ages.lab
    if(!is.null(plot.opts[["font.family"]])) fontfam = plot.opts$font.family
    if(!is.null(plot.opts[["browse"]])) browse = plot.opts$browse
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
    for(i in 1:mod$env$data$n_stocks) plot.waa(mod,"ssb", ind = i)
    #plot.waa(mod,"jan1")
    #plot.waa(mod,"totcatch")
    for(i in 1:mod$env$data$n_fleets) plot.waa(mod,"fleets", ind=i)
    for(i in 1:mod$env$data$n_indices) plot.waa(mod,"indices", ind=i)
    for(i in 1:mod$env$data$n_stocks) plot.maturity(mod, stock = i)
    dev.off()

    # PDF diagnostics -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main,"diagnostics.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    fit.summary.text.plot.fn(mod)
    plot.ll.table.fn(mod)
    for(i in 1:mod$env$data$n_fleets){
      plot.catch.4.panel(mod, use.i = i)
      plot.catch.age.comp(mod, use.i = i)
      plot.catch.age.comp.resids(mod, use.i = i)
    }
    for(i in 1:mod$env$data$n_indices){
      plot.index.4.panel(mod, use.i = i)
      plot.index.age.comp(mod, use.i = i)
      plot.index.age.comp.resids(mod, use.i = i)
    }
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov.diagnostic(mod)
    plot.NAA.4.panel(mod)
    #plot.NAA.res(mod)
    if(!is.null(mod$osa)) {
      plot.catch.age.comp.resids(mod, osa = TRUE)
      plot.index.age.comp.resids(mod, osa = TRUE)
      plot.osa.residuals(mod)
    }
    if(mod$is_sdrep) plot.all.stdresids.fn(mod)
    #plot.recruitment.devs(mod)
    dev.off()

    # PDF results -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main,"results.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    for(i in 1:mod$env$data$n_fleets) plot.sel.blocks(mod, ages.lab = mod$ages.lab, use.i=i, indices = FALSE)
    for(i in 1:mod$env$data$n_indices) plot.sel.blocks(mod, ages.lab = mod$ages.lab, use.i=i, indices = TRUE)
    if(mod$is_sdrep) plot.SSB.F.trend(mod)
    for(i in 1:mod$env$data$n_stocks) {
      plot.SSB.AA(mod, prop=FALSE, stock = i)
      plot.SSB.AA(mod, prop=TRUE, stock = i)
      for(r in 1:mod$env$data$n_regions) {
        plot.NAA(mod, prop=FALSE, stock = i, region = r)
        plot.NAA(mod, prop=TRUE, stock = i, region = r)
      }
    }
    if(mod$is_sdrep){
      for(i in 1:mod$env$data$n_stocks) {
        if(mod$env$data$recruit_model[i] %in% (3:4)) {
          plot.SR.pred.line(mod, stock = i)
        }
        plot.recr.ssb.yr(mod, loglog=FALSE, stock = i)
        plot.recr.ssb.yr(mod, loglog=TRUE, stock = i)
        plot.SARC.R.SSB(mod, stock = i)
      }
      if(mod$env$data$n_stocks>1) plot.SARC.R.SSB(mod, stock = NULL) #totals
      plot.cv(mod)
    }
    plot.fleet.F(mod)
    for(i in 1:mod$env$data$n_stocks) for(r in 1:mod$env$data$n_regions) {
      plot.M(mod, stock = i, region = r)
    }
    plot.tile.age.year(mod, type="selAA")
    plot.tile.age.year(mod, type="MAA")
    plot.tile.age.year(mod, type="NAA_devs")
    plot_q_prior_post(mod) #flag inside to plot if prior is being used. 
    plot_q(mod)
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov(mod)
    dev.off()

    # PDF reference points -----------------
    grDevices::cairo_pdf(filename=file.path(dir.main, "ref_points.pdf"), family = fontfam, height = 10, width = 10, onefile = TRUE)
    if(mod$env$data$n_stocks == 1 & mod$env$data$n_regions == 1){
      plot.SPR.table(mod, plot=TRUE)
      plot.SPR.table(mod, plot=FALSE)
      plot.annual.SPR.targets(mod)
      plot.yield.curves(mod, plot=TRUE)
      plot.yield.curves(mod, plot=FALSE)
    }
    if(mod$env$data$do_SPR_BRPs){
      if(mod$is_sdrep) plot.FXSPR.annual(mod)
    }
    if(mod$env$data$do_MSY_BRPs){
      if(mod$env$data$recruit_model %in% (3:4) & mod$is_sdrep){ # these only work if Bev-Holt S-R was fit
        plot.MSY.annual(mod)
      }
    }
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
  #take out any spaces
  sfns <- chartr(" ", "_", mod$input$stock_names)
  rfns <- chartr(" ", "_", mod$input$region_names)
  ffns <- chartr(" ", "_", mod$input$fleet_names)
  ifns <- chartr(" ", "_", mod$input$index_names)

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
      png(file.path(dir.data, paste0("catch_age_comp_", ffns[i],"_", rfns[mod$input$data$fleet_regions[i]], ".png")),
        width=10,height=10,units="in",res=res,family=fontfam)
      plot.catch.age.comp.bubbles(mod, i=i)
      dev.off()
    }
    png(file.path(dir.data,"index.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.index.input(mod)
    dev.off()
    for(i in 1:mod$env$data$n_indices){
      png(file.path(dir.data, paste0(ifns[i], "_", rfns[mod$input$data$index_regions[i]], "_age_comp.png")),
        width=10,height=10,units="in",res=res,family=fontfam)
      plot.index.age.comp.bubbles(mod, i=i)
      dev.off()
    }
    for(i in 1:mod$env$data$n_stocks) {
      png(file.path(dir.data,paste0("weight_at_age_SSB_", sfns[i],".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.waa(mod,"ssb", ind = i)
      dev.off()
    }
    # png(file.path(dir.data,"weight_at_age_Jan1.png"),width=10,height=10,units="in",res=res,family=fontfam)
    # plot.waa(mod,"jan1")
    # dev.off()
    # png(file.path(dir.data,"weight_at_age_catch.png"),width=10,height=10,units="in",res=res,family=fontfam)
    # plot.waa(mod,"totcatch")
    # dev.off()
    for(i in 1:mod$env$data$n_fleets){
      png(file.path(dir.data, paste0("weight_at_age_", ffns[i],"_fleet.png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.waa(mod,"fleets", ind=i)
      dev.off()
    }
    for(i in 1:mod$env$data$n_indices){
      png(file.path(dir.data, paste0("weight_at_age_",ifns[i],"_index.png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.waa(mod,"indices", ind=i)
      dev.off()
    }
    for(i in 1:mod$env$data$n_stocks) {
      png(file.path(dir.data,paste0("maturity_", sfns[i], ".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.maturity(mod, stock = i)
      dev.off()
    }

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
    #plot.NAA.res(mod, do.png = TRUE, fontfam=fontfam, od=dir.diag)
    if(!is.null(mod$osa)) {
      plot.catch.age.comp.resids(mod, osa = TRUE, do.png=TRUE, fontfam=fontfam, res=res, od=dir.diag)
      plot.index.age.comp.resids(mod, osa = TRUE, do.png=TRUE, fontfam=fontfam, res=res, od=dir.diag)
      plot.osa.residuals(mod, do.png=TRUE, fontfam=fontfam, res=res, od=dir.diag)
    }
    if(mod$is_sdrep) plot.all.stdresids.fn(mod, do.png=TRUE, fontfam=fontfam, res=res, od=dir.diag)

    #plot.recruitment.devs(mod)

    # PNG results -----------------
    dir.res <- file.path(dir.plots, "results")
    dir.create(dir.res, showWarnings = FALSE)
    for(i in 1:mod$env$data$n_fleets){
      plot.sel.blocks(mod, ages.lab = mod$ages.lab, indices = FALSE, do.png=TRUE, fontfam=fontfam, use.i=i, od=dir.res)
      #plot.fleet.sel.blocks(mod, ages.lab = mod$ages.lab, do.png=TRUE, fontfam=fontfam, use.i=i, od=dir.res)
    }
    for(i in 1:mod$env$data$n_indices){
      plot.sel.blocks(mod, ages.lab = mod$ages.lab, indices = TRUE, do.png=TRUE, fontfam=fontfam, use.i=i, od=dir.res)
      #plot.index.sel.blocks(mod, ages.lab = mod$ages.lab, do.png=TRUE, fontfam=fontfam, use.i=i, od=dir.res)
    }
    if(mod$is_sdrep){
      png(file.path(dir.res,"SSB_F_trend.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SSB.F.trend(mod)
      dev.off()
    }
    for(i in 1:mod$env$data$n_stocks) {
      png(file.path(dir.res,paste0("SSB_at_age_",sfns[i],".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SSB.AA(mod, prop=FALSE, stock = i)
      dev.off()
      png(file.path(dir.res,paste0("SSB_at_age_proportion_",sfns[i],".png")),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SSB.AA(mod, prop=TRUE, stock = i)
      dev.off()
      for(r in 1:mod$env$data$n_regions) {
        png(file.path(dir.res,paste0("Numbers_at_age_", sfns[i],"_", rfns[r],".png")),
          width=10,height=10,units="in",res=res,family=fontfam)
        plot.NAA(mod, prop=FALSE, stock = i, region = r)
        dev.off()
        png(file.path(dir.res,paste0("Numbers_at_age_proportion_", sfns[i],"_", rfns[r],".png")),
          width=10,height=10,units="in",res=res,family=fontfam)
        plot.NAA(mod, prop=TRUE, stock = i, region = r)
        dev.off()
      }
    }
    if(mod$is_sdrep){
      for(i in 1:mod$env$data$n_stocks) {
        if(mod$env$data$recruit_model[i] %in% (3:4)) {
          png(file.path(dir.res,paste0("SSB_Rec_", sfns[i],"_fit.png")),width=10,height=10,units="in",res=res,family=fontfam)
          plot.SR.pred.line(mod, stock = i)
          dev.off()
        }
        png(file.path(dir.res,paste0("SSB_Rec_", sfns[i],".png")),width=10,height=10,units="in",res=res,family=fontfam)
        plot.recr.ssb.yr(mod, loglog=FALSE, stock = i)
        dev.off()
        png(file.path(dir.res,paste0("SSB_Rec_loglog_", sfns[i],".png")),width=10,height=10,units="in",res=res,family=fontfam)
        plot.recr.ssb.yr(mod, loglog=TRUE, stock = i)
        dev.off()
        png(file.path(dir.res,paste0("SSB_Rec_time_", sfns[i],".png")),width=10,height=10,units="in",res=res,family=fontfam)
        plot.SARC.R.SSB(mod, stock = i)
        dev.off()
      }
      if(mod$env$data$n_stocks>1) { 
        png(file.path(dir.res,paste0("SSB_Rec_time_total.png")),width=10,height=10,units="in",res=res,family=fontfam)
        plot.SARC.R.SSB(mod, stock = NULL) #totals
        dev.off()
      }

      png(file.path(dir.res,"CV_SSB_Rec_F.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.cv(mod)
      dev.off()
    }
    png(file.path(dir.res,"F_byfleet.png"),width=10,height=10,units="in",res=res,family=fontfam)
    plot.fleet.F(mod)
    dev.off()
    for(i in 1:mod$env$data$n_stocks) for(r in 1:mod$env$data$n_regions) {
      png(file.path(dir.res,paste0("M_at_age_", sfns[i], "_", rfns[mod$input$data$stock_regions[i]],".png")),
        width=10,height=10,units="in",res=res,family=fontfam)
      plot.M(mod, stock = i, region = r)
      dev.off()
    }
    plot.tile.age.year(mod, type="selAA", do.png=TRUE, fontfam=fontfam, od=dir.res)
    plot.tile.age.year(mod, type="MAA", do.png=TRUE, fontfam=fontfam, od=dir.res)
    plot.tile.age.year(mod, type="NAA_devs", do.png=TRUE, fontfam=fontfam, od=dir.res)
    plot_q_prior_post(mod, do.png=TRUE, fontfam=fontfam, od=dir.res) #flag inside to plot if prior is being used. 
    plot_q(mod, do.png=TRUE, fontfam=fontfam, od=dir.res)
    if(!all(mod$env$data$Ecov_model == 0) & mod$is_sdrep) plot.ecov(mod, do.png=TRUE, fontfam=fontfam, od=dir.res, res=res)

    # PNG reference points -----------------
    dir.refpts <- file.path(dir.plots, "ref_points")
    dir.create(dir.refpts, showWarnings = FALSE)
    
    if(mod$env$data$n_stocks ==1 & mod$env$data$n_regions == 1) {
      png(file.path(dir.refpts,"SPR_targets_ave_plot.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SPR.table(mod, plot=TRUE)
      dev.off()
      png(file.path(dir.refpts,"SPR_targets_ave_table.png"),width=10,height=10,units="in",res=res,family=fontfam)
      plot.SPR.table(mod, plot=FALSE)
      dev.off()
      plot.annual.SPR.targets(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam)
      plot.yield.curves(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam, plot=TRUE)
      plot.yield.curves(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam, plot=FALSE)
    }
    if(mod$env$data$do_SPR_BRPs){
      if(mod$is_sdrep) plot.FXSPR.annual(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam)
    }
    if(mod$env$data$do_MSY_BRPs){
      if(mod$env$data$recruit_model %in% (3:4) & mod$is_sdrep) plot.MSY.annual(mod, od=dir.refpts, do.png=TRUE, fontfam=fontfam)
    }

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


  if(rmarkdown::pandoc_available()){
    if(table.type == "pdf"){
      cat(paste0("Making LaTeX/pdf tables: ", file.path(dir.res.tables, "wham_par_tables.pdf"), "\n"))
      par_tables_fn(mod, od = dir.res.tables, do.tex = T)
    }
    # if(table.type == "html"){
    #   cat("Making HTML tables.\n")
    #   par_tables_fn(mod, od = dir.res.tables, do.html = T)
    #   cat("Opening HTML tables in your default web-browser.\n")
    #   browseURL(file.path(dir.res.tables, "wham_par_tables.html"))
    # }
    if(out.type == "html"){
      cat("Making HTML document.\n")
      par_tables_fn(mod, od = dir.res.tables)
      make_html_figs_tables_fn(od = dir.main, do.html = T)
      html_file <- file.path(dir.main, "wham_figures_tables.html")
      if(browse) {
        cat("Opening HTML document in your default web-browser.\n")
        browseURL(html_file)
      }
      else cat(paste0("HTML document is \n", html_file, "\n"))
    }
  } else {
    #uses png output, and old html product that automatically opens in browser
    if(out.type == 'html'){
      if(browse) wham_html(dir.main = dir.main)
    }
    cat("document of parameter tables was not generated because rmarkdown could not find pandoc. Try plot_wham_output() from Rstudio.\n")
  }
}


