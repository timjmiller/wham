source("other_wham_functions.r")
source('wham_plots_tables.r')
RData.file = "mod.RData"
mod.name = "mod"
load(RData.file)
mod = get(mod.name)

wd = getwd()

plots.dir <- paste(getwd(),"/plots/", sep="")
system(paste0("mkdir ", plots.dir))

graphics.off()     # close any open windows
origpar = par()


cairo_pdf(paste0(plots.dir, "input_plots.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
#data
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

#diagnostics
cairo_pdf(paste0(plots.dir, "diagnostic_plots.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
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

#results
cairo_pdf(paste0(plots.dir, "results_plots.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
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

cairo_pdf(paste0(plots.dir, "spr_msy_plots.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
plot.SPR.table(mod)
plot.annual.SPR.targets(mod)
plot.FXSPR.annual(mod)
plot.SR.pred.line(mod)
plot.MSY.annual(mod)
plot.yield.curves(mod)
dev.off()

cairo_pdf(paste0(plots.dir, "retro_plots.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
plot.retro(mod, what = "SSB")
plot.retro(mod, what = "Fbar")
plot.retro(mod, what = "NAA")
dev.off()

cairo_pdf(paste0(plots.dir, "misc_plots.pdf"), family = "Times", height = 10, width = 10, onefile = TRUE)
plot_catch_at_age_consistency(mod)
plot_index_at_age_consistency(mod)
plot_catch_curves_for_catch(mod)
plot_catch_curves_for_index(mod)
dev.off()



